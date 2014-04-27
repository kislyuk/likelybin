#!/usr/bin/env perl

# MCUtils.pm: utility subroutine library for MCMC metagenomic binning methods
# Author: Andrey Kislyuk (kislyuk@gatech.edu)

# NB: Some of the functionality below duplicates that in the C core.

package MCUtils;
require 5.005;
my ($VERSION) = ('$Id$' =~ /,v\s+(\d+\S+)/o);

use strict;
use Exporter;
use IPC::Open2;
use Cwd;
use List::Util qw(min max sum reduce shuffle);
use File::Basename;
use File::Spec;
use File::Copy;
#use Bit::Vector;
use Class::Struct;
use AKUtils;

our @ISA = "Exporter";
our @methods = qw();
our %EXPORT_TAGS = (all => [@methods]);
Exporter::export_ok_tags('all');

local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH = $FSFind::LOG || *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

my (%kmer_vector_cache, %all_kmer_vector_cache);

1;

# Counts all k-mers of order k and below in all given sequences, gives their counts and frequencies
# Input: hash of sequences {seqname1 => seq1, seqname2 => seq2, ...}, maximum order of k-mers to count
sub getKmerCounts($$) {
	my ($seqs, $k) = @_;
	my (%counts, %freqs);
	foreach my $seqname (keys %$seqs) {
		my $nc_nucs;
		# pseudocounts (disabled)
		#foreach my $n1 (qw(A T G C)) { foreach my $n2 (qw(A T G C)) { $counts{$seqname}->{$n1.$n2}=1; }}
		foreach my $start (0..(length($$seqs{$seqname})-$k)) {
			my $mer = substr($$seqs{$seqname}, $start, $k);
			($nc_nucs++, next) unless $mer =~ /^[ATGC]+$/; # do not admit non-canonical nucleotides
			my $rc_mer = reverse($mer); $rc_mer =~ tr/ATGC/TACG/;
			$counts{$seqname}->{$mer}++;
			$counts{$seqname}->{$rc_mer}++;
		}
		my $kmer_count = sum(values %{$counts{$seqname}});
		foreach my $kmer (keys %{$counts{$seqname}}) {
			$freqs{$seqname}->{$kmer} = $counts{$seqname}->{$kmer} / $kmer_count;
		}
		for (my $n = $k-1; $n > 0; $n--) { # marginals
			foreach my $kmer (keys %{$counts{$seqname}}) {
				next if length($kmer) != $n+1;
				$counts{$seqname}->{substr($kmer, 0, $n)} += $counts{$seqname}->{$kmer};
				$freqs{$seqname}->{substr($kmer, 0, $n)} += $freqs{$seqname}->{$kmer};
			}
		}
		if ($nc_nucs/length($$seqs{$seqname}) > 0.05) {
			warn "Significant fraction of non-canonical nucleotides ($nc_nucs / ".length($$seqs{$seqname}).") in sequence $seqname";
		}
		#my $t; $t+=$counts{$seqname}->{$_} for 'A', 'T', 'G', 'C'; print "Total: $t\n";
		#foreach my $n1 (qw(A T G C)) { foreach my $n2 (qw(A T G C)) {print "$seqname\t$n1$n2\t".$counts{$seqname}->{$n1.$n2}."\n"; die if $counts{$seqname}->{$n1.$n2}==0;}}
	}
	return (\%counts, \%freqs);
}

# Returns [AA, AT, ...]
sub getKmerVector($) {
	my ($order) = @_;
	return $kmer_vector_cache{$order} if $kmer_vector_cache{$order};
	my @nucs = ('');
	for (1..$order) {
		# All combinations of A, T, G, C of length $order
		@nucs = map(($_.'A', $_.'T', $_.'G', $_.'C'), @nucs);
	}
	$kmer_vector_cache{$order} = \@nucs;
	return \@nucs;
}

# Returns [A, T, G, C, AA, AT, ...]
sub getAllKmerVector($) {
	my ($order) = @_;
	return $all_kmer_vector_cache{$order} if $all_kmer_vector_cache{$order};
	my @vector; my @nucs = ('');
	for (1..$order) {
		# All combinations of A, T, G, C of length $order
		@nucs = map(($_.'A', $_.'T', $_.'G', $_.'C'), @nucs);
		push(@vector, @nucs);
	}
	$all_kmer_vector_cache{$order} = \@vector;
	return \@vector;
}

# Computes second-order deviations from identities specified by the k-mer dimension redundancy model
sub checkModelDeviance($) {
	my ($M) = @_;
	return 0 if $$M[0]->{order} < 2;
	my $total_deviance;
	for (0..$#$M) {
		my $model = $$M[$_];
		my @deviances = ($$model{AG} - ($$model{A} - $$model{AA} - $$model{AC} - $$model{AT}),
			$$model{CC} - (1/2 - 2*$$model{A} + $$model{AA} + $$model{AC} + $$model{AT} - $$model{CA} - $$model{CG}),
			$$model{TA} - (2*$$model{AC} + $$model{AT} - 2*$$model{CA} - $$model{CG} + $$model{GC}),
			$$model{TC} - ($$model{A} - $$model{AA} - 2*$$model{AC} - $$model{AT} + $$model{CA} + $$model{CG} - $$model{GC}));
		$total_deviance += abs($_) for @deviances;
	}
	return $total_deviance;
}

# Compute aggregate model log probability
# log(P(S_i|M)) = log(P(S_i|M_1)) + log(sum_j(f_j * P(S_i|m_j) / P(S_i|m_1)))
# Model posterior probability using uninformative priors and uniform P(S):
# P(M|S) = P(S|M) = Prod_seqs(Sum_models( freq(model) * prob(seq|model)))
# With informative priors (e.g. Dirichlet), multiply this by P(model|expected distribution of models)
sub getModelLogProb($$) {
	my ($Model, $kmer_counts) = @_;
	my $logP;

	my $logModel;
	foreach my $source (0..$#$Model) {
		$$logModel[$source] = logModel($$Model[$source]);
	}

	foreach my $seqname (keys %$kmer_counts) {
		my @seq_logprobs;
		my ($max_logprob, $max_logprob_index);
		foreach my $source (0..$#$Model) {
			my $seq_logprob = getSeqLogProb($$kmer_counts{$seqname}, $$Model[$source], $$logModel[$source]);
			push(@seq_logprobs, $seq_logprob);
			$max_logprob = $seq_logprob if $max_logprob < $seq_logprob or not defined $max_logprob;
			$max_logprob_index = $source if $max_logprob == $seq_logprob;
		}
		my $prob_partsum;

		foreach my $source (0..$#$Model) {
			$prob_partsum += exp($seq_logprobs[$source] - $seq_logprobs[$max_logprob_index]);
		}
		$logP += $seq_logprobs[$max_logprob_index] + log($prob_partsum);
	}
	return $logP;
}

sub logModel($) {
	my ($model) = @_;
	my %log_model;
	foreach my $mer (@{getAllKmerVector($$model{order})}) {
		$log_model{$mer} = log($$model{$mer});
	}
	return \%log_model;
}

# log P(S_i|M_j) = log(P(AA)^count(AA)*.../P(A)^count(A)) = sum_dimers(log(P(dimer)^count(dimer) ...)) = sum(count(dimer)*log(P(dimer))...
sub getSeqLogProb($$;$) {
	my ($seq_counts, $model, $log_model) = @_;
	$log_model ||= logModel($model);

	my $log_prob;
	foreach my $mer (@{getKmerVector($$model{order})}) {
		$log_prob += $$log_model{$mer} * $$seq_counts{$mer};
	}
	if ($$model{order} > 1) { # Divide by marginals of the next highest order
		foreach my $mer (@{getKmerVector($$model{order}-1)}) {
			$log_prob -= $$log_model{$mer} * $$seq_counts{$mer};
		}
	}
	return $log_prob;
}

# Let p11, p12, ... be kmer probs for model 1, p21, p22, ... for model 2
# s11, s12, ... = standard deviations
# Define divergence D = (p11-p21)^2/(s11^2+s21^2) + (p12-p22)^2/(s12^2+s22^2) * 4 + ...
sub getDistribDivergence($$) {
	my ($frags, $settings) = @_;
	my $moments = getFragKmerMoments($frags, $settings);

	my $total_d;

	my @kmers = keys %{$$moments{(keys %$moments)[0]}};
	my @sources = keys %$moments;

	my $iter = AKUtils::combinations(2, \@sources);
	while (my ($source1, $source2) = $iter->()) {
		my $d;
		foreach my $kmer (@kmers) {
			(warn("unable to compute deviance for kmer $kmer"), next) unless $$moments{$source1}->{$kmer}->{mean} and $$moments{$source2}->{$kmer}->{mean};
			my $udiff2 = ($$moments{$source1}->{$kmer}->{mean} - $$moments{$source2}->{$kmer}->{mean}) ** 2;
			my $v2diff = $$moments{$source1}->{$kmer}->{stdev}**2 + $$moments{$source2}->{$kmer}->{stdev}**2;
			#(warn("unable to compute deviance for kmer $kmer"), next) if $v2diff == 0;
			$d += ($udiff2 / $v2diff) / (4**length($kmer));
			#print "kmer $kmer\tu1 = ".$$moments{$source1}->{$kmer}->{mean}."\tu2 = ".$$moments{$source2}->{$kmer}->{mean}."ud2=$udiff2 v2d=$v2diff sf=".(4**length($kmer))."\t\td=$d\n";
		}
#		$total_d += $d;
		# TODO: revisit the distance formula for N>2
		$total_d = $d if $total_d > $d or not defined $total_d;
#		warn "source $source1:$source2 d $d total_d $total_d\n";
	}
	return $total_d;
}

sub getTwoSeqDivergence($$) {
	my ($seqs, $settings) = @_;
	my %frags;
	die("Invalid arguments: expected exactly 2 sequences") if keys(%$seqs) != 2;
	$$settings{frag_length} = 400 if not defined $$settings{frag_length};
	die("Invalid frag_length") unless $$settings{frag_length} >= 40 and $$settings{frag_length} <= 2000;
	$$settings{frag_num} = 500 if not defined $$settings{frag_num};
	die("Invalid frag_num") unless $$settings{frag_num} >= 10 and $$settings{frag_num} <= 100000;
	$$settings{chain_order} = 3 if not defined $$settings{chain_order};
	die("Invalid chain_order") unless $$settings{chain_order} >= 1 and $$settings{chain_order} <= 6;

	foreach my $seqname (keys %$seqs) {
		foreach my $i (1..$$settings{frag_num}) {
			my $start_pos = int(rand(length($$seqs{$seqname}) - $$settings{frag_length}));
			my $frag = substr($$seqs{$seqname}, $start_pos, $$settings{frag_length});
			$frags{"$seqname:$start_pos"} = $frag;
		}
	}
	return MCUtils::getDistribDivergence(\%frags, $settings);
}

# Compute mean and standard deviation of kmer distributions across fragments within sources
sub getFragKmerMoments($$) {
	my ($frags, $settings) = @_;
	my ($counts, $freqs) = MCUtils::getKmerCounts($frags, $$settings{chain_order});
	my %moments;

	my %source_kmer_freqs;
	foreach my $frag_name (keys %$frags) {
		$frag_name =~ /^(.+):\d+$/ or warn;
		my $source_name = $1;
		foreach my $kmer (keys %{$$freqs{$frag_name}}) {
			# WARNING: assumes equal-length fragments
			push(@{$source_kmer_freqs{$source_name}->{$kmer}}, $$freqs{$frag_name}->{$kmer});
		}
	}

	foreach my $source (keys %source_kmer_freqs) {
		foreach my $kmer (keys %{$source_kmer_freqs{$source}}) {
			my $mean = sum(@{$source_kmer_freqs{$source}->{$kmer}}) / @{$source_kmer_freqs{$source}->{$kmer}};
			my $stdev = AKUtils::stdev($source_kmer_freqs{$source}->{$kmer});

			my $outlier_frac = 0.1;
			my @culled_freqs = sort @{$source_kmer_freqs{$source}->{$kmer}};
			my $outlier_bandwidth = int(@culled_freqs * $outlier_frac/2);
			@culled_freqs = @culled_freqs[$outlier_bandwidth .. @culled_freqs-$outlier_bandwidth];

			$moments{$source}->{$kmer} = {mean => $mean, stdev => $stdev,
				lo_outlier_bd => $culled_freqs[0], hi_outlier_bd => $culled_freqs[$#culled_freqs]};
		}
	}
	return \%moments;
}

sub runCompostBin($$) {
	my ($frags, $settings) = @_;
	
	$$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/".(split("::",(caller(0))[3]))[1].".$$.XXXXX", CLEANUP => !($$settings{keep}));

	my $mfa_file = "$$settings{tempdir}/cb_in.mfa";
	logmsg "Preparing to run CompostBin on $mfa_file...";

	open(OUT, '>', $mfa_file) or die;
	foreach my $seqname (shuffle keys %$frags) {
		print OUT ">$seqname\n";
		my $seq = $$frags{$seqname};
		$seq =~ s/(.{80})/$1\n/g;
		$seq .= "\n" unless $seq =~ /\n$/;
		print OUT $seq;
	}
	close OUT;

	system("cd '$$settings{tempdir}'; fvec '$mfa_file'");
#	die("Error running fvec: $!") if $?;
#cat data.mat|wc -l|perl -e'my $l=<>; print("0.5\n" x $l);' > weights.mat


#	CompostBin-0.0.1/fvec/fvec /tmp/run_compostbin.pl.490.Jig_b/markers/join.fa


}
