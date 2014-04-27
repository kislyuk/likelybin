#!/usr/bin/env perl

# TestBinner: test harness for binning applications
# Author: Andrey Kislyuk (kislyuk@gatech.edu)

package TestBinner;
my ($VERSION) = ('$Id: $' =~ /,v\s+(\d+\S+)/o);

my %settings = (
	frag_num => 500,
	frags_per_source => 'yes',
	frag_length => 800,
#	chain_order => 4,
	chain_order => 3,
	max_search_steps => 20000,
);
my %stats;

use strict;
use Getopt::Long;
use File::Temp ('tempdir');
use File::Path;
use File::Spec;
use File::Copy;
use File::Basename;
use List::Util qw(min max sum shuffle);
use Class::Struct;
use FindBin;
use lib "$FindBin::RealBin/lib";
$ENV{PATH} .= ":$FindBin::RealBin";
use AKUtils;
use MCUtils;
use Data::Dumper;
use Sys::Hostname;

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH = $FSFind::LOG || *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

exit(main());

sub main() {
	die("Usage: $0 input.mfa") if @ARGV < 1;

	my @cmd_options = ('plotfile=s', 'keep', 'frag_num=i', 'frag_length=i', 'model_var=s', 'chain_order=i', 'frags_per_source=s');
	GetOptions(\%settings, @cmd_options) or die;

	my $input_file_full = shift @ARGV;
	$input_file_full = File::Spec->rel2abs($input_file_full);
	my ($input_seq_file, $input_seq_dir) = fileparse($input_file_full);

	$settings{tempdir} = tempdir($settings{tempdir} or File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($settings{keep}));
	logmsg "$0: Temporary directory is $settings{tempdir}";

	my $seqs = AKUtils::readMfa($input_file_full);
	die("Fewer than 2 sources supplied") if keys(%$seqs) < 2;

	my $rates = {}; # relative sampling rates of the sequences
	if (-f "$input_file_full.rates") {
		die("To use rates, disable frags_per_source") if $settings{frags_per_source} eq 'yes';
		open(FH, '<', "$input_file_full.rates") or die;
		my $r = <FH>; close FH; chomp $r;
		my @rates = split /:/, $r;
		open(FH, '<', $input_file_full) or die;
		while (<FH>) {
			if (/^\>\s*(.+)/) {
				$$rates{$1} = shift @rates;
			}
		}
		close FH;
		die if @rates > 0;
	} else {
		foreach my $seqname (keys %$seqs) {
			# Sample proportionally to length
			$$rates{$seqname} = length($$seqs{$seqname});
		}
	}

	print "Rates used:\n";
	foreach my $seqname (keys %$rates) {
		print "$$rates{$seqname}\t$seqname\n";
	}

	my %frags;
	foreach my $seqname (keys %$seqs) {
		die("Sequence $seqname too short") if length($$seqs{$seqname}) < $settings{frag_length} * $settings{frag_num};
		my $num_frags_this_seq = int($settings{frag_num} * scalar(keys %$seqs) * $$rates{$seqname} / sum(values %$rates));
		$num_frags_this_seq = $settings{frag_num} if $settings{frags_per_source} eq 'yes';
		foreach my $i (1..$num_frags_this_seq) {
			my $start_pos = int(rand(length($$seqs{$seqname}) - $settings{frag_length}));
			my $frag = substr($$seqs{$seqname}, $start_pos, $settings{frag_length});
			#print "$seqname $i $start_pos ".length($frag)."\n";
			$frags{"$seqname:$start_pos"} = $frag;
		}
	}
	open(OUT, '>', "$settings{tempdir}/test.fna") or die;
	foreach my $seqname (shuffle keys %frags) {
		print OUT ">$seqname\n";
		my $seq = $frags{$seqname};
		$seq =~ s/(.{80})/$1\n/g;
		$seq .= "\n" unless $seq =~ /\n$/;
		print OUT $seq;
	}
	close OUT;

#	$ENV{PATH} .= ":$FindBin::RealBin/CompostBin-0.0.1/fvec";
#	MCUtils::runCompostBin(\%frags, \%settings);

	# Run the binner
	my $binner_opts = "-num_sources=".scalar(keys %$seqs);
	$binner_opts .= " -max_search_steps=$settings{max_search_steps}";
	$binner_opts .= " -chain_order=$settings{chain_order}";
	$binner_opts .= " -model_logfile=$settings{tempdir}/model_log" if $settings{plotfile};
	$binner_opts .= " -model_var=$settings{model_var}" if $settings{model_var};
	$binner_opts .= " -keep" if $settings{keep};
#	my $invoke_string = "perl -d:DProf ";
	my $invoke_string = "$FindBin::RealBin/mcmc.pl $settings{tempdir}/test.fna $binner_opts";
	logmsg $invoke_string;
	system($invoke_string);
	die("Error running binner") if $?;

	# Retrieve source assignments made by binner
	my %binner_assignments;
	open(IN, '<', "$settings{tempdir}/test.fna.binning") or die;
	while (<IN>) {
		if (/^(\d+)\t(.+)$/) {
			my ($source_id, $seq_name) = ($1, $2);
			push(@{$binner_assignments{$source_id}}, $seq_name);
		} else { die "Unable to parse mcmc.pl output line: \"$_\""; }
	}
	close IN;

my %binner_probs;
open(IN, '<', "$settings{tempdir}/test.fna.binning.allprobs") or die;
while (<IN>) {
	chomp;
	if (/^(.+?):(\d+)\t([^\t]+)\t([^\t]+)/) {
		my ($seq_name, $seq_pos, $prob1, $prob2) = ($1, $2, $3, $4);
		push(@{$binner_probs{$seq_name}}, [$prob1, $prob2]);
	} else { die "Unable to parse mcmc.pl output line: \"$_\""; }
}
close IN;
open(OUT, '>', "/storage2/kislyuk/tmp/cloud.".hostname());
my $i;
foreach my $seq_name (keys %binner_probs) {
	$i++;
	foreach my $prob_set (@{$binner_probs{$seq_name}}) {
		print OUT "$i\t".join("\t", @$prob_set)."\n";
	}
}
close OUT;



	# TODO: there is a less brute force way to do this
	my @seq_permutations = AKUtils::permutations(keys %$seqs);
	my ($best_error, $best_assignment);
	foreach my $assignment (@seq_permutations) {
		my $error;
		print "Computing accuracy of assignment:\n";
		foreach my $source_id (0..$#$assignment) {
			my $correct_assigns;
			foreach my $seq_name (@{$binner_assignments{$source_id}}) {
				$seq_name =~ /^(.+?):(\d+)$/ or die;
				$correct_assigns++ if $1 eq $$assignment[$source_id];
			}
			print "$source_id => $$assignment[$source_id]: $correct_assigns correct\n";
			$error += $settings{frag_num} - $correct_assigns;
		}
		$best_error = $error if $best_error > $error or not defined $best_error;
		$best_assignment = $assignment if $error == $best_error;
	}
	print "Best assignment ($best_error of ".($settings{frag_num}*scalar(keys %$seqs))." incorrect):\n";
	foreach my $source_id (0..$#$best_assignment) {
		print "$source_id => $$best_assignment[$source_id]\n";
	}
	$stats{accuracy} = sprintf("%.4f", 1 - $best_error/($settings{frag_num}*scalar(keys %$seqs)));
	print "Accuracy = $stats{accuracy}\n";

	print "Best unbiased fit accuracy (fragments only) = ".estimateMaxAccuracy(\%frags, \%frags, $rates, \%settings)." (order $settings{chain_order})\n";
	print "Best unbiased fit accuracy (entire genomes) = ".estimateMaxAccuracy($seqs, \%frags, $rates, \%settings)." (order $settings{chain_order})\n";

	my $d = MCUtils::getDistribDivergence(\%frags, \%settings);
	print "Distribution divergence = $d\n";

	plotVisuals($seqs, \%frags, $rates, "$settings{tempdir}/model_log", \%settings) if $settings{plotfile};
#	print "Analyzing profiler data...\n";
#	system("dprofpp -u");
	return 0;
}

# NB: this relies on the evaluation model in mcmc-core remaining the same as in MCUtils
sub estimateMaxAccuracy($$$$) {
	my ($source_seqs, $frags, $rates, $settings) = @_;
	my %model_seqs; # input on which to compute k-mer stats

	my ($true_counts, $true_Model) = getTrueModel($source_seqs, $rates, $settings);

	my $counts = MCUtils::getKmerCounts($frags, $$settings{chain_order});

	my $max_tps;
	foreach my $seqname (keys %$counts) {
		$seqname =~ /^(.+):\d+$/;
		my $true_source = $1;
		my ($best_sourcename, $best_source_logprob);
		foreach my $sourcename (keys %$true_Model) {
			my $logprob = MCUtils::getSeqLogProb($$counts{$seqname}, $$true_Model{$sourcename})
				+ log($$true_Model{$sourcename}->{freq});
			$best_source_logprob = $logprob if $best_source_logprob < $logprob or not defined $best_source_logprob;
			$best_sourcename = $sourcename if $best_source_logprob == $logprob;
		}
		$max_tps++ if $true_source eq $best_sourcename;
	}

	return sprintf("%.4f", $max_tps/scalar(keys %$frags));
}

sub plotVisuals($$$$$) {
	my ($seqs, $frags, $rates, $model_logfile, $settings) = @_;
	logmsg "Generating convergence plots...";
	open(IN, '<', $model_logfile) or (logmsg("Could not open convergence log\n"),return);
	my $model_log;
	{	no strict 'vars';
		local $/ = undef;
		my $d = <IN>;
		$model_log = eval $d;
		use strict 'vars';
	}
	close IN;
	my $sp_start_timestamp = $$model_log{sp_start_timestamp};
	delete $$model_log{sp_start_timestamp};

	my ($target_counts, $target_freqs) = MCUtils::getKmerCounts($seqs, $settings{chain_order});

	my @rgb_palette = ([255,51,0], [51,0,153], [204,0,102], [102,153,0], [102,0,0], [153,204,51], [255,102,102], [102,102,204]);
	my @palette; push(@palette, join(" ", ($$_[0]/255, $$_[1]/255, $$_[2]/255))) for @rgb_palette;

	logmsg "Printing to $$settings{plotfile}...";
	open(OUT, '>', $$settings{plotfile}) or (warn("Unable to open file $$settings{plotfile} for writing: $!"), return);
	print OUT "figure;\n";
#	my @mers = qw(A AA AC AT AG CA CC);
	my @mers = qw(A AA AC AT CA CG GC);
	my @legends;
#	my @timestep_vector; foreach my $timestep (sort {$a <=> $b} keys %$model_log) { push(@timestep_vector, $timestep); }
	print OUT "timesteps=[".join(" ", (sort {$a <=> $b} keys %$model_log))."];\n";

	my $max_range;
	foreach my $i (1..$#mers) { # dimers only
		my ($dimer_min, $dimer_max);
		foreach my $source (0..$#{$$model_log{0}->{models}}) {
			foreach my $timestep (sort {$a <=> $b} keys %$model_log) {
				my $model = $$model_log{$timestep}->{models}->[$source];
				$dimer_min = $$model{$mers[$i]} if $dimer_min > $$model{$mers[$i]} or not defined $dimer_min;
				$dimer_max = $$model{$mers[$i]} if $dimer_max < $$model{$mers[$i]} or not defined $dimer_max;
			}
		}
		my $range = $dimer_max - $dimer_min;
		$max_range = $range if $max_range < $range or not defined $max_range;
	}


	foreach my $i (0..$#mers) {
		@legends = ();
		my $mer = $mers[$i];
		print OUT "subplot(3, 6, [".(2*$i+1)." ".(2*$i+2)."]); hold on\n";
		print OUT "title('$mer: variance=$settings{model_var}, steps=$settings{max_search_steps}, acc=$stats{accuracy}'); xlabel('Timestep'); ylabel('Frequency')\n";
		print OUT "plot([$sp_start_timestamp $sp_start_timestamp], [0 1], ':k')\n";
		my $j=0;
		foreach my $seq (keys %$seqs) {
			print OUT "plot([1, $settings{max_search_steps}], [$$target_freqs{$seq}->{$mer}, $$target_freqs{$seq}->{$mer}], ':r')\n";
			print OUT "set(findobj('Color', 'r'), 'Color', [$palette[$j]]);\n";
			push(@legends, "'".substr($seq, 0, 15)."'");
			$j++;
		}
		my ($ymin, $ymax);
		foreach my $source (0..$#{$$model_log{0}->{models}}) {
			my (@x_vector, @y_vector);
			foreach my $timestep (sort {$a <=> $b} keys %$model_log) {
				my $model = $$model_log{$timestep}->{models}->[$source];
				push(@x_vector, $timestep);
				push(@y_vector, sprintf("%.6f", $$model{$mer}));
				$ymin = $$model{$mer} if $ymin > $$model{$mer} or not defined $ymin;
				$ymax = $$model{$mer} if $ymax < $$model{$mer} or not defined $ymax;
			}
			print OUT "plot(timesteps, [".join(" ", @y_vector)."], 'r')\n";

			print OUT "set(findobj('Color', 'r'), 'Color', [$palette[$source]]);\n";
			push(@legends, "'Model $source'");
		}
		my $yavg = ($ymin+$ymax)/2; $ymin = $yavg - $max_range/2; $ymax = $yavg + $max_range/2;
		print OUT "ylim([".sprintf("%.2f", $ymin)." ".sprintf("%.2f", $ymax)."])\n" if length($mer) == 2;
		s/_/\\_/ for @legends;
		print OUT "legend(".join(",", @legends).", 'Location', 'Best')\n" if length($mer) == 1; # first plot only
	}

	print OUT "subplot(3, 6, [15 16]); hold on;\n";
	print OUT "title('Moments of source fragment 1, 2-mer distributions');\n";

	my $moments = MCUtils::getFragKmerMoments($frags, $settings);

	my $i;
	foreach my $source (keys %$moments) {
		$i++;
		my (@means, @stdevs, @lo_outlier_bds, @hi_outlier_bds);
#		foreach my $kmer (keys %{$source_kmer_freqs{$source}}) {
		foreach my $kmer (@mers) {
=head1
			my $outlier_frac = 0.1;
			my @culled_freqs = sort @{$source_kmer_freqs{$source}->{$kmer}};
			my $outlier_bandwidth = int(@culled_freqs * $outlier_frac/2);
			@culled_freqs = @culled_freqs[$outlier_bandwidth .. @culled_freqs-$outlier_bandwidth];
=cut
			push(@means, $$moments{$source}->{$kmer}->{mean});
			push(@stdevs, $$moments{$source}->{$kmer}->{stdev});
			push(@lo_outlier_bds, $$moments{$source}->{$kmer}->{lo_outlier_bd});
			push(@hi_outlier_bds, $$moments{$source}->{$kmer}->{hi_outlier_bd});
		}
		print OUT "errorbar_x([".join(" ", @means)."], [".join(" ", 1..@means)."], [".join(" ", @lo_outlier_bds)."], [".join(" ", @hi_outlier_bds)."], 'r.');\n";
		print OUT "set(findobj('Color', 'r'), 'Color', [$palette[$i]]);\n";
	}

	print OUT "set(gca, ['y' 'TickLabel'], {'".join("', '", ('', @mers, ''))."'});\n";

	if (keys(%$moments) == 2) { # 2 sources
		print OUT "subplot(3, 6, 17); hold on;\n";
		print2SourceScatterplot(\*OUT, $frags, $rates, \@palette, $settings);
	}

	print OUT "subplot(3, 6, 18); hold on;\n";
	print OUT "title('Model log likelihood'); xlabel('Timestep'); ylabel('Log P(data | best model)')\n";
	my $j=0;
	my (@x_vector, @y_vector);
	foreach my $timestep (sort {$a <=> $b} keys %$model_log) {
		push(@x_vector, $timestep);
		push(@y_vector, sprintf("%.6f", $$model_log{$timestep}->{loglik}));
	}
	print OUT "plot(timesteps, [".join(" ", @y_vector)."], 'k')\n";
	print OUT "plot([$sp_start_timestamp $sp_start_timestamp], [$y_vector[0] $y_vector[$#y_vector]], ':k')\n";
	print OUT "xlim([0 $settings{max_search_steps}]);\n";

	close OUT;
}

# Get a model representing the actual average kmer frequencies found in the source fragments (or full-length sources)
sub getTrueModel($$$) {
	my ($source_seqs, $rates, $settings) = @_;
	my %model_seqs; # input on which to compute k-mer stats

	foreach my $seqname (keys %$source_seqs) {
		if ($seqname =~ /^(.+):\d+$/) {
			$model_seqs{$1} .= $$source_seqs{$seqname};
		} else {
			$model_seqs{$seqname} = $$source_seqs{$seqname};
		}
	}

	my ($true_counts, $true_Model) = MCUtils::getKmerCounts(\%model_seqs, $$settings{chain_order});

	foreach my $seqname (keys %model_seqs) {
		$$true_Model{$seqname}->{freq} = $$rates{$seqname} / sum(values %$rates);
		$$true_Model{$seqname}->{order} = $$settings{chain_order};
	}
	return ($true_counts, $true_Model);
}

# NB: the scatterplot only makes sense if sequences are all of the same length!
sub print2SourceScatterplot($$$$$) {
	my ($FH, $frags, $rates, $palette, $settings) = @_;
	print OUT "title('Log likelihoods of sequence fragments with true models'); xlabel('LL(seq|M1)'); ylabel('LL(seq|M2)')\n";

	my $counts = MCUtils::getKmerCounts($frags, $$settings{chain_order});
	my ($true_counts, $true_Model) = getTrueModel($frags, $rates, $settings);

	my $source1 = (keys %$true_Model)[0];
	my $source2 = (keys %$true_Model)[1];

	my (@x, @y, @colors);
	foreach my $seqname (keys %$counts) {
		$seqname =~ /^(.+):\d+$/;
		my $true_source = $1;

		my $lp1 = MCUtils::getSeqLogProb($$counts{$seqname}, $$true_Model{$source1})
			+ log($$true_Model{$source1}->{freq});
		my $lp2 = MCUtils::getSeqLogProb($$counts{$seqname}, $$true_Model{$source2})
			+ log($$true_Model{$source2}->{freq});

		push(@x, sprintf("%.6f", $lp1));
		push(@y, sprintf("%.6f", $lp2));
		push(@colors, $true_source eq $source1 ? 'c1' : 'c2');
	}
	print $FH "c1=[$$palette[0]]; c2=[$$palette[1]];\n";
	print $FH "scatter([".join(" ", @x)."], [".join(" ", @y)."], 5, [".join("; ", @colors)."]);\n";
}
