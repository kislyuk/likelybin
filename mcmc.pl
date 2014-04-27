#!/usr/bin/env perl

# mcmc_binner: metagenomic binning utility
# Author: Andrey Kislyuk (kislyuk@gatech.edu)

package MCMC;
my ($VERSION) = ('$Id: $' =~ /,v\s+(\d+\S+)/o);

my %settings = (
	logfile => '-',
	datadir => $FindBin::RealBin,
	max_search_steps => 40000, # NOTE: test_binner overrides this
	num_sources => 2, # NOTE: test_binner overrides this
	chain_order => 3, # NOTE: test_binner overrides this
	model_var => 5e-5, # NOTE: mcmc-core overrides this
#	model_var => 1e-7,
	freq_var => 5e-5,
	max_deviance => 1e-10,
	nullmat_file => "mats/nulls.pl",
	mcmc_core_exec => "mcmc-core",
);
my %stats;

use strict;
use Getopt::Long;
use File::Temp ('tempdir');
use File::Path;
use File::Spec;
use File::Copy;
use File::Basename;
use List::Util qw(min max sum);
use Class::Struct;
use FindBin;
use IO::Select;
use lib "$FindBin::RealBin/lib";
$ENV{PATH} .= ":$FindBin::RealBin/c";
use AKUtils;
use MCUtils;
#use PDL;
#use Math::BigFloat;

sub main();

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH = $FSFind::LOG || *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

my $usage = "$0 input.mfa [options]
Options can be abbreviated. Options are:
-num_sources=N	Number of source species in the dataset
-max_search_steps=N	Number of search steps to take in the MCMC simulation
-chain_order=N	K-mer order (valid values are 2 through 5)
-keep		Keep all intermediate files
";

exit(main());

sub main() {
	die("Usage: $usage") if @ARGV < 1;

	my ($input_file_full) = @ARGV;
	$input_file_full = File::Spec->rel2abs($input_file_full);
	my ($input_seq_file, $input_seq_dir) = fileparse($input_file_full);

	my @cmd_options = ('outfile=s', 'num_sources=i', 'model_logfile=s', 'max_search_steps=i', 'chain_order=i', 'keep', 'num_threads=i',
		'model_var=s');
	GetOptions(\%settings, @cmd_options) or die("Error parsing command line switches. Check the input options. Usage:".$usage);

	my @uinfo; # = getpwuid($>);
	warn "$0: Version ".$VERSION." started ".localtime()." by ".$uinfo[0]."\n";

	$settings{outfile} ||= "$input_seq_dir/$input_seq_file.binning";
	$settings{logfile} ||= "$settings{outfile}.log";
	open(FH, '>>', $settings{outfile}) or die "Unable to open file $settings{outfile} for writing: $!"; close FH;
	if ($settings{logfile} eq '-') {
		$FSFind::LOG = *STDERR;
	} else {
		open($FSFind::LOG, '>', $settings{logfile}) or die "Unable to open file $settings{logfile} for writing: $!";
		warn "Logging to $settings{logfile}. Use \"-logfile=-\" to direct logging to stderr.\n";
	}
	$settings{tempdir} = tempdir($settings{tempdir} or File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($settings{keep}));
	logmsg "Temporary directory is $settings{tempdir}";

	logmsg "Loading input sequences";
	my $seqs = AKUtils::readMfa($input_file_full);
	my ($kmer_counts, $kmer_freqs) = MCUtils::getKmerCounts($seqs, $settings{chain_order});

	my ($best_model, $best_model_ll);

	open(CORE_IN, '>', "$settings{tempdir}/core_in") or die;
	foreach my $seqname (keys %$kmer_counts) {
		print CORE_IN "$seqname\n";
		foreach my $mer (@{MCUtils::getAllKmerVector($settings{chain_order})}) {
			print(CORE_IN ($$kmer_counts{$seqname}->{$mer} or 0)." ");
		}
		print CORE_IN "\n";
	}
	close CORE_IN;

	open(HINT, '>', "$settings{tempdir}/core_hint_model") or die;
	for (1..$settings{num_sources}) {
		print(HINT (1/$settings{num_sources})."\n");
		foreach my $mer (@{MCUtils::getAllKmerVector($settings{chain_order})}) {
			# TODO: counts or freqs? (counts weigh by sequence length, freqs do not)
			my $avg_freq;
			$avg_freq += $$kmer_freqs{$_}->{$mer} for keys %$kmer_counts;
			$avg_freq /= scalar(keys %$kmer_counts);
			print HINT "\t$avg_freq\n";
		}
	}
	close HINT;

	my $invoke_string = "$settings{mcmc_core_exec} -o $settings{chain_order} -s $settings{max_search_steps} -n $settings{num_sources}";
	$invoke_string .= " -m $settings{model_var} -f $settings{freq_var}";
	if ($settings{model_logfile}) { $invoke_string .= " -l $settings{model_logfile}"; }

$invoke_string .= " -h $settings{tempdir}/core_hint_model";

#$invoke_string .= " ";# > $settings{tempdir}/core_out";
#	logmsg $invoke_string;
#	system($invoke_string);
#	die("Error running $settings{mcmc_core_exec}") if $?;
#	my $pid = open2(\*CORE_OUT, \*CORE_IN, $invoke_string);

	$settings{num_cpus} = AKUtils::getNumCPUs();
	$settings{num_threads} ||= $settings{num_cpus};

	# Launch the cores
	my $th_set = new IO::Select();
	my %thread_ids;
	foreach my $thread_id (0..$settings{num_threads}-1) {
		my $th;
		my $my_invoke_str = "$invoke_string -i $thread_id < $settings{tempdir}/core_in";
		logmsg "Launching thread $thread_id via \"$my_invoke_str\"";
		open($th, '-|', $my_invoke_str);
		$th_set->add($th);
		$thread_ids{$th} = $thread_id;
	}

	# Wait on them and process each output when it's ready
#my @lls;
	while (my @ready_rhs = $th_set->can_read) {
		foreach my $CORE_OUT (@ready_rhs) {
			# NB: must proceed to close CORE_OUT, otherwise it will come back in @ready_rhs
			my $core_reported_ll = <$CORE_OUT>;
			chomp $core_reported_ll;
#push(@lls, $core_reported_ll);
#print "core reported LL: $core_reported_ll\n";
			local $/;
			my $core_output = <$CORE_OUT>;

			my $model;
			foreach my $m (split(/\n(?=[^\t])/, $core_output)) {
				my @mer_freqs = split(/\n\t/, $m);
				my $model_freq = shift @mer_freqs;
				my %freq_hash;

				foreach my $kmer (@{MCUtils::getAllKmerVector($settings{chain_order})}) {
					$freq_hash{$kmer} = shift @mer_freqs;
					#die unless $freq_hash{$kmer} > 0; TODO: need proper checks
				}

				push(@$model, \%freq_hash);
				$$model[$#$model]->{freq} = $model_freq;
				$$model[$#$model]->{order} = $settings{chain_order};
				$$model[$#$model]->{loglik} = $core_reported_ll;
			}

			$th_set->remove($CORE_OUT);
			close $CORE_OUT;

			$best_model_ll = $core_reported_ll if $best_model_ll < $core_reported_ll or not defined $best_model_ll;
			if ($best_model_ll == $core_reported_ll) {
				$best_model = $model;
				move("$settings{model_logfile}.$thread_ids{$CORE_OUT}", $settings{model_logfile}) if $settings{model_logfile};
			} else {
				unlink("$settings{model_logfile}.$thread_ids{$CORE_OUT}") if $settings{model_logfile};
			}
			die("core reported inappropriate LL $core_reported_ll") if $core_reported_ll > -10000;
		}
	}
#print "core reported LLs: ".join(",", @lls)."\n";
=head1
	open(CORE_OUT, '<', "$settings{tempdir}/core_out");
	foreach my $m (split(/\n(?=[^\t])/, $core_output)) {
		my @mer_freqs = split(/\n\t/, $m);
		my $model_freq = shift @mer_freqs;
		my %freq_hash;

		foreach my $kmer (@{MCUtils::getAllKmerVector($settings{chain_order})}) {
			$freq_hash{$kmer} = shift @mer_freqs;
			#die unless $freq_hash{$kmer} > 0; TODO: need proper checks
		}

		push(@model, \%freq_hash);
		$model[$#model]->{freq} = $model_freq;
		$model[$#model]->{order} = $settings{chain_order};
	}
	close CORE_OUT;
=cut
	printModel($best_model, \*STDOUT);
	my $predictions = estimateSources($kmer_counts, $best_model);
	return 0;
}

sub printModelRecord($$) {
	my ($record, $settings) = @_;
	require Data::Dumper;
	open(OUT, '>', $$settings{model_logfile}) or die "Unable to open file $$settings{model_logfile} for writing: $!";
	print OUT Data::Dumper::Dumper($record);
	close OUT;
}

sub estimateSources($$) {
	my ($counts, $model) = @_;
	die("Internal error") unless defined $counts and @$model > 0;
	my %sources;
	open(OUT, '>', $settings{outfile}) or die "Unable to open file $settings{outfile} for writing: $!";
open(OUT2, '>', "$settings{outfile}.allprobs") or die;
	foreach my $seqname (keys %$counts) {
print OUT2 "$seqname";
		my ($best_model_i, $best_model_logprob);
		foreach my $i (0..$#$model) {
			my $logprob = MCUtils::getSeqLogProb($$counts{$seqname}, $$model[$i])
				+ log($$model[$i]->{freq});
print OUT2 "\t$logprob";
#warn "logprob $logprob, best $best_model_logprob\n";
			$best_model_logprob = $logprob if $best_model_logprob < $logprob or not defined $best_model_logprob;
			$best_model_i = $i if $best_model_logprob == $logprob;
		}
print OUT2 "\n";
#warn "selected best model $best_model_i ($best_model_logprob) ($seqname)\n";
		print OUT "$best_model_i\t$seqname\n";
	}
close OUT2;
	close OUT;
}

sub printModel($;$) {
	my ($model, $FH) = @_;
	$FH ||= \*STDOUT;
	print $FH "Model LL = $$model[0]->{loglik}\n";
	foreach my $i (0..$#$model) {
		print $FH "Source $i: ID $$model[$i]->{id}, order $$model[$i]->{order}, frequency ".sprintf("%.6f", $$model[$i]->{freq})."\n";
		print $FH join("\t", ("1 2>", qw(A T G C)))."\n";
		foreach my $nuc1 (qw(A T G C)) {
			print $FH "$nuc1\t";
			foreach my $nuc2 (qw(A T G C)) {
				print $FH sprintf("%.4f", $$model[$i]->{$nuc1.$nuc2})."\t";
			}
			print $FH "\n";
		}
		print $FH '(%G+C = '.($$model[$i]->{G}+$$model[$i]->{C}).")\n";
	}
}
