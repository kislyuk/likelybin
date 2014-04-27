#!/usr/bin/env perl

# AKUtils.pm: utility subroutine library
# Author: Andrey Kislyuk (kislyuk@gatech.edu)

package AKUtils;
require 5.005;
my ($VERSION) = ('$Id$' =~ /,v\s+(\d+\S+)/o);

use strict;
use Exporter;
use IPC::Open2;
use Cwd;
use List::Util qw(min max sum reduce shuffle);
use File::Basename;
use File::Spec;
use File::Temp ('tempdir');
use Storable ('dclone'); # for deep copying
#use Bit::Vector;
#use Data::Dumper; # useful for printing deep structures, e.g. to print %hash: print Data::Dumper->new([\%hash],[qw(hash)])->Indent(3)->Quotekeys(0)->Dump;
#use Fatal qw(:void open);
# Useful modules to look at:
# FindBin DynaLoader IO::Poll IO::Select IO::Socket Math::* Net::*
# algorithm::diff archive::tar authen::pam compress::* datemanip crypt::* event::* soap-lite svg term::* BerkeleyDB
# Perl tail recursion: @_=(arg1, arg2, arg3); goto &subroutine_name;

# TODO: allpairs method: all adjacent pairs in an array, also implement with checking for in-place modification
# TODO: seqwiz -> AKGenomeUtils (how to maintain processivity?)

local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };

our @ISA = "Exporter";
our @methods = qw(argmax compactArray concatFiles curSub filesDiffer flatten fullPathToExec nanmean safeGlob indexMfa getSeq getSubSeq printSeqsToFile lruQueue readPSL getFreeMem totalSize isFasta readMfa allIndexes startStopPos nanmin nanmax numEltsAtDepth alnum permuteRanges loadCSV intersectLength gaussian);
our %EXPORT_TAGS = (all => [@methods]);
Exporter::export_ok_tags('all');

our %codon_table = (
	'TTT'=>'F', 'TTC'=>'F', 'TTA'=>'L', 'TTG'=>'L',
	'CTT'=>'L', 'CTC'=>'L', 'CTA'=>'L', 'CTG'=>'L',
	'ATT'=>'I', 'ATC'=>'I', 'ATA'=>'I', 'ATG'=>'M',
	'GTT'=>'V', 'GTC'=>'V', 'GTA'=>'V', 'GTG'=>'V',
	'TCT'=>'S', 'TCC'=>'S', 'TCA'=>'S', 'TCG'=>'S',
	'CCT'=>'P', 'CCC'=>'P', 'CCA'=>'P', 'CCG'=>'P',
	'ACT'=>'T', 'ACC'=>'T', 'ACA'=>'T', 'ACG'=>'T',
	'GCT'=>'A', 'GCC'=>'A', 'GCA'=>'A', 'GCG'=>'A',
	'TAT'=>'Y', 'TAC'=>'Y', 'TAA'=>'*', 'TAG'=>'*',
	'CAT'=>'H', 'CAC'=>'H', 'CAA'=>'Q', 'CAG'=>'Q',
	'AAT'=>'N', 'AAC'=>'N', 'AAA'=>'K', 'AAG'=>'K',
	'GAT'=>'D', 'GAC'=>'D', 'GAA'=>'E', 'GAG'=>'E',
	'TGT'=>'C', 'TGC'=>'C', 'TGA'=>'*', 'TGG'=>'W',
	'CGT'=>'R', 'CGC'=>'R', 'CGA'=>'R', 'CGG'=>'R',
	'AGT'=>'S', 'AGC'=>'S', 'AGA'=>'R', 'AGG'=>'R',
	'GGT'=>'G', 'GGC'=>'G', 'GGA'=>'G', 'GGG'=>'G',
);

our %seqs;
our (%seq_cache, %cdb_handle_cache);
our (@seq_cache_queue); # need a linkedhashmap type of thing
our ($max_seq_cache_items, $max_seq_cache_size) = (1024, 200000);

our ($cur_cdb_file, $cdb_pid, $cur_getseq_file, $seq_cache_size, $totsize_warned, $devel_size_loaded, $string_approx_loaded, $quiet, $leftover_gaussian);
local (*CDB_RH, *CDB_WH, *GETSEQ_RH);

END {
	close GETSEQ_RH if defined fileno GETSEQ_RH;
#	close CDB_WH if defined fileno CDB_WH;
#	close CDB_RH if defined fileno CDB_RH;
#	waitpid($cdb_pid, 0) if defined $cdb_pid;
}

sub logmsg {my $FH = $AKUtils::LOG || *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

sub commonSetup($) {
	# Accept a settings hash with options for setup
	# Return a settings hash with tmpdir, etc

	# Load common modules, include current path in INC and PATH, set up tmpdir, export useful symbols, parse options, etc.
}

sub dna2aa($) {
	my ($dna_seq) = @_;
	my $aminoacid_seq;
	for (my $i=0; $i<length($dna_seq); $i += 3) {
		my $codon = substr($dna_seq, $i, 3);
		$aminoacid_seq .= $codon_table{$codon};
		$aminoacid_seq .= "X" unless $codon_table{$codon};
	}
	return $aminoacid_seq;
}

# Loads data from file using Data::Dumper
sub loadData($) {
	my ($file) = @_;
	open(IN, '<', $file) or die("Unable to open file $file for reading: $!");
	my $ref;
	{	no strict 'vars';
		local $/;
		my $d = <IN>;
		$ref = eval $d;
		use strict 'vars';
	}
	close IN;
	return $ref;
}

# Saves data to a file using Data::Dumper
sub saveData($$) {
	my ($ref, $file) = @_;
	require Data::Dumper;
}

# Assumes the keys are unique. TODO: keyless load
# Expects header line, keys on column 0 unless otherwise specified (by column index or name)
# TODO: sanity checks
sub loadCSV($;$) {
	my ($file, $key) = @_;
	my %table; my $key_col;
	open(FH, '<', $file) or die("Could not open file $file for reading: ".$!);
	my $header = <FH>; chomp $header;
	my @fieldnames = split(/,/, $header);
	my %f; for (@fieldnames) { die("Empty or duplicate header field") if $f{$_} or $_ eq ''; $f{$_} = 1; }
	$key = 0 if not defined $key;
	for (0..$#fieldnames) { $key_col = $_ if $key eq $fieldnames[$_]; }
	$key_col = $key if $key =~ /^\d+$/;
	while (<FH>) {
		chomp;
		my @fields = split /,/;
		for (0..$#fieldnames) {
			$table{$fields[$key_col]}->{$fieldnames[$_]} = $fields[$_];
		}
	}
	close FH;
	return(\%table, \@fieldnames);
}

# TODO: optimize this into one pass
# TODO: distinguish between unbiased vs. biased estimator
sub stdev {
	my ($list) = @_;
	return undef if @$list < 1;
	my $mean = sum(@$list) / @$list;
#	return sqrt(sum(map { ($_-$mean)**2 } @$list) / (@$list-1));
	return sqrt(sum(map { ($_-$mean)**2 } @$list) / @$list);
}

# can also use PDL
# Box-Muller transform (http://www.taygeta.com/random/gaussian.html)
sub gaussian($$) {
	my ($mean, $variance) = @_;
	if (defined $leftover_gaussian) {
		my $y = $leftover_gaussian;
		undef $leftover_gaussian;
		return $y;
	}
	my ($x1, $x2, $w, $y1, $y2);
	$variance = 1 unless defined $variance;

	do {
		$x1 = 2.0 * rand() - 1.0;
		$x2 = 2.0 * rand() - 1.0;
		$w = $x1 * $x1 + $x2 * $x2;
	} while ($w >= 1.0);

	$w = sqrt((-2.0 * log($w)) / $w);
	$y1 = ($x1 * $w) * $variance + $mean;
	$y2 = ($x2 * $w) * $variance + $mean;
	$leftover_gaussian = $y2;
	return $y1;
}

sub argmax {
	my ($hash) = @_;
	my ($best, $key, $val);

	my $bestval = -1e999;
	while (($key,$val) = each(%$hash)) {
		if ($val > $bestval) {
			$bestval = $val;
			$best = $key;
		}
	}
	return $best;
}

# http://www.perlmonks.org/?node_id=371228
# Usage:
#	my $iter = combinations( 3 => ['a' .. 'f'] );
#	while ( my @c = $iter->() ) { print "@c\n"; }
sub combinations($$) {
	my ($num, $arr) = @_;
	no strict 'refs';
	return sub { return }
		if $num == 0 or $num > @$arr;
	use strict 'refs';
	my @pick;

	return sub {
		return @$arr[ @pick = ( 0 .. $num - 1 ) ]
			unless @pick;

		my $i = $#pick;
		$i-- until $i < 0 or $pick[$i]++ < @$arr - $num + $i;
		return if $i < 0;

		@pick[$i .. $#pick] = $pick[$i] .. $#$arr;

		return @$arr[@pick];
	};
}

# TODO: optimize me
# Usage: @list_of_lists = permutations(@list)
sub permutations(@) {
	my (@set) = @_;
	return [@set] if @set < 2;
	my @perms;
	foreach my $i (0..$#set) {
		my @tail_perms = permutations(@set[0..$i-1, $i+1..$#set]);
		foreach my $subperm (@tail_perms) {
			push(@perms, [$set[$i], @$subperm]);
		}
	}
	return @perms;
}

# Input: { a => [b, c], d => [e, f], ... }
# Output: [ { a => b, d => e }, { a => b, d => f }, ... ]
sub permuteRanges($) {
	my ($ranges) = @_;
	my $permutations = [{}];
	foreach my $attrib (keys %$ranges) {
		my $new_permutations;
		foreach my $value (@{$$ranges{$attrib}}) {
			my $new_set = dclone($permutations); # deep copy
			foreach my $p (@$new_set) {
				$$p{$attrib} = $value;
			}
			push(@$new_permutations, @$new_set);
		}
		$permutations = $new_permutations;
	}
	return $permutations;
}

# Remove null, zero elements from array
# Input, output: array ref
sub compactArray {
	my $array = (@_ == 1 ? $_[0] : \@_);
	my @new_array;
	foreach my $element (@$array) { push(@new_array, $element) if $element; }
	return wantarray ? @new_array : \@new_array;
}

# Remove zero elements from array
# Input, output: array ref
sub compactArrayWithZeros {
	my $array = (@_ == 1 ? $_[0] : \@_);
	my @new_array;
	foreach my $element (@$array) { push(@new_array, $element) if defined $element; }
	return wantarray ? @new_array : \@new_array;
}

# TODO: make this use File::Temp
sub concatFiles($;$) {
	my ($files, $tmp_dir_or_file) = @_;
	die("No files supplied") unless @$files > 0;
	$tmp_dir_or_file = ($ENV{TMPDIR} or "/tmp") unless $tmp_dir_or_file;
	if (-d $tmp_dir_or_file) {
		my $rand_name;
		do { $rand_name = substr(rand(), 2, 4) } while -e "$tmp_dir_or_file/$$.$rand_name";
		$tmp_dir_or_file .= "/$$.$rand_name";
	}
	open(OUT, '>', $tmp_dir_or_file) or die("Could not open file $tmp_dir_or_file for writing: ".$!);
	foreach my $file (@$files) {
		open(FH, '<', $file) or die("Could not open file $file for reading: ".$!);
		print OUT while <FH>;
		close FH;
	}
	close OUT;
	return $tmp_dir_or_file;
}

sub curSub() {
	return (caller(1))[3];
}

sub filesDiffer($$) {
# TODO: use Algorithm::Diff
# @f1=<FH>; $diff = Algorithm::Diff->new(\@f1, \@f2); next   if  $diff->Same();
	my ($file1, $file2) = @_;
	open(FH1, '<', $file1) or die("Could not open file $file1 for reading: ".$!);
	open(FH2, '<', $file2) or die("Could not open file $file2 for reading: ".$!);
	while (1) {
		my $line1 = <FH1>;
		my $line2 = <FH2>;
		last unless defined $line1 and defined $line2;
		return 1 if $line1 ne $line2;
	}
	return 1 if <FH1> or <FH2>;
	return 0;
}

# Traverse array of arrays down to scalars and make it a flat array of those scalars.
# optimize me
sub flatten($) {
	my ($array) = @_;
	my @flattened;

	for (ref($array) eq 'ARRAY' ? @$array : ref($array) eq 'HASH' ? values %$array : $array) {
		push(@flattened, (ref($_) ? flatten($_) : $_));
	}
	return @flattened;
}

# If argument is an executable in the current path, returns the full path to it, otherwise returns undef
sub fullPathToExec($) {
	my ($executable) = @_;
	my $fullpath;
	for ("", split(/:/, $ENV{PATH})) {
		if (-x $_."/".$executable) { $fullpath = File::Spec->rel2abs($_."/".$executable); last; }
	}
	warn("Error finding full path to file ($executable)") unless -x $fullpath;
	return $fullpath;
}

sub safeGlob($) {
	my ($regexp) = @_;
	my $dir;
	if ($regexp =~ /\//) { ($dir, $regexp) = ($regexp =~ /\A(.*\/)([^\/]+)\Z/); }
	$regexp .= "\$" unless $regexp =~ /\$$/;
	$regexp = "^".$regexp unless $regexp =~ /^\^/;
	my (@files);
	local (*DIR);

	$dir ||= ".";
	$regexp ||= ".*";
	opendir(DIR, $dir) or return;
	@files = grep { /$regexp/ } readdir(DIR);
	closedir(DIR);
	foreach my $file (@files) {
		(undef($file), next) if $file eq "." or $file eq "..";
		$file = $dir.$file if $dir;
	}
	return @{compactArray(\@files)};
}

sub md5sum($) {
	my ($filename) = @_;
	my $md5sum;
	die("$filename: file not found") unless -f $filename;
	if (`which md5sum 2>/dev/null`) {
		$md5sum = `md5sum $filename`; $md5sum =~ /^([^\s]+)/; $md5sum = $1;
	} else {
		require Digest;
		open(FH, '<', $filename) or die("Could not open file $filename for reading: ".$!);
		my $md5obj = Digest->new('MD5');
		$md5obj->addfile(\*FH);
		$md5sum = $md5obj->hexdigest();
		close FH;
	}
	return $md5sum;
}

sub cdb_indexMfa($;$$) {
	require Digest;
	my ($filename, $seqs, $quiet) = @_;
	$seqs = \%seqs if @_ < 2;
	my ($seq_name, $start_fpos, $end_fpos, $seq_len, $last_line, $no_cached_index, $md5sum, $old_md5sum);

	print curSub().": Checksumming $filename... \n" unless $quiet;
	open(FH, '<', $filename) or die("Could not open file $filename for reading: ".$!);
	my $md5obj = Digest->new('MD5');
	$md5obj->addfile(\*FH);
	$md5sum = $md5obj->hexdigest();
	close FH;

	if (-f $filename.'.md5') {
		open(FH, '<', $filename.'.md5') or die("Could not open file $filename.md5 for reading: ".$!);
		$old_md5sum = <FH>;
		chomp $old_md5sum;
		close FH;
	}

	print curSub().": Loading stored index for $filename... \n" if $md5sum eq $old_md5sum and not $quiet;
	if ($md5sum ne $old_md5sum) {
		warn(curSub().": Warning: Index checksum mismatch in $filename. Re-indexing.\n") if defined $old_md5sum and not $quiet;
		print curSub().": Indexing $filename... \n" unless $quiet;
		my $invoke_string = "cdbfasta $filename"; $invoke_string .= ">/dev/null 2>&1" if $quiet;
		system($invoke_string); die("Error $? executing \"$invoke_string\". Stopped") if $?;
		open(FH, '>', $filename.'.md5') or die("Could not open file $filename.md5 for writing: ".$!);
		print FH $md5sum."\n";
		close FH;
	}

	open(IN, "cdbyank -l $filename.cidx|") or die("Could not run \"cdbyank -l $filename.cidx\": ".$!);
	while (<IN>) {
		/^([^\t]+)\t(\d+)$/ or die; # sequence name, length
		$$seqs{$1} = [$filename, $2];
	}
	close IN;
	return $seqs;
}

# Load MFA index for $filename, append it to %$seqs, return the md5sum for the MFA file from which the index was created
=head1
sub doIndexMfaFile($$) {
	my ($filename, $seqs) = @_;
	my ($seq_name, $start_fpos, $end_fpos, $seq_len);

	print curSub().": loading $filename\n";
	open(FH, "< ".$filename) or die(curSub().": Could not open file $filename for reading: ".$!);
	while (<FH>) {
		if (/^\>[\s]*([^\s]+)/) {
			my $new_seq_name = $1;
			if (defined $seq_name) {
				warn(curSub().": duplicate sequence name \"$seq_name\" in $filename") if defined $$seqs{$seq_name};
				die(curSub().": empty sequence \"$seq_name\" in $filename") unless $start_fpos < $end_fpos;
				$$seqs{$seq_name} = [$filename, $start_fpos, $end_fpos, $seq_len];
			}
			$start_fpos = tell FH;
			$seq_name = $new_seq_name; $seq_len = 0;
		} else {
			chomp; s/\s//g; #s/[^ATGCNatgcn]//g;
			$seq_len += length;
		}
		$end_fpos = tell FH;
	}
	warn(curSub().": duplicate sequence name \"$seq_name\" in $filename") if defined $$seqs{$seq_name};
	die(curSub().": empty sequence \"$seq_name\" in $filename") unless $start_fpos < $end_fpos;
	$$seqs{$seq_name} = [$filename, $start_fpos, $end_fpos, $seq_len];
	seek FH, 0, 0; # back to start
	my $md5obj = Digest->new("MD5");
	$md5obj->addfile(\*FH);
	my $md5sum = $md5obj->hexdigest();
	close FH;

#	print curSub().": done indexing $filename\n";
	if (open(FH, '>', $filename.'.faindex')) {
		foreach my $seq_name (keys %$seqs) {
			print FH "$seq_name\t${$$seqs{$seq_name}}[1]\t${$$seqs{$seq_name}}[2]\t${$$seqs{$seq_name}}[3]\n";
		}
		print FH $md5sum."\n";
		close FH;
		print curSub().": saved mfa index in $filename.faindex\n";
	} else {
		print STDERR curSub().": failed to save mfa index in $filename.faindex\n";
	}
}
=cut


# TODO: the only thing to optimize is index creation, saving and loading, pack it with pack() and gzip
# TODO: if write to cur dir fails, write to global cache dir
# TODO: corrupt or wrong faindex will pollute %seqs before being discarded

# Input: name of file with multi-fasta sequences, ref to hash to append to
# Output: hash(key->sequence name, value->[filename, seek position of first char in sequence, seek pos of last char, real length in bp])
sub indexMfa($;$$) {
	my ($filename, $seqs, $quiet) = @_;
	$seqs = \%seqs if @_ < 2;

	print curSub().": Checksumming $filename... \n" unless $quiet;
	my $md5sum = md5sum($filename);

	my $saved_md5sum;
	if (-f $filename.".faindex") {
		open(FH, '<', $filename.".faindex") or warn(curSub().": Could not open file $filename.faindex for reading: ".$!);
		while (<FH>) {
			/^([^\t]+)\t([\d]+)\t([\d]+)\t([\d]+)$/o or (/^([\w]+)$/, $saved_md5sum = $1, last) or return;
			$$seqs{$1} = [$filename, $2, $3, $4];
		}
		chomp $saved_md5sum;
	}

	if (-f $filename.".faindex" and $md5sum eq $saved_md5sum) {
		print curSub().": Loaded stored index for $filename\n" unless $quiet;
	} else {
		print curSub().": Indexing $filename... \n" unless $quiet;

		my ($seq_name, $start_fpos, $end_fpos, $seq_len);
		open(FH, "< ".$filename) or die("Could not open file $filename for reading: ".$!);
		while (<FH>) {
			if (/^\>[\s]*([^\s]+)/) {
				my $new_seq_name = $1;
				if (defined $seq_name) {
					warn(curSub().": duplicate sequence name \"$seq_name\" in $filename") if defined $$seqs{$seq_name};
					die("Empty sequence \"$seq_name\" in $filename") unless $start_fpos < $end_fpos;
					$$seqs{$seq_name} = [$filename, $start_fpos, $end_fpos, $seq_len];
				}
				$start_fpos = tell FH;
				$seq_name = $new_seq_name; $seq_len = 0;
			} else {
				chomp; s/\s//g; #s/[^ATGCNatgcn]//g;
				$seq_len += length;
			}
			$end_fpos = tell FH;
		}
		warn(curSub().": duplicate sequence name \"$seq_name\" in $filename") if defined $$seqs{$seq_name};
		die("Empty sequence \"$seq_name\" in $filename") unless $start_fpos < $end_fpos;
		$$seqs{$seq_name} = [$filename, $start_fpos, $end_fpos, $seq_len];
		seek FH, 0, 0; # back to start
		my $md5obj = Digest->new("MD5");
		$md5obj->addfile(\*FH);
		my $md5sum = $md5obj->hexdigest();
		close FH;

		if (open(FH, '>', $filename.'.faindex')) {
			foreach my $seq_name (keys %$seqs) {
				print FH "$seq_name\t${$$seqs{$seq_name}}[1]\t${$$seqs{$seq_name}}[2]\t${$$seqs{$seq_name}}[3]\n";
			}
			print FH $md5sum."\n";
			close FH;
			print curSub().": saved mfa index in $filename.faindex\n";
		} else {
			warn curSub().": failed to save mfa index in $filename.faindex\n";
		}
	}
	return $seqs;
}

sub getSeq($;$) {
	my ($seq_name, $seqs) = @_;
	die("Error: sequence name not supplied") unless $seq_name;
	$seqs ||= \%seqs;
#print "getseq: retrieving $seq_name from ".join(":",keys %$seqs)."\n";

	my $retr_seq;
	my ($filename, $start_fpos, $end_fpos) = @{$$seqs{$seq_name}};
	die("Error: cannot find sequence seek info. Stopped") unless $filename and defined $start_fpos and $end_fpos;

#	my $csize; $csize = sum(map { $csize += length $$_ } values %seq_cache);

	my $cache_line = $seq_cache{$filename.':'.$seq_name};
	return $cache_line if defined $cache_line;

# EVICTION
# estimate new sequence to be ($end_fpos - $start_fpos) bytes
# totalSize is expensive so estimate
	$seq_cache_size += ($end_fpos - $start_fpos);

	my $shrink_seq_cache = 0;
	my $cur_seq_not_added = 1;
	my $must_free_mem_for_seq = max($end_fpos - $start_fpos - getFreeMem(), 0);
# must execute once to add tag initially!

#while ($cur_seq_not_added or $seq_cache_size > max($max_seq_cache_size,
# getFreeMem())) { - WRONG - must be "size of seq to add exceeds free mem"
# evict items from queue until it's small enough...
	my $line_to_evict = lruQueue(\@seq_cache_queue, $max_seq_cache_items - $shrink_seq_cache, $filename.':'.$seq_name);
	if ($line_to_evict) {
# TODO: check that the length below equals $end_fpos - $start_fpos
		$seq_cache_size -= length(${$seq_cache{$line_to_evict}});
		undef ${$seq_cache{$line_to_evict}};
		delete $seq_cache{$line_to_evict};
	}
	die("Internal error") unless $cur_seq_not_added or $line_to_evict;
	$cur_seq_not_added = 0;
	$shrink_seq_cache++;
#}
# END EVICTION


#print "getseq: not in cache\n";
	# NB: this must be line buffered, but Perl does that by default so we're fine
	if ($cur_getseq_file ne $filename) {
#print "getseq: reopening FH for $filename\n";
		close GETSEQ_RH if defined fileno GETSEQ_RH;
		open(GETSEQ_RH, '<', $filename) or die("Could not open file $filename for reading: ".$!);
		$cur_getseq_file = $filename;
	}
#print "getseq: $seq_name in $filename rlen ".($end_fpos - $start_fpos)."\n";

	seek(GETSEQ_RH, $start_fpos, 0);
	read(GETSEQ_RH, $retr_seq, $end_fpos - $start_fpos);
	$retr_seq =~ s/[\s\n]//g;

	$seq_cache{"$filename:$seq_name"} = \$retr_seq;
#print "getseq: done retrieving $seq_name\n";
	return \$retr_seq;
}

# TODO: make this cache aware
sub getSubSeq($$$;$) {
	my ($seq_name, $start, $length, $seqs) = @_;
	my $seq = getSeq($seq_name, $seqs);
	return substr($$seq, $start, $length);
}

sub cdb_printSeqsToFile($$;$) {
	my ($seqs_to_print, $to_file, $seqs) = @_;
	die("Error: sequence name not supplied") unless @$seqs_to_print > 0;
	$seqs ||= \%seqs;
	my %seqs_by_filename;

	foreach my $seq_name (@$seqs_to_print) {
		my ($filename, $seq_len) = @{$$seqs{$seq_name}};
		die("Error: cannot find sequence seek info. Stopped") unless defined $filename and defined $seq_name;
		push(@{$seqs_by_filename{$filename}}, $seq_name);
	}

	foreach my $filename (keys %seqs_by_filename) {
		# TODO: make processive together with getseq
		open(OUT, "|cdbyank $filename.cidx >> $to_file");
		print OUT join("\n", @{$seqs_by_filename{$filename}});
		close OUT;
	}
#	system("cdbyank $filename.cidx -a '$seq_name' >> $to_file");
}

# Output: reference to a string
# TODO: implement cdb handle cache (is it necessary though, with seq cache?)
# TODO: this is insanely slow
sub cdb_getSeq($;$) {
	my ($seq_name, $seqs) = @_;
	die("Error: sequence name not supplied") unless $seq_name;
	$seqs ||= \%seqs;

#print "getseq: retrieving $seq_name\n";
	my $retr_seq;
	my ($filename, $seq_len) = @{$$seqs{$seq_name}};
	die("Error: cannot find sequence seek info. Stopped") unless defined $filename and defined $seq_name;

	my $cache_line = $seq_cache{$filename.':'.$seq_name};
	my $line_to_evict = lruQueue(\@seq_cache_queue, $max_seq_cache_items, $filename.':'.$seq_name);
	die("Internal error") if $cache_line and $line_to_evict; # if it was in the cache, there should be nothing to evict
	delete $seq_cache{$line_to_evict} if $line_to_evict;
	return $cache_line if defined $cache_line;
#print "getseq: not in cache\n";
	#open(IN, "cdbyank $filename.cidx -a '$seq_name'|"); while(<IN>) { next if /^\>/; chomp; $retr_seq .= $_; } close IN;
	if ($cur_cdb_file ne $filename) {
#print "getseq: setting up handle\n";
		close CDB_WH if defined fileno CDB_WH;
		close CDB_RH if defined fileno CDB_RH;
		waitpid($cdb_pid, 0) if defined $cdb_pid;
		$cdb_pid = open2(\*CDB_RH, \*CDB_WH, "cdbyank $filename.cidx"); # or die...
		$cur_cdb_file = $filename;
	}
#print "getseq: requesting $seq_name\n";
	# NB: this must be line buffered, but Perl does that by default so we're fine
	print CDB_WH $seq_name."\n";

	while (<CDB_RH>) {
		next if /^\>/; # TODO: sanity check if debug: see if seq_name is present in fasta header
		last if /^\n/;
		chomp;
		$retr_seq .= $_; # TODO: more sanity checks, avoid deadlock in case of error
	}

	$seq_cache{"$filename:$seq_name"} = \$retr_seq;
	return \$retr_seq;
}

# TODO: this needs more work
sub isFasta($) {
	my ($file) = @_;
	my $char;
	open(FH, "< $file") or die("Could not open file $file for reading: ".$!);
	read(FH, $char, 1); close FH;
	if ($char eq ">") { return 1; } else { return 0; }
}

# Input: queue, max items in it, label of item to add to queue
# Output: label of item to be deleted, if any
# TODO: variable strategies, dynamic queue #items, optimize with a DLL
sub lruQueue($$$) {
	my ($queue, $max_queue_items, $item_to_add, $strategy) = @_;
	die("Bad arguments") unless $queue and $max_queue_items and defined $item_to_add;
	my $item_exists_at;
	foreach my $i (0..$#$queue) {
		if ($$queue[$i] eq $item_to_add) { $item_exists_at = $i; last; }
	}
	if (defined $item_exists_at) {
		splice(@$queue, $item_exists_at, 1);
		unshift(@$queue, $item_to_add);
#my @new_queue = ($item_to_add);push(@new_queue, @$queue[0..$item_exists_at-1]) if $item_exists_at > 0;push(@new_queue, @$queue[$item_exists_at+1..$#$queue]) if $item_exists_at < $#$queue;@$queue = @new_queue;
	} else {
		unshift(@$queue, $item_to_add);
	}

	if (@$queue > $max_queue_items) {
		die("Internal error") if @$queue > $max_queue_items+1;
		return pop(@$queue);
	} else {
		return undef;
	}
}

=head1
   1. matches - Number of bases that match that aren't repeats
   2. misMatches - Number of bases that don't match
   3. repMatches - Number of bases that match but are part of repeats
   4. nCount - Number of 'N' bases
   5. qNumInsert - Number of inserts in query
   6. qBaseInsert - Number of bases inserted in query
   7. tNumInsert - Number of inserts in target
   8. tBaseInsert - Number of bases inserted in target
   9. strand - '+' or '-' for query strand. For translated alignments, second '+'or '-' is for genomic strand
  10. qName - Query sequence name
  11. qSize - Query sequence size
  12. qStart - Alignment start position in query
  13. qEnd - Alignment end position in query
  14. tName - Target sequence name
  15. tSize - Target sequence size
  16. tStart - Alignment start position in target
  17. tEnd - Alignment end position in target
  18. blockCount - Number of blocks in the alignment (a block contains no gaps)
  19. blockSizes - Comma-separated list of sizes of each block
  20. qStarts - Comma-separated list of starting positions of each block in query
  21. tStarts - Comma-separated list of starting positions of each block in target
=cut
# $clusters = {key: db_seq, value: [[align1, ... alignN], [align1, ... alignN], ...]
sub readPSL($;$) {
	my ($filename, $clusters) = @_;
	$clusters ||= {};
	open(FH, '<', $filename) or die("Could not open file $filename for reading: ".$!);
	$_ = <FH>;
	warn(curSub().": File $filename doesn't look like a PSL file") unless /^psLayout version 3/;
	while (<FH>) {
		/^(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\+|-)\t([^\t]+)\t(\d+)\t(\d+)\t(\d+)\t([^\t]+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t([\d\,]+)\t([\d\,]+)\t([\d\,]+)$/
			or next;
		my ($n_match, $n_mism, $repm, $n_count, $q_nins, $q_bpins, $t_nins, $t_bpins, $strand, $q_name, $q_size, $q_start, $q_end, $t_name, $t_size, $t_start, $t_end, $block_cnt, $block_sizes, $q_starts, $t_starts)
			= ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21);

		push(@{$$clusters{$t_name}}, []);
		my $cur_cluster = ${$$clusters{$t_name}}[$#{$$clusters{$t_name}}];

		my ($i1, $i2, $i3);
		for (1..$block_cnt) {
			my ($next_i1, $next_i2, $next_i3) = (index($block_sizes, ',', $i1), index($q_starts, ',', $i2), index($t_starts, ',', $i3));
			my ($cur_blocksize, $cur_qstart, $cur_tstart)
				= (substr($block_sizes, $i1, $next_i1-$i1), substr($q_starts, $i2, $next_i2-$i2), substr($t_starts, $i3, $next_i3-$i3));
			($i1, $i2, $i3) = ($next_i1+1, $next_i2+1, $next_i3+1);

push(@$cur_cluster, [$q_name, $q_start, $q_end]);
		}
	}
	close FH;
	return $clusters;
}

# Usage: @indexes = allIndexes(string, substring)
# or @fuzzy_indexes = allIndexes(string, substring, #insertions, #deletions, #substitutions)
# WARNING: string::approx can produce wrong results in many situations, especially for short strings
# TODO: case insensitive, warn on min length/complexity w/amatch
# TODO: matches start at only certain position ranges or modulos
sub allIndexes($$;$$$) {
	my ($string, $substring, $in, $del, $sub) = @_;
	my @matches;
	if ($in or $del or $sub) {
		unless ($string_approx_loaded) {
			foreach my $dir (@INC) {
				if (-f "$dir/String/Approx.pm") {
					require String::Approx;
					$string_approx_loaded = 1;
				}
			}
			die("Unable to load String::Approx") unless $string_approx_loaded;
		}
		# string::approx is weird and can get stuck, limit the number of matches
		my ($nummatches, $maxmatches) = (0, 1000);
		my $pm = -length($substring);
		while (($pm = String::Approx::aindex($substring, ['I'.$in.'D'.$del.'S'.$sub.' initial_position='.($pm+length($substring))], $string)) != -1) {
			$nummatches++;
			(warn("Warning: Possible runaway match engine, indexes truncated to $maxmatches"), last) if $nummatches > $maxmatches;
			push(@matches, $pm);
		}
	} else {
		my $pm = -1;
		push(@matches, $pm) while ($pm = index($string, $substring, $pm+1)) != -1;
	}
	return wantarray ? @matches : \@matches;
}

# returns total size of structure and all its refs, in bytes, or undef if Devel::Size not installed
sub totalSize($) {
	my ($ptr) = @_;
	return undef unless $ptr;
	return Devel::Size::total_size($ptr) if $devel_size_loaded;
	foreach my $dir (@INC) {
		if (-f "$dir/Devel/Size.pm") {
			require Devel::Size;
			$devel_size_loaded = 1;
			return Devel::Size::total_size($ptr);
		}
	}
	warn(curSub().": Warning: Devel::Size not found. Unable to determine memory consumption\n") unless $totsize_warned or $quiet;
	$totsize_warned = 1;
	return undef;
}

# on Linux, returns estimated amount of free memory (MemFree + Cached), in bytes, or undef if unknown
sub getFreeMem() {
	my $free_mem;
	open(FH, '<', '/proc/meminfo') or return undef;
	while (<FH>) {
		if (/^MemFree.+?(\d+)\s+kB/) {
			$free_mem = $1;
		} elsif (/^Cached.+?(\d+)\s+kB/) {
			$free_mem += $1; last;
		}
	}
	close FH;
	return ($free_mem * 1024 or undef);
}

sub nanmean {
	my $array = ((@_ == 1 and ref($_[0]) eq 'ARRAY') ? $_[0] : \@_);
	my $size = (defined $$array[0] ? 1 : 0);
	my $sum = reduce { $size++ if defined $b; $a + $b } @$array;
	return undef if $size == 0;
	return $sum / $size;
}

sub nanmin {
	my $array = ((@_ == 1 and ref($_[0]) eq 'ARRAY') ? $_[0] : \@_);
	return undef unless $array;
	shift @$array while @$array and not defined $$array[0];
	my $min = reduce { defined $b ? $a < $b ? $a : $b : $a } @$array;
	return $min;
}

sub nanmax {
	my $array = ((@_ == 1 and ref($_[0]) eq 'ARRAY') ? $_[0] : \@_);
	return undef unless $array;
	shift @$array while @$array and not defined $$array[0];
	my $max = reduce { defined $b ? $a > $b ? $a : $b : $a } @$array;
	return $max;
}

# Compute the number of items/pairs at given traverse depth of a tree of arrays or hashes, e.g. numEltsAtDepth([[1,2],[3,4]], 1) = 4
sub numEltsAtDepth($$) {
	my ($ref, $depth) = @_;
	my @subrefs = ($ref);
	my $total;
	@subrefs = map { ref($_) eq 'ARRAY' ? @$_ : ref($_) eq 'HASH' ? values %$_ : die "Scan depth exceeds structure depth" } @subrefs while $depth--;
	map { $total += (ref($_) eq 'ARRAY' ? @$_ : ref($_) eq 'HASH' ? values %$_ : die "Scan depth exceeds structure depth") } @subrefs;
	return $total;
}

# TODO: error checking
# Usage: $seqs = readMfa($seqfile); foreach $seqname (keys %$seqs) { $seq = $$seqs{$seqname}; }
# Usage: (values %{readMfa($seqfile)})[0]
sub readMfa($;$) {
	my ($mfa_file, $settings) = @_;
	open(FH, '<', $mfa_file) or die("Could not open file $mfa_file for reading: ".$!);
	my %seqs;
	my ($cur_seq_hdr, $seq);
	while (<FH>) {
		if (/^\>\s*(.+)/) {
			$seqs{$cur_seq_hdr} = $seq if $cur_seq_hdr;
			undef $seq;
			$cur_seq_hdr = ($$settings{first_word_only} ? (split /\s/, $1)[0] : $1);
		} else {
			chomp; s/\s//g; $seq .= uc $_;
		}
	}
	close FH;
	$seqs{$cur_seq_hdr} = $seq if $cur_seq_hdr;
	die("Error: No sequences found in $mfa_file") unless %seqs;
	return \%seqs;
}

sub readGff($) {
	my ($gff_file) = @_;
	my @features;
	open(GFF, '<', $gff_file) or die "Error opening file $gff_file for reading: $!";
	while (<GFF>) {
		next if /^\#/;
		chomp;
		my @l = split /\t/;
		my %feature;
		foreach my $field (qw(seqid source type start end score strand phase attributes)) {
			my $value = shift @l;
			$feature{$field} = $value if $value ne '.';
		}
		# NB: Though the standard says attributes are in the form "attr1=value1;attr2=value2",
		# many files contain "attr1 value1;attr2 value2"
		foreach my $attrib (split /;/, $feature{attributes}) {
			$attrib =~ /^(.+)=(.+)$/ or next;
			foreach my $tag (qw(ID Name Alias Parent Target Gap Derives_from Note Dbxref Ontology_term)) {
				# NB: these can be separated by commas
				$feature{$1} = $2 if $1 eq $tag;
			}
		}
		push(@features, \%feature);
	}
	return \@features;
}

# TODO: this is a low-level parser, break the ugly out of readPSL and make that a high-level parser
sub readPSL2($) {
	my ($psl_file) = @_;
	my %hits;
	open(PSL, '<', $psl_file) or die "Error opening file $psl_file for reading: $!";
	<PSL> for 1..5; # skip header
	while (<PSL>) {
		chomp;
		my @l = split /\t/;
		my %hit;
		foreach my $field (qw(matches misMatches repMatches nCount qNumInsert qBaseInsert tNumInsert tBaseInsert strand qName qSize qStart qEnd tName tSize tStart tEnd blockCount blockSizes qStarts tStarts)) {
			my $value = shift @l;
			$hit{$field} = $value;
		}
		push(@{$hits{$hit{tName}}}, \%hit);
	}
	return \%hits;
}

# TODO: GC% dependence
# Sequence: Single sequence to mutate
# Indels: Target # indels to introduce
# Substs (TODO)
# Rate: Indel rate per Kchar (either rate or indels must be defined)
# Subregions: A list of subregions in the form [{start=>, end=>}, {start=>, end=>}, ...] to insert errors into
# TODO: record coordinates must be adjusted when producing indels upstream of them
sub mutateSequence($$$;$$$) {
	my ($seq, $indels, $substs, $rate, $subregions, $whiteout) = @_;
	my @record;
	die("Cannot accept both rate and fixed target") if $rate and ($indels or $substs);
	die("invalid rate") if defined $rate and ($rate >= 1 or $rate <= 0);
	die("invalid target") if defined $indels and (int($indels) != $indels or $indels < 0);
	die("invalid target") if defined $substs and (int($substs) != $substs or $substs < 0);
	die("invalid subregions") if defined $subregions and ref($subregions) ne 'ARRAY';
	die("invalid whiteout") if defined $whiteout and ($whiteout < 1 or $whiteout > 1000);

	$rate = 1/10000 unless $rate or $indels or $substs;

	if ($rate and not ($indels or $substs)) {
		$indels = int(length($seq) * $rate);
		die("Internal error computing mutations from rate") if $indels > 10000;
		return ($seq, \@record) if $indels < 1;
	}

	$subregions ||= [[1, length($seq)-1]];

	MUT: while (1) {
		my $pos = int(rand(length($seq)));
		CDS: foreach my $region (@$subregions) {
			my ($start, $end) = ($$region{lo}, $$region{hi});
			die("Error reading subregions") unless $start and $end;
#			die("Subregion coordinate $end exceeds sequence length ".(length($seq)+) if $end > length($seq); # this check needs adjustment
			if ($pos > $start+$whiteout and $pos < $end-$whiteout) {
				#print "position $pos admitted in REGION $start $end\n";
				if (rand() < 0.5) { # insert
					substr($seq, $pos, 1) .= substr("ATGC", int(rand(length("ATGC"))), 1);
					push(@record, [$pos, 'I']);
				} else { # delete
					substr($seq, $pos, 1) = "";
					push(@record, [$pos, 'D']);
				}
				$indels--; last MUT if $indels == 0;
				last CDS;
			}
		}
	}
	return ($seq, \@record);
}

sub alnum {
	my ($i);
	my ($len1, $len2) = (length($a), length($b));
	for ($i = 0; ($i < $len1) && ($i < $len2); ++$i) {
		my $c1 = substr($a, $i, 1);
		my $c2 = substr($b, $i, 1);
		($c1 =~ /^\d/o) || ($c2 =~ /^\d/o) || ($c1 ne $c2) and last;
	}
	my $a_r = ($i < $len1) ? substr($a, $i) : "";
	my $b_r = ($i < $len2) ? substr($b, $i) : "";
	my ($a_n, $a_s) = ($a_r =~ /^(\d+)(.*)$/o);
	my ($b_n, $b_s) = ($b_r =~ /^(\d+)(.*)$/o);
	return (defined($a_n) && defined($b_n)) ?
	(($a_n <=> $b_n) || ($a_s cmp $b_s)) : ($a cmp $b);
}

sub intersectLength($$$$) {
	my ($start1, $end1, $start2, $end2) = @_;
	my $ilen;
	die("Internal error") unless $start1 and $end1 and $start2 and $end2;
	die("Internal error") unless $start1 <= $end1 and $start2 <= $end2;
#	($start1, $end1) = ($end1, $start1) if $start1 > $end1;
#	($start2, $end2) = ($end2, $start2) if $start2 > $end2;
	if ($start1 < $start2) {
		if ($end1 < $end2) {
			$ilen = $end1 - $start2 + 1;
		} else {
			$ilen = $end2 - $start2 + 1;
		}
	} else {
		if ($end2 < $end1) {
			$ilen = $end2 - $start1 + 1;
		} else {
			$ilen = $end1 - $start1 + 1;
		}
	}
	$ilen = 0 if $ilen < 0;
	return $ilen;
}

# Returns tRNAscanSE predictions as {seqname => [{lo =>, hi =>, strand =>, trna_type =>, anticodon =>}, ...], ...}
# Settings: optional: tse_exec: executable name for tRNAscanSE, tempdir: path to temporary directory
sub gettRNAscanSEPredictions($;$) {
	my ($input_file_full, $settings) = @_;
	die("Internal error: no input supplied") unless $input_file_full;
	$input_file_full = File::Spec->rel2abs($input_file_full);
	my ($input_file, $input_dir) = fileparse($input_file_full);

	logmsg "Running tRNAscanSE on $input_file_full...";
	$$settings{tse_exec} = AKUtils::fullPathToExec($$settings{tse_exec} or 'tRNAscan-SE');

	my $tse_opts = "-q";
	$tse_opts .= "" if $$settings{trnascan_modeltype} eq 'eukaryotic';
	$tse_opts .= " -G" if $$settings{trnascan_modeltype} eq 'general';
	$tse_opts .= " -O" if $$settings{trnascan_modeltype} eq 'organellar';
	$tse_opts .= " -B" if $$settings{trnascan_modeltype} eq 'bacteria';
	$tse_opts .= " -A" if $$settings{trnascan_modeltype} eq 'archaea';

	$$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/".(split("::",(caller(0))[3]))[1].".$$.XXXXX", CLEANUP => !($$settings{keep}));

	system("$$settings{tse_exec} $tse_opts -o '$$settings{tempdir}/trnascan-se.out' '$input_file_full'");
	die("Error running $$settings{tse_exec}: $!") if $?;

	return loadtRNAscanSEPredictions("$$settings{tempdir}/trnascan-se.out");
}

sub loadtRNAscanSEPredictions($) {
	my ($predfile) = @_;
	open(PRED, '<', $predfile) or die("Could not open file $predfile for reading: ".$!);
	my $header .= <PRED> for 1..3;

	my %predictions;
	while (<PRED>) {
		chomp;
		my ($seqname, $index, $begin, $end, $trna_type, $anticodon, $intron_begin, $intron_end, $cove_score) = split /\t/;
		$begin = int($begin); $end = int($end);
		my $strand = ($begin < $end ? '+' : '-');
		push(@{$predictions{$seqname}}, {seqname => $seqname, lo => min($begin, $end), hi => max($begin, $end), strand => $strand,
			trna_type => $trna_type, anticodon => $anticodon,
			intron_begin => $intron_begin, intron_end => $intron_end, cove_score => $cove_score});
	}
	close PRED;
	return \%predictions;
}

# TODO: this might still be wrong - compare with bitmap-based
# TODO: check partial (end-of-sequence) orfs
# TODO: add options to discard ORFs with excessive Ns, require start (Met), accept partial (end-of-seq) orfs
sub findOrfs($$) {
	my ($seqs, $min_len) = @_;
	my ($start_codons, $stop_codons) = ("ATG|GTG", "TAA|TAG|TGA"); # FIXME: use startCodons()
	my ($rs_start_codons, $rs_stop_codons) = (scalar(reverse($start_codons)), scalar(reverse($stop_codons)));
	$rs_start_codons =~ tr/ATGC/TACG/; $rs_stop_codons =~ tr/ATGC/TACG/;
	my %orfs;

	while (my ($seqname, $seq) = each(%$seqs)) {
		warn("Sequence $seqname contains spaces or newlines") if $seq =~ /(\n|\s)/;
		warn("Sequence $seqname contains characters other than ATGCN") unless $seq =~ /^[ATGCN]+$/;

		foreach my $i (0, 1, 2) { # frame
			pos($seq) = $i;
			while ($seq =~ /\G((?:.{3})+?(?:$stop_codons))/go) { # forward strand
				my $orf = $1;
	#			if ($orf =~ /($start_codons)((?:.{3})+)$/o and length($2)+3 >= $min_len) { # contains at least one in-frame start
				if (length($orf)+3 >= $min_len) {
					$orfs{$seqname}->{$i}->{pos($seq)} = pos($seq)-length($orf)+1;
				}
			}
		}

		my $rc_seq = reverse($seq); $rc_seq =~ tr/ATGC/TACG/;
		foreach my $i (0, 1, 2) { # frame
			pos($rc_seq) = ((length($rc_seq) % 3) - $i) % 3;
			while ($rc_seq =~ /\G((?:.{3})+?(?:$stop_codons))/go) {
				my $orf = $1;
	#			if ($orf =~ /($start_codons)((?:.{3})+)$/o and length($2)+3 >= $min_len) { # contains at least one in-frame start
				if (length($orf)+3 >= $min_len) {
					$orfs{$seqname}->{$i.'R'}->{length($rc_seq)-pos($rc_seq)+1} = length($rc_seq) - pos($rc_seq) + length($orf);
				}
			}
		}
	}
	return \%orfs;
}

=head1
Extract prokaryotic (uninterrupted) ORFs from a nucleotide sequence.
Output: {seq1name => {frame => {stop_coord => {start=>N, stop=>N, seq=>..., aa_seq=>...}}}}
Settings:
	min_len: minimum length of ORFs to report. Default is 90 nt.
	need_seq: output ORF nucleotide sequence in the "seq" field. Default is false.
	need_aa_seq: output ORF amino acid sequence in the "aa_seq" field (starting at the reported start coordinate). Default is false.
x	confirm_start: Require a Met codon to be in the ORF, and modify the meaning of min_len to mean
x		the distance from the most upstream Met codon to the end of the ORF.
	get_truncated_orfs: Report ORFs truncated by the start or end of the sequence. The confirm_start setting is
		ignored for these ORFs; the min_len setting is in effect unless min_partorf_len is set. The coordinate
		closest to the truncated end is reported as the first nucleotide of a full in-frame codon on that end
		(i.e. for sequence start, it's 1, 2, or 3). The fields truncated_lo or truncated_hi are set accordingly.
x	min_partorf_len: Require truncated ORFs to be at least this length. Default is equal to min_len.
	orftype: if set to start2stop, 
The output fields are:
	ustop: the first nucleotide downstream of the stop codon of the preceding in-frame ORF
	start: the most upstream nucleotide of the first Met codon in the ORF
	stop: the most downstream nucleotide of the stop codon in the ORF
	truncated_upstream: set to true if the ORF may be truncated upstream (no in-frame stop codon precedes the ORF)
	truncated_downstream: set to true if the ORF is truncated downstream (sequence terminates before the ORF's stop codon)
	seq: nucleotide sequence from start to stop
	aa_seq: nucleotide sequence from start to stop
=cut
sub findOrfs2($;$) {
	my ($seqs, $settings) = @_;
	if (ref($seqs) ne 'HASH') { $seqs = {'seq1' => $seqs}; }
	logmsg "Extracting ORFs from ".keys(%$seqs)." sequences...";
	my ($start_codons, $stop_codons) = ("ATG|GTG", "TAA|TAG|TGA"); # FIXME: use startCodons()
	my ($rs_start_codons, $rs_stop_codons) = (scalar(reverse($start_codons)), scalar(reverse($stop_codons)));
	$rs_start_codons =~ tr/ATGC/TACG/; $rs_stop_codons =~ tr/ATGC/TACG/;
	$$settings{min_len} = 90 unless defined $$settings{min_len};
	my %orfs;
	while (my ($seqname, $seq) = each(%$seqs)) {
		die("Sequence $seqname contains spaces or newlines") if $seq =~ /(\n|\s)/;
		warn("Sequence $seqname contains characters other than ATGCN") unless $seq =~ /^[ATGCN]+$/;

		foreach my $frame (0, 1, 2) {
			my $cur_orf_start = $frame;
			my $cur_pos;
			for ($cur_pos = $frame; $cur_pos < length($seq)-2; $cur_pos += 3) { # walk codon-by-codon
				my $codon = substr($seq, $cur_pos, 3);
				if ($codon =~ /^(?:$stop_codons)$/ or $cur_pos+3 >= length($seq)-2) { # codon is a stop codon or this is the end of sequence
					my $orf_len = $cur_pos - $cur_orf_start + 3;
					if ($orf_len < $$settings{min_len}) {
						$cur_orf_start = $cur_pos+3; next;
					} elsif (not($$settings{get_truncated_orfs}) and $cur_pos+3 >= length($seq)-2 and $codon !~ /^(?:$stop_codons)$/) {
						$cur_orf_start = $cur_pos+3; next;
					}
					
					my $truncated_upstream = 1 if $cur_orf_start == $frame;
					my $truncated_downstream = 1 if $cur_pos+3 >= length($seq)-2 and $codon !~ /^(?:$stop_codons)$/;
					my $orf_nt_seq = substr($seq, $cur_orf_start, $orf_len);
					
					my $start_offset = findStartInORF($orf_nt_seq, $start_codons);
					if ($start_offset < 0 or $orf_len-$start_offset < $$settings{min_len}) {
						$cur_orf_start = $cur_pos+3; next;
					}
					if ($$settings{orftype} eq "start2stop") {
						$orf_nt_seq = substr($seq, $cur_orf_start+$start_offset, $orf_len-$start_offset);
					}

					$orfs{$seqname}->{$frame}->{$cur_pos} = {
						ustop => $cur_orf_start+1,
						start => $cur_orf_start+1+$start_offset,
						stop => $cur_pos+3,
						length => $orf_len,
						strand => '+'};
					$orfs{$seqname}->{$frame}->{$cur_pos}->{truncated_upstream} = 1 if $truncated_upstream;
					$orfs{$seqname}->{$frame}->{$cur_pos}->{truncated_downstream} = 1 if $truncated_downstream;
					$orfs{$seqname}->{$frame}->{$cur_pos}->{seq} = $orf_nt_seq if $$settings{need_seq};
					$orfs{$seqname}->{$frame}->{$cur_pos}->{aa_seq} = dna2aa($orf_nt_seq) if $$settings{need_aa_seq};
					$cur_orf_start = $cur_pos+3;
				}
			}
		}
		
		foreach my $frame (0, 1, 2) { # Reverse frame
			my $cur_orf_start = length($seq) - ((length($seq)-$frame) % 3) - 1;
			my $cur_pos;
			for ($cur_pos = $cur_orf_start; $cur_pos >= 2; $cur_pos -= 3) { # walk codon-by-codon
				my $codon = substr($seq, $cur_pos-2, 3);
				if ($codon =~ /^(?:$rs_stop_codons)$/ or $cur_pos-3 < 0) { # codon is a stop codon
					my $orf_len = $cur_orf_start - $cur_pos + 3;
					if ($orf_len < $$settings{min_len}) {
						$cur_orf_start = $cur_pos-3; next;
					} elsif (not($$settings{get_truncated_orfs}) and $cur_pos-3 < 0 and $codon !~ /^(?:$rs_stop_codons)$/) {
						$cur_orf_start = $cur_pos-3; next;
					}

					my $truncated_upstream = 1 if $cur_orf_start == length($seq) - ((length($seq)-$frame) % 3) - 1;
					my $truncated_downstream = 1 if $cur_pos-3 < 0 and $codon !~ /^(?:$rs_stop_codons)$/;

					my $orf_nt_seq = substr($seq, $cur_pos-2, $orf_len);
die unless length($orf_nt_seq) % 3 == 0;
#print "$frame R:$cur_pos: $orf_nt_seq\n";
					$orf_nt_seq = reverse($orf_nt_seq); $orf_nt_seq =~ tr/ATGC/TACG/;
#print "$frame R rc:$cur_pos: $orf_nt_seq\n";

					my $start_offset = findStartInORF($orf_nt_seq, $start_codons);
					if ($start_offset < 0 or $orf_len-$start_offset < $$settings{min_len}) {
						$cur_orf_start = $cur_pos-3; next;
					}
					if ($$settings{orftype} eq "start2stop") {
						$orf_nt_seq = substr($seq, $cur_pos-2, $orf_len-$start_offset);
						$orf_nt_seq = reverse($orf_nt_seq); $orf_nt_seq =~ tr/ATGC/TACG/;
					}
#print "$frame R [start2stop]:$cur_pos: $orf_nt_seq\n";
die(length($orf_nt_seq)) unless length($orf_nt_seq) % 3 == 0;
					$orfs{$seqname}->{$frame.'R'}->{$cur_pos-2} = {
						ustop => $cur_orf_start+1,
						start => $cur_orf_start+1-$start_offset,
						stop => $cur_pos-1,
						length => $orf_len,
						strand => '-',
					};
					$orfs{$seqname}->{$frame.'R'}->{$cur_pos-2}->{truncated_upstream} = 1 if $truncated_upstream;
					$orfs{$seqname}->{$frame.'R'}->{$cur_pos-2}->{truncated_downstream} = 1 if $truncated_downstream;
					$orfs{$seqname}->{$frame.'R'}->{$cur_pos-2}->{seq} = $orf_nt_seq if $$settings{need_seq};
					$orfs{$seqname}->{$frame.'R'}->{$cur_pos-2}->{aa_seq} = dna2aa($orf_nt_seq) if $$settings{need_aa_seq};
					$cur_orf_start = $cur_pos-3;
				}
			}
		}
	}
	return \%orfs;
}

# Internal for findOrfs
sub findStartInORF($$) {
	my ($seq, $start_codons) = @_;
	for (my $cur_pos = 0; $cur_pos < length($seq)-2; $cur_pos += 3) {
		return $cur_pos if substr($seq, $cur_pos, 3) =~ /^(?:$start_codons)$/;
	}
	return -1;
}


sub getGlimmer3Predictions($;$$) {
	my ($input_file_full, $seqs, $settings) = @_;
	die("Internal error: no input supplied") unless $input_file_full;
	$input_file_full = File::Spec->rel2abs($input_file_full);
	my ($input_file, $input_dir) = fileparse($input_file_full);

	logmsg "Preparing to run glimmer3 on $input_file_full...";

	$$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/".(split("::",(caller(0))[3]))[1].".$$.XXXXX", CLEANUP => !($$settings{keep}));

	$seqs ||= readMfa($input_file_full);

	my $longorfs_infile = $input_file_full;
	if (scalar(keys %$seqs) > 1) {
		$longorfs_infile = "$$settings{tempdir}/glimmer3_longorfs_in.fasta";
		open(FH, '>', $longorfs_infile) or die "Could not open file $longorfs_infile for writing: ".$!;
		# TODO: CHECK - spacer might bias model
		my $seqs = join('TTAGTTAGTTAG', values %$seqs); # sequences separated by all-frame stop spacer
		$seqs =~ s/(.{80})/$1\n/g;
		$seqs .= "\n" unless $seqs =~ /\n$/;
		print FH ">glimmer3_longorfs_in\n$seqs";
		close FH;
	}

	my ($longorfs_opts, $glimmer_opts);
	# Sequences are assumed to be non-circular if there's more than one of them
	$longorfs_opts = "--linear" if $$settings{linear_genome} or scalar(keys %$seqs) > 1;
	$glimmer_opts = " --linear" if $$settings{linear_genome} or scalar(keys %$seqs) > 1;
	$glimmer_opts .= " --rbs_pwm $$settings{gl3_rbs_pwm_file}" if $$settings{gl3_rbs_pwm_file};
	$glimmer_opts .= " --gene_len $$settings{gl3_min_gene_len}" if $$settings{gl3_min_gene_len};
	$glimmer_opts .= " --threshold $$settings{gl3_calling_threshold}" if $$settings{gl3_calling_threshold};
	$glimmer_opts .= " --max_olap $$settings{gl3_max_overlap}" if $$settings{gl3_max_overlap};
	$glimmer_opts .= " --extend" unless $$settings{gl3_no_extended_predicts};

# FIXME: don't spam the input file location
	system("long-orfs $longorfs_opts --no_header --cutoff 1.15 '$longorfs_infile' '$input_file.longorfs' 2>'$$settings{tempdir}/glimmer3.log'");
	die("Error running long-orfs: $!") if $?;
	system("extract -t '$longorfs_infile' '$input_file.longorfs' > '$$settings{tempdir}/$input_file.train'");
	die("Error running extract: $!") if $?;
	system("build-icm -r '$$settings{tempdir}/$input_file.icm' < '$$settings{tempdir}/$input_file.train'");
	die("Error running build-icm: $!") if $?;
	logmsg "Running glimmer3 on $input_file_full...";
	system("glimmer3 $glimmer_opts '$input_file_full' '$$settings{tempdir}/$input_file.icm' '$$settings{tempdir}/$input_file.glimmer3' 2>>'$$settings{tempdir}/glimmer3.log'");
	die("Error running glimmer3: $!") if $?;

	return loadGlimmer3Predictions("$$settings{tempdir}/$input_file.glimmer3.predict");
}

# NB: On extended hits, glimmer3 can predict coordinates 1 or 2 nt beyond the edge of the sequence to preserve length %3 == 0
# This behavior is unlike gmhmmp, which predicts those as "<0" or ">[end]" from which </> are stripped
sub loadGlimmer3Predictions($) {
	my ($predfile) = @_;
	open(PRED, '<', $predfile) or die("Could not open file $predfile for reading: ".$!);
	my %predictions;
	my $cur_seq = 'unknown';
	while (<PRED>) {
		if (/^\>\s*(.+?)\s*$/) {
			$cur_seq = $1;
		} else {
			my @line = split /\s+/;
			my ($gene_id, $start, $end, $frame, $raw_score) = @line;
			my $strand = ($frame =~ /\+/ ? '+' : '-');
			push(@{$predictions{$cur_seq}}, { seqname => $cur_seq, lo => min($start, $end), hi => max($start, $end), strand => $strand,
				gl3_score => $raw_score, type => 'CDS', start => $start, stop => $end, predictor => 'Glimmer3'});
		}
	}
	close PRED;
	return \%predictions;
}

=head1 getGenemarkPredictions
Accepts filename with mfa data
Outputs GMHMMP predictions (see loadGMHMMPredictions for format)
Settings:
	gm_predictor - name of or full path to gmhmmp
	gm_trainer - name of or full path to gmsn.pl
	gm_predictor_xopts - extra options to pass to gmhmmp
	gm_trainer_xopts - extra options to pass to gmsn.pl
	gmhmm_model -  (If not specified, trains the model using GeneMarkS)
=cut
sub getGenemarkPredictions($;$$) {
	my ($input_file_full, $seqs, $settings) = @_;
	die("Internal error: no input supplied") unless $input_file_full;
	$input_file_full = File::Spec->rel2abs($input_file_full);
	my ($input_file, $input_dir) = fileparse($input_file_full);

	$$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/".(split("::",(caller(0))[3]))[1].".$$.XXXXX", CLEANUP => !($$settings{keep}));
	$$settings{gm_trainer} ||= AKUtils::fullPathToExec('gmsn.pl') or die;
	$$settings{gm_predictor} ||= AKUtils::fullPathToExec('gmhmmp') or die;

	$seqs ||= readMfa($input_file_full);

	if (defined $$settings{gmhmm_model}) {
		die("Model file $$settings{gmhmm_model} not found") unless -f $$settings{gmhmm_model};
	} else {
		logmsg "Running $$settings{gm_trainer} on $input_file_full...";
		my $gm_trainer_opts = "--combine --gm $$settings{gm_trainer_xopts}";
		my $invoke_string = "cd $$settings{tempdir}; $$settings{gm_trainer} $gm_trainer_opts '$input_file_full' >/dev/null 2>&1";
		system($invoke_string);
		die("Error running $$settings{gm_trainer}: $!") if $?;
		$$settings{gmhmm_model} = "$$settings{tempdir}/GeneMark_hmm_combined.mod";
		die("Error running $$settings{gm_trainer}: no model generated") unless -f $$settings{gmhmm_model};
	}

	logmsg "Running $$settings{gm_predictor} on $input_file_full...";

	my %predictions;
	my $lst_file = "$$settings{tempdir}/gm_out.lst";
	unlink $lst_file;
	while (my ($seqname, $seq) = each(%$seqs)) {
		my $temp_fna = "$$settings{tempdir}/temp.fna";
		open(OUT, '>', $temp_fna) or die("Unable to open file $temp_fna for writing: $!");
		print OUT ">$seqname\n$seq";
		close OUT;

		my $gm_predictor_opts = "-r -k -m $$settings{gmhmm_model} -o $temp_fna.lst $$settings{gm_predictor_xopts}";
		system("cd $$settings{tempdir}; $$settings{gm_predictor} $gm_predictor_opts $temp_fna");
		die("Error running $$settings{gm_predictor}: $!") if $?;
		
		open(IN, '<', "$temp_fna.lst") or die("Internal error");
		open(LST, '>>', $lst_file);
		print LST ">$seqname\n";
		print LST while <IN>;
		close LST;
		close IN;

		unlink $temp_fna;
		unlink "$temp_fna.lst";
	}
	return loadGMHMMPredictions($lst_file);
}

=head1
Load a gmhmmp output file. Return a hash:
TODO
=cut
sub loadGMHMMPredictions($) {
	my ($predfile) = @_;
	open(PRED, '<', $predfile) or die("Could not open file $predfile for reading: ".$!);
	my %predictions;
	my $cur_seq = 'unknown';
	while (<PRED>) {
		if (/^\>\s*(.+?)\s*$/) {
			$cur_seq = $1;
		} else {
			my @line = split /\s+/;
			shift(@line) if $line[0] eq '';
			next unless @line > 5 and $line[0] =~ /^\d+$/;
			my ($Gene, $Strand, $LeftEnd, $RightEnd, $GeneLength, $Class, $Spacer, $RBS_score) = @line;
			my $extended_lo = 1 if $LeftEnd =~ /[<>]/;
			my $extended_hi = 1 if $RightEnd =~ /[<>]/;
			$LeftEnd =~ s/[<>]//g; $RightEnd =~ s/[<>]//g;
			my ($start, $end) = ($Strand eq '+' ? ($LeftEnd, $RightEnd) : ($RightEnd, $LeftEnd));
			push(@{$predictions{$cur_seq}}, { seqname => $cur_seq, lo => $LeftEnd, hi => $RightEnd, strand => $Strand,
				pg_class => $Class, rbs_spacer => $Spacer, rbs_score => $RBS_score,
				extended_lo => $extended_lo, extended_hi => $extended_hi, 
				type => 'CDS', start => $start, stop => $end, predictor => 'gmhmmp'});
		}
	}
	close PRED;
	return \%predictions;
}

# DEPRECATED
sub getGenePredictions($;$$$$) {
	my ($input_file, $gm_trainer, $gm_predictor, $seqs, $model) = @_;
	return getGenemarkPredictions($input_file, $seqs, {gm_trainer => $gm_trainer, gm_predictor => $gm_predictor, gmhmm_model => $model});
}

sub printSeqsToFile($$) {
	my ($seqs, $outfile) = @_;
	open(OUT, '>', $outfile) or die "Unable to open file $outfile for writing: $!";
	foreach my $seqname (keys %$seqs) {
		my $seq = $$seqs{$seqname};
		$seq =~ s/(.{80})/$1\n/g;
		$seq .= "\n" unless $seq =~ /\n$/;
		print OUT ">$seqname\n$seq";
	}
	close OUT;
	return $outfile;
}

sub getNumCPUs() {
	my $num_cpus;
	open(IN, '<', '/proc/cpuinfo'); while (<IN>) { /processor\s*\:\s*\d+/ or next; $num_cpus++; } close IN;
	return $num_cpus || 1;
}

sub accessions2gb($$;$) {
	my ($accessions, $gb_filename, $settings) = @_;

	$$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/".(split("::",(caller(0))[3]))[1].".$$.XXXXX", CLEANUP => !($$settings{keep}));

	open(LST, '>', "$$settings{tempdir}/acc_lst") or die "Unable to open $$settings{tempdir}/acc_lst for writing: $!";
	for (my $i=0; $i<@$accessions; $i+=128) {
		my $qs = "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nucleotide&sendto=t&list_uids=";
		for (my $j=$i; $j<min($i+128, scalar(@$accessions)); $j++) {
			$qs .= $$accessions[$j].",";
		}
		print LST "$qs\n";
	}
	close LST;
	system("wget -i $$settings{tempdir}/acc_lst -O $gb_filename");
	die("Error running wget: $!") if $?;
	return $gb_filename;
}

=head1
Create a BLAST database. Return the database location.
Settings:
	formatdb_protein_db: set to true if protein sequence, false/undef if nucleotide.
	formatdb_in_place: create database in the given location, not in temporary space
=cut
sub formatBLASTdb($;$) {
	my ($input_fasta_file, $settings) = @_;
	my $input_file_full = File::Spec->rel2abs($input_fasta_file);
	my ($input_file, $input_dir) = fileparse($input_file_full);

	my $formatdb_opts = ($$settings{formatdb_protein_db} ? "-p T" : "-p F");
	$formatdb_opts .= " -t $$settings{formatdb_db_title}" if defined $$settings{formatdb_db_title};
	$formatdb_opts .= " -n $$settings{formatdb_db_name}" if defined $$settings{formatdb_db_name};
	my ($invoke_string, $db_loc);
	if ($$settings{formatdb_in_place}) {
		$invoke_string = "formatdb $formatdb_opts -i $input_file_full";
		$db_loc = $input_file_full;
	} else {
		$$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/".(split("::",(caller(0))[3]))[1].".$$.XXXXX", CLEANUP => !($$settings{keep}));
		symlink($input_file_full, "$$settings{tempdir}/$input_file");
		$invoke_string = "cd $$settings{tempdir}; formatdb $formatdb_opts -i $input_file";
		$db_loc = "$$settings{tempdir}/$input_file"
	}
	system($invoke_string); die("Error running formatdb: $!") if $?;
	return $db_loc;
}

# BLAST all sequences in %$seqs against a supplied database. Return results in tabulated format.
# mode = blastp|blastn|blastx|psitblastn|tblastn|tblastx (required)
# settings required: blast_db
# settings optional: blastseqs_xopts, num_cpus, tempdir
sub blastSeqs($$$) {
	my ($seqs, $mode, $settings) = @_;

	die("Internal error: invalid mode") if $mode !~ /^(blastp|blastn|blastx|psitblastn|tblastn|tblastx)$/;
	die("No BLAST database name supplied") unless $$settings{blast_db};

	$$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/".(split("::",(caller(0))[3]))[1].".$$.XXXXX", CLEANUP => !($$settings{keep}));

	$$settings{num_cpus} ||= getNumCPUs();

	my $blast_infile = "$$settings{tempdir}/blastseqs.in";
	my $blast_outfile = "$$settings{tempdir}/blastseqs.out";
	printSeqsToFile($seqs, $blast_infile);

	my $blast_qs = "blastall -p $mode -d $$settings{blast_db} -m 8 -a $$settings{num_cpus} -i $blast_infile -o $blast_outfile ";
	$blast_qs .= $$settings{blast_xopts};

	logmsg "Running $blast_qs";
	system($blast_qs);
	if ($?) {
		my $er_str = "Error running \"$blast_qs\": $!";
		$$settings{ignore_blast_errors} ? warn($er_str) : die($er_str);
	}
#	open(BLAST_OUT, "$blast_qs |") or die "Unable to run \"$blast_qs\": $!";
	return loadBLAST8($blast_outfile);
}

sub loadBLAST8($) {
	my ($blast_outfile) = @_;
	my @hits;
	open(BLAST_OUT, '<', $blast_outfile) or die "Unable to open $blast_outfile for reading: $!";
	while (<BLAST_OUT>) {
		chomp;
		my @line = split /\t/;
		my %hit;
		for (qw(name1 name2 percent_id al_len mismatch_bp gap_openings start1 end1 start2 end2 Evalue bitscore)) {
			$hit{$_} = shift @line;
		}
		push(@hits, \%hit);
	}
	close BLAST_OUT;

	return \@hits;
}

__END__

# + strand only for now
# seq coords start at 1
sub startStopPos($) {
	my ($seq) = @_;
	my ($start_codons, $stop_codons) = (['ATG','GTG'], ['TAA','TAG','TGA']);
	my (%start_codon_pos, %stop_codon_pos);

	foreach my $pos (0..length($seq)/3) {
		FRAME: foreach my $frame (0, 1, 2) {
			my $codon = substr($seq, ($pos*3)+$frame, 3);
			foreach my $startc (@$start_codons) {
				if ($codon eq $startc) {
					push(@{$start_codon_pos{$frame}}, ($pos*3)+$frame+1);
					next FRAME;
				}
			}
			foreach my $stopc (@$stop_codons) {
				if ($codon eq $stopc) {
					push(@{$stop_codon_pos{$frame}}, ($pos*3)+$frame+1);
					next FRAME;
				}
			}
		}
	}
	return (\%start_codon_pos, \%stop_codon_pos);
}
