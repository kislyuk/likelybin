#!/usr/bin/env perl

use strict;

#use PDL;
#use PDL::LinearAlgebra;
use Data::Dumper;
$Data::Dumper::Indent=0;

my $mfile = "mats.m";
my $ofile = "mats.octave";
my $max_ord = 5;

open(OUT, '>', $mfile) or die;
foreach my $o (1..$max_ord) {
	my $m = getMat($o);
	print OUT "A$o=[";
	foreach my $r (@$m) {
		$_ ||= 0 for @$r;
		print OUT "[".join(" ", @$r)."];\n";
	}
	print OUT "];\n";
	warn "$o\t".scalar(@{$$m[0]})."\n";
}
print OUT "N$_ = null(A$_)';\n" for 1..$max_ord;
print OUT "save -text $ofile ".join(' ', "N1".."N$max_ord")."\n";
close OUT;
system("octave $mfile >/dev/null"); die if $?;

open(IN, '<', $ofile) or die;
my $name; my %mats;
while (<IN>) {
	chomp;
	if (/name: (.+)/) {
		$name = $1; next;
	} elsif (/\#/) {
		next;
	} else {
		my @l=split /\s+/;
		shift @l;
		push(@{$mats{$name}}, \@l);
	}
}
close IN;

if (1) { # c
	my $separators = 0;
	foreach my $name (sort keys %mats) {
#		print "const double ${name}[".@{$mats{$name}}."][".@{$mats{$name}->[0]}."] = {\n";
		print "const double ${name}[] = {";
		foreach my $row (@{$mats{$name}}) {
			print "\t{" if $separators;
			my $prow = join(", ", @$row);
			$prow =~ s/(.{80}, )/$1\n\t/g;
			print "$prow, \n";
			print "}," if $separators;
		}
		print "};\n";
	}
	print "\nconst double * Nmats[] = {".join(", ", sort keys %mats)."};\n";
	print "\nconst int Ndims[".keys(%mats)."][2] = {";
	foreach my $name (sort keys %mats) {
		print "{".@{$mats{$name}}.", ".@{$mats{$name}->[0]}."}, ";
	}
	print "};\n";
} else { # perl
	print Dumper(\%mats);
}

unlink $mfile; unlink $ofile;

sub getMat($) {
	my ($order) = @_;
	my $row_len; $row_len += 4 ** $_ for 1..$order;
	my @matrix;
	foreach my $ord (1..$order) {
		my $start_index; $start_index += 4 ** $_ for 1..$ord-1; $start_index++;
		my $stop_index; $stop_index += 4 ** $_ for 1..$ord;

		# summation row
		my @sum_row; $sum_row[$row_len-1] = undef;
		$sum_row[$_-1] = 1 for $start_index..$stop_index;
		push(@matrix, \@sum_row);

		# complementarity rows
		foreach my $i ($start_index..$stop_index) {
			next unless nucmer2index(revcomp(index2nucmer($i))) > $i;
			my @comp_row; $comp_row[$row_len-1] = undef;
			$comp_row[$i-1] = 1; $comp_row[nucmer2index(revcomp(index2nucmer($i)))-1] = -1;
			push(@matrix, \@comp_row);
		}

		next if $ord == 1;
		# marginal constraint rows
		my @marginals = qw(A T G C);
		# All combinations of A, T, G, C of length $ord-1
		@marginals = map(($_.'A', $_.'T', $_.'G', $_.'C'), @marginals) for 1..$ord-2;

		foreach my $m (@marginals) {
			my @marginal_row1; $marginal_row1[$row_len-1] = undef;
			my @marginal_row2; $marginal_row2[$row_len-1] = undef;
			$marginal_row1[nucmer2index($m)-1] = 1;
			$marginal_row2[nucmer2index($m)-1] = 1;
			foreach my $n (qw(A T G C)) {
#				next if $m =~ /^$n+$/;
				$marginal_row1[nucmer2index($n.$m)-1] = -1;
				$marginal_row2[nucmer2index($m.$n)-1] = -1;
			}
			push(@matrix, (\@marginal_row1, \@marginal_row2));
		}

	}
	return \@matrix;
}

sub revcomp($) {
	my ($s) = @_;
	$s =~ tr/ATGC/TACG/;
	$s = reverse($s);
	return $s;
}

# Performs the following mapping:
# 1 => A, 2 => T, 3 => G, 4 => C,
# 5 => AA, 6 => AT, 7 => AG, 8 => AC, 9 => TA, ... 20 => CC,
# 21 => AAA, ...
sub index2nucmer($) {
	my ($i) = @_;
	die if $i < 1;
	my ($length, $nucmer);
	while (1) {
		my $upper_bd; $upper_bd += 4**$_ for 1..$length;
		$i-1 < $upper_bd ? last : $length++;
	}
	my $base_i=0; $base_i += 4 ** $_ for 1..$length-1; # indices taken up by lower orders
	$i -= $base_i + 1; # remove offset for lower orders; index at 0 instead of 1
	for (my $b=$length; $b>0; $b--) {
		my $stride = 4**($b-1);
		$nucmer .= int($i/$stride);
		$i %= $stride;
	}
	$nucmer =~ tr/0123/ATGC/;
#print "length $length\tbaseline $base_i\tnucmer $nucmer\n";
	return $nucmer;
}

sub nucmer2index($) {
	my ($nucmer) = @_;
	die if $nucmer =~ /[^ATGC]/;

	my $i;
	$i += 4**$_ for 1..length($nucmer)-1;
	foreach my $pos (0..length($nucmer)-1) {
		my $c = substr($nucmer, $pos, 1);
		$c =~ tr/ATGC/0123/;
		my $stride = 4**(length($nucmer)-$pos-1);
		$i += $c*$stride;
	}
	return $i+1;
}

__END__
# Don't use on negative numbers
#use Math::BaseCalc;
#my $calc = new Math::BaseCalc(digits => [0..3]);
sub inbase($$) {
	my ($n, $k) = @_;
	my $d = int($n / $k);
	return (!!$d && inbase($d, $k)) . $n % $k;
}
