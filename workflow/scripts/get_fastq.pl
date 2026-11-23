#!/usr/bin/env perl

# Writes the reads of one or more FASTQ files (a pair, or one file of single-end reads) as one FASTA file.

use strict;
use warnings;
use Getopt::Long;

my $out;

GetOptions(
    "out=s" => \$out,
) or die "Usage: $0 -out output.fasta reads.fastq [reads_2.fastq ...]\n";

die "Error: -out and at least one FASTQ file are required.\n"
    unless defined $out && @ARGV;

open(my $OUT, ">", $out) or die "Cannot write $out: $!\n";

foreach my $file (@ARGV) {
    open(my $IN, "<", $file) or die "Cannot open $file: $!\n";

    while (my $header = <$IN>) {
        my $seq  = <$IN>;
        my $plus = <$IN>;
        my $qual = <$IN>;

        chomp($header, $seq);

        $header =~ s/^@/>/;

        print $OUT "$header\n$seq\n";
    }

    close($IN);
}

close($OUT);
