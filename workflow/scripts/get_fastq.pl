#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my ($r1, $r2, $out);

GetOptions(
    "1=s"   => \$r1,
    "2=s"   => \$r2,
    "out=s" => \$out,
) or die "Usage: $0 -1 reads_1.fastq -2 reads_2.fastq -out output.fasta\n";

die "Error: -1, -2 and -out are required.\n"
    unless defined $r1 && defined $r2 && defined $out;

open(my $OUT, ">", $out) or die "Cannot write $out: $!\n";

foreach my $file ($r1, $r2) {
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
