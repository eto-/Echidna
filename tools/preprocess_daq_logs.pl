#!/usr/bin/perl

use Time::ParseDate;

my $ifile = shift;

open (IFILE, $ifile) or die "can not open $ifile file: $!\n";

while (defined (my $line = <IFILE>)) {
  chomp $line;
  if ($line =~ /Syncing/) {
    die "unknown line format for $line\n" unless ($line =~ /^(.*)\ ppc(\d+).*board\ (\d+)/);
    my ($ppc, $board) = ($2, $3);
    my $timet = parsedate ($1);
    next unless (defined $timet);
    print "$timet\t1\t$ppc\t$board\t0\n";
  } elsif ($line =~ /Too\ many/) {
    die "unknown line format for $line\n" unless ($line =~ /^(.*)\ ppc(\d+).*many\ (\d+).*board\ (\d+)/);
    my ($ppc, $board, $hits) = ($2, $4, $3);
    my $timet = parsedate ($1);
    next unless (defined $timet);
    print "$timet\t2\t$ppc\t$board\t$hits\n";
  } elsif ($line =~ /DPR\ FULL/) {
    die "unknown line format for $line\n" unless ($line =~ /^(.*)\ ppc(\d+).*board (\d+)/);
    my ($ppc, $board) = ($2, $3);
    my $timet = parsedate ($1);
    next unless (defined $timet);
    print "$timet\t3\t$ppc\t$board\t0\n";
  } elsif ($line =~ /mismatch,\ disabling/) {
    die "unknown line format for $line\n" unless ($line =~ /^(.*)\ ppc(\d+).*board\ (\d+)/);
    my ($ppc, $board) = ($2, $3);
    my $timet = parsedate ($1);
    next unless (defined $timet);
    print "$timet\t4\t$ppc\t$board\t0\n";
  } elsif ($line =~ /Disabling/) {
    die "unknown line format for $line\n" unless ($line =~ /^(.*)\ ppc(\d+).*board=(\d+)/);
    my ($ppc, $board) = ($2, $3);
    my $timet = parsedate ($1);
    next unless (defined $timet);
    print "$timet\t10\t$ppc\t$board\t0\n";
  } elsif ($line =~ /old\ trigger/) {
    die "unknown line format for $line\n" unless ($line =~ /^(.*)\ ppc(\d+).*fifo\ (\d+)/);
    my ($ppc, $board) = ($2, $3);
    my $timet = parsedate ($1);
    next unless (defined $timet);
    print "$timet\t12\t$ppc\t$board\t0\n";
  } else {
    #print "unknown $line\n";
  }
}

