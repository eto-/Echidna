#!/usr/bin/perl
# Author: Stefano Davini <Stefano.Davini@ge.infn.it>
# 
# 
# gen_bx_process_cmd.pl: generate bx_process.pl commands over stepped pair of runs
# Arguments are 
# first run: the run number of the firs run to be processed 
# last run : the run number of the last run to be processed
# step: the step between pair of runs to be processed 
# output file (optional): the name of the output file with bx_process.pl commands 
#

use warnings;
use strict;


my $outfile = 'process_normal_runs.sh';
if (scalar(@ARGV)<3){
	die "gen_bx_process_cmd.pl [first run] [last run] [step] [output file (default $outfile)] \n";
}

my $first = int($ARGV[0]);
my $last  = int($ARGV[1]); 
my $step  = int($ARGV[2]);

if (!$first){
	die "Error: $ARGV[0] is not a valid argument for 'first run number' \n";
}
if (!$last){
	die "Error: $ARGV[1] is not a valid argument for 'last run number' \n";
}
if (!$step){
	die "Error: $ARGV[2] is not a valid argument for 'step' \n";
}
if ($first>=$last){
	die "Error: the last run number must be > than the first run number \n";
}

if ($ARGV[3]){
	$outfile = $ARGV[3];
}

my $OUT;
open $OUT, ">", $outfile or die "Error opening $outfile\n";

for (my $i=$first; $i<$last; $i+=$step){
	my $thisfirst = $i;
	my $thislast  = ($i + $step-1<$last) ? ($i + $step-1) : $last;

	my $cmd = "./bx_process.pl =reprocess =rewrite_files =valid $thisfirst-$thislast pbs_echidna.sh \n";
	print $OUT $cmd;

}

close $OUT;
print "Reprocess commands written on $outfile\n";
