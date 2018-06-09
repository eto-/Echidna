#!/usr/bin/perl
# Author: Stefano Davini <stefano.davini@ge.infn.it>
#
# reprocess_submitter.pl: call this script via nohup
# i.e. $ nohup ./reproces_submitter.pl &
# This script checks the number of reprocess jobs pending/running;
# if the number of jobs  under a certain threshold,
# a new reprocess command is executed.
# The input file containing bx_process.pl command is $input;
# this input file is usually generated via the script GenBxProcessCmd.pl
#

use warnings;
use strict;

my $name ='production';
my $here = $ENV{PWD};

my $input = "$here/process_normal_runs.sh";
my $INP;

open $INP, "<", $input or die "Error: can't open $input\n";
my @submit_cmds = <$INP>;
close $INP;

my $sleep_time = 180;

while ((scalar(@submit_cmds))>0){

    my $run  = int(`qstat | grep $name | grep reprocess | grep -c R `);
    my $qued = int(`qstat | grep $name | grep reprocess | grep -c Q `);

    if (($run<60) && ($qued<30)){
	
	my $newprocess = shift (@submit_cmds);
	print "Submitting $newprocess";
	system("$newprocess");
	sleep 1;
    }
    else{
	my $todo = scalar(@submit_cmds);
	print "R: $run, Q: $qued; $todo commands to submit; Sleeping $sleep_time seconds... \n";
	sleep $sleep_time;
    }
}
