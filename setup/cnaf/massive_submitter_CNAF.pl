#!/usr/bin/perl
#
# Call this script via nohup:    i.e. $ nohup ./DST_massive_submitter.pl &
# This script checks the number of reprocess jobs pending/running;
# If the number of jobs  under a certain threshold, a new process is submitted.
# Based on S. Davini script.
# Revised: Alessandra Re (August 2013)
# 

use warnings;
use strict;

my $process = shift;
die "\n\t Usage: ./massive_submitter_CNAF.pl input_file.sh\n\n" if(!$process);

open INPUT, $process or die $!;
my @commands = <INPUT>;
close INPUT;

my $sleep_time = 100;

while ($#commands >= 0){
	my $running = int(`bjobs | grep -c RUN`);
	my $pending = int(`bjobs | grep -c PEND`);

    	if(($running < 50) && ($pending < 10)) {
		my $new_process = shift(@commands);
		print "\t Submitting $new_process";
		system("$new_process");
		sleep 10;
    	} else{
		my $still = $#commands+1;
		print "\t RUNNING: $running, PENDING: $pending; still $still process to submit.\n";
		sleep $sleep_time;
    	}
}
