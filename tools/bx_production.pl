#!/usr/bin/perl -w
# Author: Stefano Davini <stefano.davini@ge.infn.it>
#
# Production utility: this script is used to run the Echidna data production procedure:
# precalibrations, electronics_calibration, echidna (production).
# Several checks are done before starting production.

use POSIX qw(strftime _exit);
use Pg;
use Cwd;
use Time::ParseDate;
use Time::Local;
use Net::FTP;
use File::Temp;
use strict;
use warnings;
use threads;

my $user = $ENV{USER};
my $pwd  = $ENV{PWD};
my $home = $ENV{HOME};

# check if we are in a HOME subfolder
unless ($pwd =~ /^${home}/){
	die " ERROR: bx_production works only in subforlders of home directory ($home}), while you are in $pwd\n";
}

# check if echidna is executable
unless (-x "echidna"){
        die " ERROR: Echidna is not installed\n";      
}

# check if qsub is executable
unless (-x "/usr/bin/qsub"){
        die " ERROR: qsub is not executable\n";      
}

# check if bx_repository.pl is present
unless (-e 'bx_repository.pl'){
	die " ERROR: bx_repository.pl is not here ($pwd)\n";
}

# gobal variables, default values
my $check_raw				= 1; # bx_repository.pl check (rawdata is ok)
my $check_type				= 1; # check run type  (normal/calibration)
my $check_gzip				= 0; # check zipped file integrity
my $run_precalibrations			= 1; # run precalibrations
my $run_electronics_calibrations	= 1; # run electronics calibrations
my $run_laser_calibrations		= 0; # run laser calibrations
my $run_laser_calibrations_validate	= 0; # run laser calibrations validate
my $run_echidna 			= 1; # run echidna (production)
my $run_fadc_reader 			= 1; # run fadc_reader 
my $run_validation 			= 0; # run validation (to be implemented)
my $run_number				= 0; # run number to process
my $num_events				= 0; # number of event to process; if 0, process all events 
my $config				= ''; # Echidna configuration flag; 'source' for source runs
my $type				= 'normal'; # run type (normal/calibration)
my $queue				= ''; # batch queue; default becomes '-q production'
my $sleep				= 0; # sleep before checking rawdata (used wait rawfile to be copied)
my $retry				= 0; # retry checking of rawdata (used to wait rawfile to be copied)
my $db_host 				= "bxdb"; # database host
my $db_port 				= 5432; # database port
my $db_name				= 'daq_config';
my $db_user 				= 'borex_guest';
my $sendmail				= 0; # send mail to shifters and RC when the process is finished
my $toaddr				= ''; # will be filled with shifters e-mail


# read arguments
foreach my $arg (@ARGV){
        if (($arg =~ /^help$/i) || ($arg =~/^usage$/i)){
		usage();
		exit;
	}
	elsif ($arg =~ /^check_raw\=(\d+)/){
		$check_raw = $1;
	}
	elsif ($arg =~ /^check_type\=(\d+)/){
		$check_type = $1;
	}
	elsif ($arg =~ /^check_gzip\=(\d+)/){
		$check_gzip = $1;
	}
	elsif ($arg =~ /^run_precalibrations\=(\d+)/){
		$run_precalibrations = $1;
	}
	elsif ($arg =~ /^run_electronics_calibrations\=(\d+)/){
		$run_electronics_calibrations = $1;
	}
	elsif ($arg =~ /^run_laser_calibrations\=(\d+)/){
		$run_laser_calibrations = $1;
	}
	elsif ($arg =~ /^run_laser_calibrations_validate\=(\d+)/){
		$run_laser_calibrations_validate = $1;
	}
	elsif ($arg =~ /^run_echidna\=(\d+)/){
		$run_echidna = $1;
	}
	elsif ($arg =~ /^run_validation\=(\d+)/){
		$run_validation = $1;
	}
	elsif ($arg =~ /^run_fadc_reader\=(\d+)/){
		$run_fadc_reader = $1;
	}
	elsif ($arg =~ /^run_number\=(\d+)/){
		$run_number = $1;
	}
	elsif ($arg =~ /^num_events\=(\d+)/){
		$num_events = $1;
	}
	elsif ($arg =~ /^config\=(\w+)/){
		$config = $1;
	}
	elsif ($arg =~ /^type\=(\w+)/){
		$type = $1;
	}
	elsif ($arg =~ /^queue\=(\w+)/){
		$queue = ' -q '.$1;
	}
	elsif ($arg =~ /^sleep\=(\d+)/){
		$sleep = $1;
	}
	elsif ($arg =~ /^retry\=(\d+)/){
		$retry = $1;
	}
	elsif ($arg =~ /^db_port\=(\d+)/){
		$db_port = $1;
	}
	elsif ($arg =~ /^db_host\=(\w+)/){
		$db_host = $1;
	}
	elsif ($arg =~ /^sendmail\=(\d+)/){
		$sendmail = $1;
	}
	else {
		die " ERROR: unknown option $arg \n";
	}
}

# check arguments
if (!($run_number)){
	warn "ERROR: run number is not specified;\n";
	usage();
	exit();
}

# check run Type (normal/calibration)
if ($check_type){
	my $db_connection = Pg::connectdb("host=$db_host port=$db_port dbname=$db_name user=$db_user");
	die "Unable to connect to the database\n" if ($db_connection->status != PGRES_CONNECTION_OK);
	#print "Connected do database\n";
	my $query = "SELECT \"Type\" FROM \"Run\" WHERE \"RunNumber\"=$run_number  ";
	my $result = $db_connection->exec($query);
	die "Error in query: " . $db_connection->errorMessage . "\n" unless ($result->resultStatus == PGRES_TUPLES_OK);
	my $ntuples = $result->ntuples;
	my $fields  = $result->nfields;
	my $type_   = $result->fnumber ("\"Type\"");
	$type    = $result->getvalue(0, $type_);
	unless (($type eq 'normal') || ($type eq 'calibration')){
		die "ERROR: unknow type $type for run $run_number\n"; 
	}
	if ($type eq "calibration"){
		$run_precalibrations = 1;
		$run_electronics_calibrations = 1;
		$run_laser_calibrations = 1;
		$run_laser_calibrations_validate = 1;
		$run_echidna = 0;
		$run_fadc_reader = 0;
		$run_validation = 0;
	}
}

# check the config option flag
die unless (($config eq '') || ($config eq 'source') || ($config eq 'normal'));
if ($config eq 'normal'){
	$config = '';   # the normal configuration of Echidna has no flags
}

# if Echidna rootfile exists on $PWD, do not proceed
my $zeros  = ($run_number>9999) ? '0' : '00';
my $suffix = ($config eq 'source') ? '_source_c14.root' : '_c14.root';
my $run_string = 'Run'.$zeros.$run_number.$suffix;
if (-e $run_string){
	warn "WARNING: $run_string is already present here\n";
 	warn "If you want to re-process this run, please rename it, or remove it, or move it in another folder\n";
 	exit();
}
 

# check the raw data file
if ($check_raw)	{
	my $is_ok = 0;
	while ($retry>=0){ # retry = 0 on default, do this cicle at least once
		# sleep in order for the rawdatafile to be copied on storage
		sleep($sleep);
		# bx_repository.pl check
		print "Checking raw data files for run $run_number: ";
		my $ret = `./bx_repository.pl check $run_number 2>&1`;
		# bx_repository prints on  standard error, 2&1 redirect on std output
		# the output of bx_repository should be i.e '5003 ok' on normal conditions
		if ($ret =~ /(\d+)(\s+)(\w+)/){ 
			my $ch_run_number = $1;
			my $ch_ok = $3;
			print "$ch_run_number $ch_ok\n";
			if (($ch_run_number eq $run_number) && ($ch_ok eq 'ok')){
				# the rawdata file is present on bxstorage
				$is_ok = 1;
				last;
			}
			else {
				$retry--;
			}
		}
		else{
			$retry--;
		}
	} # end of while

	die "ERROR: Rawdata file for run number $run_number is not on bxstorage or may have problems\n" unless ($is_ok eq 1);

	if ($check_gzip){
		# now check the compressed (g-zipped) rawfile integrity
		print "Checking the compresed raw data file integrity\n";
		# I need the start time of the run to get year and week to create path
		my $start_time = get_run_start_time ($run_number);
		die "ERROR: undefined start $start_time time for run $run_number\n" unless (defined $start_time);
		my ($year_index, $week_index) = get_run_start_times ($start_time);
		my $rawfile = '/bxstorage/rawdata/'.$year_index.'/'.$week_index.'/Run'.$zeros.$run_number.'_0?.out.gz';
		$is_ok = 0;
		$retry = 5;
		while ($retry--){
			my $ret = `gzip -t $rawfile 2>&1`;
			if ($ret eq '') {
				$is_ok = 1;
				last;
			}
			sleep($sleep);
		}
		die "ERROR: $rawfile integrity check failed\n" unless ($is_ok eq 1);
	}
}

# set the mail address (shifters + run coordinator)
if ($sendmail){
	# i need the start date of this run in order to know who are the shifters 
	my $start_time = get_run_start_time ($run_number);
	warn "WARNING: undefined start $start_time time for run $run_number\n" unless (defined $start_time);
	my ($year_index, $week_index, $month_index, $day_index) = get_run_start_times ($start_time);
	my $dayzero = ($day_index<10) ? '0' : ''; 
	my $monzero = ($month_index<10) ? '0' : ''; 
	my $date = $year_index.'-'.$monzero.$month_index.'-'.$dayzero.$day_index;

	# i need to know the user names of captain 1 and 2
	my ($cap1, $cap2) = get_shifters ($date);

	# i need the email address of the shifters
	my $eml1 = get_email ($cap1);
	my $eml2 = get_email ($cap2);

	# here is my address list
	$toaddr= $toaddr.$eml1.' '.$eml2;
}


# set the batch queue
if (($queue eq '') && ($user eq 'production')){
	$queue = '-q production';
}

# set job name
my $jobname = 'production_'.$run_number;

# temp folder for pbs script file
my $tempfile = File::Temp::mktemp( "/tmp/production.XXXX" ); 

# Write pbs script
my $FILE;
open  $FILE, ">", $tempfile or die " ERROR: can't open $tempfile \n";
print $FILE '#!/bin/sh'."\n";
print $FILE '#set -x'."\n";
#if [ -n "$root" ]; then choose_root $root; fi
print $FILE 'hostname'."\n";
print $FILE 'date >&2'."\n";

if ($sendmail){
	print $FILE 'echo " " | mailx -s "Run '.$run_number.': production started" bxruncoord@lngs.infn.it '.$toaddr."\n";
}
if ($run_precalibrations){
	print $FILE './echidna -f run://'.$run_number.' precalibrations'."\n";
	print $FILE 'date >&2'."\n";
}
if ($run_electronics_calibrations){
	print $FILE './echidna -f run://'.$run_number.' electronics_calibrations'."\n";
	print $FILE 'date >&2'."\n";
}
if ($run_laser_calibrations){
	print $FILE './echidna -f run://'.$run_number.' laser_calibrations'."\n";
	print $FILE 'date >&2'."\n";
}
if ($run_laser_calibrations_validate){
	print $FILE './echidna -f run://'.$run_number.' laser_calibrations_validate'."\n";
	print $FILE 'date >&2'."\n";
}
if ($run_echidna){
	my $ev = ($num_events) ? " -e $num_events " : '';
	print $FILE './echidna -f run://'.$run_number." $ev $config \n";
	print $FILE 'date >&2'."\n";
}
if ($run_fadc_reader){
	print $FILE './fadc_reader '.$run_number." \n";
	print $FILE 'date >&2'."\n";
}
if ($run_validation){
	# to be implemented
}
if ($sendmail){
	print $FILE 'echo " " | mailx -s "Run '.$run_number.': production finished" bxruncoord@lngs.infn.it '.$toaddr."\n";
}
close $FILE;

print "Submitting:\n";
if ($run_precalibrations) { print "- precalibrations\n";}
if ($run_electronics_calibrations) { print "- electronics_calibrations\n";}
if ($run_laser_calibrations) { print "- laser_calibrations\n";}
if ($run_laser_calibrations_validate) { print "- laser_calibrations_validate\n";}
if ($run_echidna) { print "- echidna\n";}
if ($run_fadc_reader) { print "- fadc_reader\n";}
if ($run_validation) { print "- run_validation\n"};
print "for run number $run_number ($type run)";
if ($num_events) { print " for $num_events events"};
if ($config) { print "\nconfiguration: $config"};
print "\n";

# Submit the job
system ("qsub -N $jobname -d $pwd $queue $tempfile");

sub usage{
	print "bx_production.pl [global options]\n";
	print "Global options:\n";
	print "\trun_number=[run]\t\t\t Set the run number to process\n";
	print "\tcheck_raw=[0/1]\t\t\t\t Optional: check if the rawdata file is on bxstorage before doing the production (default 1)\n";
	print "\tcheck_type=[0/1]\t\t\t Optional: check the run type, force arguments for calibration runs (default 1)\n";
	print "\trun_precalibrations=[0/1]\t\t Optional: run precalibrations (default 1)\n";
	print "\trun_electronics_calibrations=[0/1]\t Optional: run electronics_calibrations (default 1)\n";
	print "\trun_laser_calibrations=[0/1]\t\t Optional: run laser_calibrations (default 0)\n";
	print "\trun_laser_calibrations_validate=[0/1]\t Optional: run laser_calibrations_validate (default 0)\n";
	print "\trun_echidna=[0/1]\t\t\t Optional: run echidna (default 1)\n";
	print "\trun_fadc_reader=[0/1]\t\t\t Optional: run fadc_reader (default 1)\n";
	print "\tnum_events=[num]\t\t\t Optional: set the number of event to process in echidna; if not specified, process all events\n";
	print "\tconfig=[config]\t\t\t\t Optional: set the configuration flag of echidna; set config=source for source runs; (default '')\n";
	print "\ttype=[num]\t\t\t\t Optional: set the run type (normal/calibration); if check_type=1 the run will be check via db\n";
	print "\tqueue=[queue]\t\t\t\t Optional: set the batch queue for the job; if not specified, the job is submitted in the 'production' queue\n";

}

sub get_run_start_times {
	my $start_time = shift;
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday, $isdst) = localtime $start_time;
	my $week_time = $start_time - (($wday * 24 + $hour) * 60 + $min) * 60 + 4800;
	my @tmp = localtime $week_time;
	my $week_index = strftime ("%b_%d", @tmp); 
	my $year_index = $tmp[5] + 1900;
	return ($year_index, $week_index, $mon+1, $mday);
}

sub get_run_start_time {
	my $db_connection = Pg::connectdb("host=$db_host port=$db_port dbname=$db_name user=$db_user");
	die "Unable to connect to the database\n" if ($db_connection->status != PGRES_CONNECTION_OK);
	my $run_id = shift;
	my $query = "SELECT * FROM \"Run\" WHERE \"RunNumber\"=$run_id;";
	my $result = $db_connection->exec($query);
	die "Error in query: " . $db_connection->errorMessage . "\n" unless ($result->resultStatus == PGRES_TUPLES_OK);
	return undef unless ($result->ntuples);
	my $start_date_pos = $result->fnumber ("\"StartTime\"");
	my $start_date = $result->getvalue (0, $start_date_pos);

	$query = "SELECT EXTRACT(EPOCH FROM TIMESTAMP '$start_date')";
	$result = $db_connection->exec($query);
	die "Error in query: " . $db_connection->errorMessage . "\n" unless ($result->resultStatus == PGRES_TUPLES_OK);
	return undef unless ($result->ntuples == 1);
	my $start_time_pos = $result->fnumber ("date_part");
	my $start_time = $result->getvalue (0, $start_time_pos);

	return $start_time;
}

sub get_shifters {
	my $db_connection = Pg::connectdb("host=$db_host port=$db_port dbname=bx_shift user=$db_user");
	die "Unable to connect to the database\n" if ($db_connection->status != PGRES_CONNECTION_OK);
	#print "Connected do database\n";
	my $date = shift;
	my $query = "SELECT \"cap1\",\"cap2\" FROM shifts WHERE '$date' BETWEEN beg AND fin";
	my $result = $db_connection->exec($query);
	die "Error in query: " . $db_connection->errorMessage . "\n" unless ($result->resultStatus == PGRES_TUPLES_OK);
	my $ntuples = $result->ntuples;
	my $fields  = $result->nfields;
	my $cap1_   = $result->fnumber ("\"cap1\"");
	my $cap2_   = $result->fnumber ("\"cap2\"");
	my $cap1    = $result->getvalue(0, $cap1_);
	my $cap2    = $result->getvalue(0, $cap2_);

	return ($cap1, $cap2);
}

sub get_email {
	my $db_connection = Pg::connectdb("host=$db_host port=$db_port dbname=bx_shift user=$db_user");
	die "Unable to connect to the database\n" if ($db_connection->status != PGRES_CONNECTION_OK);
	my $shifter = shift;
	my $query = "SELECT \"email\" FROM users WHERE \"login\"='$shifter'";
	my $result = $db_connection->exec($query);
	warn "Error in query: " . $db_connection->errorMessage . "\n" unless ($result->resultStatus == PGRES_TUPLES_OK);
	return undef unless ($result->ntuples);
	my $email_pos = $result->fnumber ("\"email\"");
	my $email = $result->getvalue (0, $email_pos);

	return $email;
}
