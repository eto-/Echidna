#!/usr/bin/perl -w

# Author: Alessandra Re <alessandra.re@mi.infn.it>
# Last revision: 2015-03-24

use DBI;
use POSIX qw(strftime _exit);
use Time::ParseDate;
use Time::Local;
use strict;

# Routines
sub get_proper_cycle;
sub get_proper_place;
sub get_valid_events;
sub is_source_run;

# Run to be analyzed
my $run_id = shift;
die "\n\tUsage:\t$0 <run_number>\n\n" if(!$run_id);

# Database parameters
my $db_connection = 0;
my $db_host = "bxdb.lngs.infn.it";
my $db_port = 5432;
my $db_user = "borex_guest";
my $db_passwd = "";

my $source_flag = is_source_run($run_id);
my $cycle = get_proper_cycle();
my ($number_of_events, $duration) = get_valid_events($run_id, $source_flag);
if(($source_flag) && ($number_of_events)) {
	print "\n\tThe number of valid events for calibration source run $run_id ($duration s) is $number_of_events\n";
} elsif(!$number_of_events) {
	die "\n\tNumber of validated events not avaliable.\n\n";
} else {
	print "\n\tThe number of valid events for run $run_id ($duration s) is $number_of_events\n";
}
	
my ($yyyy,$mm,$day, $year_index, $week_index, $hour, $min, $sec) = get_proper_place($run_id);
my $RunTime = $hour.':'.$min.':'.$sec;
my $RunDate = $yyyy."-".$mm."-".$day;
my $RunPath = "/storage/gpfs_data/borexino/rootfiles/cycle_".$cycle."/".$year_index."/".$week_index."/";
$RunPath = $RunPath."ancillary/" if($source_flag);

$run_id = '0'.$run_id if(int($run_id) < 10000);
my $RunName = $RunPath.'Run0'.$run_id.'_c'.$cycle.'.root';
$RunName = $RunPath.'Run0'.$run_id.'_source_c'.$cycle.'.root' if ($source_flag);

print "\n\tRun date: $RunDate\t Run start time: $RunTime\n\t-->  Position: $RunName\n\n";


#########################################################
sub get_valid_events {
	my $run_id = shift;
	my $source_flag = shift;
	my $dbname_id = 'bx_runvalidation';
  my $db_data = "SELECT \"LastValidEvent\", \"Duration\" FROM \"ValidRuns\" WHERE \"RunNumber\"=$run_id;";
	$db_data = "SELECT \"LastValidEvent\", \"Duration\" FROM \"SourceValidRuns\" WHERE \"RunNumber\"=$run_id;" if($source_flag);
	$db_connection = DBI->connect("DBI:Pg:dbname=$dbname_id;host=$db_host;port=$db_port", "$db_user", "$db_passwd");
	die "\n\t Error! Database connection absent: $DBI::errstr\n\n" if(!$db_connection);
	my $query = $db_connection->prepare($db_data);
	$query->execute();
	die "\n\tError in query: $DBI::errstr\n\n" if(!$query);
	my ($number_of_events, $duration) = $query->fetchrow_array;
	$query->finish;
	$db_connection->disconnect;
	return ($number_of_events, $duration);
}

sub get_proper_place {
	my $run_id = shift;
	$db_connection = DBI->connect("DBI:Pg:dbname=daq_config;host=$db_host;port=$db_port", "$db_user", "$db_passwd");
 	die "\n\t Error! Database connection absent: $DBI::errstr\n\n" if(!$db_connection);
	my $db_data = "SELECT \"StartTime\" FROM \"Run\" WHERE \"RunNumber\"=$run_id;";
	my $query = $db_connection->prepare($db_data);      
	$query->execute();
 	die "\n\tError in query: $DBI::errstr\n\n" if(!$query);
	
	my @abbr = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );

	my $RunDate = $query->fetchrow_array;
	# EXTRACT YEAR_WEEK INDEX
	my ($yyyy,$mm,$day) = $RunDate =~ /(\d+)-(\d+)-(\d+)/;
	my ($hour,$min,$sec) = $RunDate =~ /(\d+):(\d+):(\d+)/;
	my $month = $mm-1;
	my $year = $yyyy-1900;
	my $start_time = timelocal($sec,$min,$hour,$day,$month,$year);
	my $wday = (localtime($start_time))[6];
	my $week_time = $start_time - (($wday * 24 + $hour) * 60 + $min) * 60 + 4800;
	my @tmp = localtime $week_time;
	my $right_day = sprintf("%.2d",$tmp[3]);
	my $week_index = "$abbr[$tmp[4]]\_$right_day";
	my $year_index = $tmp[5]+1900;
	
	$query->finish;
	$db_connection->disconnect;
	return($yyyy,$mm,$day, $year_index, $week_index, $hour, $min, $sec);
}

sub get_proper_cycle {
	$db_connection = DBI->connect("DBI:Pg:dbname=bx_runvalidation;host=$db_host;port=$db_port", "$db_user", "$db_passwd");
	die "\n\t Error! Database connection absent: $DBI::errstr\n\n" if(!$db_connection);
	my $db_data = "SELECT \"ProductionCycle\" FROM \"ValidRuns\" ORDER BY \"ValidationDate\" DESC limit 1;";
	my $query = $db_connection->prepare($db_data);
	$query->execute();
	die "\n\tError in query: $DBI::errstr\n\n" if(!$query);
	my ($cycle) = $query->fetchrow_array;
	$query->finish;
	$db_connection->disconnect;
	return ($cycle);
}

sub is_source_run {
	my $run_id = shift;
	my $source_flag;
	$db_connection = DBI->connect("DBI:Pg:dbname=bx_runvalidation;host=$db_host;port=$db_port", "$db_user", "$db_passwd");
	die "\n\t Error! Database connection absent: $DBI::errstr\n\n" if(!$db_connection); 
  my $db_data1 = "SELECT \"LastValidEvent\" FROM \"ValidRuns\" WHERE \"RunNumber\"=$run_id;";
  my $query_check1 = $db_connection->prepare($db_data1);
  $query_check1->execute();
  my $check1 = $query_check1->fetchrow_array;
  $query_check1->finish;
  if ($check1) {
		$source_flag = 0;
	} else {
    my $db_data2 = "SELECT \"LastValidEvent\" FROM \"SourceValidRuns\" WHERE \"RunNumber\"=$run_id;";
    my $query_check2 = $db_connection->prepare($db_data2);
    $query_check2->execute();
    my $check2 = $query_check2->fetchrow_array;
    $query_check2->finish;
    if ($check2) {
      $source_flag = 1;
    } else {
			$db_connection->disconnect;
      die "\n\t Error! Run $run_id is NOT a valid normal or calibration source run!\n\n";
    }
  }
	$db_connection->disconnect;
	return $source_flag;
}
