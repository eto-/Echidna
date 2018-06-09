#!/usr/bin/perl -w
#
# Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
#
# $Id: bx_process.pl,v 1.20 2011/03/10 14:14:06 razeto Exp $
#
# bx_process.pl: run a specified command over a run list (using only real
# runs)
#
use POSIX qw(strftime _exit);
use Pg;
use Time::ParseDate;
use Time::Local;
use strict;

# Connection handlers
my %db_connections;

# Defaults
my $do_nothing = 0;
my $process_no_files = 0;
my $run_type_selection = undef;
my $valid = 0;
my $reprocess = 0;
my $rewrite_files = "";
my $queue = undef;

my $db_host = "bxdb";
my $db_port = 5432;
my %db_names = ( 
  daq_config => "daq_config",
  bx_runvalidation => "bx_runvalidation"
);
my $db_user = "borex_guest";
my $db_password = "";
my $output_folder = undef;
my $sleep_sec = 2;
my $max_events = undef;

while (defined $ARGV[0] and $ARGV[0] =~ /^=/) {
  $_ = shift;
  last if /^==/;
  if (/^=do_nothing$/) {
    $do_nothing = 1;
  } elsif (/^=process_no_files$/) {
   $process_no_files = 1;
  } elsif (/^=reprocess$/) {
    $reprocess = 1;
  } elsif (/^=no_reprocess$/) {
    $reprocess = -1;
  } elsif (/^=rewrite_files$/) {
    $rewrite_files = "=rewrite_files";
  } elsif (/^=valid$/) {
    $valid = 1;
    usage () if (defined $run_type_selection and $run_type_selection !~ /^normal$/);
    $run_type_selection = "normal";
  } elsif (/^=queue$/) {
    $queue = shift;
  } elsif (/^=max_events$/) {
    $max_events = shift;
    die "=max_events is not an integer" unless ($max_events == int($max_events));
    test_opt ("=max_events", $max_events);
  } elsif (/^=output$/) {
    $output_folder = shift;
    die "$output_folder is not a directory\n" unless (-d $output_folder);
    test_opt ("=output_folder", $output_folder);
  } elsif (/^=sleep$/) {
    $sleep_sec = shift;
    test_opt ("=sleep", $sleep_sec);
    die "sleep argument ($sleep_sec) must be an integer\n" if ($sleep_sec =~ /\D/)
  } elsif (/^=run_type$/) {
    $run_type_selection = shift;
    test_opt ("=run_type", $run_type_selection);
    usage () if ($valid and $run_type_selection !~ /^normal$/);
  } elsif (/^=db_host$/) {
    $db_host = shift;
    test_opt ("=db_host", $db_host);
#  } elsif (/^=db_name$/) {
#    $db_name = shift;
#    test_opt ("=db_name", $db_name);
  } elsif (/^=db_user$/) {
    $db_user = shift;
    test_opt ("=db_user", $db_user);
  } elsif (/^=db_password$/) {
    $db_password = shift;
    test_opt ("=db_password", $db_password);
  } else {
    usage ($ARGV[0]);
  } 
}

$reprocess = 1 if ($valid and $reprocess == 0);
$output_folder = "/scratch" if ($reprocess > 0 and not $output_folder);

my $run_list_pattern = "";
while (defined $ARGV[0]) {
  $_ = shift;
  last if (/[^0-9\-,]/);
  $run_list_pattern .= $_;
}

unless ($run_list_pattern) {
  my $run_file = $_;
  if (defined $run_file and -r $run_file) {
    $run_list_pattern = read_run_file ($run_file);
    $_ = shift;
  }
}

usage ($_) unless ($run_list_pattern);
if (/^echidna$/ or /^pbs_echidna.sh$/) {
  process ($run_list_pattern, $_, join (" ", @ARGV));
} else {
  usage ($_);
}


sub process {
  my $run_list_pattern = shift;
  my $program = shift;
  my $program_argument = shift;
  die "echidna syntax: run_list echidna echidna_options\n" unless (defined $run_list_pattern); 

  connect_db ();

  my @run_list = extract_unique_run_list ($run_list_pattern);
  foreach my $run (@run_list) {
    process_run ($run, $program, $program_argument);
  } 
}

sub process_run {
  my $run = shift;
  my $program = shift;
  my $program_argument = shift;

  my $start_time = get_run_start_time ($run);
  return unless (defined $start_time);

  my $slices = get_run_slices ($run);
  return unless ($slices > 0 or $process_no_files);

  my $run_type = get_run_type ($run);
  return if (defined ($run_type_selection) && $run_type !~ /^$run_type_selection$/);

  my $nevents = 0;
  $nevents = get_valid_events ($run) if ($valid);
  $nevents = $max_events if (defined $max_events and ($nevents <= 0 or $max_events < $nevents));
  return if ($nevents < 0);

  my $specifier;
  if ($program =~ /echidna/) {
    die "program argument ($program_argument) contains reserved options [fO]\n" if ($program_argument =~ /-[fO]/);
    $specifier = "-f run://" . $run . " ";
    if ($nevents) {
      $specifier .= "-e " . $nevents . " ";
      die "program argument ($program_argument) contains reserved options e when using =valid\n" if ($program_argument =~ /-e/);
    }
    $specifier .= "-O " . $output_folder . " " if ($output_folder);
  } else {
    die "unknow run specifier syntax";
  }

  $specifier = "=reprocess " . $rewrite_files . " " . $specifier if ($program =~ /pbs_echidna.sh/ and $reprocess > 0);
  $specifier = "=queue " . $queue . " " . $specifier if ($program =~ /pbs_echidna.sh/ and $queue);

  my $cmd = "./" . $program . " " . $specifier . $program_argument;

  $cmd =~ s/%r/$run/g;
  $cmd =~ s/%t/$run_type/g;

  if ($do_nothing) { print "executing " . $cmd . "\n"; } 
  else { system $cmd; }
  sleep $sleep_sec if ($sleep_sec > 0); # avoid submitting jobs too fast
}



sub test_opt {
  my $opt_name = shift;
  my $opt_value = shift;
  die "$opt_name requires an argument ($opt_value)\n" if (!defined $opt_value || $opt_value =~ /^=/);
}

sub connect_db {
  if (scalar (%db_connections)) {
    warn "Already connected to db\n";
    return;
  }
  $db_connections{daq_config} = Pg::connectdb("host=$db_host port=$db_port dbname=$db_names{daq_config} user=$db_user");
  die "Unable to connect to the database\n" if ($db_connections{daq_config}->status != PGRES_CONNECTION_OK);
  $db_connections{bx_runvalidation} = Pg::connectdb("host=$db_host port=$db_port dbname=$db_names{bx_runvalidation} user=$db_user");
  die "Unable to connect to the database\n" if ($db_connections{daq_config}->status != PGRES_CONNECTION_OK);
}

sub get_last_run {
  my $query = "SELECT \"RunNumber\" FROM \"Run\" ORDER by \"RunNumber\" DESC;";
  my $result = $db_connections{daq_config}->exec($query);
  die "Error in query: " . $db_connections{daq_config}->errorMessage . "\n" unless ($result->resultStatus == PGRES_TUPLES_OK);
  return undef unless ($result->ntuples);
  my $run_number_pos = $result->fnumber ("\"RunNumber\"");
  my $run_number = $result->getvalue (0, $run_number_pos);

  return $run_number;
}

sub get_run_type {
  my $run_id = shift;
  my $query = "SELECT \"Type\" FROM \"Run\" WHERE \"RunNumber\"=$run_id;";
  my $result = $db_connections{daq_config}->exec($query);
  die "Error in query: " . $db_connections{daq_config}->errorMessage . "\n" unless ($result->resultStatus == PGRES_TUPLES_OK);
  return undef unless ($result->ntuples);
  my $type_pos = $result->fnumber ("\"Type\"");
  my $type = $result->getvalue (0, $type_pos);
  return $type;
}

sub get_run_slices {
  my $run_id = shift;
  my $query = "SELECT \"NumberOfFiles\" FROM \"Run\" WHERE \"RunNumber\"=$run_id;";
  my $result = $db_connections{daq_config}->exec($query);
  die "Error in query: " . $db_connections{daq_config}->errorMessage . "\n" unless ($result->resultStatus == PGRES_TUPLES_OK);
  return undef unless ($result->ntuples);
  my $number_of_files_pos = $result->fnumber ("\"NumberOfFiles\"");
  my $number_of_files = $result->getvalue (0, $number_of_files_pos);
  return $number_of_files;
}

sub get_run_start_time {
  my $run_id = shift;
  my $query = "SELECT * FROM \"Run\" WHERE \"RunNumber\"=$run_id;";
  my $result = $db_connections{daq_config}->exec($query);
  die "Error in query: " . $db_connections{daq_config}->errorMessage . "\n" unless ($result->resultStatus == PGRES_TUPLES_OK);
  return undef unless ($result->ntuples);
  my $start_date_pos = $result->fnumber ("\"StartTime\"");
  my $start_date = $result->getvalue (0, $start_date_pos);

  $query = "SELECT EXTRACT(EPOCH FROM TIMESTAMP '$start_date')";
  $result = $db_connections{daq_config}->exec($query);
  die "Error in query: " . $db_connections{daq_config}->errorMessage . "\n" unless ($result->resultStatus == PGRES_TUPLES_OK);
  return undef unless ($result->ntuples == 1);
  my $start_time_pos = $result->fnumber ("date_part");
  my $start_time = $result->getvalue (0, $start_time_pos);
  
  return $start_time;
}

sub get_run_start_times {
  my $start_time = shift;
  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday) = localtime $start_time;
  my $week_time = $start_time - ((($wday - 1) * 24 + $hour) * 60 + $min) * 60;
  my @tmp = localtime $week_time;
  my $week_index = strftime ("%b_%d", @tmp); 
  my $year_index = $tmp[5] + 1900;
  return ($year_index, $week_index);
}

sub extract_unique_run_list {
  my $run_list_pattern = shift;
  my @run_list;
  my @fields = split (/,/, $run_list_pattern);
  foreach my $field (@fields) {
    if ($field =~ /^(\d+)?-(\d+)?$/) {
      my $low_limit = $1;
      my $high_limit = $2;
      $low_limit = 1200 unless (defined $low_limit);
      $high_limit = get_last_run () unless (defined $high_limit);
      warn "Ignoring range $field ($high_limit < $low_limit)\n" if ($high_limit < $low_limit);
      for (my $i = $low_limit; $i <= $high_limit; $i++) { push @run_list, $i; }
    } elsif ($field =~ /^(\d+)$/) {
      push @run_list, $1;
    } else {
      die "wrong range format for field \"$field\"\n";
    }
  }

  my @tmp_list = sort @run_list;
  @run_list = ();
  my $old_field = -1;
  foreach my $field (@tmp_list) {
    next if ($field == $old_field);
    push @run_list, $field;
    $old_field = $field;
  }
  
  return @run_list;
}

sub get_valid_events {
  my $run_id = shift;
  my $query = "SELECT \"LastValidEvent\" FROM \"ValidRuns\" WHERE \"RunNumber\"=$run_id;";
  my $result = $db_connections{bx_runvalidation}->exec($query);
  die "Error in query: " . $db_connections{bx_runvalidation}->errorMessage . "\n" unless ($result->resultStatus == PGRES_TUPLES_OK);
  return -1 unless ($result->ntuples);
  my $number_of_events_pos = $result->fnumber ("\"LastValidEvent\"");
  my $number_of_events = $result->getvalue (0, $number_of_events_pos);
  return $number_of_events;
}

sub read_run_file {
  my $run_file = $_;
  warn "file $run_file is executable\n" if (-x $run_file); 
  if (-B $run_file) {
    warn "file $run_file is no ascii\n";
    return undef;
  }
  open (RF, $run_file) or die "Can not open file $run_file\n";
  my $ret = undef;
  while (<RF>) {
    chomp;
    s/\s//g;			# remove spaces
    next if (/^$/ or /^#/);	# ignore empty or comment lines
    $ret .= "," if (defined $ret);
    $ret .= $_;
  } 
  close RF;
  return $ret;
}

sub usage {
  print "bx_process.pl [global_options] run_list command <command flags>\n";
  print "Available commands:\n";
  print "\techidna:\t\t\t  run echidna with the specified command flags\n";
  print "\tpbs_echidna.sh:\t\t\t  run batch echidna with the specified command flags\n";
  print "Command flags substitution:\n";
  print "\t%r:\t\t\t  current run number\n";
  print "\t%t:\t\t\t  current run type\n";
  print "Run List format:\n";
  print "\tnum:\t\ta simple run number\n";
  print "\tnum1,num2:\ttwo runs (more can be added)\n";
  print "\tnum1-num2:\tthe range from run1 to run2\n";
  print "\tfile:\t\tascii file with one run per line\n";
  print "\tExample: \"1610, 1800-\" gets run 1610 and every run after 1800 (included)\n";
  print "Global options:\n";
  print "\t=do_nothing:\t\t\t  does only a simulation\n";
  print "\t=process_no_files:\t\t  ignore run with no files\n";
  print "\t=no_reprocess:\t\t\t  do not call pbs_echidna with =reprocess and -O sractch flags (meaningful only with =valid)\n";
  print "\t=reprocess:\t\t\t  force calling pbs_echidna with =reprocess flag (implies =output /scatch)\n";
  print "\t=rewrite_files:\t\t\t  call pbs_echidna with =rewrite_files flag (meaningful only with =valid)\n";
  print "\t=valid:\t\t\t\t  select only run with good validation with the right number of events (implies \"=run_type normal =reprocess\")\n";
  print "\t=queue queue:\t\t\t  set the queue to run on (valid only with pbs_echidna.sh)\n";
  print "\t=max_events num:\t\t  set a maximum number of event to process x run (requires =valid)\n";
  print "\t=run_type type:\t\t\t  select only run with the specified type (laser, pulser, normal, test, random, source, no_file)\n";
  print "\t=output folder:\t\t\t  specify a new default output folder\n";
  print "\t=sleep sec:\t\t\t  sleep for sec seconds after each jobs (0 mean no sleep)\n";
  print "\t=db_host host:\t\t\t  specifies the name/IP of the database machine\n";
  print "\t=db_name name:\t\t\t  specifies the db name to connect\n";
  print "\t=db_user user:\t\t\t  specifies the db user\n";
  print "\t=db_password pwd:\t\t  specifies the db password\n";
  exit (0);
}
