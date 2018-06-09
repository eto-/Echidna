#!/usr/bin/perl -w
#
# Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
#
# $Id: bx_repository.pl,v 1.56 2014/03/05 16:53:38 misiaszek Exp $
#
# Repository management utility: copy local run to a repository (via ftp) or
# copy run files from a remote repository to a local repository tree (which 
# can be used in echidna specifying the bx_event_reader.repository_url to
# file://local_dir)
#
#
use POSIX qw(strftime _exit);
use Pg;
use Time::ParseDate;
use Time::Local;
use Net::FTP;
use threads;
use strict;

$ENV{PATH} = "$ENV{PATH}:/usr/sbin/";

# Connection handlers
my $db_connection = 0;
my $ftp_connection = 0;
my $http_connection = 0;

# Defaults
my $do_nothing = 0;
my $rewrite_files = 0;
my $keep_days = 7;
my $parallel = 0;
my $running_jobs = 0;
my $delete_uploaded = 0;

my $db_host = "bxdb";
my $db_port = 5432;
my $db_name = "daq_config";
my $db_user = "borex_guest";
my $db_password = "";

my $remote_repository_host = "bxmaster-data.lngs.infn.it";
my $remote_repository_ftp_path = "/bxstorage";
my $remote_repository_ftp_user = "storage";
my $remote_repository_ftp_password = "Data";
my $remote_repository_http_path = "/bxstorage";
my $remote_repository_http_user = "borex";
my $remote_repository_http_password = "b0rexin0";
my $my_base = $ENV{"PWD"};
my $lockfile = "/var/lock/bx_repository";

END {
  unlink $lockfile or die "can not remove lock file: $!\n" if (-f $lockfile);
  if ($parallel) {
    foreach my $thr (threads->list ()) {
      $thr->join if ($thr->tid && $thr != threads->self)
    }
  }
}

$SIG{'TERM'} = $SIG{'INT'} = sub { exit (1); };

while (defined $ARGV[0] and $ARGV[0] =~ /^-/) {
  $_ = shift;
  last if /^--/;
  if (/^-do_nothing$/) {
    $do_nothing = 1;
  } elsif (/^-delete_uploaded$/) {
    $delete_uploaded = 1;
  } elsif (/^-rewrite_files$/) {
    $rewrite_files = 1;
  } elsif (/^-parallel$/) {
    $parallel = shift;
    test_opt ("-parallel", $parallel);
  } elsif (/^-keep_days$/) {
    $keep_days = shift;
    test_opt ("-keep_days", $keep_days);
  } elsif (/^-db_host$/) {
    $db_host = shift;
    test_opt ("-db_host", $db_host);
  } elsif (/^-db_name$/) {
    $db_name = shift;
    test_opt ("-db_name", $db_name);
  } elsif (/^-db_user$/) {
    $db_user = shift;
    test_opt ("-db_user", $db_user);
  } elsif (/^-db_password$/) {
    $db_password = shift;
    test_opt ("-db_password", $db_password);
  } elsif (/^-remote_repository_host$/) {
    $remote_repository_host = shift;
    test_opt ("-remote_repository_host", $remote_repository_host);
  } elsif (/^-remote_repository_ftp_path$/) {
    $remote_repository_ftp_path = shift;
    test_opt ("-remote_repository_ftp_path", $remote_repository_ftp_path);
  } elsif (/^-remote_repository_ftp_user$/) {
    $remote_repository_ftp_user = shift;
    test_opt ("-remote_repository_ftp_user", $remote_repository_ftp_user);
  } elsif (/^-remote_repository_ftp_password$/) {
    $remote_repository_ftp_password = shift;
    test_opt ("-remote_repository_ftp_password", $remote_repository_ftp_password);
  } elsif (/^-remote_repository_http_path$/) {
    $remote_repository_http_path = shift;
    test_opt ("-remote_repository_http_path", $remote_repository_http_path);
  } elsif (/^-remote_repository_http_user$/) {
    $remote_repository_http_user = shift;
    test_opt ("-remote_repository_http_user", $remote_repository_http_user);
  } elsif (/^-remote_repository_http_password$/) {
    $remote_repository_http_password = shift;
    test_opt ("-remote_repository_http_password", $remote_repository_http_password);
  } else {
    usage ($ARGV[0]);
  } 
}

my $mode = shift;
usage () if (!defined $mode || $mode =~ /^-/);

if ($mode =~ /^upload$/) {
  upload (@ARGV);
} elsif ($mode =~ /^download$/) {
  download (@ARGV);
} elsif ($mode =~ /^reorganize$/) {
  reorganize (@ARGV);
} elsif ($mode =~ /^check$/) {
  check (@ARGV);
} elsif ($mode =~ /^url$/) {
  url (@ARGV);
} elsif ($mode =~ /^trim$/) {
  trim (@ARGV);
} else { usage ($ARGV[0]); }



sub upload {
  my @files = @_;
  die "upload syntax: upload file/directory list\n" unless (scalar @files); 
  my @file_list;

  foreach my $file (@files) {
    if (-d $file) {
      push @file_list, glob ("$file/Run*.out*");
      push @file_list, glob ("$file/Run*.root");
      push @file_list, glob ("$file/Run*.log");
      push @file_list, glob ("$file/dst_*.root");
      warn "directory $file does not contain data files\n" unless (scalar @file_list);
    } elsif (-f $file) {
      push @file_list, $file;
    } else {
      warn "unknown or missing file type for $file\n";
    }
  }

  die "no valid file found\n" unless (scalar @file_list);

  my %seen = ();
  @file_list = grep { ! $seen{$_} ++ } @file_list;
  foreach my $file (keys %seen) { warn "found duplicated entry $file with $seen{$file} entries\n" if ($seen{$file} > 1); }

  create_lock ();
  connect_db ();
  connect_ftp ();

  foreach my $file (@file_list) {
    unless (-s $file) {
      warn "ignoring empty file $file\n";
      next;
    }
    if ($file =~ /(Run(\d{6})_\d{2}.out(.gz)?)$/) {
      my $run_id = int $2;
      my $file_name = $1;
      my $gz = $3;
      unless (defined $gz) {
        my $lsout = `lsof -t $file`;
	chomp $lsout;
	unless ($lsout =~ /^$/) {
          warn "ignoring open file $file\n";
          next;
	}
        warn "file $file is not compressed\n";
	unless ($do_nothing) {
          `gzip --fast $file`;
          push @file_list, "$file.gz";
	}
        next;
      }
      my $start_time = get_run_start_time ($run_id);
      next unless (defined $start_time);
      my ($year_index, $week_index) = get_run_start_times ($start_time);
      copy_file_over_ftp ("$remote_repository_ftp_path/rawdata/$year_index/$week_index/", $file, $file_name);
    } elsif ($file =~ /(Run(\d{6})([[:alpha:]_]+)?(_fadc_c(\d+))?.(root|log))/) {
      my $file_name = $1;
      my $run_id = int $2;
      my $cycle = $5;
      my $start_time = get_run_start_time ($run_id);
      next unless (defined $start_time);
      my ($year_index, $week_index) = get_run_start_times ($start_time);
      copy_file_over_ftp ( "$remote_repository_ftp_path/fadc/rootfiles/cycle_$cycle/$year_index/$week_index/", $file, $file_name); 
    } elsif ($file =~ /(Run(\d{6})([[:alpha:]_]+)?(_c(\d+))?.(root|log))/) {
      my $file_name = $1;
      my $run_id = int $2;
      my $file_tag = $3;
      my $cycle_tag = $4;
      my $cycle = $5;
      unless (defined $cycle_tag) {
	warn "ignoring file ($file_name) without cycle tag\n";
	next;
      }
      if ($file_tag and $file_tag !~ /^_(precalibrations|electronics_calibrations|online|laser_calibrations|laser_calibrations_validate|water|source|cngs)$/) {
        warn "ignoring unknown file type ($file_tag) for file $file_name\n";
        next;
      }
      die "cycle tag not defined for file $file" unless ($cycle);
      my $start_time = get_run_start_time ($run_id);
      next unless (defined $start_time);
      my ($year_index, $week_index) = get_run_start_times ($start_time);
      if ($file_tag) { copy_file_over_ftp ("$remote_repository_ftp_path/rootfiles/cycle_$cycle/$year_index/$week_index/ancillary", $file, $file_name); }
      else { copy_file_over_ftp ("$remote_repository_ftp_path/rootfiles/cycle_$cycle/$year_index/$week_index/", $file, $file_name); }
    } elsif ($file =~ /(Run(\d{6})_(m(\d+)c(\d+))(_[[:alpha:]_]+)?.(root|log))/) {
      my $file_name = $1;
      my $run_id = int $2;
      my $cycle_tag = $3;
      my $m4 = $4;
      my $cycle = $5;
      my $file_tag = $6;
      unless (defined $cycle_tag) {
	warn "ignoring file ($file_name) without cycle tag\n";
	next;
      }
      my $start_time = get_run_start_time ($run_id);
      next unless (defined $start_time);
      my $ancillary = "";
      $ancillary = "ancillary/" if ($file_tag);
      copy_file_over_ftp ("$remote_repository_ftp_path/rootfiles/m4_${m4}-c${cycle}/$ancillary", $file, $file_name);
    } elsif ($file =~ /(Calib_.*_Run_(\d+)_(m(\d+)c(\d+)).root)/) {
      my $file_name = $1;
      my $run_id = int $2;
      my $cycle_tag = $3;
      my $m4 = $4;
      my $cycle = $5;
      unless (defined $cycle_tag) {
	warn "ignoring file ($file_name) without cycle tag\n";
	next;
      }
      my $start_time = get_run_start_time ($run_id);
      next unless (defined $start_time);
      copy_file_over_ftp ("$remote_repository_ftp_path/rootfiles/m4_${m4}-c${cycle}/", $file, $file_name); 
    } elsif ($file =~ /(dst_(\d+)_[[:alpha:]]{3}_\d{2}_c(\d+).root)/) {
      my $file_name = $1;
      my $year = $2;
      my $cycle = $3;
      copy_file_over_ftp ("$remote_repository_ftp_path/dst/cycle_$cycle/$year/", $file, $file_name);
    } elsif ($file =~ /(dst_(\d+)_[[:alpha:]]{3}_\d{2}_(c(\d+)m(\d+)).(root|log))/) {
      my $file_name = $1;
      my $year = $2;
      my $cycle = $4;
      my $m4 = $5;
      copy_file_over_ftp ("$remote_repository_ftp_path/dst/cycle_$cycle-m$m4/$year/", $file, $file_name);
    } elsif ($file =~ /(dst-pp_(\d+)_[[:alpha:]]{3}_\d{2}_c(\d+).root)/) {
      my $file_name = $1;
      my $year = $2;
      my $cycle = $3;
      copy_file_over_ftp ("$remote_repository_ftp_path/dst-pp/cycle_$cycle/$year/", $file, $file_name);
    } else { 
      warn "ignoring $file\n";
    }
  }
}

sub download {
  my $local_repository = shift;
  my $run_list_pattern = join ("", @_);
  die "download syntax: download local_repository run_list\n" unless (defined $local_repository && defined $run_list_pattern); 

  unless (-d $local_repository) {
    mkdir $local_repository or die "Cannot create local repository $!\n";
  }
  unless (-d "$local_repository/rawdata") {
    mkdir "$local_repository/rawdata" or die "Cannot create $local_repository/rawdata $!\n";
  }

  chdir "$local_repository/rawdata" or die "Cannot chdir $local_repository/rawdata $!\n";
  connect_db ();

  my @run_list = extract_unique_run_list ($run_list_pattern);
  
  foreach my $run (@run_list) {
    my $slices = get_run_slices ($run);
    my $start_time = get_run_start_time ($run);
    next unless (defined $start_time);
    my ($year_index, $week_index) = get_run_start_times ($start_time);
    for (my $i = 0; $i < $slices; $i++) {
      my $file_name = sprintf ("Run%06d_%02d.out.gz", $run, $i + 1);
      my $file_url = "http://$remote_repository_host$remote_repository_http_path/rawdata/$year_index/$week_index/$file_name";
      unless (-d "$year_index") {
        mkdir "$year_index" or die "Cannot create $year_index $!\n"; 
      }
      unless (-d "$year_index/$week_index") {
        mkdir "$year_index/$week_index" or die "Cannot create $year_index/$week_index $!\n"; 
      }
      `wget -nv -O $year_index/$week_index/$file_name $file_url` if (!-f "$year_index/$week_index/$file_name" or (-f "$year_index/$week_index/$file_name" and $rewrite_files));
    }
  } 
}

sub reorganize {
  my $file = shift;
  my $local_repository = shift;
  die "reorganize syntax: reorganize file/directory local_reporitory\n" unless (defined $file && defined $local_repository); 

  die "$file not found\n" unless (-e $file);
  unless (-d $local_repository) {
    mkdir $local_repository or die "Cannot create local repository $!\n";
  }
  unless (-d "$local_repository/rawdata") {
    mkdir "$local_repository/rawdata" or die "Cannot create $local_repository/rawdata $!\n";
  }

  chdir "$local_repository/rawdata" or die "Cannot chdir $local_repository/rawdata $!\n";

  my @file_list;
  if (-d $file) {
    @file_list = glob ("$file/Run*.out*");
    die "directory $file does not contain raw data files\n" unless (scalar @file_list);
  } elsif (-f $file) {
    push @file_list, $file;
  } else {
    die "unknown file type for $file\n";
  }

  create_lock ();
  connect_db ();

  foreach $file (@file_list) {
    if ($file =~ /(Run(\d{6})_\d{2}.out(.gz)?)/) {
      my $run_id = int $2;
      my $file_name = $1;
      unless (defined $3) {
        warn "file $file_name is not compressed\n";
        `gzip $file` unless ($do_nothing);
        push @file_list, "$file.gz";
        next;
      }
      my $start_time = get_run_start_time ($run_id);
      next unless (defined $start_time);
      my ($year_index, $week_index) = get_run_start_times ($start_time);
      unless ($do_nothing) {
        mkdir "$year_index" or die "Cannot create $year_index $!\n" unless (-d "$year_index");
        mkdir "$year_index/$week_index" or die "Cannot create $year_index/$week_index $!\n" unless (-d "$year_index/$week_index"); 
        my $dest_file = "$year_index/$week_index/$file_name";
        `mv $file $dest_file` if (!-f $dest_file  or (-f $dest_file and $rewrite_files));
      } 
    }
  }
}

sub check {
  my $repository = shift;
  if ($repository =~ /^[0-9\-,]+$/) {
    unshift @_, $repository;
    $repository = "http://$remote_repository_host$remote_repository_http_path";
  }
  my $run_list_pattern = join ("", @_);
  die "check syntax: check [repository] run_list\n" unless (defined $repository && defined $run_list_pattern); 

  connect_db ();

  my @run_list = extract_unique_run_list ($run_list_pattern);
  
  foreach my $run (@run_list) {
    my $ret = check_raw ($run, $repository);
    warn "$run ok\n" if (!$ret and scalar (@run_list) < 10 and get_run_slices ($run) > 0);
  } 
}

sub url {
  my $repository = shift;
  if ($repository =~ /^[0-9\-,]+$/) {
    unshift @_, $repository;
    $repository = "http://$remote_repository_host$remote_repository_http_path";
  }
  my $cycle_tag = pop;
  if ($cycle_tag =~ /^[0-9\-,]+$/) {
    push @_, $cycle_tag;
    $cycle_tag = undef;
  } elsif ($cycle_tag =~ /^c(ycle_)?([0-9]+)$/) {
    $cycle_tag = $2;
  } else {
    die "unknow cycle tag $cycle_tag\n";
  }

  my $run_list_pattern = join ("", @_);
  die "check syntax: check [repository] run_list [url]\n" unless (defined $repository && defined $run_list_pattern); 

  connect_db ();

  my @run_list = extract_unique_run_list ($run_list_pattern);
  my $url_base;
  my $filename_tail;
  if ($cycle_tag) { 
    $url_base = "$repository/rootfiles/cycle_$cycle_tag"; 
    $filename_tail = "_c$cycle_tag.root";
  } else { 
    $url_base = "$repository/rawdata"; 
    $filename_tail = "_01.out.gz";
  }
  foreach my $run (@run_list) {
    my $start_time = get_run_start_time ($run);
    next unless (defined $start_time);
    my ($year_index, $week_index) = get_run_start_times ($start_time);
    my $file_name = sprintf ("Run%06d%s", $run, $filename_tail);
    print "$url_base/$year_index/$week_index/$file_name\n";
  } 
}

sub trim {
  my $dir = shift;
  die "trim syntax: trim directory [repository]\n" unless (defined $dir or -d $dir); 
  my $repository = shift;
  $repository = "http://$remote_repository_host$remote_repository_http_path" unless (defined $repository);

  my @file_list = glob ("$dir/Run*.out*");;
  push @file_list, glob ("$dir/Run*.root");
  push @file_list, glob ("$dir/Run*.log");
  die "directory $dir does not contain data files\n" unless (scalar @file_list);

  create_lock ();
  connect_db ();

  my $time_limit = time - $keep_days * 24 * 3600;
  my $free = 0;
  my $count = 0;
  foreach my $file (@file_list) {
    if ($file =~ /Run(\d{6})_\d{2}.out(.gz)?$/) {
      my $run_id = int $1;
      unless (-s $file or defined $2) {
        warn "triming empty file $file\n";
	unlink $file or warn "removing $file: $!\n" unless ($do_nothing);
	next;
      }
      my $start_time = get_run_start_time ($run_id);
      next unless (defined $start_time);
      next unless ($start_time < $time_limit);
      unless (check_raw ($run_id, $repository, 1)) {
	warn "triming $file\n";
	$free += -s $file; $count ++;
	unlink $file or warn "removing $file: $!\n" unless ($do_nothing);
      }
    } elsif ($file =~ /(Run(\d{6})([[:alpha:]_]+)?(_c(\d+))?.(root|log))/) {
      my $file_name = $1;
      my $run_id = int $2;
      my $file_tag = $3;
      my $cycle_tag = $4;
      my $cycle = $5;
      unless (defined $cycle_tag) {
	warn "ignoring file ($file_name) without cycle tag\n";
	next;
      }
      if ($file_tag and $file_tag !~ /^_(precalibrations|electronics_calibrations|online|laser_calibrations|laser_calibrations_validate|water|source)$/) {
        warn "ignoring unknown file type ($file_tag) for file $file_name\n";
        next;
      }
      my $start_time = get_run_start_time ($run_id);
      next unless (defined $start_time);
      next unless ($start_time < $time_limit);
      $file_name = "ancillary/$file_name" if ($file);
      unless (check_root ($run_id, $repository, $file_name, $cycle)) {
	$free += -s $file; $count ++;
	unlink $file or warn "removing $file: $!\n" unless ($do_nothing);
      }
    } else {
      warn "ignoring $file\n";
    }
  }
  $free = int ($free / 1024);
  print "removed $count files for $free kbytes\n" if ($count);
}

sub check_raw {
  my $run = shift;
  my $repository = shift;
  my $second_test = shift;

  my $slices = get_run_slices ($run);
  my $start_time = get_run_start_time ($run);
  return -1 unless (defined $start_time);
  my ($year_index, $week_index) = get_run_start_times ($start_time);
  my $found_slices = 0;
  for (my $i = 0; $i < 100; $i++) {
    my $file_name = sprintf ("Run%06d_%02d.out.gz", $run, $i + 1);
    my $file_url = "$repository/rawdata/$year_index/$week_index/$file_name";
    if ($repository =~ /^http/) {
      my $ret = `wget -nv --spider $file_url 2>&1`;
      last unless ($ret =~ /200.OK/);
    } else {
      last unless (-f $file_url);
      if (system ("gzip -t $file_url 2>/dev/null")) {
	warn "run $run does not match gzip, transfer incomplete\n";
	return 2;
      }
      last if ($?);
    }
    $found_slices++;
  }
  if ($found_slices != $slices) {
    unless ($repository =~ /^http/ or defined ($second_test) or getpwuid ($<) !~ /^daqman$/) {
      #warn "run $run found $found_slices slices on storage while $slices in db, trying to fix\n";
      my @fix_db = ("$my_base/fix_db", "$run");
      chdir "$repository/rawdata/$year_index/$week_index/";
      system (@fix_db);
      return check_raw ($run, $repository, 1);
    } else {
      warn "run $run found $found_slices slices on storage while $slices in db, please fix db\n";
      return 2;
    }
  } 
  #elsif ($found_slices < $slices) {
  #  warn "run $run found $found_slices slices on storage while $slices in db, data missing\n";
  #  return 1;
  #}
  return 0;
}

sub check_root {
  my $run = shift;
  my $repository = shift;
  my $file = shift;
  my $cycle = shift;

  my $start_time = get_run_start_time ($run);
  return -1 unless (defined $start_time);
  my ($year_index, $week_index) = get_run_start_times ($start_time);
  my $file_url = "$repository/rootfiles/cycle_$cycle/$year_index/$week_index/$file";
  if ($repository =~ /^http/) {
    my $ret = `wget -nv --spider $file_url 2>&1`;
    return 0 if ($ret =~ /200.OK/);
  } else {
    return 0 if (-f $file_url);
  }
  return 1;
}

sub test_opt {
  my $opt_name = shift;
  my $opt_value = shift;
  die "$opt_name requires an argument ($opt_value)\n" if (!defined $opt_value || $opt_value =~ /^-/);
}

sub connect_db {
  if ($db_connection) {
    warn "Already connected to db\n";
    return;
  }
  $db_connection = Pg::connectdb("host=$db_host port=$db_port dbname=$db_name user=$db_user");
  die "Unable to connect to the database\n" if ($db_connection->status != PGRES_CONNECTION_OK);
}

sub get_last_run {
  my $query = "SELECT \"RunNumber\" FROM \"Run\" ORDER by \"RunNumber\" DESC;";
  my $result = $db_connection->exec($query);
  die "Error in query: " . $db_connection->errorMessage . "\n" unless ($result->resultStatus == PGRES_TUPLES_OK);
  return undef unless ($result->ntuples);
  my $run_number_pos = $result->fnumber ("\"RunNumber\"");
  my $run_number = $result->getvalue (0, $run_number_pos);

  return $run_number;
}

sub get_run_slices {
  my $run_id = shift;
  my $query = "SELECT \"NumberOfFiles\" FROM \"Run\" WHERE \"RunNumber\"=$run_id;";
  my $result = $db_connection->exec($query);
  die "Error in query: " . $db_connection->errorMessage . "\n" unless ($result->resultStatus == PGRES_TUPLES_OK);
  return undef unless ($result->ntuples);
  my $number_of_files_pos = $result->fnumber ("\"NumberOfFiles\"");
  my $number_of_files = $result->getvalue (0, $number_of_files_pos);
  return $number_of_files;
}

sub get_run_start_time {
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

sub get_run_start_times {
  my $start_time = shift;
  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday) = localtime $start_time;
  my $week_time = $start_time - (($wday * 24 + $hour) * 60 + $min) * 60 + 4800;
  my @tmp = localtime $week_time;
  my $week_index = strftime ("%b_%d", @tmp); 
  my $year_index = $tmp[5] + 1900;
  return ($year_index, $week_index);
}

sub connect_ftp {
  $ftp_connection = Net::FTP->new ("$remote_repository_host", Debug => 0, Timeout => 250) or die "Cannot connect to $remote_repository_host: $@\n";
  $ftp_connection->login ($remote_repository_ftp_user, $remote_repository_ftp_password) or die "Cannot login: " . $ftp_connection->message . "\n";
  warn "connection done\n" unless ($parallel);
}

sub disconnect_ftp {
  $ftp_connection->quit ();
}

sub copy_file_over_ftp {
  my $dir = shift;
  my $file = shift;
  my $remote_file = shift;
  my $local_size = -s $file;

  connect_ftp () if ($parallel);

  $ftp_connection->cwd ("/") or die "Cannot chdir / " . $ftp_connection->message . "\n";
  unless ($do_nothing) {
    unless ($ftp_connection->cwd ($dir)) {
      $ftp_connection->mkdir ($dir, 1) or $ftp_connection->mkdir ($dir, 1) or die "Cannot mkdir $dir: " . $ftp_connection->message . "\n";
      $ftp_connection->cwd ($dir) or die "Cannot chdir $dir " . $ftp_connection->message . "\n";
    }
    my $size = $ftp_connection->size ($remote_file);
    if (!$size or ($size and $rewrite_files)) {
      $ftp_connection->binary ();
      $ftp_connection->put ($file) or die "Cannot put $file " . $ftp_connection->message . "\n";
      $size = $ftp_connection->size ($remote_file);
      unless ($size == $local_size) {
	warn "transfer size mismatch for file $file ($local_size != $size)\n";
      } else {
	print "copied $file in ftp://$remote_repository_host$dir\n";
	unlink $file if ($delete_uploaded);
      }
    } else {
      unless ($size == $local_size) {
	warn "file $file has different size than ftp://$remote_repository_host$dir/$remote_file ($local_size != $size)\n";
      } else {
	print "ignored $file already present at ftp://$remote_repository_host$dir$remote_file\n";
      }
    }
  }

  disconnect_ftp () if ($parallel);
}

sub parallel_copy_file_over_ftp {
  unless ($parallel) {
    copy_file_over_ftp (@_);
    return;
  }


  while ($running_jobs >= $parallel) {
    foreach my $thr (threads->list ()) {
      next unless ($thr->is_joinable());
      $thr->join ();
      $running_jobs --;
    }
    sleep 1;
  }

  threads->create (\&copy_file_over_ftp, @_);
  $running_jobs ++;
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

sub usage {
  print "bx_repository.pl [global_options] command <command_flags>\n";
  print "Available commands:\n";
  print "\tupload file:\t\t\tupload file (or directory contents) to the main repository\n";
  print "\tdownload local_rep run_list:\tdownload run specified in run list to local repository folder\n";
  print "\tcheck [repository] run_list:\tcheck the repository for missing files\n";
  print "\turl [repository] run_list [cycle]:\tdump the url of a file, if cycle is missing raw data file will be dumped\n";
  print "\treorganize file local_rep:\tput file (or dir contents) in the local repository\n";
  print "\ttrim local_dir [repository]:\tremove old files in local directory if they are present in the repository\n";
  print "Run List format:\n";
  print "\tnum:\t\ta simple run number\n";
  print "\tnum1,num2:\ttwo runs (more can be added)\n";
  print "\tnum1-num2:\tthe range from run1 to run2\n";
  print "\tExample: \"1610, 1800-\" gets run 1610 and every run after 1800 (included)\n";
  print "Global options:\n";
  print "\t-do_nothing:\t\t\t\t  does only a simulation\n";
  print "\t-rewrite_files:\t\t\t\t  permits file overrite (normally disabled)\n";
  print "\t-delete_uploaded:\t\t\t  deleted sucessully uploaded files (normally disabled)\n";
  print "\t-keep_days:\t\t\t\t  how old a file should be to deleted by trim (default 7)\n";
  print "\t-db_host host:\t\t\t\t  specifies the name/IP of the database machine\n";
  print "\t-db_name name:\t\t\t\t  specifies the db name to connect\n";
  print "\t-db_user user:\t\t\t\t  specifies the db user\n";
  print "\t-db_password pwd:\t\t\t  specifies the db password\n";
  print "\t-remote_repository_host host:\t\t  specifies the name/IP of the repository machine\n";
  print "\t-remote_repository_ftp_path path:\t  specifies the path for the ftp repository\n";
  print "\t-remote_repository_ftp_user user:\t  specifies the user for the ftp repository\n";
  print "\t-remote_repository_ftp_password pwd:\t  specifies the password for the ftp repository\n";
  print "\t-remote_repository_http_path path:\t  specifies the path for the http repository\n";
  print "\t-remote_repository_http_user user:\t  specifies the user for the http repository\n";
  print "\t-remote_repository_http_password pwd:\t  specifies the password for the http repository\n";
  exit (0);
}

sub create_lock {
  die "lock file is already present, maybe $0 is already running (else remove $lockfile)\n" if (-f $lockfile);
  open (LOCK, ">$lockfile") or die "can not create lock file: $!\n";
  close (LOCK);
}

