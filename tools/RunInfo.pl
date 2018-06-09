#!/usr/bin/perl -w
# Mapping Utility
# Access the DB and provides infos about a Borexino RUN
# M.Pallavicini created - 20-12-2002
# 

use Pg;
use POSIX;
use Time::Local;


#
# check if we are running on the local cluster in Hall C or not
# and set data base server name accordingly
# 
$_=`hostname`;
if (/^bxweb/ || /^bxdb/ || /^bxoff1/ || /^bxbuild/ || /^bxterm1/  ) {
  $dbhost = "bxdb";
} else {
  $dbhost = "bxdb.lngs.infn.it";
}

#
# set default value to variables
# 
my $html=0;
my $search=0;                    # search string in comments
my $d_start="2000-01-01 00:00";  # start date
my $d_end="2050-12-31 00:00";    # end date
my $UseColors=1;                 # color or white printout
my $Search=0;                    
my $ProfileKey=0;                # select profile
my $Type=0;                      # trigger type
my $Sort="RunNumber";            # sorting field
my $MinEvents = 0;               # minimum number of events
my $Short = 0;                   # short printout
my $VeryShort = 0;               # run numbers only
my $Valid = 0;                   # show valid runs only
my $DutyCycle = 0;               # compute duty cycle in run interval

#
# parse parameters
# 

$r_start = 1;
$r_end = 99999;

foreach $token (@ARGV) {
   
   # select number
   $i = index($token,"run=");
   if ($i>=0) {
     my $r = substr($token,4);
     my $j = index($r,"-");
     if ($j>=0) {
        $r_start = substr($r,0,$j);
	$r_end = substr($r,$j+1,length($r) - $j);
	for($k=$r_start; $k<$r_end; $k++) {
	  @runs = (@runs,$k);
	}
     } else {
       @runs = split(',',$r);
     }
     next;
   }
  
   # select profile name
   $i = index($token,"profile=");
   if ($i>=0) {
     $_ = substr($token,8);
     if ( /[^0-9]/ ) {
       $Digit=0;
     } else {
       print "$_\n";
       $Digit=1;
     }
     $ProfileKey=$_;
     next;
   }

   # select minimum event size
   $i = index($token,"events=");
   if ($i>=0) {
	   $_ = substr($token,7);
	   if ( /[^0-9]/ ) {
	   	die "Wrong events number\n";
	   }
	   $MinEvents=$_;
	   next;
   }
   
   # select start date
   $i = index($token,"after=");
   if ($i>=0) {
     $d_start = substr($token,6);
     $d_start = $d_start." 00:00";
     next;
   }
   
   # select end date
   $i = index($token,"before=");
   if ($i>=0) {
     $d_end = substr($token,7);
     $d_end = $d_end." 23:59";
     next;
   }
   
   # search string in comment
   $i = index($token,"search=");
   if ($i>=0) {
     $Search = substr($token,7);
     next;
   }
   
   # select run type
   $i = index($token,"type=");
   if ($i>=0) {
     $Type = substr($token,5);
     next;
   }
   
   # select sort field
   $i = index($token,"sort=");
   if ($i>=0) {
     $Sort = substr($token,5);
     next;
   }


   #help
   $i = index($token,"help");
   if ($i>=0) {
     &Usage();
     exit;
   }
   
   #do not use colors on screen
   $i = index($token,"nocol");
   if ($i>=0) {
     $UseColors=0;
     next;
   }
   
   #produce a compact output (1 line per run)
   $i = index($token,"short");
   if ($i>=0) {
     $Short=1;
     next;
   }
   
   #produce the list of run only
   $i = index($token,"list");
   if ($i>=0) {
     $VeryShort=1;
     next;
   }
   
   #show valid runs only
   $i = index($token,"valid");
   if ($i>=0) {
     $Valid=1;
     next;
   }
   
   #compute duty cycle
   $i = index($token,"dutycycle");
   if ($i>=0) {
     $DutyCycle=1;
     next;
   }
   
   # option to write output as HTML output
   $i = index($token,"html");
   if ($i>=0) {
     $html=1;
     print "Option key html not implemented yet\n";
     exit;
     next;
   }

   printf("Unknown option $token. Ignored\n");
}

@data=();
%ValidRuns=();

# get run max e min
my $rmin = 9999999;
my $rmax = 0;
for $r (@runs) {
 	if ( $r > $rmax ) {
		$rmax = $r;
	}
	if ($r < $rmin ) {
		$rmin = $r;
	}
}
$r_start = $rmin;
$r_end = $rmax;

for($i=$r_start; $i<=$r_end; $i++) {
	$ValidRuns{$i}=0;
}

if ($Valid==1) {
	&GetValidRuns($dbhost,$r_start,$r_end);
}

my $query = &PrepareQuery();
my $recs = &GetData($query,$dbhost);
if (!$html && (!($Short || $VeryShort)) ) {
  &TextDump($recs);
} 
if (!$html && $Short ) {
  &ShortTextDump($recs);
}
if (!$html && $VeryShort ) {
  &VeryShortTextDump($recs);
}
if ($html) {
  &HtmlDump();
}
if ($DutyCycle) {
	$duty = &ComputeDutyCycle($recs);
	printf("\nDutycycle for runs = %d - %d: %4.3f \n\n",$r_start,$r_end,$duty);
}
exit(0);

###
###

sub ComputeDutyCycle() {
	my $recs = $_[0];
	$tot = 0;
	$duty = 1;
	for($i=0; $i<$recs; $i++) {
		($a1,$b1) = split(' ',$data[$i*7+3]);
		($y1,$m1,$d1) = split('-',$a1);
		($h1,$mm1,$s1) = split(':',$b1);
		$y1 -= 1900;
		$m1 -= 1;

		$t1 = timelocal($s1,$mm1,$h1,$d1,$m1,$y1);
		if ($i==0) {$t_begin = $t1;}

		($a2,$b2) = split(' ',$data[$i*7+4]);
		($y2,$m2,$d2) = split('-',$a2);
		($h2,$mm2,$s2) = split(':',$b2);
		$y2 -= 1900;
		$m2 -= 1;
		$t2 = timelocal($s2,$mm2,$h2,$d2,$m2,$y2);
		if ($i==($recs-1)) {$t_end = $t2;}
		
		$deltat = $t2 - $t1;
		if ( $Valid==0 || ($ValidRuns{ $data[$i*7+0] }>100 )) {
			$tot += $deltat;
		}
	}
	if ( $t_end != $t_begin && $recs != 0 ) {
		$duty = $tot / ($t_end-$t_begin);
	} 
	return $duty;
}

## prepare query with all options required
sub PrepareQuery() {
  my $query = "SELECT \"RunNumber\",\"ProfileID\",\"Type\",
                      \"StartTime\",\"StopTime\",\"Events\",\"Comments\"";
  $query = $query . " FROM \"Run\" WHERE ( \"StartTime\" \> \'$d_start\' AND \"StopTime\" \< \'$d_end\' AND
	              \"Events\" >= $MinEvents  AND \"RunNumber\"<=$r_end AND \"RunNumber\">=$r_start ) ";
  if ( $ProfileKey ) {
    my $id = $ProfileKey;
    if ( !$Digit ) {
      $id = &GetProfileByName($ProfileKey);
    }
    if ( $id >= 0 ) {
      $query = $query . " AND (\"ProfileID\"=$id)";
    }
  }
  if ($Search) {
    my $ss = "\%$Search\%";
    $query = $query . " AND (\"Comments\" like \'$ss\' ) ";
  }
  if ($Type) {
    $query = $query . " AND (\"Type\" like \'$Type\' ) ";
   }
  $query = $query . " ORDER BY \"$Sort\""; 

  return $query;
}

## get profile number knowing name
sub GetProfileByName() {
  my $name = $_[0];
  my $conn = Pg::connectdb("host=$dbhost port=5432 user=borex_guest dbname=daq_config");
  $conn->status && error("Error connecting to Data Base daq_config");
  my $qq = "SELECT \"ProfileID\" FROM \"Profiles\" WHERE (\"Name\" like \'$name\')";
  my $result = $conn -> exec( $qq );
  ##$result->resultStatus == PGRES_TUPLES_OK || error("Error reading from table \"Profiles\"");
  if ($result->ntuples == 0) {
    print "Profile $name does not exist. Ignored \n";
    return -1;
  }
  return $result->getvalue(0,0);
}

## connect to DB and get data
sub GetData() {
  my $q = $_[0];
  my $h = $_[1];
  my $conn = Pg::connectdb("host=$h port=5432 user=borex_guest dbname=daq_config");
  $conn->status && error("Error connecting to Data Base daq_config");
  my $result = $conn -> exec( $q );
  $result->resultStatus == PGRES_TUPLES_OK || error("Error reading from table \"Run\"");
  if ($result->ntuples == 0) {
    print "No data available. \n";
    print "Input parameters not consistent. \n";
    exit;
  }
  my $records = $result->ntuples;
  if ( $result->nfields != 7 ) {
    print "Error reading DB.\n";
    print "Wrong number of fields\n";
    exit;
  }
  # save into data array
  @data=();
  my $i=0;
  my $f=0;
  for ($i = 0; $i < $records; $i++) {
    for ($f = 0; $f < 7; $f++) {
      $data[$i*7+$f] = $result->getvalue($i,$f);
    }
  }
  return $records;
}


## dump on screen (very compact version, run numbers only)
sub VeryShortTextDump() {
  my $recs = $_[0];
  for($i=0; $i<$recs; $i++) {
    printf("%04d\n",$data[$i*7+0]);
  }
}


## dump on screen (compact version)
sub ShortTextDump() {
	my $recs = $_[0];
	for($i=0; $i<$recs; $i++) {
		if ($Valid==0 || $ValidRuns{$data[$i*7+0]}>0) {
			my $dd=0;
			printf("\033[01;31m") if ($UseColors); #magenta
				printf("Run=%04d",$data[$i*7+0]);
			printf("\033[01;33m") if ($UseColors); #yellow
				printf(" Type:");
			printf("\033[01;32m") if ($UseColors); #green
				printf("%s",$data[$i*7+2]);
			printf("\033[01;33m") if ($UseColors); #yellow
				printf(" Events:");
			printf("\033[01;32m") if ($UseColors); #green
				printf("%s",$data[$i*7+5]);
			$ll=length($data[$i*7+5]);
			for($j=$ll;$j<8;$j++) {
				print " ";
			}
			printf("\033[01;33m") if ($UseColors); #yellow
				printf(" Last Valid:");
			printf("\033[01;32m") if ($UseColors); #green
				printf("%s",$ValidRuns{$data[$i*7+0]});
			$ll=length($ValidRuns{$data[$i*7+0]});
			for($j=$ll;$j<8;$j++) {
				print " ";
			}
			
			printf("\033[01;33m") if ($UseColors); #yellow
				printf(" Start:");    
			printf("\033[01;32m") if ($UseColors); #green
				printf("%s",$data[$i*7+3]); 
			printf("\033[01;33m") if ($UseColors); #yellow
				printf(" End:"); 
			printf("\033[01;32m") if ($UseColors); #green
				printf("%s\n",$data[$i*7+4]);
		}
	}
	printf("\033[0m") if ($UseColors); # back to white
}


## dump on screen (verbose version)
sub TextDump() {
	my $recs = $_[0];
	for($i=0; $i<$recs; $i++) {
		if ($Valid==0 || $ValidRuns{$data[$i*7+0]}>0) {
			my $dd=0;
			printf("\033[01;31m") if ($UseColors); #magenta
				printf("Run = %04d\n",$data[$i*7+0]);
			printf("\033[01;33m") if ($UseColors); #yellow
				printf(" Profile = "); 
			printf("\033[01;32m") if ($UseColors); #green
				printf("%02d",$data[$i*7+1]);
			printf("\033[01;33m") if ($UseColors); #yellow
				printf("   Type: ");
			printf("\033[01;32m") if ($UseColors); #green
				printf("%s",$data[$i*7+2]);
			printf("\033[01;33m") if ($UseColors); #yellow
				printf("   Events: ");
			printf("\033[01;32m") if ($UseColors); #green
				printf("%s",$data[$i*7+5]);
			printf("\033[01;33m") if ($UseColors); #yellow
				printf("   Last Valid: ");
			printf("\033[01;32m") if ($UseColors); #green
				printf("%s\n",$ValidRuns{$data[$i*7+0]});
			printf("\033[01;33m") if ($UseColors); #yellow
				printf(" Start Date: ");    
			printf("\033[01;32m") if ($UseColors); #green
				printf("%s",$data[$i*7+3]); 
			printf("\033[01;33m") if ($UseColors); #yellow
				printf("   End Date: "); 
			printf("\033[01;32m") if ($UseColors); #green
				printf("%s\n",$data[$i*7+4]);
			printf("\033[01;33m") if ($UseColors); #yellow
				$ll = length($data[$i*7+6]);
			$off=0;
			printf("\033[01;36m") if ($UseColors); #cyan
				while($ll>0) {
					if ($ll>59) {
						$l1 = 59;
					} else {
						$l1 = $ll;
					}
					$ss = substr($data[$i*7+6],$off,$l1);
					printf("\t%s\n",$ss);
					$ll -= $l1;
					$off += $l1;
				}
		}
	}
	printf("\033[0m") if ($UseColors); # back to white
}

## dump in html format
sub HtmlDump() {
}

## Print out Error Message and exit
sub error{
  print ("@_\n");
  exit 1;
}   

## get list of valid runs 
sub GetValidRuns() {
	my $h = $_[0];
	my $r1 = $_[1];
	my $r2 = $_[2];
	my $conn = Pg::connectdb("host=$h port=5432 user=borex_guest dbname=bx_runvalidation");
	$conn->status && error("Error connecting to Data Base daq_config");
	my $q = "SELECT \"RunNumber\",\"LastValidEvent\" FROM \"ValidRuns\" WHERE \"RunNumber\">=$r1 AND \"RunNumber\"<=$r2;";
	my $result = $conn -> exec( $q );
	$result->resultStatus == PGRES_TUPLES_OK || error("Error reading from table \"Run\"");
	my $records = $result->ntuples;
	if ($records==0) {
		return;
	}
	my $i=0;
	for ($i = 0; $i < $records; $i++) {
		my $r = $result->getvalue($i,0);
		$ValidRuns{$r} = $result->getvalue($i,1);
	}
	return;
}

## help infos
sub Usage() {
  print "\nBorexino RunInfo Utility.\n";
  print "Dump useful informations about data run\n";
  print "\n";
  print "Usage: RunInfo [field1=val] [field2=val] [key] [sort=sort_f]\n";
  print " field1,2..n can be any one among the following names:\n";
  print "  run=r1,r2,r3              runs r1,r2,r3\n";
  print "  run=r1-r2                 all runs from r1 to r2\n";
  print "  events=min                at least min events are required\n";
  print "  after=yyyy-mm-dd          runs after date\n";
  print "  before=yyyy-mm-dd         runs before date\n";
  print "  profile=prof_num          runs of profile prof_num only (digit)\n";
  print "  profile=prof_name         runs of profile prof_name only (name)\n";
  print "  search=string             search for string in comment\n";
  print "  type=type_name            runs with Type = type_name\n";
  print "  sort=sort_field           sort output by sort_field\n";
  print " keys can be:\n";
  print "  help                      dump this screen\n";
  print "  valid                     show valid runs only\n";
  print "  nocol                     do not use colors in print out\n";
  print "  html                      output in html format\n";
  print "  short                     one line of output per run \n";
  print "  list                      run list only\n";
  print "  dutycycle                 compute duty cycle for run-list \n";
  print " sort=sort_f                sort_f name of sorting field. Default RunNumber\n";
  print "\n";
}

