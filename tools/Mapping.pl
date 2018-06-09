#!/usr/bin/perl -w
# Mapping Utility
# Access the DB and provides mapping infos
#  

use Pg;
use POSIX;

#
# set default value to variables
# 
my $crate=-1;
my $channel=-1;
my $feb=-1;
my $fec=-1;
my $lbnb=-1;
my $lbnc=-1;
my $hvb=-1;
my $hvc=-1;
my $hole=0;
my $cluster=-1;
my $profileID=1;
my $sort="";
my $polar=0;
my $html=0;
my $all=0;

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
# hash table to associate internal names to data base names
#
%names = ( "crate" => "\"RackID\"" ,
           "channel" => " \"RackID\", \"LabenBoard\", \"LabenChannel\"",
           "feb" => "\"FEBoard\"",
           "fec" => "\"FEChannel\"",
           "lbnb" => "\"LabenBoard\"",
           "lbnc" => "\"LabenChannel\"",
           "hvb" => "\"HVBoard\"",
           "hvc" => "\"HVChannel\"",
           "hole" => "\"HoleLabel\"",
           "cluster" => "\"FiberBundle\""
);

#
# parse parameters
# 
foreach $token (@ARGV) {
   $i = index($token,"crate=");
   if ($i>=0) {
     $crate = substr($token,6);
     if ($crate<=0 || $crate>14) {
       print "Invalid crate number. Must be 1..14\n";
       exit;
     }
     next;
   }   
   $i = index($token,"rack=");
   if ($i>=0) {
     $crate = substr($token,5);
     if ($crate<=0 || $crate>14) {
       print "Invalid crate number. Must be 1..14\n";
       exit;
     }
     next;
   }

   $i = index($token,"channel=");
   if ($i>=0) {
     $channel = substr($token,8);
     if ($channel<=0 || $channel>2240) {
       print "Invalid channel number. Must be 1..2240\n";
       exit;
     }
     next;
   }
   $i = index($token,"feb=");
   if ($i>=0) {
     $feb = substr($token,4);
     if ($feb<=0 || $feb>14) {
       print "Invalid front end board number. Must be 1..14\n";
       exit;
     }
     next;
   }
   $i = index($token,"fec=");
   if ($i>=0) {
     $fec = substr($token,4);
     if ($fec<=0 || $fec>12) {
       print "Invalid front end channel number. Must be 1..12\n";
      exit;
     }
      next;
   }
   $i = index($token,"lbnb=");
   if ($i>=0) {
     $lbnb = substr($token,5);
     if ($lbnb<=0 || $lbnb>20) {
       print "Invalid digital board number. Must be 1..20\n";
       exit;
     }
     next;
   }
   $i = index($token,"lbnc=");
   if ($i>=0) {
     $lbnc = substr($token,5);
     if ($lbnc<=0 || $lbnc>8) {
       print "Invalid digital channel number. Must be 1..8\n";
       exit;
     }
     next;
   }
   $i = index($token,"hvb=");
   if ($i>=0) {
     $hvb = substr($token,4);
     if ($hvb<=0 || $hvb>7) {
       print "Invalid hv board number. Must be 1..7\n";
       exit;
     }
     next;
   }
   $i = index($token,"hvc=");
   if ($i>=0) {
     $hvc = substr($token,4);
     if ($hvc<=0 || $hvc>24) {
       print "Invalid hv channel number. Must be 1..24\n";
       exit;
     }
     next;
   }
   $i = index($token,"hole=");
   if ($i>=0) {
     $hole = substr($token,5);
     next;
   }
   $i = index($token,"cluster=");
   if ($i>=0) {
     $cluster = substr($token,8);
     next;
   }   
   $i = index($token,"fiber=");
   if ($i>=0) {
     $cluster = substr($token,6);
     next;
   }

   $i = index($token,"profileID="); 
##   $_ = $token;
##   if (/^profileID=(\d+)$/ && $1 >= 0)
   if ($i>=0) {
     $profileID = substr($token,10);
     next;
   }

   #help
   $i = index($token,"help");
   if ($i>=0) {
     &Usage();
     exit;
   }
   
   # sort field
   $i = index($token,"sort=");
   if ($i>=0) {
     $sort=substr($token,5);
     next;
   }
   
   # option to write pmt polar coordinates instead of cartesian
   $i = index($token,"polar");
   if ($i>=0) {
     $polar=1;
     next;
   }
   
   # option to count all channels, including not mapped one in this profile
   $i = index($token,"all");
   if ($i>=0) {
     $all=1;
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

if ( $channel >= 0 ) {
  my $cc = POSIX::floor(($channel-1)/160) + 1;
  if ($crate != -1 && $crate!=$cc) {
    print "Crate number $crate not consistent with logical channel $channel. Changed to $cc\n";
  }
  $crate = $cc;
  my $bb = POSIX::floor(($channel-1-($crate-1)*160)/8) + 1;
  if ($lbnb != -1 && $lbnb!=$bb) {
    print "Laben Board number $lbnb not consistent with logical channel $channel. Changed to $bb\n";
  }
  $lbnb = $bb;
  my $ch = ($channel-($crate-1)*160-($lbnb-1)*8);
  if ($lbnc != -1 && $lbnc!=$ch) {
    print "Laben channel number $lbnc not consistent with logical channel $channel. Changed to $ch\n";
  }
  $lbnc = $ch;
}

@data_list=();
my $recs=0;
if ($cluster<0) {
  $recs = &GetMap1($profileID,$hole,$crate,$feb,$fec,$lbnb,$lbnc,$hvb,$hvc,$sort);
  if ($recs>0) {
    &GetMap2($profileID,$cluster,$recs);
  }
} else {
  print "Fiber cluster cannot be a search key so far. Sorry!\n";
  exit;
}

# output 
if ($html==0) {
  &DumpTextData($recs,$polar,$all);
} else {
  &HtmlOutput($recs,$all);
}

exit;

###
### End of main program ###
###

#
# Print out Error Message and exit
#
sub error{
  print ("@_\n");
  exit 1;
}   
#
# Dump Help Message
#
sub Usage() {
  print "\nUtility Mapping.pl\n";
  print "Dumping Borexino Data Base Mapping Infos\n";
  print "\nUsage :  Mapping [field1=val] [field2=val] [key] [sort=sort_f]\n";
  print "\n field1,2..n can be any one among the following names: \n";
  print "    crate       (crate number from 1 to 14)\n";
  print "    rack        (same as crate)\n";
  print "    channel     (logical channel number)\n";
  print "    feb         (front end board number 1..14\n";
  print "    fec         (front end channel number 1..12\n";
  print "    feb         (front end board number 1..14\n";
  print "    lbnb        (laben board number 1..20\n";
  print "    lbnc        (laben board number 1..8\n";
  print "    hvb         (high voltage board number 1..7\n";
  print "    hvc         (high voltage number 1..24\n";
  print "    hole        (SSS hole ID\n";
  print "    cluster     (Fiber cluster number\n";
  print "    fiber       (same as cluster\n";
  print "    profileID   (Profile number\n";
  print "\n key can be any of the follwing:\n";
  print "    help        (to get this message)\n";
  print "    html        (to get output in html form)\n";
  print "    polar       (to get PMT positions as Theta,Phi angles, not X,Y,Z)\n";
  print "    all         (print all channels including not cabled one)\n";
  print "\n sort=sort_f specifies the output sort field\n";
}

#
# Dump data on standard output
#
sub DumpTextData() {
  my $r=$_[0];
  my $pol=$_[1];
  my $i=0;
  if ($pol==0) {
    print "Chan  Hole Fiber    (  X,  Y,  Z)     Crate   FE   Laben    HV\n";
  } else {
    print "Chan  Hole Fiber   (Theta,Phi)  Crate   FE   Laben    HV\n";
  }
  for($i=0; $i<$r; $i++) {
    if ( ( $data_list[$i*16+8] > 0 ) || ( $all==1) ) {
      printf "%04d ",$data_list[$i*16+8];
      printf "%+05d ",$data_list[$i*16+0];
      printf "  %02d ", $data_list[$i*16+15];
      if ($pol==0) {
        printf " (%+4.2f,%+4.2f,%+4.2f) ",$data_list[$i*16+9],$data_list[$i*16+10],
                                          $data_list[$i*16+11];
      } else {
        printf " (%+06.1f,%+06.1f)",$data_list[$i*16+12],$data_list[$i*16+13];
      }
      printf "  %2d",$data_list[$i*16+1];
      printf "   %02d.%02d",$data_list[$i*16+2],$data_list[$i*16+3];
      printf "  %02d.%02d",$data_list[$i*16+4],$data_list[$i*16+5];
      printf "  %02d.%02d",$data_list[$i*16+6],$data_list[$i*16+7];
      print "\n";
    }
  }
}

# check if a hole id is in the valid range
sub ValidHole() {
  $_ = shift;
  return 1 if ( /^-?(\d+)$/ && $1>100 && $1<=2618 );
  return 0;
}

# read cable mapping table
sub GetMap1() {
  my $p  = $_[0];
  my $ho = $_[1];
  my $r  = $_[2];
  my $fe = $_[3];
  my $fc = $_[4];
  my $lb = $_[5];
  my $lc = $_[6];
  my $hb = $_[7];
  my $hc = $_[8];
  my $srt = $_[9];

  my $conn = Pg::connectdb("host=$dbhost port=5432 user=borex_guest dbname=bx_geometry");
  $conn->status && error("Error connecting to Data Base bx_geometry");
  my $query = "SELECT \"HoleLabel\",\"RackID\",\"FEBoard\",\"FEChannel\",\"LabenBoard\",\"LabenChannel\",\"HVBoard\",\"HVChannel\"  FROM \"CableMapping\" WHERE \"ProfileID\"=$p";
  if ($ho!=0) {
    $query = $query." AND \"HoleLabel\"=$ho";
  }
  if ($r>0) {
    $query = $query." AND \"RackID\"=$r";
  }
  if ($fe>0) {
    $query = $query." AND \"FEBoard\"=$fe";
  }
  if ($fc>0) {
    $query = $query." AND \"FEChannel\"=$fc";
  }
  if ($lb>0) {
    $query = $query." AND \"LabenBoard\"=$lb";
  }
  if ($lc>0) {
    $query = $query." AND \"LabenChannel\"=$lc";
  }
  if ($hb>0) {
    $query = $query." AND \"HVBoard\"=$hb";
  }
  if ($hc>0) {
    $query = $query." AND \"HVChannel\"=$hc";
  }
  if ($srt) {
    $query = $query." ORDER BY $names{$srt}"; 
  }
  
  my $result = $conn -> exec( $query );
  $result->resultStatus==PGRES_TUPLES_OK || error("Error reading from table \"Profiles\"");
  if ($result->ntuples == 0) {
    print "No data available. \n";
    print "Input parameters not consistent. \n";
    exit;
  }
  my $records = $result->ntuples;
  if ( $result->nfields != 8 ) {
    print "Error reading DB.\n";
    print "Wrong number of fields\n";
    exit;
  }
  # save into data array
  # not all fields are filled by this query
  @data_list=();
  my $i=0;
  my $f=0;
  for ($i = 0; $i < $records; $i++) {
    for ($f = 0; $f < 8; $f++) {
      $data_list[$i*16+$f] = $result->getvalue($i,$f);
    }
    for($f = 8; $f<16; $f++) {
      $data_list[$i*16+$f]=0;
    }
  }
  return $records;
}

sub GetMap2() {
  my $p  = $_[0];
  my $clust = $_[1];
  my $records = $_[2];
  
  printf("Estimated time to complete: %d seconds\n",$records/6);

  my $conn = Pg::connectdb("host=$dbhost port=5432 user=borex_guest dbname=bx_geometry");
  $conn->status && error("Error connecting to Data Base bx_geometry");
  my $query = "SELECT \"ChannelID\",\"X\",\"Y\",\"Z\",\"Theta\",\"Phi\",\"Conc\",\"FiberBundle\"  FROM \"HolesMapping\" WHERE \"ProfileID\"=$p";
  if ( $clust >= 0 ) {
    $query = $query." AND \"FiberBundle\"=$clust";
  }
  my @result=();
  for($i=0; $i<$records; $i++) {
    my $hole = $data_list[$i*16];
    $result[$i]=0;
    if ( &ValidHole($hole) ) {
      my $query1 = $query." AND \"HoleLabel\"=$hole";
      $result[$i] = $conn->exec( $query1 );
      $result[$i]->resultStatus==PGRES_TUPLES_OK || die("Error!!!");
      if ($result[$i]->ntuples == 1 ) {
        if ( $result[$i]->nfields != 8 ) {
          print "Error reading DB.\n";
          print "Wrong number of fields\n";
          exit;
        } 
        $fields = $result[$i]->nfields;
        for($f=0; $f < $fields; $f++) {
          $data_list[$i*16+8+$f] = $result[$i]->getvalue(0,$f);
        }
        ##$data_list[$i*16+15] = compute_cluster($data_list[$i*16])
      }
    }
    else {
      $data_list[$i*16]=0;
    }
  }
}

sub compute_cluster(){
   $ho = $_[0];
   $a = abs($ho);
   if ($a<100) {
     print "Invalid SSS Hole identifier\n";
     exit;
   }
   if ($a>999) {
     $ring=substr($a,0,2);
     $ind=substr($a,2,2);
   } else {
     $ring=substr($a,0,1);
     $ind=substr($a,1,2);
   }
   if ($ring>21 || $ring<1 || $ind>72 || $ind<1 || $a>2118) {
     print "Invalid SSS Hole identifier\n";
     exit;   
   }
   my $clus;
   if ($ring<=12) {
     if ($ind>72 || $ind<1) {
       print "Invalid SSS Hole identifier\n";
       exit;    
     }
     $clus = int(($ind-1)/6) + 3;
   }
   if ($ring>12 && $ring<=15) {
     if ($ind>36 || $ind<1) {
       print "Invalid SSS Hole identifier\n";
       exit;    
     }   
     $clus = int(($ind-1)/3) + 3;
   }
   if ($ring>15 && $ring<=18) {
     if ($ind>36 || $ind<1) {
       print "Invalid SSS Hole identifier\n";
       exit;    
     }   
     $clus = int(($ind-1)/18) + 1;
   }
   if ($ring>18) {
     if ($ind>18 || $ind<1) {
       print "Invalid SSS Hole identifier\n";
       exit;    
     }   
     $clus = int(($ind-1)/9) + 1;  
   }
   if ($ho<0) {
     if ($clus>2) {
       $clus += 12;
     }
     if ($clus==1) {
       $clus = 27;
     }
     if ($clus==2) {
       $clus = 28;
     }
   }
   return $clus;
}

