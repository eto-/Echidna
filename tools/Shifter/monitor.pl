#!/usr/bin/perl

# mass normalization
$tons = 270.2 ;

$run = $ARGV[0];

$filtering = 1 ;
if($ARGV[1]) {$filtering = 0;}

# check if loop or single mode
if(!$run) {$loop = 1;}

# single mode
if(!$loop) {
  $file = 'pippo';
#  open(LS,'ls  /bxstorage/rootfiles/cycle_13/2011/*/*_c13.root |');   # c13->c14
#  open(LS,'ls  /bxstorage/rootfiles/cycle_14/2011/*/*_c14.root |');
  open(LS,'ls  /bxstorage/rootfiles/cycle_14/*/*/*_c14.root |');
  while(<LS>) {
#    if(/Run0${run}\_c13\.root/){   # c13->c14
    if(/Run0${run}\_c14\.root/){
      chomp($_);
      $file = $_;
    }
  } close(LS);
  if($file =~ 'pippo') { print "file not found!\n"; exit ;} 
  exe();
  
# loop mode
} else {
  print "\nloop mode\n\n";
  open(INF,'inputfiles.dat');
  while(<INF>) {
#    if(/Run0(\d+)\_c13\.root/) {   # c13->c14
    if(/Run0(\d+)\_c14\.root/) {
      $run = $1 ;
      $file = $_;
      chomp($file);
      print "$run $file \n";
      exe();
    }
  } close(INF);

}

  
sub exe() {
  print "Starting monitor.pl for run $run\n";
  #remove files and date corresponding to $run
  system("perl clean.pl $run");   
  print "Cleaning..... done!\n";
  # create all needed files
  $list = "l${run}.list";
  open(IN,">$list");
  print IN "http://bxmaster-data/$file";
  close(IN);
  $outfile = "rootfiles/monitor\_run${run}\.root";
  $hook = "rootfiles/monitor\_run${run}\.hook";


#print "Mo te creo un root file: $outfile \n";

  # start bxfilter
  if($filtering) {system("./main $list $outfile");}

  # remove useless files
  system("rm -f $list; rm -f $hook"); 

  # look for the livetime
  open(LT,'/home/production/run_validation_out.txt');
  while(<LT>) {
    if(/${run}\s+\d+\s+\d+\s+\d+\s+\d+\s+(\d+)/) { $LT = $1; }
  } close(LT);

  print "Livetime = $LT seconds\n"; 

  # write macro for monitor.C
  $macro = "macro/macro\_run${run}.C";
  open(MAC,">${macro}");
  print MAC<<EOF;
{ 
gROOT->LoadMacro("monitor.C");
monitor(\"$outfile\", $LT);
}

EOF

  close(MAC);
  print "write macro...\n";
  # run monitor.C and  take the data results
  open(ROO,"root.exe -q -b $macro |");
  while(<ROO>) {

    if(/Po210\s+(\d+\.*\d*)\s+(\d+\.*\d*)/) {
      $po210 = $1/$LT*86400/$tons;
      $po210_e = $2/$LT*86400/$tons;
      print "Po210 count rate = $po210 +- $po210_e cpd/t\n";    
    }

    if(/NPo210N\s+(\d+\.*\d*)\s+(\d+\.*\d*)/) {
      $Npo210N = $1/$LT*86400/145.4;
      $Npo210_eN = $2/$LT*86400/145.4;
      print "Po210 count rate in the North = $Npo210N +- $Npo210_eN cpd/t\n";    
    }

    if(/SPo210S\s+(\d+\.*\d*)\s+(\d+\.*\d*)/) {
      $Spo210S = $1/$LT*86400/124.96;
      $Spo210_eS = $2/$LT*86400/124.96;
      print "Po210 count rate in the South = $Spo210S +- $Spo210_eS cpd/t\n";    
    }

    if(/Radon\s+(\d+)/) {
      $radon = $1/$LT*86400/$tons; 
      $radon_e = sqrt($1)/$LT*86400/$tons;
      print "Radon count rate = $radon +- $radon_e cpd/t\n";    
    }  
    
    if(/Thoron\s+(\d+)/) {
      $thoron = $1/$LT*86400; 
      $thoron_e = sqrt($1)/$LT*86400;
      print "Thoron count rate = $thoron +- $thoron_e cpd\n";    
    }  

    if(/Day\s+(\-*\d+\.*\d*)/) {
      $day = $1; 
      print "Day of run $run = $day\n";        
    }  

    if(/Count\s+(\d+)/) {
      $count = $1/$LT*86400/$tons;
      $count_e = sqrt($1)/$LT*86400/$tons;
      print "Count rate = $count +- $count_e cpd/t\n";    
    }

    if(/CountZ0\s+(\d+)/) {
      $countZ0 = $1/$LT*86400/145.2; #/calotta(0); 145.2 tons from CCD cameras 
      $countZ0_e = sqrt($1)/$LT*86400/145.2; #/calotta(0); 145.2 tons from CCD cameras /calotta(0);
      print "CountZ0 rate = $countZ0 +- $countZ0_e cpd/t\n";    
    }

    if(/CountZ1\s+(\d+)/) {
      $countZ1 = $1/$LT*86400/97.6; #/calotta(0); 97.6 tons from CCD cameras
      $countZ1_e = sqrt($1)/$LT*86400/97.6; #/calotta(0); 97.6 tons from CCD cameras
      print "CountZ1 rate = $countZ1 +- $countZ1_e cpd/t\n";    
    }

    if(/CountZ2\s+(\d+)/) {
      $countZ2 = $1/$LT*86400/54.0; #/calotta(0); 54.0 tons from CCD cameras
      $countZ2_e = sqrt($1)/$LT*86400/54.0; #/calotta(0); 54.0 tons from CCD cameras
      print "CountZ2 rate = $countZ2 +- $countZ2_e cpd/t\n";    
    }

    if(/CountZ3\s+(\d+)/) {
      $countZ3 = $1/$LT*86400/20.3; #/calotta(0); 20.3 tons from CCD camera
      $countZ3_e = sqrt($1)/$LT*86400/20.3; #/calotta(0); 20.3 tons from CCD camera
      print "CountZ3 rate = $countZ3 +- $countZ3_e cpd/t\n";    
    }

    if(/Peak\s+(\d+\.*\d*)\s+(\d+\.*\d*)/) {
      $peak   = $1;
      $peak_e = $2;
      print "Peak = $peak +- $peak_e \n";    
    }
  } close(ROO);

  $LTday = $LT/86400 ;
  open(OUT,">>monitor_results.dat");
  print OUT "$run $day $LTday $count $count_e  $countZ0 $countZ0_e $countZ1 $countZ1_e $countZ2 $countZ2_e $countZ3 $countZ3_e $radon $radon_e $po210 $po210_e $peak $peak_e $Npo210N $Npo210_eN $Spo210S $Spo210_eS\n";
  close(OUT);
  open(OUT,">>thoron_results.dat");
  print OUT "$run $day $LTday $thoron $thoron_e\n";
  close(OUT);
  
  
  # sort  monitor_results.dat
  system("perl detector.pl $run");
  system('perl sort.pl monitor_results.dat');
  system('perl sort.pl info_detector.dat');
  
  # clean corrupted data
  system("perl automatic_clean_data.pl");

  # move gifs to the right folders
  system('mv /home/production/www/gif/xz_*gif /home/production/www/xz');
  system('mv /home/production/www/gif/yz_*gif /home/production/www/yz');
  system('mv /home/production/www/gif/xy_*gif /home/production/www/xy');
  system('mv /home/production/www/gif/rz_*gif /home/production/www/rz');
  system('mv /home/production/www/gif/spectrum_*gif /home/production/www/spectra');
    
  # produce the animated gifs  
  system("perl converter.pl /home/production/www/xz");
  system("perl converter.pl /home/production/www/yz");
  system("perl converter.pl /home/production/www/xy");
  system("perl converter.pl /home/production/www/rz");
  
  #system("convert -delay 50 -loop 0  /home/production/www/xz/*.gif /home/production/www/xz_animated.gif");
  #system("convert -delay 50 -loop 0  /home/production/www/yz/*.gif /home/production/www/yz_animated.gif");
  #system("convert -delay 50 -loop 0  /home/production/www/xy/*.gif /home/production/www/xy_animated.gif");
  #system("convert -delay 50 -loop 0  /home/production/www/rz/*.gif /home/production/www/rz_animated.gif");
  #system("convert -delay 50 -loop 0  /home/production/www/spectra/*.gif /home/production/www/spectra_animated.gif");
  
  # sort radon data
  system('perl sortradon.pl');
  system('perl sortthoron.pl');
  system('perl cleanradon.pl');
  system('perl cleanthoron.pl');

  print "About to lounch the radon.exe programm. \n";
  
  system("/home/production/monitor/bxfilter/radon.exe");

  print "The radon.exe programm is done. \n";
  
  #system("convert -delay 50 -loop 0  /home/production/www/radon/xz/*.gif /home/production/www/radon/xz_rn_animated.gif");
  #system("convert -delay 50 -loop 0  /home/production/www/radon/yz/*.gif /home/production/www/radon/yz_rn_animated.gif");
  #system("convert -delay 50 -loop 0  /home/production/www/radon/xy/*.gif /home/production/www/radon/xy_rn_animated.gif");
  #system("convert -delay 50 -loop 0  /home/production/www/radon/rz/*.gif /home/production/www/radon/rz_rn_animated.gif");
  
  system("perl converter.pl /home/production/www/radon/xz");
  system("perl converter.pl /home/production/www/radon/yz");
  system("perl converter.pl /home/production/www/radon/xy");
  system("perl converter.pl /home/production/www/radon/rz");
  
  
  # produce the summary plot
  system("/home/production/monitor/bxfilter/graphzoom.exe 0 0");
  system("/home/production/monitor/bxfilter/graphthoron.exe 0 0");
  system("/home/production/monitor/bxfilter/graphzoomdetector.exe 0 0");
  
  
  # update the web site
  system('perl index.pl');
  
  system('cp monitor_results.dat www/');
  system('cp thoron_results.dat www/cgi-bin/');
  system('cp thoron_results.dat www/');
  system('cp thoron.dat www/cgi-bin/');
  system('cp radon.dat www/');
  
  print "DONE\n";
  
}

 
sub calotta {
  my($z) = shift;
  $r = 4.25;
  $h = $r - $z;
  $ratio = ($h*$h)*($r - $h/3.)/(4./3.*$r*$r*$r);
  return $ratio ;  
}
