#!/usr/bin/perl

@files = glob('/home/room3/production/Echidna/logs/*_calib_electronincs.log');
$output     = 'general.dat';
$outchannel = 'channel.dat';
$num_of_files = @files;
if($ARGV[0] =~ '--help' || $ARGV[0] =~ '-h') {
  print "Usage: ./electronics.pl [channel number] [root]\n\n";
  print "    -h or --help:	   print this message \n";
  print "    no input:  	   look the summary plot with root\n";
  print "    channel number:	   insert the channel number; produce a gif for the single channel number\n";
  print "    root:		   in the channel case, show the plot with root\n";
  print "\n ...by davide.franco\@mi.infn.it\n";
  exit(0); 
}
if(!$ARGV[0]) {
  &general ;
} else {
  &channel_chooser($ARGV[0]);
}
#general
sub general {
  open(OUT,">$output");
  foreach $file (@files) {
    if(/(\d+)_calib_electronincs.log/) {
      $run = $1
    }
    @Status = ("DEAD","Hot", "Low_eff", "Retriggering");
    open(IN,$file);
    while(<IN>) {
      if(/bx_calib_laben_electronics/) {
	if(/Lg (\d+):* (\w+)_in_\s*(\w+)\s+/) {
	  $type = $3; # neutrino laser pulser
	  chomp($type);
	  $channel = $1;  
	  $status = $2; # DEAD Hot Low_eff
	  if($status =~ "DEAD") {
	    $status = 0;
	  } elsif ($status =~ "Hot") {
	    $status = 1;
	  } elsif ($status =~ "Low_eff") {
	    $status = 2;
	  } elsif ($status =~ "Retriggering") {
	    $status = 3;
	  }	  
	  if($type =~ "neutrino") {
	    if($neutrino[$status][$channel]) {
	      $neutrino[$status][$channel] += 1;
	    } else {
	      $neutrino[$status][$channel] = 1;
	    }
	  } elsif($type =~ "laser") {
	    if($laser[$status][$channel]) {
	      $laser[$status][$channel] += 1;
	    } else {
	      $laser[$status][$channel] = 1;
	    }
	  } elsif($type =~ "pulser") {
	    if($pulser[$status][$channel]) {
	      $pulser[$status][$channel] += 1;
	    } else {
	      $pulser[$status][$channel] = 1;
	    }	  
	  }
	}
      }
    } close(IN);
  }
  
  for($i=1; $i<=2240;$i++) {
    print     "$i ";
    print OUT "$i ";
    for ($k = 0; $k< 4 ; $k++) {
    
      $value_neutrino = $neutrino[$k][$i];
      $value_pulser = $pulser[$k][$i];
      $value_laser = $laser[$k][$i];
      if(!$neutrino[$k][$i]) {$value_neutrino = 0;}
      if(!$pulser[$k][$i]) {$value_pulser = 0;}
      if(!$laser[$k][$i]) {$value_laser = 0;}
      
      print     " $value_neutrino $value_pulser $value_laser";
      print OUT " $value_neutrino $value_pulser $value_laser";

    }
    print OUT "\n";
    print     "\n";
  }
    #system("root.exe electronics.C");
}


sub channel_chooser(){
  open(OUT,">$outchannel");
  my $channel = shift ;
  chomp($ARGV[1]);
  print OUT "$channel $num_of_files\n";
  foreach $file (@files) {
    if($file =~ m/.*\/(\d+)_calib/) {
      $run = $1;
      push(@runs,$run);
    }
    open(IN,$file);
    while(<IN>) {
      if(/bx_calib_laben_electronics/) {
	if(/Lg $channel:* (\w+)_in_\s*(\w+)\s+/) {
	  $type = $2; # neutrino laser pulser
	  chomp($type);
	  $status = $1; # DEAD Hot Low_eff
	  if($status =~ "DEAD") {
	    $status = 0;
	  } elsif ($status =~ "Hot") {
	    $status = 1;
	  } elsif ($status =~ "Low_eff") {
	    $status = 2;
	  } elsif ($status =~ "Retriggering") {
	    $status = 3;
	  }	  
	  if($type =~ "neutrino") {
	    if($neutrino[$status][$run]) {
	      $neutrino[$status][$run] += 1;
	    } else {
	      $neutrino[$status][$run] = 1;
	    }
	  } elsif($type =~ "laser") {
	    if($laser[$status][$run]) {
	      $laser[$status][$run] += 1;
	    } else {
	      $laser[$status][$run] = 1;
	    }
	  } elsif($type =~ "pulser") {
	    if($pulser[$status][$run]) {
	      $pulser[$status][$run] += 1;
	    } else {
	      $pulser[$status][$run] = 1;
	    }	  
	  }
	}
      }

    } close(IN);
  }
  foreach (@runs) {
    print  "$_ ";
    print OUT "$_ ";
    for ($k = 0; $k< 4 ; $k++) {
    
      $value_neutrino = $neutrino[$k][$_];
      $value_pulser = $pulser[$k][$_];
      $value_laser = $laser[$k][$_];
      if(!$neutrino[$k][$_]) {$value_neutrino = 0;}
      if(!$pulser[$k][$_]) {$value_pulser = 0;}
      if(!$laser[$k][$_]) {$value_laser = 0;}
      
      print OUT " $value_neutrino $value_pulser $value_laser   ";
      print     " $value_neutrino $value_pulser $value_laser   ";

    }
    print OUT "\n";
    print     "\n";
 
  }  
  close(OUT);
  if($ARGV[1] =~ "root") { 
    system("root.exe electronics_channel.C");
  } else {
    system("root.exe -b -q electronics_channel.C"); 
  }


}










