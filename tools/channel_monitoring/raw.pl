#!/usr/bin/perl

@files = glob('/home/room3/production/Echidna/logs/*_calib_electronincs.log');
$output     = 'general.dat';
$outchannel = 'channel.dat';
if($ARGV[0] =~ '--help' || $ARGV[0] =~ '-h') {
  print "Usage: ./raw.pl [channel number] [root]\n\n";
  print "    -h or --help:	   print this message \n";
  print "    no input:  	   look the summary plot with root\n";
  print "    channel number:	   insert the channel number; produce a gif for the single channel number\n";
  print "    root:		   in the channel case, show the plot with root\n";
  print "\n ...by davide.franco\@mi.infn.it\n";
  exit(0); 
}
  &general ;

#channel_history
sub general() {
  foreach $file (@files) {
    if($file =~ m/.*\/(\d+)_calib/) {
      $run = $1;
      push(@runs,$run);
    }
    open(IN,$file);
    while(<IN>) {
      if(/bx_calib_laben_electronics: Lg (\d+): fifo (\w+)\s+/) {
        $channel = $1;
        $type = $2;
	#print "$run $channel $type\n";
	if($type =~ "full") {
	  $full{$channel}++;
	  if($ARGV[0] == $channel) {$fullchannel{$run}++;}
	}
	if($type =~ "empty") {
	  $empty{$channel}++;
	  if($ARGV[0] == $channel) {$emptychannel{$run}++;}
	} 
      }
      if(/bx_calib_laben_electronics: Lg (\d+) only (\d+\.*\d*)\s+\%.*\s+(\w+)\s+triggers/) {
        $trigger = $3;
	$value = $2;
	$channel = $1;
	if($trigger =~ "neutrino") {
	  $neutrino{$channel}++;
	  if($ARGV[0] == $channel) {$neutrinochannel{$run} = $value;}
	} elsif ($trigger =~ "pulser") {
	  $pulser{$channel}++;	
	  if($ARGV[0] == $channel) {$pulserchannel{$run} = $value;}
	} elsif ($trigger =~ "laser") {
	  $laser{$channel}++;
	  if($ARGV[0] == $channel) {$laserchannel{$run} = $value;}
	}
	#print "$run $1 $2 $3\n";
      }
    }
  }
  if(!$ARGV[0]) {
    open(OUT,">$output");
    for ($i =  1; $i<=2240; $i++) {
      if(!$full{$i}) {$full{$i} = 0;}
      if(!$empty{$i}) {$empty{$i} = 0;}
      if(!$neutrino{$i}) {$neutrino{$i} = 0;}
      if(!$laser{$i}) {$laser{$i} = 0;}
      if(!$pulser{$i}) {$pulser{$i} = 0;}
      print "$i $full{$i} $empty{$i} $neutrino{$i} $laser{$i} $pulser{$i}\n";
      print OUT "$i $full{$i} $empty{$i} $neutrino{$i} $laser{$i} $pulser{$i}\n";
    }
  } else {
    open(OUT,">$outchannel");
    print OUT "$ARGV[0]\n";
    foreach $run (@runs)  {
      if(!$fullchannel{$run}) {$fullchannel{$run} = 0;}
      if(!$emptychannel{$run}) {$emptychannel{$run} = 0;}
      if(!$neutrinochannel{$run}) {$neutrinochannel{$run} = -1;}
      if(!$laserchannel{$run}) {$laserchannel{$run} = -1;}
      if(!$pulserchannel{$run}) {$pulserchannel{$run} = -1;}
      print "$run $fullchannel{$run} $emptychannel{$run} $neutrinochannel{$run} $laserchannel{$run} $pulserchannel{$run}\n";
      print OUT "$run $fullchannel{$run} $emptychannel{$run} $neutrinochannel{$run} $laserchannel{$run} $pulserchannel{$run}\n";
    }
  }
  close(OUT);
  if($ARGV[0]) {
    if($ARGV[1] =~ "root") { system("root.exe raw_channel.C"); } else { system("root.exe  -b -q raw_channel.C") ;}
  } else {
    system("root.exe raw.C");
  }
}
