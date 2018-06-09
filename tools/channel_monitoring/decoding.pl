#!/usr/bin/perl

@files = glob('/home/room3/production/Echidna/logs/*_calib_electronincs.log');
$output     = 'general.dat';
$outchannel = 'channel.dat';
$num_of_files = @files;
if($ARGV[0] =~ '--help' || $ARGV[0] =~ '-h') {
  print "Usage: ./decoding.pl [channel number] [root]\n\n";
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
  &channel_history($ARGV[0]);
}

#channel_history
sub channel_history() {
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
      if(/bx_calib_laben_decoding/) {
	if(/channel $channel has less hits/) {
	  if($nohits{$run}) {
	    $nohits{$run}++;  
	  } else {
	    $nohits{$run} = 1;
	  }
	}

	if(/channel $channel too many/) {
          if($toohits{$run}) {
	    $toohits{$run}++;  
	  } else {
	    $toohits{$run} = 1;
	  }
	}

	if(/channel $channel , peak_mean is (.*) and peak_rms is (.*)/) {
	    $peak{$run} = $1;  
	    $rms{$run}  = $2;
	}
      }
    } close(IN);
  }
  foreach (@runs) {
    if(!$nohits{$_})  {$nohits{$_} = 0;}
    if(!$toohits{$_}) {$toohits{$_} = 0;}
    if(!$peak{$_})    {$peak{$_} = 1000.;}
    if(!$rms{$_})     {$rms{$_} = 1000.;}
    print     "$_ $nohits{$_} $toohits{$_} $peak{$_} $rms{$_} \n";
    print OUT "$_ $nohits{$_} $toohits{$_} $peak{$_} $rms{$_} \n";
  }  
  close(OUT);
  if($ARGV[1] =~ "root") { 
    system("root.exe decoding_channel.C");
  } else {
    system("root.exe -b -q decoding_channel.C"); 
  }
}

#general
sub general {
  open(OUT,">$output");
  foreach $file (@files) {
    if(/(\d+)_calib_electronincs.log/) {
      $run = $1
    }
    open(IN,$file);
    while(<IN>) {
      if(/bx_calib_laben_decoding/) {
	if(/channel (\d+) has less hits/) {
          if($nohits{$1}) {
	   $nohits{$1}++;  
	  } else {
	   $nohits{$1} = 1;
	  }
	}

	if(/channel (\d+) too many/) {
          if($toohits{$1}) {
	   $toohits{$1}++;  
	  } else {
	   $toohits{$1} = 1;
	  }
	}

	if(/channel (\d+) , peak/) {
          if($peak{$1}) {
	   $peak{$1}++;  
	  } else {
	   $peak{$1} = 1;
	  }
	}
      }
    } close(IN);

  }

  for ($i = 1; $i<=2240; $i++) {
    if(!$nohits{$i})  { $nohits{$i} = 0; }
    if(!$toohits{$i}) { $toohits{$i} = 0; }
    if(!$peak{$i})    { $peak{$i} = 0; }
    print     "$i $nohits{$i} $toohits{$i} $peak{$i}\n";
    print OUT "$i $nohits{$i} $toohits{$i} $peak{$i}\n";
  }
  close(OUT);  
  
  system("root.exe decoding.C");

}













