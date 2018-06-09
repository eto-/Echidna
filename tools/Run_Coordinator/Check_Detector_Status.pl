#!/usr/bin/perl -w
use Pg;

my @time = localtime(time);

my $month = $time[4];
my $year  = $time[5]+1900;
if($month==0) {
  $month=12;
  $year--;
}
my $help  = 0;
my $check_channels = 0;
my $check_spectrum = 0;

foreach (@ARGV){
  if(/month\=(.*)/){
    $month=$1;
  } 
  elsif(/year\=(.*)/){
    $year=$1;
  }
  elsif(/help/){
    $help = 1;
  }
  elsif(/check_channels/){
    $check_channels = 1;
  }
  elsif(/check_spectrum/){
    $check_spectrum = 1;
  }
}

if($help){
    print "\nUsage: ./Check_Detector_Status.pl [options]

Options:
  month=[month num]     Selection by month number (if month number X < 10, write 0X. Default: previous month).
  year=[year]           Self explaining. Optional.
  check_channels        Monitor the channels that are dead in neutrino, hot in neutrino or failing precalibrations
  check_spectrum        Extrapolate the LY from the 14C end-point and determine the 210Po peak position

";
  exit;
}

my @abbr = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
my $mm = $month-1;

print "\n";
print "Analyzing month: $abbr[$mm] $year.\n";
if($month<10) { $month = "0$month"; }
#print "Analyzing month: $month $year.\n";
#######################################
#### Calculation of the duty cycle ####
#######################################

# The month is considered to start at the end of the last run of the previous month
# It is also considered to finish at the end of the last run of the month
# The gap between the real beginning of the month and the end of the last run of the previous month is named $delay
# The gap between the real end of the month and the end of the last run of the month is named $delay_tail
# Validated runs, their start time and duration are taken from /home/production/run_validation_out.txt

my $file = 'temp.txt';
system "grep \"$year-$month\" /home/production/run_validation_out.txt > $file";
open IN,"<$file";

my @lines = <IN>;
close(IN);
system "rm $file";

@lines = sort@lines;

# Calculation of $delay_tail
my @last_line = split(" ",$lines[@lines-1]);
$last_line[20] =~ /^(.{2}):(.{2}):(.{2})$/;
my $delay_tail = $1*60*60 + $2*60 + $3 + $last_line[5] - 86400;
if( $delay_tail < 0 ) { $delay_tail = 0; }


my $bool = 0;
my $delay = 0;
my $livetime = 0;
my $last_day = 0;

my @run_list;
my $num_runs = 0;

for(my $i = 0; $i < @lines; $i++){

  my @line = split(" ",$lines[$i]);

  #Calculation of $delay
  if( $line[19] =~ /^$year-$month-.{2}$/ && $bool == 0 ) {
    my @prev_line = split(" ",$lines[$i-1]);
    $prev_line[20] =~ /^(.{2}):(.{2}):(.{2})$/;
    $delay = $1*60*60 + $2*60 + $3 + $prev_line[5] - 86400;    
    if( $delay < 0 ) { $delay = 0; }
    $bool = 1;
  }

  if( $line[19] =~ /^$year-$month-(.{2})$/ ) {
    if( $1 > $last_day ) { $last_day = $1; }
    $livetime += $line[5];
    $run_list[$num_runs] = $line[0];
    $num_runs++;
  }
}

my $total_time = $last_day*86400 - $delay + $delay_tail;
my $duty_cycle = sprintf("%3.2f",$livetime/$total_time);

print "\n";
print "***************************\n";
print "**** DUTY CYCLE:  $duty_cycle ****\n";
print "***************************\n";
print "\n";



if( !($check_channels) && !($check_spectrum) ){
  #print "If you want to monitor the \"bad\" channels, type the command \'./Check_Detector_Status.pl month=XX check_channels\'.\n";
  #print "Else, if you want to monitor the energy spectrum, type the command\'./Check_Detector_Status.pl month=XX check_spectrum\'.\n\n";
  print "Type the command \'./Check_Detector_Status.pl help\' to know how to monitor the \"bad\" channels or the spectrum.\n\n";
}



my $prod_list = 'production_log.list';
my $precalib_list = 'precalibrations_log.list';

if( $check_channels || $check_spectrum ){

  #####################################################################
  #### Producing the list of runs for the detector stability check ####
  #####################################################################

  print "Selecting the runs of the month...\n\n";

  for($i = 0; $i < @run_list; $i++){
    system "find /bxstorage/rootfiles/cycle_14/ -name \"Run0$run_list[$i]_c14.log\" >> $prod_list";
    if($check_channels) { system "find /bxstorage/rootfiles/cycle_14/ -name \"Run0$run_list[$i]_precalibrations_c14.log\" >> $precalib_list"; }
  }
}


if($check_channels){
# Channels failing precalibrations 
  print "Looking for channels failing precalibrations...\n\n";

  my $failing_file_name = "failing_precalib_vs_run_$abbr[$mm]_$year.txt";

  open FIN,"<$precalib_list";
  open FOUT1,">$failing_file_name";
  open FOUT2,">channels_failing_precalib_$abbr[$mm]_$year.txt";

  my @lg_fail;
  my @lg_fail_tot;

  for(my $i = 0; $i < 4000; $i++) { $lg_fail_tot[$i] = 1; }


  while(<FIN>){
    chomp $_;
    $_ =~ /^.*\/Run0(.{5})_precalibrations_c14.log$/;
    my $run = $1;
    open IN_TEMP,"<$_";

    for($i = 0; $i < 4000; $i++) { $lg_fail[$i] = 0; }

    my $failing_precalib = 0;
    while(<IN_TEMP>){
      if($_ =~ /^.*ordinary\slogical\schannel\s(.{1,4})\shas.*$/) { $lg_fail[$1] = 1; }
      if($_ =~ /^.*bad\selectronic\schannels\s(.{1,4})\,.*$/){ $failing_precalib += $1; }
      if($_ =~ /^.*off\selectronic\schannels\s(.{1,4})\,.*$/){ $failing_precalib += $1; }
    }
    close(IN_TEMP);

    for($i = 0; $i < @lg_fail_tot; $i++) { $lg_fail_tot[$i] = $lg_fail_tot[$i] && $lg_fail[$i]; }

    print FOUT1 "$run   $failing_precalib\n";
  }

  for($i = 1; $i < @lg_fail_tot; $i++) { if( $lg_fail_tot[$i] == 1 ) { print FOUT2 "$i\n"; } } 

  close(FIN);
  close(FOUT1);
  close(FOUT2);


# Channels hot_in_neutrino or dead_in_neutrino
  print "Looking for channels hot_in_neutrino or dead_in_neutrino...\n\n";

  my $dead_file_name = "dead_in_nu_vs_run_$abbr[$mm]_$year.txt";
  my $dead_od_file_name = "dead_in_nu_od_vs_run_$abbr[$mm]_$year.txt";
  my $hot_file_name  = "hot_in_nu_vs_run_$abbr[$mm]_$year.txt";

  open FIN,"<$prod_list";
  open FOUT,">$dead_file_name";
  open FOUT1,">$hot_file_name";
  open FOUT2,">channels_dead_in_neutrino_$abbr[$mm]_$year.txt";
  open FOUT3,">channels_hot_in_neutrino_$abbr[$mm]_$year.txt";
  open FOUT4,">$dead_od_file_name";

  my @channels;
  my @channels_od;
  my @channels_tot;
  my @lg_hot_in_nu;
  for(my $i = 0; $i < 4000; $i++) { $channels_tot[$i] = 1; $lg_hot_in_nu[$i] = 0; }

  while(<FIN>){
    chomp $_;
    $_ =~ /^.*\/Run0(.{5})_c14.log$/;
    my $run = $1;
    open IN_TEMP,"<$_";

    for($i = 0; $i < 4000; $i++) { $channels[$i] = 0; $channels_od[$i] = 0; }

    my $hot_in_nu = 0;
    my $dead_in_nu = 0;
    my $dead_in_nu_od = 0;
    my $dead_in_nu_tot = 0;
    while(<IN_TEMP>){
      if($_ =~ /^.*channel\s(.*)\sbecause\sof\shot_in_neutrino.*$/){
        $lg_hot_in_nu[$1] = 1;
        $hot_in_nu++;
      }
      if($_ =~ /^.*channel\s(.*)\sbecause\sof\sdead_in_neutrino.*$/){
        $dead_in_nu_tot++;
        if($1 < 3000) { 
          $dead_in_nu++; 
          $channels[$1] = 1;
        }
        elsif($1 > 2999){
          $dead_in_nu_od++;
          $channels_od[$1] = 1;
        }
      }
    }
    close(IN_TEMP);

    for($i = 0; $i < @channels_tot; $i++) { 
      if( $i < 3000 ) { $channels_tot[$i] = $channels_tot[$i] && $channels[$i]; }
      else            { $channels_tot[$i] = $channels_tot[$i] && $channels_od[$i]; }
    }

    print FOUT1 "$run   $hot_in_nu\n";
    print FOUT4 "$run   $dead_in_nu_od\n";
    print FOUT  "$run   $dead_in_nu\n";
  }

  for($i = 1; $i < @channels_tot; $i++) { 
    if( $channels_tot[$i] == 1 ) { print FOUT2 "$i\n"; } 
    if( $lg_hot_in_nu[$i] == 1 ) { print FOUT3 "$i\n"; }
  }

  close(FIN);
  close(FOUT);
  close(FOUT1);
  close(FOUT2);
  close(FOUT3);
  close(FOUT4);


  #Saving the plot of the previous variables as a function of run
  print "Preparing the plot of the previous variables as a function of run...\n\n";

  my $root_istro = "stability_plot_maker_$abbr[$mm]_$year.C";
  open OUT,">$root_istro";
  print OUT<<END;
{
  gROOT->LoadMacro("/home/production/Run_Coordinator/Run_Coordinator/SetAliases.C");
  a = plotxy("$dead_file_name");
  a->SetMarkerColor(kRed);
  a->SetMarkerStyle(9);
  TCanvas c1;
  a->GetYaxis()->SetRangeUser(0,85);
  a->GetXaxis()->SetTitle("Run Number");
  a->Draw("ap");

  b = plotxy("$dead_od_file_name");
  b->SetMarkerColor(kRed+2);
  b->SetMarkerStyle(9);
  b->Draw("psame");

  c = plotxy("$hot_file_name");
  c->SetMarkerColor(kOrange);
  c->SetMarkerStyle(9);
  c->Draw("psame");

  d = plotxy("$failing_file_name");
  d->SetMarkerColor(kBlue);
  d->SetMarkerStyle(5);
  d->Draw("psame");

  TLegend* leg = new TLegend(0.5,0.72,0.85,0.91);
  leg->AddEntry(a,"Dead in Neutrino ID","P");
  leg->AddEntry(b,"Dead in Neutrino OD","P");
  leg->AddEntry(c,"Hot in Neutrino","P");
  leg->AddEntry(d,"Failing Precalibration","P");
  leg->SetShadowColor(kWhite);
  leg->SetFillColor(kWhite);
  leg->Draw();

  c1->SaveAs("stability_plot_channels_dead_hot_and_failing_precalib_$abbr[$mm]_$year.pdf");
}
END

  system "root -l -b -q $root_istro > a.txt; rm a.txt";
  system "rm $prod_list $precalib_list";
}



if($check_spectrum){
  open LIST,"<$prod_list";
  my $rootfile = "root_instructions.C";

  while(<LIST>){
    $_ =~ /^\/bxstorage\/rootfiles\/cycle_14\/(.{4})\/(.{6})\/Run0(.{5})_c14.log$/;
    print "Determining 14C-LY and 210Po peak for run $3\n";

    open ROOTCMD,">$rootfile";
    print ROOTCMD<<END;
{
  gSystem->Load("/home/production/Run_Coordinator/Run_Coordinator/c14_spectrum_C.so");
  C14_Fitter("$1","$abbr[$mm]","$2","$3",50,90,500);
}
END
    close(ROOTCMD);

    system "root -l -q -b $rootfile >> fit.log";
  }

  my $root_istro = "root_instructions.C";
  open OUT,">$root_istro";
  print OUT<<END;
{
  gROOT->LoadMacro("/home/production/Run_Coordinator/Run_Coordinator/SetAliases.C");
  a = plotxyerror("c14_ly_vs_run_$abbr[$mm].txt");
  a->SetMarkerColor(kRed);
  a->SetMarkerStyle(9);
  TCanvas c1;
  a->GetYaxis()->SetRangeUser(150,560);
  a->GetXaxis()->SetTitle("Run Number");
  a->Draw("ap");

  b = plotxyerror("po210_peak_vs_run_$abbr[$mm].txt");
  b->SetMarkerColor(kBlue);
  b->SetMarkerStyle(9);
  b->Draw("psame");

  TLegend* leg = new TLegend(0.5,0.72,0.85,0.91);
  leg->AddEntry(a,"14C Light Yield","P");
  leg->AddEntry(b,"210Po Peak","P");
  leg->SetShadowColor(kWhite);
  leg->SetFillColor(kWhite);
  leg->Draw();

  c1->SaveAs("stability_plot_c14_ly_po210_peak_$abbr[$mm]_$year.pdf");
}
END

  system "root root_instructions.C";
  system "rm $prod_list $precalib_list";
}
