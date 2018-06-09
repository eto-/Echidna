#!/usr/bin/perl -w
use strict;

###########################################################################
# This script computes the duty-cycle of the week.							 		  	  #
# Data are taken from:							  																		#
# - http://bxmaster.lngs.infn.it//home/production/run_validation_out.txt  #
# - http://bxmaster.lngs.infn.it//home/production/problematic_runs.txt    #
# - Borexino "daq_config" DataBase		  			  													#
#		      						  	  																							#
# December 2012 - alessandra.re@mi.infn.it			  	  										#
# July     2014 - A. Re, general review + check on ECHIDNA & MOE	  			#
# 			  						 	 valid events		  																#	
# October  2014 - A. Re, new first day time evaluation										#
###########################################################################

use Date::Calc qw(:all);
use Date::Parse;
use Term::ANSIColor qw(:constants);
use DBI;

sub DataBase_connection ; 
sub RunCheckList ;
sub DutyCycle ;

#---------------------------------------
#	STANDARD INPUTS	FILES
#---------------------------------------
my $file_inRV = '/home/production/run_validation_out.txt';
my $file_inPR = '/home/production/problematic_runs.txt';
my $file_tmpDB = 'tempDB.txt';
my $file_tmpRV = 'tempRV.txt';

#---------------------------------------
#	GLOBAL SETTINGS
#---------------------------------------
my $days_of_shift = 7;
my $searched_dow = 2;						# 2 = Tuesday
my $searched_day = Day_of_Week_to_Text($searched_dow);

my @today = Today();
my $current_dow = Day_of_Week(@today);
my $today_dow = Day_of_Week_to_Text($current_dow);
my $today_date = sprintf("%d-%.2d-%.2d", @today);

#---------------------------------------
#	INPUTS REQUEST
#---------------------------------------
print BOLD, CYAN, "\n\t**************************************************\n";
print "\t*   	   BOREXINO END-SHIFT PROCEDURE          *\n";
print "\t*\t\t\t\t\t\t *\n";
print "\t*   - Duty-cycle computation   		         *\n";
print "\t*   - Echidna/MOE validated events cross-check   *\n";
print "\t**************************************************\n", RESET;
print "\n\tHello! Today is $today_dow, $today_date.\n"; 
print "\tPlease, enter your shift starting date (yyyy-mm-dd)\n";
print "\tor press ENTER for default (last-shift week) evaluation.\n\t";
chomp(my $given_week = <STDIN>);

#---------------------------------------
#	OPTION 1 - THERE IS AN INPUT
#---------------------------------------
if($given_week) {
	my @start = split(/-/,$given_week);
	my $dow = Day_of_Week($start[0],$start[1],$start[2]);
	if($dow != $searched_dow) {
		print RED "\n\tError! No shift started on $given_week.\n\n", RESET;
		exit();
	}	

	my @end = Add_Delta_Days(@start,+6);
	my $SSD = sprintf("%d-%.2d-%.2d", @start);		# SSD = Shift Starting Date
	my $start_shift = str2time($SSD);
	my $SED = sprintf("%d-%.2d-%.2d", @end);			# SED = Shift Ending Date
	my $end_shift = str2time($SED);
	print "\n\tYou are analyzing the shift started on $searched_day, $SSD.\n\n";
  my $firstrun_time = DataBase_connection($start_shift,$end_shift);
	my ($check_flag,$num_runs,$calib_time) = RunCheckList($start[0],$start[1],$start[2]);
	if($check_flag==0) {
		DutyCycle($num_runs,$calib_time,$firstrun_time);
	} else {
		print RED "\n\n\tBefore computing the dutycycle you must address all missing infos.\n\n", RESET;
		system("rm $file_tmpRV $file_tmpDB");
		exit();
	}

#---------------------------------------
#	OPTION 2 - NO GIVEN INPUT: GO DEFAULT
#---------------------------------------
} else {
	my (@prev, @prev2, @next);
	if($searched_dow == $current_dow) {
		@prev = Add_Delta_Days(@today,-7);
		@prev2 = Add_Delta_Days(@today,-14);
	} elsif($searched_dow > $current_dow) {
		@next = Add_Delta_Days(@today, $searched_dow - $current_dow);
		@prev = Add_Delta_Days(@next,-7);
		@prev2 = Add_Delta_Days(@next,-14);
	} else {
		@prev = Add_Delta_Days(@today, $searched_dow - $current_dow);
		@prev2 = Add_Delta_Days(@today, ($searched_dow - $current_dow) - 7);
		@next = Add_Delta_Days(@prev,+7);
	}

	my @start = @prev2;
	my @end = Add_Delta_Days(@start,+6);
	my $SSD = sprintf("%d-%.2d-%.2d", @start);		# SSD = Shift Starting Date
	my $start_shift = str2time($SSD);
	my $SED = sprintf("%d-%.2d-%.2d", @end);			# SED = Shift Ending Date
	my $end_shift = str2time($SED);
	print "\tYou are analyzing the shift started on $searched_day, $SSD.\n\n";
  my $firstrun_time = DataBase_connection($start_shift,$end_shift);
	my ($check_flag,$num_runs,$calib_time) = RunCheckList($start[0],$start[1],$start[2]);
	if($check_flag==0) {
		DutyCycle($num_runs,$calib_time,$firstrun_time);
	} else {
		print RED "\n\n\tBefore computing the dutycycle you must address all missing infos.\n\n", RESET;
		system("rm $file_tmpRV $file_tmpDB");
		exit();	
	}
}


#---------------------------------------
# SUBROUTINES
#---------------------------------------

#	CONNECTION TO DATA-BASE
sub DataBase_connection() {
	my $start_shift = shift;
	my $end_shift = shift;

	my $firstrun_flag = 0;
	my $firstrun_time;

	open TEMPDB,"> $file_tmpDB" or die $!;
 
	my $username = "borex_guest";
	my $password = "xyz";
	my $dbhost = "bxdb.lngs.infn.it";
	my $dbname = "daq_config";
	my $db = DBI->connect("DBI:Pg:dbname=$dbname;host=$dbhost", "$username", "$password");
	die RED, "\tDatabase connection absent: $DBI::errstr\n", RESET if(!$db);

	my $dbdata = "SELECT \"RunNumber\", \"Type\", \"StartTime\", \"Events\", \"Duration\" FROM \"Run\" ORDER BY \"RunNumber\"";
	my $query = $db->prepare($dbdata);	
	$query->execute();
	die RED, "$DBI::errstr\n", RESET if(!$query);
		
	while(my @array = $query->fetchrow_array) {
		my ($RunNumber, $Type, $StartTime, $Events, $Duration) = @array;
		my @db_info = split(/ +/,$StartTime);	
		#	$db_info[0]	-->	Start date (yyyy-mm-dd)
		#	$db_info[1]	-->	Start time (hh:mm:ss)
		my $date_test = str2time($db_info[0]);
		if(($date_test == $start_shift) && ($firstrun_flag == 0)) {
			$firstrun_time = $db_info[1];
			$firstrun_flag = 1;
		}	
		print TEMPDB "$RunNumber $Duration $Type\n" if(($date_test >= $start_shift) && ($date_test <= $end_shift) && ($Events > 20000));
	}
	$query->finish;
	$db->disconnect; 
	close TEMPDB;
	return($firstrun_time);
}


#	CROSS-CHECK BETWEEN DATABASE, RUN_VALIDATION_OUT.txt AND PROBLEMATIC_RUNS.txt
# CHECK ON MOE/ECHIDNA PRESENCE AND GOODNESS
sub RunCheckList() {
	local $SIG{__WARN__} = sub{warn @_ unless $_[0] =~ m(^.* too (?:big|small));};	
	# Previous line suppresses warnings due to an old version of perl LWP::UserAgent on BxMaster.
	
	my @start = (shift, shift, shift);

	# Open and read run_validation_out.txt
	open VALID, $file_inRV or die $!;
	my @rowsRV = <VALID>;
	close VALID;

	# Open and read data_base file
	open DBDATA, $file_tmpDB or die $!;
	my @rowsDB = <DBDATA>;		
	close DBDATA;
	
	# Open tmp output file
	open TEMPRV,"> $file_tmpRV" or die $!;

	my $Nrow = 0;
	LINE:	foreach my $rowRV(@rowsRV) {
		$Nrow++;
		next LINE unless($Nrow >= 206); 		# Row 206 corresponds to the first Phase I validated run (5007)

		my @info_rowRV = split(/ +/,$rowRV);
		for(my $i = 0; $i < $days_of_shift; $i++) {
			my ($run_number, $ech_cycle, $valid_ev, $livetime, $start_date, $start_time, $file_path);
			if($Nrow < 2460) {							# Row 2460 corresponds to the first 4 digits run number (10003)
				$run_number = $info_rowRV[1];
				$ech_cycle  = $info_rowRV[2];
				$valid_ev   = $info_rowRV[5];
				$livetime   = $info_rowRV[6]; 
				$start_date = $info_rowRV[20];
				$start_time = $info_rowRV[21];
				$file_path  = $info_rowRV[26];
			} else {
				$run_number = $info_rowRV[0];
				$ech_cycle  = $info_rowRV[1];
				$valid_ev   = $info_rowRV[4];
				$livetime   = $info_rowRV[5]; 
				$start_date = $info_rowRV[19];
				$start_time = $info_rowRV[20];
				$file_path  = $info_rowRV[25];
			}
			my @tmp = Add_Delta_Days(@start,+$i);
			my $date = sprintf("%d-%.2d-%.2d", @tmp);
			my $date_tmp = str2time($date);
			my $date_to_test = str2time($start_date);
			$file_path =~ /^http:\/\/bxmaster-data.lngs.infn.it\/(.*)$/;
			my $ECH_filepath = $1;
			print TEMPRV "$info_rowRV[0] $info_rowRV[1] $info_rowRV[4] $info_rowRV[5] $info_rowRV[19] $info_rowRV[20] $ECH_filepath\n" if($date_tmp == $date_to_test);
		}
	}
	close TEMPRV;

	my $num_runs = 0;
	my $check_flag = 0;
	my $calib_time = 0;
	foreach my $rowDB(@rowsDB) {
		my @info_rowDB = split(/ +/,$rowDB);
		#	$info_rowDB[0]	-->	Run number 	
		#	$info_rowDB[1]	-->	Duration
		#	$info_rowDB[2]  -->	Type of run
		if($info_rowDB[2] =~ /^normal$/) {
			my $checkPR = int(`grep -c \"$info_rowDB[0]\" $file_inPR`);
			my $checkRV = int(`grep -c \"$info_rowDB[0]\" $file_tmpRV`);
			if(($checkPR == 0) && ($checkRV == 0)) {
				print MAGENTA "\tMissing infos on Run $info_rowDB[0].\n";
				print "\tPlease edit $file_inPR or validate this run.\n", RESET;
				$check_flag = 1;
			} elsif((($checkRV == 1) && ($checkPR == 0)) || (($checkRV == 1) && ($checkPR == 1))) {
				$num_runs++;
			}
		} elsif($info_rowDB[2] =~ /^calibration$/) {
			$calib_time += $info_rowDB[1];
		}	
	}

	open TEMPRV, $file_tmpRV or die $!;
	my @rowsRV2 = <TEMPRV>;		
	close TEMPRV;
	@rowsRV2 = sort(@rowsRV2);						# Reorder lines according to run number
	my $previous_run = 0;

	foreach my $rowRV(@rowsRV2) {
		my @info_rowRV = split(/ +/,$rowRV);

		# Check on double entries in run_validation_out.txt
		if($info_rowRV[0] == $previous_run) {
			system("rm $file_tmpRV $file_tmpDB");
			print RED "\tError! Run $info_rowRV[0] is listed twice in $file_inRV!\n\n", RESET;
			exit();
		}

		# Check on MOE presence and goodness
		my $MOE_file = "/bxstorage/rootfiles/m4_98-c".$info_rowRV[1]."/Run00".$info_rowRV[0]."_m98c".$info_rowRV[1].".root";
		$MOE_file = "/bxstorage/rootfiles/m4_98-c".$info_rowRV[1]."/Run0".$info_rowRV[0]."_m98c".$info_rowRV[1].".root" if($info_rowRV[0]>9999);
		unless(-s $MOE_file) {
			print MAGENTA "\tWARNING! MOE rootfile of run $info_rowRV[0] does not exist!\n", RESET;
			$check_flag = 1;
			goto ECH_LINE;
		}

  	my $MOEroot_instruction = "MOEroot_instructions.exe";
		open ROOTMOE,"> $MOEroot_instruction" or die $!;
  	print ROOTMOE <<END;
{
TFile *fMOE=TFile::Open("$MOE_file");
if(fMOE->IsZombie()) {
  std::cerr << " File $MOE_file not found\n";
  fMOE->Close();
  return 1;
}
fMOE->cd();
TTree *triggers=(TTree*)fMOE->Get("triggers");
Double_t last_good_event;
triggers->SetBranchAddress("last_good_event", &last_good_event);
triggers->GetEntry(0);
cout << last_good_event << endl;
fMOE->Close(); 
return 0;
}                       
END
 		close ROOTMOE;
  	system("root -l -q -b $MOEroot_instruction > MOE.log");
		my $MOE_events = int(`grep -v "Processing MOEroot_instructions.exe..." MOE.log  | grep -v ">>> c16 validation & stability libraries loaded."`);
  	my $MOE_diff = abs($info_rowRV[2] - $MOE_events);
		print MAGENTA "\tWARNING! Run $info_rowRV[0]: MOE events do not correspond to validated events in $file_inRV.\n", RESET if($MOE_diff>10); 
  	system("rm MOE.log $MOEroot_instruction");

		# Check on ECHIDNA presence and goodness	
		ECH_LINE: $info_rowRV[6] =~ /(.*)root/;
		my $ECH_file = $1."root";

		unless(-s $ECH_file) {
			print MAGENTA "\tWARNING! Echidna rootfile of run $info_rowRV[0] does not exist!\n", RESET;
			$check_flag = 1;
			goto CHECK_EXIT;
		}

  	my $ECHroot_instruction = "ECHroot_instructions.exe";

		open ROOTECH,"> $ECHroot_instruction" or die $!;
  	print ROOTECH <<END;
{
TFile *fECH=TFile::Open("$ECH_file");
if(fECH->IsZombie()) {
  std::cerr << " File $ECH_file not found\n";
  fECH->Close();
  return 1;
}
fECH->cd();
TTree *bxtree=(TTree*)fECH->Get("bxtree");
Int_t ECHentries = bxtree->GetEntries();
fECH->Close(); 
cout << ECHentries << endl;
return 0;
}                       
END
 		close ROOTECH;
  	system("root -l -q -b $ECHroot_instruction > ECH.log");
		my $ECH_events = int(`grep -v "Processing ECHroot_instructions.exe..." ECH.log  | grep -v ">>> c16 validation & stability libraries loaded."`);
  	my $ECH_diff = abs($info_rowRV[2] - $ECH_events);
    print MAGENTA "\tWARNING! Run $info_rowRV[0]: ECHIDNA events do not correspond to validated events in $file_inRV.\n", RESET if($ECH_diff>10);
		system("rm ECH.log $ECHroot_instruction");
		
		CHECK_EXIT: $previous_run = $info_rowRV[0];
	}
	return($check_flag,$num_runs,$calib_time);
}


#	DUTY-CYCLE COMPUTATION
sub DutyCycle() {
	my $num_runs = shift;
	my $calib_time = shift;
	my $firstrun_time = shift;

	print "\n\tNumber of validated runs: $num_runs\n";

	open TEMPRV, $file_tmpRV or die $!;
	my @rowsRV = <TEMPRV>;		
	close TEMPRV;

	my $run_counter = 0;
	my $Live_time = 0;
	my ($FirstDayTime, $LastDayTime);
	@rowsRV = sort(@rowsRV);						# Reorder lines according to run number
	foreach my $rowRV(@rowsRV) {
		my @info_rowRV = split(/ +/,$rowRV);
		#	$info_rowRV[0] --> Run number 	
		#	$info_rowRV[1] --> Echidna cycle
		#	$info_rowRV[2] --> Validated events
		#	$info_rowRV[3] --> Livetime (s)	
		#	$info_rowRV[4] --> Run start date (yyyy-mm-dd)
		#	$info_rowRV[5] --> Run start time (hh:mm:ss)
		#	$info_rowRV[6] --> Echidna file path

		#	COMPUTATION OF FIRST DAY TOTAL TIME
		$firstrun_time =~ /^(.{2}):(.{2}):(.{2})$/;
		$FirstDayTime = 86400 - ($1*60*60 + $2*60 + $3);

		# COMPUTATION OF LAST DAY TOTAL TIME AND LIVE TIME	#
		if($run_counter == ($num_runs-1)) {
			$info_rowRV[5] =~ /^(.{2}):(.{2}):(.{2})$/;
		  $LastDayTime = $1*60*60 + $2*60 + $3 + $info_rowRV[3];
		  $Live_time += $info_rowRV[3];

		#	COMPUTATION OF LIVE-TIME DURING THE OTHER 5 DAY
		} else {
			$Live_time += $info_rowRV[3];
			$run_counter++;
		}
	}
	my $FullDays = (5*86400);
	my $Total_time = $FirstDayTime + $FullDays + $LastDayTime - $calib_time;	
	my $Duty_cycle = sprintf("%2.2f %%",($Live_time/$Total_time)*100);

	print BOLD, GREEN, "\n\t****************************\n";
	print "\t*   DUTY CYCLE:  $Duty_cycle   *\n";
	print "\t****************************\n\n", RESET;
	system("rm $file_tmpRV $file_tmpDB");
}
