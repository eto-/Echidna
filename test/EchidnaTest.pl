#!/usr/bin/perl -w
#
# Test Script for Echidna CVS code
# $Id: EchidnaTest.pl,v 1.3 2004/07/16 12:55:44 pallas Exp $
# Purpose: 
#    check out latest version or specific branch,
#    compile it, link it, run on several test files,
#    send automatic error report to developers
# Author: M. Pallavicini
# Created: 15/7/04
# History:
#      15/7/4: first version, compile and link only
#

use POSIX qw(strftime);

$domail=1;
$workdir="/tmp/";

foreach $token (@ARGV) {
   
   #flag to select specific branch instead of CVS top
   $i = index($token,"--branch=");
   $j = index($token,"-b=");
   if ($i>=0) {
     $branch = substr($token,10);
     next;
   }
   if ($j>=0) {
	 $branch = substr($token,3);
     next;
   }

   #flag to disable mail message
   $i = index($token,"--nomail");
   $j = index($token,"-n");
   if ($i>=0 || $j>=0) {
     $domail=0;
     next;
   }
   
   #flag to select work dir other than /tmp/
   $i = index($token,"--directory=");
   $j = index($token,"-d=");
   if ($i>=0) {
     $workdir = substr($token,12);
     next;
   }
   if ($j>=0) {
	 $workdir = substr($token,3);
     next;
   }

   #help 
   $i = index($token,"--help");
   $j = index($token,"-h");
   if ($i>=0 || $j>=0) {
     print "\n*** Utility for Echidna Tests\n";
	 print "*** Usage: EchidnaTest.pl [options]\n";
	 print "*** Valid options:\n";
	 print "      --help, -h : dump this message and exit\n";
	 print "      --branch=, -b=: select specific branch instead of taking CVS top\n";
	 print "      --nomail, -n: do not send mail message to developers\n";
	 print "      --directory=, -d=: specify a working directory other than /tmp/\n"; 
	 print "\n\n";
     exit;
   }

   printf("Unknown option $token. Ignored\n");
}

@packages = ("event", "framework", "interface", "root", "modules");

#hash table with correspondance between package and developer
%package_responsible = ( 
		"event" => "Davide",
		"framework" => "Alessandro",
		"interface" => "Alessandro",
		"root" => "Davide",
		"modules" => "Daniela"
		);

#mail address of developers
%mail_address = ( 
		"Davide" => "ddangelo\@lngs.infn.it",
		"Alessandro" => "razeto\@lngs.infn.it",
		"Daniela" => "dmanuzio\@ge.infn.it"
		);

#library name of packages
%libfile = (
		"event" => "libevent.a",
		"framework" => "libframework.a",
		"interface" => "libinterface.a",
		"root" => "libbxroot.a",
		"modules" => "libmodules.a"
		);

#mail address of testers group
$test_group = "pallas\@ge.infn.it";

#create a work directory
$tm = time();
if (!(-d $workdir)) {
	die "Directory $workdir does not exist!\n";
}
$workdir = "$workdir/EchidnaTestSession_".$tm;
mkdir $workdir;
print "Directory $workdir created\n";

#main 
&CheckOutCode($branch);
&CompileCode();
&CheckCompilation();
&LinkCode();
&CheckLink();

print "Echidna Test was completed.\n";

exit;
###################
###################

## Subs

# check out code from CVS
# if $branch is not 0, use it, otherwise take CVS top
#
sub CheckOutCode() {
    my $b=$_[0];

	print "Code CheckOut starting\n";

    # check CVSROOT
	$cvs = $ENV{"CVSROOT"};
	if (!$cvs) {
		die "Error: You have to set CVSROOT environment variable";
	}	
	chdir $workdir;
	if ( !$b ) {
		print "Check out top CVS repository\n";
        `cvs co Echidna`;
	} 
	else {
		print "Check out branch ".$branch."\n";
        `cvs co -j$branch Echidna`;
	}
}

#
# compile code: compile all packages and save messages in files
#
sub CompileCode() {
	chdir "$workdir/Echidna";
	foreach $pkg (@packages) {
		print "Start compiling $pkg ...";
		`make $pkg >& $workdir/$pkg.log`;
		print " ... DONE!\n";
	}
}

#
# link code
#
sub LinkCode() {
	chdir "$workdir/Echidna";
	print "Start linking echidna...";
    `make >& $workdir/link.log`;
	print " ... DONE!\n";
}


#
# check if compilation was fine and look for errors
#
sub CheckCompilation() {
	$go_on=1;
	chdir "$workdir/Echidna";
	foreach $pkg (@packages) {
		$libf = $libfile{$pkg};
		$pkgdir = "$workdir/Echidna/$pkg";
		if (!( -f "$pkgdir/$libf" )) {
			print "Error occurred while compiling package $pkg\n";
			print "Sending report to $package_responsible{$pkg} at $mail_address{$package_responsible{$pkg}}\n";
			open (MSG,">$workdir/msg.txt");
			print MSG "Dear $package_responsible{$pkg}, \n";
			print MSG "this message was generated by Echidna Automatic Test Utilities.\n";
			print MSG "Check the report below and correct the errors as soon as possible.\n";
			$now = strftime "%a %b %e %H:%M:%S %Y", localtime;
			print MSG "Test date: $now\n";
			print MSG "\n\nMakefile output for Package: $pkg\n\n\n";
			close MSG;
			`cat $workdir/$pkg.log >> $workdir/msg.txt`;
			if ( $domail ) {
			    `cat $workdir/msg.txt | mail -n -N -s "Echidna: ERROR REPORT for package $pkg" $mail_address{$package_responsible{$pkg}}`; 
			}
			$go_on=0;
		}
	}
	if (!$go_on) {
		print "Errors occurred in compilation. Test finished.\n";
		exit;
	}
}

#
# check link
#
sub CheckLink() {
	chdir "$workdir/Echidna";
	if ( -f "./echidna" ) {
		print "Compilation and Linking were successfull!\n";
		return;
	} else {
		print "Missing executable. Something went wrong.\n";
        open (MSG,">$workdir/link.txt");
		print MSG "Dear package responsibles, \n";
		print MSG "this message was generated by Echidna Automatic Test Utilities.\n";
		print MSG "A linking problem was found. The makefile output is given below.\n";
		$now = strftime "%a %b %e %H:%M:%S %Y", localtime;
		print MSG "Test date: $now\n";
		print MSG "\n\nMakefile output for Linking\n\n\n";
		close MSG;
		`cat $workdir/link.log >> $workdir/link.txt`;
		if ( $domail ) {
			print "Sending message to dmanuzio\@ge.infn.it\n";
		    `cat $workdir/link.txt | mail -n -N -s "Echidna: ERROR REPORT for LINKING" dmanuzio\@ge.infn.it`; 
		}
	}
}
