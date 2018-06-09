#!/bin/sh

sed -e 1,/--CUT1--/\!d -e/--CUT1--/d mach4/Mach4Event.hh-scheme > event/Mach4Event.hh
awk '{ gsub (/[.]/, "__", $1); print "    double get_" $1 " () const { return " $1 "; }" }' mach4/M4_variables.txt >> event/Mach4Event.hh
sed -e /--CUT1--/,/--CUT2--/\!d -e/--CUT[12]--/d mach4/Mach4Event.hh-scheme >> event/Mach4Event.hh
sed -e s/^/\ \ \ \ double\ / -e s/[.]/__/g -e s/\ *$/\;/ mach4/M4_variables.txt >> event/Mach4Event.hh
sed -e /--CUT2--/,$\!d -e/--CUT2--/d mach4/Mach4Event.hh-scheme >> event/Mach4Event.hh


sed -e 1,/--CUT1--/\!d -e/--CUT1--/d mach4/Mach4Event.cc-scheme > event/Mach4Event.cc
#sed -e s/[.]/__/g mach4/M4_variables.txt | awk 'BEGIN { l = 0; c = 0; } END { if (c != 0) { l = l " = " } print l "-99;" } { if (c == 0) { l = "  " $1 ; c++ } else { l = l " = " $1; c++ } if (c == 5) { print l " = \\"; c = 0; } };' >> event/Mach4Event.cc
#awk 'BEGIN { l = 0; c = 0; } END { print l " \\" } { gsub (/[.]/, "__", $1); if (c == 0) { l = $1 "(-99)" ; c++ } else { l = l ", " $1 "(-99)"; c++ } if (c == 5) { print l ", \\"; c = 0; } };' mach4/M4_variables.txt >> event/Mach4Event.cc
awk 'BEGIN { l = 0; c = 0; bof=""} END { print l " \\" } { gsub (/[.]/, "__", $1); if (c == 0) { l = $1 "(-99)" ; c++ } else { l = l ", " $1 "(-99)"; c++ } if (c == 5) { printf("%s%s", bof, l); c = 0; l=""; bof=",\\\n" } };' mach4/M4_variables.txt >> event/Mach4Event.cc 
sed -e /--CUT1--/,/--CUT2--/\!d -e/--CUT[12]--/d mach4/Mach4Event.cc-scheme >> event/Mach4Event.cc
awk '{ s1=$1; gsub (/[.]/, "__", $1); print "  tree->SetBranchAddress (\"" s1 "\", &" $1 ");" }' mach4/M4_variables.txt | grep -v laben.muontrack.spd >> event/Mach4Event.cc
sed -e /--CUT2--/,$\!d -e/--CUT[12]--/d mach4/Mach4Event.cc-scheme >> event/Mach4Event.cc

