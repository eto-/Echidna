#!/bin/sh

ALLFILES=storage_runs_all.txt
GOODFILES=storage_runs_good.txt
MISSINGFILES=missing_runs.txt
cat /dev/null > $GOODFILES

if [ ! -L rootechidna.so ]; then
  echo "Launch this script in the echidna folder"
  exit
fi

cycle=`readlink -f rootechidna.so |sed 's/.*so.\(cycle_[0-9]*\)_.*/\1/'`
find /bxstorage/rootfiles/${cycle}/ -name \*.root |grep -v ancillary > $ALLFILES;
if [ "$1" ]; then
  TMPFILE=`mktemp /tmp/XXXXXX`
  for i in `cat $1`; do
    grep $i $ALLFILES >> $TMPFILE
  done
  mv $TMPFILE $ALLFILES
fi

root=`which root | cut -f3 -d\/ | grep ^root_v`
cat << EOF | qsub -d $PWD
if [ -n "$root" ]; then choose_root $root; fi
for i in \`cat $ALLFILES\`; do
  run=\`echo \$i | sed "s/.*Run0*\([0-9]*\)_.*root/\1/"\`;
  events=\`psql -h bxdb.lngs.infn.it -At -c 'SELECT "LastValidEvent" FROM "ValidRuns" WHERE "RunNumber"='\$run';' bx_runvalidation borex_guest\`
  fail=\`echo -e '.L rootechidna.so \n .L tools/CheckFiles.C \n cout << ",,,:" <<  Check("http://bxmaster-data'\$i'", '\$events') << endl;\n' |root -b -l 2>/dev/null |grep ,,,: | cut -f2 -d:\`; 
  if [ -z "\$fail" ]; then 
    echo \$i failed
  elif [ "\$fail" == "0" ]; then
    echo \$run >> $GOODFILES
  else
    echo \$i failed
  fi
done

if [ -n "$1" ]; then
  blm $1 -$GOODFILES  > $MISSINGFILES
else
  TMPFILE=\`mktemp /tmp/valid_runs.XXXXXX\`
  psql -h bxdb.lngs.infn.it -At -c 'SELECT "RunNumber" FROM "ValidRuns"' bx_runvalidation borex_guest > \$TMPFILE
  blm \$TMPFILE -$GOODFILES  > $MISSINGFILES
  rm \$TMPFILE
fi
EOF
