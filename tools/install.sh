#!/bin/sh

if [ -z "$1" ]; then
  echo "usage: `basename $0` directory (/borexbin/Echidna/ or /home/production/Echidna)"
  echo "  echidna is searched in $PWD, $PWD/../ and `dirname $0`/../"
  exit 0
fi

if [ ! -d $1 ]; then
  echo "error: $1 is not a directory"
  exit 1
fi

target=$1
source=$PWD
if [ ! -f main.cc ]; then
  pwd | grep -q tools\$
  if [ $? -eq 0 -a -f "../main.cc" ]; then 
    source="$PWD/../"
  elif [ -f `dirname $0`/../main.cc ]; then
    target=$PWD
    source="`dirname $0`/../"
  else
    echo "error: echidna source not found"
    exit 1
  fi
fi

cd $source
cp -dp echidna pbs_echidna.sh libechidna* rootechidna.* echidna.cfg user.cfg $target

mkdir -p $target/particleid/data
cp -dp particleid/data/[^C][^V][^S]* $target/particleid/data/

[ -e $target/tools ] || ln -sf ./ $target/tools
cp -dp tools/RunInfo.pl tools/bx_repository.pl tools/bx_process.pl tools/CheckFiles.C tools/check_storage.sh tools/gen_bx_process_cmd.pl tools/reprocess_submitter.pl tools/process_source_runs.sh tools/bx_production.pl tools/Shift_analyzer.pl $target

mkdir -p $target/event
cp -dp event/BxEvent.hh event/Mach4Event.hh $target/event


if [ -x tools/omon_read ]; then
  cp -dp tools/omon_read $target
  cp -dp tools/online_echidna.sh $target
elif `whoami | grep -q daqman && hostname | grep -q bxbuildi`; then
  echo "WARNING online stuff not installed"
fi

if [ -x tools/fix_db ]; then
  cp -dp tools/fix_db $target
elif `whoami | grep -q production`; then
  echo "WARNING fix_db not installed"
fi

echo "installed from $source to $target"

