#!/bin/sh

echo "`basename $0` $1 at `date`"
(dirname $0 | grep -q tools\$) && cd `dirname $0`/../

function do_stop () {
  killall omon_read 2>/dev/null
  ipcrm shm `ipcs | grep 4e1f | cut -f2 -d" "` 2>/dev/null >/dev/null
  [ -x check_online_echidna ] && killall -9 check_online_echidna
  [ -x messagebox ] && killall -9 messagebox
  sleep 1
  pidof omon_read echidna
  if [ $? -eq 0 ]; then
    sleep 5
    killall -9 -q omon_read echidna
    echo "KILLED"
  fi
}

case "$1" in
  start)
    do_stop
    if [ -n "$2" ]; then run=$2;
    else run=0000;
    fi
    [ -x messagebox ] && ./messagebox &
    [ -x check_online_echidna ] && ./check_online_echidna &
    sleep 1
    cycle=`readlink -f rootechidna.so |sed 's/.*so.cycle_\([0-9]*\)_.*/_c\1/'`
    ./tools/omon_read | ./echidna -l /online/Run00${run}_online${cycle}.log -f /dev/stdin -vV online >/bxwww/data/omon.log 2>/tmp/echidna.err &
  ;;
  stop)
    do_stop
  ;;
  *)
    echo "usage `basename $0` {start|stop}"
    exit 1
  ;;
esac
