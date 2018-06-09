#!/bin/sh

# Predefined variables
qwait="-W"
interactive=""
mail=""
streams="-k oe"
#queue="-q production"
options=""

if test $# = 0 -o "$1" = "-h" -o "$1" = "--help" ; then
  echo "Sample use: ./Vshape.sh year month week, e.g. ./Vshape.sh 2011 Jan 16 "
  exit 0
fi

my_arg0=`basename "$1"`
my_arg1=`basename "$2"`
my_arg2=`basename "$3"`

until shift 1 ; do 
  options="$options $1"
done

if [ ! -e "./VShape_DAQ" ]; then
  echo "ERROR: file "VShape_DAQ" not found in $PWD"
  exit 1
fi

if [ ! -x "./VShape_DAQ" ]; then
  echo "ERROR: file $PWD/VShape_DAQ is not executable."
  exit 2
fi

if [ ! -e "/home/runcoord/VesselShapes/fit_rt" ]; then
  echo "ERROR: file "fit_rt" not found in $PWD"
  exit 1
fi

if [ ! -e "/home/runcoord/VesselShapes/dst_list.lst" ]; then
  echo "ERROR: file "dst_list.lst" not found in $PWD"
  exit 1
fi

if [ ! -x "/home/runcoord/VesselShapes/fit_rt" ]; then
  echo "ERROR: file $PWD/fit_rt is not executable."
  exit 2
fi

if [ ! -e "/home/runcoord/VesselShapes/fit_vol_ellip" ]; then
  echo "ERROR: file "fit_rt" not found in $PWD"
  exit 1
fi

if [ ! -e "/home/runcoord/VesselShapes/dst_list.lst" ]; then
  echo "ERROR: file "dst_list.lst" not found in $PWD"
  exit 1
fi

if [ ! -x "/home/runcoord/VesselShapes/fit_vol_ellip" ]; then
  echo "ERROR: file $PWD/fit_rt is not executable."
  exit 2
fi

# Create pbs script file
TMPFILE=`mktemp /tmp/VShape_DAQ.XXXX`
if [ $? -ne 0 ]; then
  echo "Could not create temporary script file"
  exit 1
fi
cat << EOF > $TMPFILE
#!/bin/sh
. /home/cluster_environment.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/lib
cd $PWD
./VShape_DAQ $my_arg0 $my_arg1 $my_arg2
EOF

# Launch the first job
#QSTATI=`eval qsub $queue $interactive $mail $streams $TMPFILE`
QSTATI=`eval qsub $queue $interactive $mail $TMPFILE`
echo "$QSTATI"

WAIT="depend=afterok:"
FIRST=$WAIT$QSTATI

# Create pbs script file
TMPFILEI=`mktemp /tmp/fit_rt.XXXXXX`
if [ $? -ne 0 ]; then
  echo "Could not create temporary script file"
  exit 1
fi
cat << EOF > $TMPFILEI
#!/bin/sh
. /home/cluster_environment.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/lib
cd $PWD
/home/runcoord/VesselShapes/fit_rt
EOF

# Launch the second job on hold
#QSTATII=`eval qsub $qwait $FIRST $queue $streams $TMPFILEI`
QSTATII=`eval qsub $qwait $FIRST $queue $TMPFILEI`
echo "$QSTATII"

SECOND=$WAIT$QSTATII

# Create pbs script file
TMPFILEII=`mktemp /tmp/fit_vol_ellip.XXXXXX`
if [ $? -ne 0 ]; then
  echo "Could not create temporary script file"
  exit 1
fi
cat << EOF > $TMPFILEII
#!/bin/sh
. /home/cluster_environment.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/lib
cd $PWD
/home/runcoord/VesselShapes/fit_vol_ellip
EOF

# Launch the third job on hold w.r. to the second
#QSTATIII=`eval qsub $qwait $SECOND $queue $streams $TMPFILEII`
QSTATIII=`eval qsub $qwait $SECOND $queue $TMPFILEII`
echo "$QSTATIII"


