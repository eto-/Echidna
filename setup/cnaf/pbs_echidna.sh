#!/bin/sh
# Predefined variables
interactive=""
mail=""
streams=" "
queue=" "
echidna_options=""
reprocess=0
rewrite_files=""
filename=" "
n_events=" "
#root=`which root | cut -f3 -d\/ | grep ^root_v`

# Internal state
eoo=0

# Some function handlers
function check_opt () {
  case "$1" in
    =*)
      return 1
    ;;
#    -*)
#      return 2
#    ;;
  esac
  return 0
}

function check_opt_order () {
  check_opt $1
  if [ $? -eq 1 -a $eoo -eq 1 ]; then
    echo "script options must precede every echidna options ($1)"
    exit 1
  fi
}

function check_opt_value () {
  check_opt $1
  if [ $? -ne 0 ]; then
    echo "script options =$2 takes an argument ($1)"
    exit 1
  fi
}


# PWD CHECK
echo $PWD | grep -q gpfs
if [ $? -ne 0 ]; then
  echo "pbs_echidna works only in subforlders of gpfs, while you are in $PWD"
  exit 1
fi

# Main options loop
for i in $*; do
  check_opt_order $i			# First check opt order
  if [ -z "$queue" ]; then 		# Second assign script option values (to which opt who need)
    check_opt_value $i "queue"
    queue="-q $i"
  else
    case "$i" in			# Then search the right option
      =interactive)
        interactive="-I"
        streams=""
      ;;
      =mail)
        mail="-N"
      ;;
      =streams)
        streams=""
      ;;
      =queue)
        unset queue
      ;;
      =reprocess)
        reprocess=1
      ;;
      =rewrite_files)
        rewrite_files="-rewrite_files"
      ;;
      =*)
        echo "unknown script option"
        exit 1
      ;;
      *)				# If here this must be an echidna option; chenge state
        eoo=1
        echidna_options="$echidna_options $i"
        [ -z "$filename" ] && filename="-f $i";
        [ -z "$n_events" ] && n_events="-e $i";
        case "$i" in			# search for filename
          -f)				
	    unset filename
          ;;
          -f*)
	    filename=$i
	    ;;
          -e)				
	    unset n_events
          ;;
          -e*)
	    n_events=$i
	    ;;
        esac
      ;;
    esac
  fi
done

# Check options consistency
if [ "$interactive" ]; then
  [ "$mail" ] && echo "warning mailing enabled in interactive mode is deprecated"
fi

# Check the queue
if [ "$queue" = " " -a `whoami` = "production" ]; then
  queue="-q production"
  [ $reprocess -eq 1 ] && queue="-q reprocess";
fi

if [ "$queue" = " " ]; then
  queue="-q borexino_physics"
fi

# parse the filename to show it in qstat and event number
jobname="echidna"
if [ "$filename" -a "$filename" != " " ]; then
  f=`echo $filename | sed "s/.*-f[\ \t]*\(.*\)[\ \t]*.*/\1/"`
  r=`basename $f`
  [ "$r" ] && jobname="${jobname}_$r"
fi
if [ "$n_events" -a "$n_events" != " " ]; then
  n_events=`echo $n_events | sed "s/.*-e[\ \t]*\(.*\)[\ \t]*.*/\1/"`
else
  n_events=0
fi
if [ $reprocess -eq 1 -a $n_events -le 0 ]; then
  echo "$0 error: =reprocess flag requires a definite number of events"
  exit 1
fi

# set standard error and output file names
if [ "$streams" = " " ]; then
  streams="-e ${jobname}.err -o ${jobname}.out"
fi

# Create pbs script file
TMPFILE=`mktemp echidna.XXXXXX` # created in local folder; to be improved ...
if [ $? -ne 0 ]; then
  echo "Can not create temporary script file"
  exit 1
fi

if [ $reprocess -eq 1 ]; then
  retries=5
  FAILEDRUNS=$PWD/failed_runs.txt
  PROCESSEDRUNS=$PWD/processed_runs.txt
  echidna_options="-p bx_writer.print_level info $echidna_options"
else
  FAILEDRUNS=/dev/null
  PROCESSEDRUNS=/dev/null
  retries=1
fi
cat << EOF > $TMPFILE
#!/bin/sh
#set -x
.  /opt/exp_software/borexino/root32/v5-20-00/bin/thisroot.sh 
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD:/opt/exp_software/borexino/postgresql/usr/lib/:/opt/exp_software/borexino/fftw/lib/
hostname
TMPFILE=\`mktemp tempfile.XXXXXX\`
retries=0
fail=1
start_t=\`date +%s\` 
date >&2
while [ \$retries -lt $retries -a \$fail -ne 0 ]; do
  let retries++;
  ./echidna $echidna_options > \$TMPFILE
  date >&2
  if [ $reprocess -eq 1 ]; then
    log_file=\`grep "^Running.*on logfile" \$TMPFILE | cut -f7 -d" "\`
    run_number=\`echo \$log_file | sed -e s/_.*$//\`
    if [ -n "\$log_file" ]; then
      critics=\`grep  "critic" \$log_file |wc -l\`
      grep -q "\(Echidna NORMAL END\)\|\(Echidna ENDED\)" \$log_file
      if [ \$? -eq 0 -a \$critics -eq 0 ]; then
      root_file=\`grep "info: bx_writer: root.file.*opened" \$log_file |cut -f6 -d" "\`
      if [ -n "\$root_file" ]; then
	  check_events=$n_events
	  grep -q "root tree successfully created" \$log_file
	  if [ $? -ne 0 ]; then check_events=0; fi
          root_command=".L rootechidna.so \n .L tools/CheckFiles.C \n cout << \",,,:\" <<  Check(\"\$root_file\", \$check_events) << endl; \n "
	  fail=\`echo -e "\$root_command" | root -b -l | grep ,,,: | cut -f2 -d:\`
	  if [ -z "\$fail" ]; then fail=1; fi
	fi
      fi
    fi
  fi
done
now_t=\`date +%s\`
dt=\`expr \$now_t - \$start_t\`
if [ $reprocess -eq 1 -a \$fail -ne 0 ]; then
  echo -e "\$run_number failed\t\$retries retries\t\$dt seconds\t\$HOSTNAME\t\$now_t time_t" >> $PROCESSEDRUNS
  rm \$root_file
elif [ $reprocess -eq 1 ]; then
  root_file_uploaded=0
  upload_retries=0
  while [ \$root_file_uploaded -eq 0 -a \$upload_retries -lt 5 ]; do
    root_file_uploaded=\`./tools/bx_repository.pl $rewrite_files upload \$root_file | grep "^copied" | wc -l\`
    let upload_retries++; 
    [ \$root_file_uploaded -eq 0 ] && sleep 200;
  done
  log_file_uploaded=0
  upload_retries=0
  while [ \$log_file_uploaded -eq 0 -a \$upload_retries -lt 5 ]; do
    log_file_uploaded=\`./tools/bx_repository.pl $rewrite_files upload \$log_file | grep "^copied" | wc -l\`
    let upload_retries++;
    [ \$log_file_uploaded -eq 0 ] && sleep 200;
  done
  if [ \$root_file_uploaded -eq 1 -a \$log_file_uploaded -eq 1 ]; then
    echo -e "\$run_number succeeded\t\$retries retries\t\$dt seconds\t\$HOSTNAME\t\$now_t time_t" >> $PROCESSEDRUNS
    rm \$root_file
    rm \$log_file
  else
     echo -e "\$run_number upload_failed\t\$retries retries\t\$dt seconds\t\$HOSTNAME\t\$now_t time_t" >> $PROCESSEDRUNS
     mv \$root_file $PWD
  fi
fi
cat \$TMPFILE
rm \$TMPFILE
EOF


#echo script name is $TMPFILE

# Make script file executable
chmod u+x $TMPFILE

# Launch the job
bsub -J $jobname  $queue $interactive $mail $streams ./$TMPFILE
