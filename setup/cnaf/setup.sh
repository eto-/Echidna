# check that we are under setup/cnaf folder

HERE=`basename $PWD`
DIR=`dirname $PWD`
UPDIR=`basename $DIR`

if [ "$HERE" != "cnaf" ] || [ "$UPDIR" != "setup" ]
then
	echo "Error: this script only works under setup/cnaf";
	exit 1
fi

# check that Makefiles exist

SETUP_FILES_PRESENT=1;
for SETUP_FILE in Makefile_echidna
do
	if [ ! -f $SETUP_FILE ]
	then
		SETUP_FILES_PRESENT=0;
		echo -e "Error: $SETUP_FILE not found";
	fi
done

if [ $SETUP_FILES_PRESENT -ne 1 ]
then
	exit 2
fi

# copying files
echo "Placing the CNAF massive submitter in Echidna/tools/"
cp massive_submitter_CNAF.pl ../../tools/

echo "Replacing Echidna Makefile"
cp Makefile_echidna ../../Makefile

echo "Replacing pbs_echidna.sh"
cp pbs_echidna.sh ../../

# replacing bxdb with bxdb.lngs.infn.it in echidna.cfg
echo "Replacing bxdb with bxdb.lngs.infn.it in echidna.cfg"
sed -i 's/SERVER1 bxdb[[:blank:]]/SERVER1 bxdb.lngs.infn.it /g' ../../echidna.cfg

# Adjusting bx_reader method in echidna.cfg
echo "Adjusting bx_reader method in echidna.cfg"
sed -i s?'://bxmaster-data.lngs.infn.it/bxstorage'?':///storage/gpfs_data/borexino'?g ../../echidna.cfg
sed -i s?' wget -nv --timeout=60 -O - '?' ./tools/reader_CNAF.sh '?g ../../echidna.cfg

# root setup reminder 
echo "Do not forget to setup root before running make. Use the command:"
echo ". /opt/exp_software/borexino/root32/v5-34-24/bin/thisroot.sh"
echo "You can now compile and use Echidna."
