echo -e "\n*** THIS SCRIPT MUST RUN ON IA32_SL4 MACHINES ***\n"

export ORI_LD=$LD_LIBRARY_PATH

export PINST=/afs/in2p3.fr/group/nusol/Collab/setup/files/
export P432=/afs/in2p3.fr/system/ia32_sl4/
export P564=/afs/in2p3.fr/system/amd64_sl5/

unset ROOT_ENV_SET
export LD_LIBRARY_PATH=$ORI_LD:$P432/usr/local/pgsql/8.1/lib/:$P432/usr/local/lib/
export PGLIB=$P432/usr/local/pgsql/8.1/lib/
export ROOTSYS=$P432/usr/local/root/root_v5.20.00/root/
source $P432/usr/local/shared/bin/root_env.sh
export NOXML=1

