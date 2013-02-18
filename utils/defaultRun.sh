#! /bin/bash

model=" -m GAMMA "
# FT="-am ft-enable-mpi"
FT=""

exec=./exaML

if [ $# != 4  ]; then
    echo "./script default|debug <numProc> <numthread> <dataset>"
    exit
fi

default=$1
num=$2
numThread=$3
dataset=$4

configArgs="--disable-productive" # 
if [ "$default" == "default" ]; then    
    gdb=""
elif [ "$default" == "debug" ]; then 
    configArgs="$configArgs --enable-debug"
    gdb="urxvt -e gdb -ex run --args "
else  
    echo "first argument must be either 'debug' or 'default'"
    exit 
fi

if [ ! -d testData/$dataset ]; then
    echo "could not find dataset testData/$dataset"
    exit
fi


aln="-s ./testData/$dataset/aln"
tree="-t ./testData/$dataset/tree"


execArgs=" $aln $model $tree -n tmp  "
export PATH="/lhome/labererae/lib/ompi/bin:$PATH"
export LD_LIBRARY_PATH="/lhome/labererae/lib/ompi/lib:$LD_LIBRARY_PATH"


status="$(./config.status --config | tr -d "'" )"

if  [ "$(echo $status)"  == "$(echo $configArgs)" ]; then 
    echo "no need to re-configure / re-build"
else 
    echo "calling ./configure $configArgs" 
    ./configure $configArgs  
    make clean 
fi 

rm -f exaML
make -j 4 

basecall="mpirun $FT -np $num $gdb $exec $execArgs -T $numThread"
if [ -f ./exaML ]; then 
    echo "calling exaML as $basecall"
    wait 
    $basecall    
fi 
