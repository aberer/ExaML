#! /bin/bash


if [ $# != 3  ]; then
    echo "./script <debug> <numProc> <numThread>"
    exit
fi

debugTarget=""
execOutput="./examl-*"
if [ $1 == 1 ]; then
    gdb=" urxvt -e  gdb -ex run  --args "
    debugTarget=" debug"
    execOutput="./debug-examl-*"
fi

num=$2
numThread=$3


aln="-s ./testData/49.binary"
tree="-t ./testData/49.tree"

# exec="./raxmlLight-MPI-SSE3"
exec=$(ls -tr $execOutput | tail -n 1 )
model=PSR


args=" $aln -m $model $tree -n tmp "
export PATH="/lhome/labererae/lib/ompi/bin:$PATH"
export LD_LIBRARY_PATH="/lhome/labererae/lib/ompi/lib:$LD_LIBRARY_PATH"

make clean && make $debugTarget 
rm -f  *.tmp ; 

#FT= -am ft-enable-mpi
FT=""
mpirun $FT -np $num $gdb $exec $args -T $numThread
