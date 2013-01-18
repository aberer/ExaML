#! /bin/bash


if [ $# != 2  ]; then
    echo "./script <debug> <numProc>"
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



aln="-s ./testData/49.binary"
tree="-t ./testData/49.tree"

# exec="./raxmlLight-MPI-SSE3"
exec=$(ls -tr $execOutput | tail -n 1 )
model=PSR


args=" $aln -m $model $tree -n tmp -Q "
export PATH="/lhome/labererae/lib/ompi/bin:$PATH"
export LD_LIBRARY_PATH="/lhome/labererae/lib/ompi/lib:$LD_LIBRARY_PATH"

make clean && make $debugTarget 
rm -f  *.tmp ; 

mpirun -am ft-enable-mpi -np $num $gdb $exec $args
