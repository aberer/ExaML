#! /bin/bash


if [ $# != 3  ]; then
    echo "./script <debug> <numProc> <numthread>"
    exit
fi



debugTarget=""
execOutput="./examl-*"
if [ $1 == 1 ]; then
    gdb=" urxvt -e gdb -ex run --args "
    debugTarget=" debug"
    execOutput="./debug-examl-*"
elif [ $1 == 2 ]; then 		# valgrind  
    gdb=" valgrind --tool=memcheck "
    debugTarget=" debug"
    execOutput="./debug-examl-*"    
fi


num=$2
numThread=$3


aln="-s ./testData/small-dna/49.binary"
tree="-t ./testData/small-dna/49.tree"

# exec="./raxmlLight-MPI-SSE3"
exec=$(ls -tr $execOutput | tail -n 1 )
model=GAMMA


args=" $aln -m $model $tree -n tmp "
export PATH="/lhome/labererae/lib/ompi/bin:$PATH"
export LD_LIBRARY_PATH="/lhome/labererae/lib/ompi/lib:$LD_LIBRARY_PATH"

make clean && make $debugTarget 
rm -f  *.tmp ; 

# FT="-am ft-enable-mpi"
FT=""


sleep 1 
mpirun $FT -np $num $gdb $exec $args -T $numThread
