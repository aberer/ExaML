#! /bin/bash


numThread=2

if [ $# != 2  ]; then
    echo "./script <debug> <numProc>"
    exit
fi

if [ $1 == 1 ]; then
    gdb="urxvt -e  gdb --args"
fi

numThread=$2

# aln="-s ../testdata/small.dna.binary"
# tree="-t ../testdata/small.startingTree.dna.tree"

aln="-s ../testdata/medium.dna.binary"
tree="-t ../testdata/medium.startingTree.dna.tree"

exec="./raxmlLight-MPI-SSE3"
model=GAMMA


make clean && make 

rm -f  *.tmp ; 
mpirun -np $numThread $gdb $exec $aln -m $model $tree -n tmp




