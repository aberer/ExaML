#! /bin/bash

exec="/lhome/labererae/proj/ExaML/examl/examl"
aln="/lhome/labererae/proj/my-ExaML/testData/medium-dna/aln.binary"
tree="/lhome/labererae/proj/my-ExaML/testData/medium-dna/tree"
for model in  GAMMA PSR 
do 
    for multiSched in  1 0 
    do 
	for branch in 1  0
	do 
	    args="-m $model"
	    name="$model"
	    if [ "$multiSched" == "1" ]; then
		name="$name-Q"
		args="$args -Q"
	    else
		name="$name-noQ"
	    fi

	    if [ "$branch" == "1" ]; then
		name="${name}-M"
		args="$args -M"
	    else
		name="${name}-noM"
	    fi
	    
	    mpirun -np 3 $exec -s $aln -t $tree $args -n $name 
	done  	
    done 
done 
