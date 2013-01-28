#! /bin/bash


numProc=2
numThreads=2

# :todo: make a parameter  
theDataset="small-dna"

exec=$( ls -rt $(dirname $0)/../examl-* | tail -n 1  )
aln="/lhome/labererae/proj/my-ExaML/testData/small-dna/49.binary"
tree="/lhome/labererae/proj/my-ExaML/testData/small-dna/49.tree"

for model in  PSR GAMMA 
do 
    for multiSched in  0 1 
    do 
	for branch in 0 1
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
	    call="mpirun -np $numProc $exec -s $aln -t $tree $args -n $name -T $numThreads "
	    $call > /dev/null
	    
	    # verify 	    
	    relevantFile=$(dirname $0)/../testRun/$theDataset/ExaML_info.$name 
	    if [ ! -f $relevantFile ]; then
		echo "could not find file $relevantFile"
		continue  
	    fi
	    lnlCorrect=$(grep "tree:"  $relevantFile |  cut -f 2 -d ':' | tr -d ' ') 
	    lnlHere=$(grep "tree:"  ExaML_info.$name |  cut -f 2 -d ':' | tr -d ' ')

	    if [ "$lnlCorrect" == "$lnlHere" ]; then
		echo -e  "SUCCESS: $name\t$lnlCorrect"
		rm ExaML_* 
	    else
		echo -e  "FAILURE: $name\t$lnlCorrect\t!=\t$lnlHere"
		echo -e "call was: $call"
	    fi
	    
	done 
    done 
done 

