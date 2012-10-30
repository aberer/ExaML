#! /bin/bash

latest_binary=$(ls -ltr examl-*   | tail -n 1 | sed "s/.*\(examl*\)/\1/")
echo $latest_binary
if [ "$?" -ne "0" ]; then
    echo "could not find executable"
    exit -1
fi

dir=$(realpath ./testRun/)
id=tmp


cmd="mpirun -np 2 $latest_binary -T 2  -s testData/49.binary -t testData/49.tree -m PSR   -n $id  -w $dir"
echo $cmd
rm $dir/*.$id ;  $cmd
