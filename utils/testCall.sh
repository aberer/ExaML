#! /bin/bash

latest_binary=$(ls -ltr examl-*   | tail -n 1 | sed "s/.*\(examl*\)/\1/")
dir=$(realpath ./testRun/)
id=tmp

rm $dir/*.$id ;  mpirun -np 2 $latest_binary -s testData/49.binary -t testData/49.tree -m PSR   -n $id  -w $dir
