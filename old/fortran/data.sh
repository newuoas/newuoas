#!/bin/bash

SUFFIX=.F1
rm -f PROB
for i in `ls *$SUFFIX`
do
    echo `basename $i $SUFFIX`>>PROB
done

SOLVERNUM=`cat SOLVERNUM`
for (( SOLVER=1; SOLVER <= $SOLVERNUM; SOLVER++ ))
do
    for STAT in N F #FOPTREC  
    do
        echo "Processing data for $STAT$SOLVER..."
        rm -f $STAT$SOLVER
    
        for PROBLEM in `cat PROB`
        do
            i=$PROBLEM.$STAT$SOLVER
            sed -i 's/^[ \t]*//' $i
            sed -i 's/[ \t]*$//' $i
            cat $i>>$STAT$SOLVER
        done
    done
done

for STAT in DIM
do
    echo "Processing data for $STAT..."
    rm -f $STAT

    for PROBLEM in `cat PROB`
    do
        i=$PROBLEM.$STAT
        sed -i 's/^[ \t]*//' $i
        sed -i 's/[ \t]*$//' $i
        cat $i>>$STAT
    done
done

rm -f FSTART
echo "Processing data for FSTART..."
for PROBLEM in `cat PROB`
do 
    cut -d" " -f1 < $PROBLEM.FOPTREC1 >> FSTART
done 


matlab -nodesktop -nosplash -r "cd $PWD; iternum; priternum; exit" 
