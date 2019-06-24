args=$@

IDIR=/home/gluo/nas/teaching/pdc2019/s1600012760/ch9-1.1/inputs
ODIR=/home/gluo/nas/teaching/pdc2019/s1600012760/ch9-1.1/results

g++ delta_step_openmp.cpp -o delta_step_openmp -std=c++11 -fopenmp

for TYPE in d
do
    OFILE=$ODIR/dstep-openmp-s-USA-road-$TYPE.ss.res
    cat /dev/null > $OFILE
    for CITY in NY
    do
        GRFILE=$IDIR/USA-road-$TYPE/USA-road-$TYPE.$CITY.gr
        SSFILE=$IDIR/USA-road-$TYPE/USA-road-$TYPE.$CITY.ss
        echo "f $GRFILE $SSFILE" >> $OFILE
        ./delta_step_openmp $GRFILE $SSFILE $OFILE 
    done
done
