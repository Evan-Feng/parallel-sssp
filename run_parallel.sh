args=$@

IDIR=/home/gluo/nas/teaching/pdc2019/s1600012760/ch9-1.1/inputs
ODIR=/home/gluo/nas/teaching/pdc2019/s1600012760/ch9-1.1/results

g++ parallel.cpp -O3 -o parallel -std=c++11 -fopenmp

for TYPE in d t
do
    OFILE=$ODIR/parallel-USA-road-$TYPE.ss.res
    cat /dev/null > $OFILE

    CITY=NY
    DELTA=20000
    GRFILE=$IDIR/USA-road-$TYPE/USA-road-$TYPE.$CITY.gr
    SSFILE=$IDIR/USA-road-$TYPE/USA-road-$TYPE.$CITY.ss
    echo "f $GRFILE $SSFILE" >> $OFILE
    ./parallel $GRFILE $SSFILE $OFILE $DELTA

    CITY=NE
    DELTA=20000
    GRFILE=$IDIR/USA-road-$TYPE/USA-road-$TYPE.$CITY.gr
    SSFILE=$IDIR/USA-road-$TYPE/USA-road-$TYPE.$CITY.ss
    echo "f $GRFILE $SSFILE" >> $OFILE
    ./parallel $GRFILE $SSFILE $OFILE $DELTA

    CITY=CTR
    DELTA=20000
    GRFILE=$IDIR/USA-road-$TYPE/USA-road-$TYPE.$CITY.gr
    SSFILE=$IDIR/USA-road-$TYPE/USA-road-$TYPE.$CITY.ss
    echo "f $GRFILE $SSFILE" >> $OFILE
    ./parallel $GRFILE $SSFILE $OFILE $DELTA
done
