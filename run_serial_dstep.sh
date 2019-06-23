IDIR=/home/gluo/nas/teaching/pdc2019/s1600012760/ch9-1.1/inputs
ODIR=/home/gluo/nas/teaching/pdc2019/s1600012760/ch9-1.1/results

g++ serial_delta_stepping.cpp -o serial_delta_stepping -std=c++11

for TYPE in d t
do
    OFILE=$ODIR/serial-dstep-USA-road-$TYPE.ss.res
    cat /dev/null > $OFILE
    for CITY in NY NE CTR
    do
        GRFILE=$IDIR/USA-road-$TYPE/USA-road-$TYPE.$CITY.gr
        SSFILE=$IDIR/USA-road-$TYPE/USA-road-$TYPE.$CITY.ss
        echo "f $GRFILE $SSFILE" >> $OFILE
        ./serial_delta_stepping $GRFILE $SSFILE $OFILE
    done
done
