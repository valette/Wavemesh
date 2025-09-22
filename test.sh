
for f in armadillo fandisk cow horse rocker-arm bimba homer; do
    file=$f.obj
    wget https://github.com/alecjacobson/common-3d-test-models/raw/master/data/$file
    echo $file
    for r in {-1..2}; do
        echo "$r"
        for a in {0..1}; do
            echo "****************** $f $a $r*******************"
            bin/wavemesh c $file -a $a -r $r
            cat repport.txt
            bin/wavemesh d out.ddd
            rm out.ddd
        done
    done
done
