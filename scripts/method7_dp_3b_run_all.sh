
n_bp=5
min_len=5

for f in "data/"*; do
    f=${f:(9)}
    f=${f%.*}
    echo $f
    
    python scripts/method7_dp_3b.py -p "data/" -d $f -b $n_bp -k $min_len -l

done