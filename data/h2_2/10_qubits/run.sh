for A in 1 2 3 4 5;
do
    python las-qpe_so.py --an 10 --shots 10240 --log h4_${A} --npy results_${A} &> out_${A}.dat
    echo "$A finished."
done
