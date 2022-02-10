for A in 1 2 4 6 8 10;
do
    python las-qpe_so.py --an ${A} &> out_${A}_1k.dat
    python las-qpe_so.py --an ${A} --shots 10240 &> out_${A}_10k.dat
    echo "$A finished."
done
