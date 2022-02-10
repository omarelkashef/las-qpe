for A in 1024 2048 4096 6144 8192 10240;
do
    python las-qpe_so.py --an 10 --shots ${A} &> out_${A}_10q.dat
    echo "$A finished."
done
