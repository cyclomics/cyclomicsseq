for i in $(ls cycas_results); do
    medaka smolecule --qualities --depth $1 --model $2 --chunk_len $3 --chunk_ovlp $4 --method $5 --length $6 --threads $7 --batch_size $8 /medaka_results/$i/ cycas_results/$i
    if [ $? -eq 0 ] 
    then 
        echo "succesfully processed $i"
    else 
        echo "Could not call consensus on $i" >&2 
    fi
done
exit 0
