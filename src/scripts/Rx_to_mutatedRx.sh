rm *.backup
for i in `ls R*`
	do
	read_size=`echo $i | cut -d '_' -f4 | tr -d 'r\.fasta'`
	count=`grep -c '>' $i`
	gx=`echo $i | cut -d '_' -f1,2,3 | sed s/^R/G/g `".fasta"
        #if read_size > Gx length, then reduce read_size to this length
        gx_size=`awk 'NR==2' ../Gx/$gx | tr -d '\n' | wc -m`
        if [ $gx_size -gt $read_size ]
        then
            gx_size=$read_size
        fi
        #launch mason read simulator
	echo "exp:"$i" gx:"$gx" read_size:"$read_size" gx_size:"$gx_size" count:"$count
	cp $i $i.backup
	../../mason-0.1.2-Linux-x86_64/bin/mason illumina -s 1 -N $count -o ./test/$i.temp.fasta --haplotype-snp-rate 0.001 --haplotype-indel-rate 0.001 --read-length $gx_size --include-read-information ../Gx/$gx
	#need to reuse original read name structure
        sed "s/^>.*contig=\([^ ]\+\) .*length=\([^ ]\+\) .*orig_begin=\([^ ]\+\) .*orig_end=\([^ ]\+\).*$/>\1_r"$read_size"_\3_\4/g" ./test/$i.temp.fasta > ./test/$i
        rm ./test/$i.temp.fasta
        rm ./test/$i.temp.fasta.sam
done
