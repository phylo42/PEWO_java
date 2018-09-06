#default values
echo "usage: Rx_to_mutatedRx.sh snp_rate[float] indel_rate[float]"
snp_rate=$1
indel_rate=$2

rm *.backup
mkdir test
for i in `ls R*`
	do
	read_size=`echo $i | cut -d '_' -f4 | tr -d 'r\.fasta'`
	count=`grep -c '>' $i`
	gx=`echo $i | cut -d '_' -f1,2,3 | sed s/^R/G/g `".fasta"
        #gx_size=`awk 'NR==2' ../Gx/$gx | tr -d '\n' | wc -m`
	gx_size=$(grep -v '>' ../Gx/$gx | awk '(NR==1||length<shortest){shortest=length} END {print shortest}' )
	#if read_size > Gx length, then reduce read_size to Gx length 5% (to let space for simulated indels)
	read_size_minor=$read_size
	#take read size + XX% as the limit to take into sizes imposed by simulated indels
	equation="(("$gx_size"+("$gx_size"*0.5))+0.5)/1"
	read_size_major=$(echo $equation | bc)
	#if Gx length< read_size_major
        if [ $gx_size -lt $read_size_major ]
        then
	    #the +0.5)/1  addition is necessary so that bc command rounds the result
	    equation="(("$gx_size"-("$gx_size"*0.5))+0.5)/1"
            echo "$equation"
            read_size_minor=$(echo $equation | bc)
        fi
        #launch mason read simulator
	echo "exp:"$i" gx:"$gx" read_size:"$read_size" gx_size:"$gx_size" read_size_minor:"$read_size_minor" count:"$count
	cp $i $i.backup
	echo "../../mason-0.1.2-Linux-x86_64/bin/mason illumina --forward-only -s 1 -N $count -o ./test/$i.temp.fasta --haplotype-snp-rate $snp_rate --haplotype-indel-rate $indel_rate --read-length $read_size_minor --include-read-information ../Gx/$gx --allow-N-from-genome"
	../../mason-0.1.2-Linux-x86_64/bin/mason illumina --forward-only -s 1 -N $count -o ./test/$i.temp.fasta --haplotype-snp-rate $snp_rate --haplotype-indel-rate $indel_rate --read-length $read_size_minor --include-read-information ../Gx/$gx --allow-N-from-genome
	#need to reuse original read name structure
        sed "s/^>.*contig=\([^ ]\+\) .*length=\([^ ]\+\) .*orig_begin=\([^ ]\+\) .*orig_end=\([^ ]\+\).*$/>\1_r"$read_size"_\3_\4/g" "./test/"$i.temp.fasta > "./test/"$i
        rm "./test/"$i".temp.fasta"
        rm "./test/"$i".temp.fasta.sam"
done
