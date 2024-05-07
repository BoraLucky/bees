##################################fastp########################################	
mkdir -p 01_293_clean_data/{AB,BM,Central,HN,ML,NE,QH,TW}
for i in {AB,BM,Central,HN,ML,NE,QH,TW}
do
	cat /home/00_Apis_cerana/02_all_data_ln/293_data/${i}.pop.DataPath |xargs -n 2  |while read line
	do 
	array=($line)
	read1=${array[0]}
	read2=${array[1]}
	echo  $read1
	echo $read2
	name=`basename  $read1 |sed 's/_1.fq.gz//g'`	
	echo  $name
	fastp -i $read1 -I $read2 -q 20  -w 15 --detect_adapter_for_pe --trim_poly_g -h 01_293_clean_data/${i}/${name}.report.html -j 01_293_clean_data/${i}/${name}.report.json -o 01_293_clean_data/${i}/${name}_1_clean.fq.gz -O 01_293_clean_data/${i}/${name}_2_clean.fq.gz 2>01_293_clean_data/${i}/${name}.fastp.log 
	echo ${name}.fastp.job.end at `date`	
	done	
done

      ###############################mapping############################################

mkdir -p 02_mapping/{AB,BM,Central,HN,ML,NE,QH,TW}
for i in {AB,BM,Central,HN,ML,NE,QH,TW}
do
	find ./01_kun_paper_293_clean_data/${i} -name *.fq.gz |xargs -n 2  |while read line
	do
	array=($line)
	read1=${array[0]}
	read2=${array[1]}
	echo  $read1
	echo $read2
	name=`basename $read1 |sed 's/_1_clean.fq.gz//g'`	
	echo  $name
	/home/00_biosoft/03_map_align_soft/bwa-mem2_2.0/bwa-mem2  mem -M -t 30  -R "@RG\tID:${name}\tLB:${name}\tPL:ILLUMINA\tSM:${name}"   B1952_final_genome.fasta  $read1 $read2  2>02_mapping/${i}/${name}.bwa-mem2.log | samblaster -M  -r 2>02_mapping/${i}/${name}_markdup.log | sambamba view -S -f bam -t 20 /dev/stdin | sambamba sort -t 20 /dev/stdin -o 02_mapping/${i}/${name}.markdup.sorted.bam
	sambamba  flagstat  -t 20 02_mapping/${i}/${name}.markdup.sorted.bam  >02_mapping/${i}/${name}.markdup.sorted.bam.flagstat
	done
done
