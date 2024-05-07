echo job started at `date`
##################################population file ##################################
for i in {AB,BM,Central,HN,ML,NE,QH,TW}
do 
	find /home/00_Apis_cerana/06_SV/00_NGS_mapping/02_mapping/${i} -name *.bam |sort  >${i}.list
	for sample in `cat ${i}.list`
	do 
	echo $sample |cut -d'/'  -f9-  |sed 's/.markdup.sorted.bam//g' |awk -F '/' '{print $2"\t"$1}' >>${i}.poulation.file
	done
done

cat {AB,BM,HN,ML,NE,QH,TW,Central}.poulation.file >293sample_pop.file  && rm *.poulation.file *.list



###############################10K####################################################
for i in `grep '>' B1952_final_genome.fasta |sed 's/>//g'`
do
	pixy --stats pi fst dxy  \
             --vcf 293_hard_filter_population_filter_pass.allsite.vcf.gz \
             --populations 293sample_pop.file \
             --window_size 10000 \
             --n_cores 5 \
             --chromosomes $i \
	     --output_folder pixy_10K_window \
	     --output_prefix ${i}_pixy	
done



echo job started at `date`

