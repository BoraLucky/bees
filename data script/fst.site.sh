echo job started at `date`

for i in {AB,BM,Central,HN,ML,NE,QH,TW}
do 
	find /home/00_Apis_cerana/06_SV/00_NGS_mapping/02_mapping/${i} -name *.bam |sort  >${i}.list
	for sample in `cat ${i}.list`
	do 
	echo $sample |cut -d'/'  -f10  |sed 's/.markdup.sorted.bam//g'  >>${i}.poplist
	done
done



#################################fst,site#################################################
for i in {AB,BM,HN,ML,NE,QH,TW}
do
	vcftools --gzvcf 293_hard_filter_population_filter_pass.allsite.vcf.gz \
		 --weir-fst-pop ${i}.poplist \
		 --weir-fst-pop Central.poplist \
		 --out ${i}_Central.per_site
done

for i in {AB,BM,HN,ML,NE,QH,TW}
do
  	Rscript Fst_destribution.R ${i}_Central.10K_win.windowed.weir.fst ${i}.windowed.pi ${i}_Central.per_site.weir.fst ${i}
done

Rscript ks_test.R  
cat *_summary_statistics.txt >all.pop.summary_statistics.txt

####################################parieFst##########################
Rscript FST_matrix_plot.R 293_hard_filter_population_filter_pass.allsite.vcf.gz 293sample_pop.file
echo job end at `date`

