echo job started at `date`

for i in {AB,BM,Central,HN,NE,QH,TW}
do 
	find /home/00_Apis_cerana/00_NGS_mapping/02_mapping/${i} -name *.bam |sort  >${i}.list
	for sample in `cat ${i}.list`
	do 
	echo $sample |cut -d'/'  -f10  |sed 's/.markdup.sorted.bam//g'  >>${i}.poplist
	done
done


###############################################################selscan,xpnsl##########################################

for scaffold in `cut -f1 B1952_final.fasta.fai | grep -vE 'scaffold18802_MitoZ|WJBF0'  `
do
	vcftools --gzvcf ${scaffold}.phase.vcf.gz --recode --recode-INFO-all --stdout --keep Central.poplist |bgzip -c  > Central.${scaffold}.phase.vcf.gz

done


for pop in {AB,BM,HN,NE,QH,TW}
do
	for scaffold_round1 in `cut -f1 B1952_final.fasta.fai | grep -vE 'scaffold18802_MitoZ|WJBF0' `
	do
	vcftools --gzvcf ${scaffold_round1}.phase.vcf.gz --recode --recode-INFO-all --stdout --keep ${pop}.poplist |bgzip -c > ${pop}.${scaffold_round1}.phase.vcf.gz
	selscan --xpnsl --vcf ${pop}.${scaffold_round1}.phase.vcf.gz --vcf-ref Central.${scaffold_round1}.phase.vcf.gz --threads 10 --out ${pop}.${scaffold_round1}
	rm ${pop}.${scaffold_round1}.xpnsl.log
	done
	
	norm --xpnsl --files ${pop}.*.xpnsl.out --bp-win  --winsize 10000  --qbins 10 --min-snps 10 

	sed -n '25,34p' ${pop}.norm.logfile |cut -d' ' -f1 |sed '1i\10' |sed '1i\bin_boundary' > ${pop}.bin_boundary.txt  

	for scaffold_round2 in `cut -f1 B1952_final.fasta.fai | grep -vE 'scaffold18802_MitoZ|WJBF0' `
	do
	tail -n +2 ${pop}.${scaffold_round2}.xpnsl.out.norm | awk '{$1="'''${scaffold_round2}'''";print $0}' |sed 's/ /\t/g' > ${pop}.${scaffold_round2}.xpnsl.out.norm.temp 
	mv ${pop}.${scaffold_round2}.xpnsl.out.norm.temp ${pop}.${scaffold_round2}.xpnsl.out.norm

	awk '{print "'''${scaffold_round2}'''""\t"$0}' ${pop}.${scaffold_round2}.xpnsl.out.norm.10kb.windows >${pop}.${scaffold_round2}.xpnsl.out.norm.10kb.windows.temp 
	mv ${pop}.${scaffold_round2}.xpnsl.out.norm.10kb.windows.temp ${pop}.${scaffold_round2}.xpnsl.out.norm.10kb.windows
	done

	cat ${pop}.*.xpnsl.out.norm >${pop}.xpnsl.out.norm
	sed -i '1i\Chrom\tpos\tgpos\tp1\tsL1\tp2\tsL2\txpnsl\tnormxpehh\tcrit' ${pop}.xpnsl.out.norm 
	chmod 544 ${pop}.xpnsl.out.norm
	rm ${pop}.*.xpnsl.out.norm
	tail -n +2 ${pop}.xpnsl.out.norm |awk '$9>=2 {print}' |awk '{print "'''${pop}'''""\t"$0}' |sed  '1i\Pop\t\Chrom\tpos\tgpos\tp1\tsL1\tp2\tsL2\txpnsl\tnormxpehh\tcrit' > ${pop}.xpnsl.out.score2.norm 

	cat ${pop}.*.xpnsl.out.norm.10kb.windows >${pop}.xpnsl.out.norm.10kb.windows
	awk '{print "'''${pop}'''""\t"$0}' ${pop}.xpnsl.out.norm.10kb.windows |sed  '1i\Pop\t\Chrom\t\win_start\twin_end\tsites_in_win\tfrac_scores_gt_threshold\tfrac_scores_lt_threshold\tapprox_percentile_for_gt_threshold_wins\tapprox_percentile_for_lt_threshold_wins\tmax_score\tmin_score' >${pop}.xpnsl.out.norm.10kb.windows.temp  
	mv ${pop}.xpnsl.out.norm.10kb.windows.temp ${pop}.xpnsl.out.norm.10kb.windows
	chmod 544 ${pop}.xpnsl.out.norm.10kb.windows
	rm ${pop}.*.xpnsl.out.norm.10kb.windows 

	Rscript State.xpnsl.out.norm.10kb.windows.heterogeneity.R ${pop}.xpnsl.out.norm.10kb.windows ${pop}.bin_boundary.txt
	mv Winsite_Proportion2.pdf ${pop}.Winsite_Proportion2.pdf
	mv Bin_Boundaries_window.table ${pop}.Bin_Boundaries_window.table
	
	head -1 all_Window10K_Site10_background_outlier.table >header
	tail -n +2 all_Window10K_Site10_background_outlier.table |sort -k2,2 -k3n,3  >body
	cat header body >${pop}.Window10K_Site10_background_outlier.table
	rm header body all_Window10K_Site10_background_outlier.table
	

	grep outlier ${pop}.Window10K_Site10_background_outlier.table |awk '$0="Window"NR"_"$0' |cut -f1-4,10 |sed 's/.1_RagTag//g' |sed 's/CM02200//g' |sed 's/CM0220//g' > ${pop}.outlier_window.PreMagma_geneLOC_simulate.txt 
	cut -f2,5 ${pop}.outlier_window.PreMagma_geneLOC_simulate.txt |awk '$0=$0"\t[[:space:]]"' |awk '{print $1$3$2}' > ${pop}.outlier.window.Chrom.maxSnp.score  
	cut -f1,2,3,10 ${pop}.xpnsl.out.score2.norm |sed '1d' |awk '$0="Snp"NR"_"$0' |sed 's/.1_RagTag//g' |sed 's/CM02200//g' |sed 's/CM0220//g' |awk '{print $1"\t"$2"\t"$4"\t"$3}' |sed '1i\snp_id\tChrom\tnormxpehh\tpos' > ${pop}.xpnsl.out.score2.norm.simple
	for i in `cat ${pop}.outlier.window.Chrom.maxSnp.score`;do grep -w $i ${pop}.xpnsl.out.score2.norm.simple >> ${pop}.outlier.window.maxSnp.position;done 
	awk '{print $1"\t"$2"\t"$4"\t"$3}' ${pop}.outlier.window.maxSnp.position  > ${pop}.outlier.window.maxSnp.position.PreMagma_simulate.txt 
	magma --annotate nonhuman window=0,0 --snp-loc  ${pop}.outlier.window.maxSnp.position.PreMagma_simulate.txt --gene-loc ${pop}.outlier_window.PreMagma_geneLOC_simulate.txt --out ${pop}.maxSnp_in_outlier.window 
	for i in `cut -f1 ${pop}.outlier.window.maxSnp.position.PreMagma_simulate.txt`;do grep -w $i ${pop}.maxSnp_in_outlier.window.genes.annot && echo $i "exist in outlier window" || echo $i "not exist in outlierwindow" >>need_remove_snp.id;done  
	mv ${pop}.outlier.window.maxSnp.position.PreMagma_simulate.txt ${pop}.outlier.window.maxSnp.position.PreMagma_simulate.filter
	for i in `cut -d' ' -f1 need_remove_snp.id`;do grep -v -w $i ${pop}.outlier.window.maxSnp.position.PreMagma_simulate.filter > ${pop}.outlier.window.maxSnp.position.filter.temp; mv ${pop}.outlier.window.maxSnp.position.filter.temp ${pop}.outlier.window.maxSnp.position.PreMagma_simulate.filter;done 
	rm need_remove_snp.id  ${pop}.outlier.window.Chrom.maxSnp.score ${pop}.outlier.window.maxSnp.position *.log
	magma --annotate nonhuman window=0,0 --snp-loc  ${pop}.outlier.window.maxSnp.position.PreMagma_simulate.filter --gene-loc maker.geneLOC.PreMagma --out ${pop}.maxSnp_in_maker 
	grep -v '#' ${pop}.maxSnp_in_outlier.window.genes.annot |cut -f3- > ${pop}.window_max_snp
	paste ${pop}.outlier_window.PreMagma_geneLOC_simulate.txt ${pop}.window_max_snp > ${pop}.outlier_window.max_snp.table
	rm ${pop}.maxSnp_in_outlier.window.genes.annot  ${pop}.window_max_snp  ${pop}.outlier_window.PreMagma_geneLOC_simulate.txt ${pop}.maxSnp_in_maker.log 
	mkdir ${pop}
	mv *.pdf *.table  ${pop}.*   ${pop}
done


                      #############################dxy,pi,Fst ,selscan outlier##############
tail -n +2 all_pop6.Window10K_Site10_background_outlier.table |sed 's/.1_RagTag//g' |cut -f1,2,3 >selscan.id
tail -n +2 Allpop6_Central.background_outlier.table  |cut -f1,2,3  >dxy_pi_fst.id
sort dxy_pi_fst.id selscan.id  | uniq -d > dxy_pi_fst.selscan.overlap.id 
sort -k1,1  -k2,2 -k3,3n dxy_pi_fst.selscan.overlap.id >dxy_pi_fst.selscan.overlap.sorted.id
sed 's/\t/[[:space:]]/g' dxy_pi_fst.selscan.overlap.sorted.id > grep.overlap.dxy_pi_fst.id
awk '{print $1,$2".1_RagTag",$3}' dxy_pi_fst.selscan.overlap.sorted.id |sed 's/ /\t/g' > grep.overlap.selscan.id
sed -i 's/\t/[[:space:]]/g'  grep.overlap.selscan.id 
LC_ALL=C grep -wf grep.overlap.dxy_pi_fst.id Allpop6_Central.background_outlier.table |cut -f1-8  > dxy_pi_fst.overlap.table
LC_ALL=C grep -wf grep.overlap.selscan.id all_pop6.Window10K_Site10_background_outlier.table  |cut -f12 >selscan.overlap.table
sort -k1,1  -k2,2 -k3,3n selscan.overlap.table > selscan.overlap.sorted.table
sort -k1,1  -k2,2 -k3,3n dxy_pi_fst.overlap.table > dxy_pi_fst.overlap.sorted.table
cut -f12 selscan.overlap.sorted.table |paste dxy_pi_fst.overlap.sorted.table -  >selscan.dxy_pi_fst.heterogeneity.table
sed -i '1i\population\tchromosome\twindow_pos_1\twindow_pos_2\tavg_dxy\tavg_wc_fst\tavg_pi.pop\tavg_pi.central\tselscan_heterogeneity' selscan.dxy_pi_fst.heterogeneity.table
Rscript Selscan_Fst_pi_dxy_window.background_outlier_boxplot.R selscan.dxy_pi_fst.heterogeneity.table Dxy_Selscan_background_outlier_boxplot.pdf Fst_Selscan_background_outlier_boxplot.pdf Pop_Pi_Selscan_background_outlier_boxplot.pdf Central_Pi_Selscan_background_outlier_boxplot.pdf mean_pi.pdf

          ########################Rho selscan outlier################################################
grep -vE 'ML|Central' all_pop_FastEPRR.result |tail -n +2  |cut -f1,2,3  >FastEPRR.id
sort FastEPRR.id selscan.id  | uniq -d > FastEPRR.selscan.overlap.id 
sort -k1,1  -k2,2 -k3,3n FastEPRR.selscan.overlap.id >FastEPRR.selscan.overlap.sorted.id 
sed 's/\t/[[:space:]]/g' FastEPRR.selscan.overlap.sorted.id > grep.overlap.FastEPRR.id
awk '{print $1,$2".1_RagTag",$3}' FastEPRR.selscan.overlap.sorted.id |sed 's/ /\t/g' > grep.overlap1.selscan.id
LC_ALL=C grep -wf grep.overlap.FastEPRR.id all_pop_FastEPRR.result |cut -f1-5  > Rho.overlap.table
sort -k1,1  -k2,2 -k3,3n Rho.overlap.table > Rho.overlap.sorted.table 
LC_ALL=C grep -wf grep.overlap1.selscan.id all_pop6.Window10K_Site10_background_outlier.table > selscan.overlap1.table
sort -k1,1  -k2,2 -k3,3n selscan.overlap1.table > selscan.overlap1.sorted.table
cut -f12 selscan.overlap1.sorted.table|paste Rho.overlap.sorted.table -  > selscan.Rho.heterogeneity.table

sed -i '1i\population\tchromosome\twindow_pos_1\twindow_pos_2\tRho\tselscan_heterogeneity' selscan.Rho.heterogeneity.table
Rscript Selscan_Rho_window.background_outlier_boxplot.R selscan.Rho.heterogeneity.table Selscan_Rho_background_outlier_boxplot.pdf

echo job end at `data`
