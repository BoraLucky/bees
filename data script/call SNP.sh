mkdir -p 01_GATK_per_Sample/{AB,BM,Central,HN,LN_mixture,ML,NE,QH,TW}
for i in {AB,BM,Central,HN,LN_mixture,ML,NE,QH,TW}
do
	cd  01_GATK_per_Sample/$i
	find /home/00_Apis_cerana/00_NGS_mapping/02_mapping/$i -name *.bam >${i}.bam.list
	split -10  -d ${i}.bam.list  ${i}.bam_list.split
	for split_part in `ls ${i}.bam_list.split*`
	do 
	mkdir ${split_part}.dir
	mv $split_part ${split_part}.dir
	cd ${split_part}.dir
	echo "echo GATK.job.started at \`date\`
ln -s /home/00_Apis_cerana/04_assembly/B1952/B1952_final_genome.tosingle.upper.fasta B1952_final_genome.fasta
samtools faidx B1952_final_genome.fasta
gatk CreateSequenceDictionary -R B1952_final_genome.fasta  -O B1952_final_genome.dict
for sample_bams in \`cat $split_part\` " >>${split_part}.gatk.sh
echo 'do
	sample_name=`basename $sample_bams |cut -d"." -f1`
	gatk --java-options "-Xmx8g -XX:ParallelGCThreads=8" HaplotypeCaller --verbosity ERROR -R B1952_final_genome.fasta -I ${sample_bams} -O  ${sample_name}.g.vcf.gz -ERC GVCF --exclude-intervals scaffold18802_MitoZ
	echo ${sample_name}.GATK.job.end at `date`	
done'	 >>${split_part}.gatk.sh
echo "echo GATK.job.end at \`date\`" >>${split_part}.gatk.sh
	cd ..
	done
	cd ../..
done

echo job_started at `date`
mkdir 02_GenomicsDBImport && cd 02_GenomicsDBImport/



find /home/00_Apis_cerana/07_island/01_GATK_per_Sample/ -name *.g.vcf.gz >293samples_gvcf_Path_list

for i in `cat 293samples_gvcf_Path_list`;do basename $i |sed 's/.g.vcf.gz//g' >>smples_name;done
paste smples_name 293samples_gvcf_Path_list >293sample_mapfile

for scaffold in `grep '>' B1952_final_genome.fasta |grep -v scaffold18802_MitoZ |sed  's/>//g' `
do
	mkdir ${scaffold} && cp 293sample_mapfile ${scaffold}
	echo "
echo  job_stared at \`date\`
gatk --java-options "-Xmx8g"  GenomicsDBImport --sample-name-map 293sample_mapfile  --genomicsdb-workspace-path ${scaffold}_GenomicsDB --intervals ${scaffold}  --reader-threads 5
echo job_end  at \`date\` " >${scaffold}/${scaffold}.sh
done


date
mkdir 03_genetypeGVCFs_allSites && cd 03_genetypeGVCFs_allSites
ln -s /home/00_Apis_cerana/04_assembly/B1952/B1952_final_genome.tosingle.upper.fasta B1952_final_genome.fasta

for scaffold in `grep '>' B1952_final_genome.fasta |grep -v scaffold18802_MitoZ |sed  's/>//g' `
do
	mkdir ${scaffold}
	echo "
echo  job_stared at \`date\`
ln -s  /home/00_Apis_cerana/04_assembly/B1952/B1952_final_genome.tosingle.upper.fasta B1952_final_genome.fasta
samtools faidx B1952_final_genome.fasta
gatk CreateSequenceDictionary -R B1952_final_genome.fasta  -O B1952_final_genome.dict
ln -s ../../02_GenomicsDBImport/${scaffold}/${scaffold}_GenomicsDB ./
gatk --java-options "-Xmx8g" GenotypeGVCFs -L ${scaffold} -all-sites -R B1952_final_genome.fasta -V gendb://${scaffold}_GenomicsDB -O ${scaffold}.GenotypeGVCFs.allSites.vcf.gz
echo job_end  at \`date\` " >${scaffold}/${scaffold}.GenotypeGVCFs.allSites.sh
done


echo job started at `date`
mkdir 04_merge_filter/
find `pwd`/03_genetypeGVCFs_allSites -name *.GenotypeGVCFs.allSites.vcf.gz >04_merge_filter/293_genetypeGVCFs_list
cd 04_merge_filter/

gatk MergeVcfs -I 293_genetypeGVCFs_list  -O 293_MergeVcfs.raw.allSites.vcf.gz

ln -s /home/00_Apis_cerana/04_assembly/B1952/B1952_final_genome.tosingle.upper.fasta B1952_final_genome.fasta
samtools faidx  B1952_final_genome.fasta
gatk CreateSequenceDictionary -R B1952_final_genome.fasta  -O B1952_final_genome.dict


gatk SelectVariants -R B1952_final_genome.fasta -V 293_MergeVcfs.raw.allSites.vcf.gz --select-type-to-include SNP  -O 293_MergeVcfs.snp.vcf.gz 
gatk SelectVariants -R B1952_final_genome.fasta -V 293_MergeVcfs.raw.allSites.vcf.gz --select-type-to-include INDEL  -O 293_MergeVcf.INDEL.vcf.gz 
gatk SelectVariants -R B1952_final_genome.fasta -V 293_MergeVcfs.raw.allSites.vcf.gz --select-type-to-include NO_VARIATION -O 293_MergeVcf.invariant_sites.vcf.gz  
gatk VariantFiltration -R B1952_final_genome.fasta -V 293_MergeVcfs.snp.vcf.gz --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "Hard_Filter"  -O 293_MergeVcfs.snp.hard_filter.vcf.gz 
gatk VariantFiltration -R B1952_final_genome.fasta -V 293_MergeVcf.INDEL.vcf.gz --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "Hard_Filter" -O 293_MergeVcf.INDEL.hard_filter.vcf.gz 


java -jar /home/00_biosoft/04_file_process_soft/snpEff/SnpSift.jar filter "( na FILTER ) | (FILTER = 'PASS')" 293_MergeVcfs.snp.hard_filter.vcf.gz >293_MergeVcfs.snp.hard_filter_pass.vcf 
bgzip 293_MergeVcfs.snp.hard_filter_pass.vcf
tabix 293_MergeVcfs.snp.hard_filter_pass.vcf.gz
java -jar /home/00_biosoft/04_file_process_soft/snpEff/SnpSift.jar filter "( na FILTER ) | (FILTER = 'PASS')"  293_MergeVcf.INDEL.hard_filter.vcf.gz > 293_MergeVcf.INDEL.hard_filter_pass.vcf  
bgzip 293_MergeVcf.INDEL.hard_filter_pass.vcf
tabix 293_MergeVcf.INDEL.hard_filter_pass.vcf.gz
rm -rf 293_MergeVcfs.snp.vcf.gz 293_MergeVcf.INDEL.vcf.gz 293_MergeVcfs.snp.hard_filter.vcf.gz 293_MergeVcf.INDEL.hard_filter.vcf.gz 293_MergeVcfs.snp.vcf.gz.tbi 293_MergeVcf.INDEL.vcf.gz.tbi 293_MergeVcfs.snp.hard_filter.vcf.gz.tbi 293_MergeVcf.INDEL.hard_filter.vcf.gz.tbi  
 
gatk MergeVcfs -I 293_MergeVcfs.snp.hard_filter_pass.vcf.gz -I 293_MergeVcf.INDEL.hard_filter_pass.vcf.gz  -O 293_MergeVcf.hard_filter_pass.variant_sites.vcf.gz # 
bcftools filter -e 'AC==0 || AC==AN' --SnpGap 3 293_MergeVcf.hard_filter_pass.variant_sites.vcf.gz |bcftools view -m2 -M2 -v snps -O z -o 293_MergeVcf.hard_filter_pass.bialleic_snp.vcf.gz 
tabix 293_MergeVcf.hard_filter_pass.bialleic_snp.vcf.gz
rm -rf 293_MergeVcfs.snp.hard_filter_pass.vcf.gz 293_MergeVcf.INDEL.hard_filter_pass.vcf.gz 293_MergeVcfs.snp.hard_filter_pass.vcf.gz.tbi 293_MergeVcf.INDEL.hard_filter_pass.vcf.gz.tbi


vcftools --gzvcf 293_MergeVcf.hard_filter_pass.bialleic_snp.vcf.gz \
	 --remove-indels --maf 0.01  --minQ 20 \
	 --min-meanDP 5 --max-meanDP 50 \
	 --max-missing 0.75 \
	 --minDP 3 --minGQ 20 \
	 --recode --recode-INFO-all --stdout |bgzip -c > 293_MergeVcf.hard_filter_pass.population_filter.bialleic_snp.vcf.gz
tabix 293_MergeVcf.hard_filter_pass.population_filter.bialleic_snp.vcf.gz
rm -rf 293_MergeVcf.hard_filter_pass.bialleic_snp.vcf.gz 293_MergeVcf.hard_filter_pass.bialleic_snp.vcf.gz.tbi

vcftools --gzvcf 293_MergeVcf.invariant_sites.vcf.gz \
	 --remove-indels --minQ 20 --max-missing 0.75 --min-meanDP 5 --max-meanDP 50 \
	 --recode  --recode-INFO-all --stdout |bgzip -c > 293_MergeVcf.invariant_sites.population_filter.vcf.gz
tabix 293_MergeVcf.invariant_sites.population_filter.vcf.gz 
rm -rf 293_MergeVcf.invariant_sites.vcf.gz 293_MergeVcf.invariant_sites.vcf.gz.tbi


gatk MergeVcfs -I 293_MergeVcf.hard_filter_pass.population_filter.bialleic_snp.vcf.gz -I 293_MergeVcf.invariant_sites.population_filter.vcf.gz  -O 293_hard_filter_population_filter_pass.allsite.vcf.gz
rm -rf 293_MergeVcf.hard_filter_pass.population_filter.bialleic_snp.vcf.gz 293_MergeVcf.invariant_sites.population_filter.vcf.gz 293_MergeVcf.hard_filter_pass.population_filter.bialleic_snp.vcf.gz.tbi 293_MergeVcf.invariant_sites.population_filter.vcf.gz.tbi
echo job end at `date`
