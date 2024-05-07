date

samtools faidx  B1952_final_genome.fasta
gatk CreateSequenceDictionary -R B1952_final_genome.fasta  -O B1952_final_genome.dict

for i in `grep '>' B1952_final_genome.fasta |sed 's/>//g' |grep -v scaffold18802_MitoZ`
do
     gatk SelectVariants \
          -R B1952_final_genome.fasta \
          -V 293_MergeVcf.hard_filter_pass.population_filter.bialleic_snp.vcf.gz \
          -L $i \
          -O ${i}.vcf.gz
     java -Xmx300g  -jar /home/00_biosoft/06_Population_soft/Beagle5.2/beagle.21Apr21.304.jar nthreads=20  gt=${i}.vcf.gz  out=${i}.bialleic_snp.phase
done


################################Rho###############################################################
cp /home/00_Apis_cerana/07_island/05_fst/293sample_pop.file ./
mkdir {AB,BM,HN,ML,NE,QH,TW,Central}
for pop in {AB,BM,HN,ML,NE,QH,TW,Central}
do
	grep $pop 293sample_pop.file | cut -f1 >${pop}.sample
	awk '$2>10000' B1952_final_genome.fasta.fai  |while  read line
	do
		array=($line)
		scaffold=${array[0]}
		end_pos=${array[1]}
	        vcftools --gzvcf ${scaffold}.bialleic_snp.phase.vcf.gz --recode --recode-INFO-all --stdout  --keep ${pop}.sample |bgzip -c  > ${pop}/${pop}_${scaffold}.bialleic_snp.phase.out.vcf.gz
        ############################FastEPRR_VCF_step1######################
		echo "
#!/usr/bin/R 
rm(list=ls()) 
suppressPackageStartupMessages(library(FastEPRR))
suppressPackageStartupMessages(library(mboost))
FastEPRR_VCF_step1(vcfFilePath=\"`pwd`/${pop}/${pop}_${scaffold}.bialleic_snp.phase.out.vcf.gz\",erStart=1,erEnd=${end_pos},winLength=10000,srcOutputFilePath=\"`pwd`/${pop}/step1/${scaffold}\")" >${pop}/${pop}_${scaffold}.FastEPRR.step1.R
       	done
	
done		
##########################################FastEPRR_VCF_step1#######################################
for pop in {AB,BM,HN,ML,NE,QH,TW,Central}
do
	cd $pop
	mkdir step1	
	for i in `ls *FastEPRR.step1.R`
	do
		/usr/bin/Rscript $i
	done
	cd .. 
done
#########################################FastEPRR_VCF_step2 and step3#############################################
for pop in {BM,HN,ML,NE,QH,TW,Central}
do
	cd $pop
	mkdir step2 step3
	echo "#!/usr/bin/R 
rm(list=ls()) 
suppressPackageStartupMessages(library(FastEPRR))
suppressPackageStartupMessages(library(mboost))

FastEPRR_VCF_step2(srcFolderPath=\"/home/00_Apis_cerana/07_island/10_FastEPRR/${pop}/step1\",DXOutputFolderPath=\"/home/00_Apis_cerana/07_island/10_FastEPRR/${pop}/step2\")

FastEPRR_VCF_step3(srcFolderPath=\"/home/00_Apis_cerana/07_island/10_FastEPRR/${pop}/step1\",DXFolderPath=\"/home/00_Apis_cerana/07_island/10_FastEPRR/${pop}/step2\",finalOutputFolderPath=\"/home/00_Apis_cerana/07_island/10_FastEPRR/${pop}/step3\")" >${pop}.FastEPRR.step2_3.R
      cd ..
done
COMMENT
########################################all result####################################################################
mkdir all_pop_FastEPRR_result 
cd all_pop_FastEPRR_result

for Popname in {AB,BM,Central,HN,ML,NE,QH,TW}
do
	for scaffold in  `grep CM0220 B1952_final_genome.fasta.fai|cut -f1`
	do
		awk '$0="'''$Popname''' '''$scaffold''' "$0'  /home/00_Apis_cerana/07_island/10_FastEPRR/${Popname}/step3/chr${scaffold}    >> all_pop_FastEPRR.result.temp
	done
done

grep -v Start  all_pop_FastEPRR.result.temp |sed  '1i population chromosome window_pos_1 window_pos_2 Rho CIL CIR' |sed 's/ /\t/g' |sed -i 's/.1_RagTag//g'  >all_pop_FastEPRR.result

date

