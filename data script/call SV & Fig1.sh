    ##################################delly#######################################
mkdir -p 03_delly/{AB,BM,Central,HN,LN_mixture,ML,NE,QH,TW}
for i in {AB,BM,Central,HN,LN_mixture,ML,NE,QH,TW}
do 
	find 02_mapping/${i} -name *.markdup.sorted.bam |while read line
	do
	bamfile=($line)
	name=`basename $bamfile |sed 's/.markdup.sorted.bam//g'`
	delly call -t all -s 15  -q 20 -g  B1952_final_genome.fasta  -o 03_delly/${i}/${name}.delly.bcf $bamfile 1>03_delly/${i}/${name}.delly.str.out  2>03_delly/${i}/${name}.delly.str.error
	bcftools convert  -O v --threads 8  -o 03_delly/${i}/${name}.delly.vcf  03_delly/${i}/${name}.delly.bcf
	vcf-sort -c 03_delly/${i}/${name}.delly.vcf >03_delly/${i}/${name}.delly.sorted.vcf 
	done
done
     ##################################menta (py2)######################################

mkdir -p 04_manta/{AB,BM,Central,HN,LN_mixture,ML,NE,QH,TW}
for i in {AB,BM,Central,HN,LN_mixture,ML,NE,QH,TW}
do
	find 02_mapping/${i} -name *.markdup.sorted.bam |while read line
	do
	bamfile=($line)
	name=`basename $bamfile |sed 's/.markdup.sorted.bam//g'`
	/home/00_biosoft/02_SV_soft/manta-1.6.0/bin/configManta.py --bam $bamfile  --referenceFasta B1952_final_genome.fasta --runDir 04_manta/${i}/${name} 
	04_manta/${i}/${name}/runWorkflow.py -m local -j 20  -m local -g 100  1>04_manta/${i}/${name}/${name}.str.out 2>04_manta/${i}/${name}/${name}.str.err
                  
	gunzip  --stdout 04_manta/${i}/${name}/results/variants/diploidSV.vcf.gz >04_manta/${i}/${name}/results/variants/${name}.manta.vcf
	gunzip  --stdout 04_manta/${i}/${name}/results/variants/candidateSV.vcf.gz | grep -v '^#' >> 04_manta/${i}/${name}/results/variants/${name}.manta.vcf
	gunzip  --stdout 04_manta/${i}/${name}/results/variants/candidateSmallIndels.vcf.gz | grep -v '^#' >>04_manta/${i}/${name}/results/variants/${name}.manta.vcf
	                               
 	python /home/00_biosoft/02_SV_soft/manta-1.6.0/libexec/convertInversion.py /home/00_biosoft/04_file_process_soft/Samtools-1.10/bin/samtools B1952_final_genome.fasta 04_manta/${i}/${name}/results/variants/${name}.manta.vcf >04_manta/${i}/${name}/results/variants/${name}.mantaConvertINV.vcf 
                 
	mkdir -p 04_manta/Menta_allsample_out/${i}
	vcf-sort -c 04_manta/${i}/${name}/results/variants/${name}.mantaConvertINV.vcf >04_manta/Menta_allsample_out/${i}/${name}.Menta.ConvertINV.sorted.vcf
	done
done

                  ##########################smoov########################     

mkdir -p 05_smoove/{AB,BM,Central,HN,LN_mixture,ML,NE,QH,TW}
for i in {AB,BM,Central,HN,LN_mixture,ML,NE,QH,TW}

do 
	find 02_mapping/${i} -name *.markdup.sorted.bam |while read line
	do
	bamfile=($line)
	name=`basename $bamfile |sed 's/.markdup.sorted.bam//g'`
	smoove call --outdir 05_smoove/${i}/${name}  --name ${name} --fasta B1952_final_genome.fasta -p 3 --genotype $bamfile
	zcat 05_smoove/${i}/${name}/${name}-smoove.genotyped.vcf.gz |vcf-sort -c >05_smoove/${i}/${name}/${name}.smoove.sorted.vcf
	done
done

                   ########################merge#################################
mkdir -p 06_each_sample_merge_sv/{AB,BM,Central,HN,LN_mixture,ML,NE,QH,TW}

for i in {AB,BM,Central,HN,LN_mixture,ML,NE,QH,TW}
do 
	find 02_mapping/${i} -name *.markdup.sorted.bam |while read line
	do
	bamfile=($line)
	name=`basename $bamfile |sed 's/.markdup.sorted.bam//g'`
	find `pwd` -name ${name}.*sorted.vcf | grep -v survivor_ant  >06_each_sample_merge_sv/${i}/${name}.3vcf.list
	SURVIVOR merge 06_each_sample_merge_sv/${i}/${name}.3vcf.list 1000 2 1 0 0 30  06_each_sample_merge_sv/${i}/${name}.merge.vcf
	done
done
                #############################merge_again_filter_anno#####################
mkdir 07_all_sample_merge_sv
ls 06_each_sample_merge_sv/*/*.merge.vcf >07_all_sample_merge_sv/293_vcf.list 
SURVIVOR merge 07_all_sample_merge_sv/293_vcf.list 1000 1 1 0 0 30 07_all_sample_merge_sv/293merge.vcf
bcftools filter  -i 'SVLEN>-200000 & SVLEN<200000' -Ov --threads 10 -o 07_all_sample_merge_sv/293merge.onlyChrom.filter200K.vcf 07_all_sample_merge_sv/293merge.onlyChrom.vcf
faToTwoBit -long B1952_final_genome.fasta 07_all_sample_merge_sv/B1952_final_genome.fasta.2bit   
twoBitInfo   -nBed 07_all_sample_merge_sv/B1952_final_genome.fasta.2bit 07_all_sample_merge_sv/B1952_final_genome.gap 
vcftools  --vcf 07_all_sample_merge_sv/293merge.onlyChrom.filter200K.vcf --out 07_all_sample_merge_sv/293merge.onlyChrom.filter200K.noGap --exclude-bed 07_all_sample_merge_sv/B1952_final_genome.gap --recode-INFO-all --recode
survivor_ant -o 07_all_sample_merge_sv/293merge.onlyChrom.filter200K.noGap.recode.survivor_ant.vcf -g B1952.makers.gene.gff -i 07_all_sample_merge_sv/293merge.onlyChrom.filter200K.noGap.recode.vcf
vcf-sort -c 07_all_sample_merge_sv/293merge.onlyChrom.filter200K.noGap.recode.survivor_ant.vcf >07_all_sample_merge_sv/293merge.final_temp.vcf
sed -i 's/overlapped_VCF=0//g' 07_all_sample_merge_sv/293merge.temp.vcf

               #########################################Svtyper ################################
find `pwd`/02_mapping -name *.bam >07_all_sample_merge_sv/Sample_bam.list
cd 07_all_sample_merge_sv
cp  293merge.temp.vcf 293merge.gt.vcf

                              
for Bam_Path in `cat Sample_bam.list`
do  
	name=`basename $Bam_Path |sed 's/.markdup.sorted.bam//g'`
	svtyper -B $Bam_Path  -i 293merge.gt.vcf >293merge.gt.${name}.vcf
	echo "${name} job done"
	mv 293merge.gt.${name}.vcf 293merge.gt.vcf
done

vcftools --vcf 293merge.gt.vcf --extract-FORMAT-info GT --stdout >293merge.gt.vcf.gt      
                      
#########################################################293Sample_DEL###############################################
                          ####################################PCA##############################
bcftools view -H 293.final.simplify.DEL.vcf |cut -f 1 |uniq |awk '{print $0"\t"$0}' >filename.chrom-map.txt 
vcftools --vcf 293.final.simplify.DEL.vcf --plink --chrom-map filename.chrom-map.txt --out 293.DEL
plink --file 293.DEL  --recode --out 293.DEL.PCA  --allow-extra-chr --make-bed --pca 20
Rscript tSNE.R 293.DEL.newID.PCA.eigenvec


