
             ##############################LDprune#################
zcat 293_MergeVcf.hard_filter_pass.population_filter.bialleic_snp.vcf.gz |bcftools view -H | cut -f 1 | uniq | awk '{print $0"\t"$0}' >filename.chrom-map.txt  
vcftools --gzvcf 293_MergeVcf.hard_filter_pass.population_filter.bialleic_snp.vcf.gz  --plink --chrom-map filename.chrom-map.txt --out 293_MergeVcf 
plink --file 293_MergeVcf --double-id --set-missing-var-ids @:#  --indep-pairwise 50 10 0.2 --allow-extra-chr --out 293_MergeVcf.ld 
sed -i 's/:/\t/g' 293_MergeVcf.ld.prune.in
vcftools --gzvcf  293_MergeVcf.hard_filter_pass.population_filter.bialleic_snp.vcf.gz --positions 293_MergeVcf.ld.prune.in  --stdout --recode  |bgzip -c > 293_MergeVcf.hard_filter_pass.population_filter.bialleic.LDpruned.vcf.gz 
tabix 293_MergeVcf.hard_filter_pass.population_filter.bialleic.LDpruned.vcf.gz

            #######################easySFS#############################
for i in {AB,BM,HN,NE,QH,TW}
do
	grep -E "$i|CT" down.sample.poplist > ${i}_Central.downsample
	vcftools --gzvcf 293_MergeVcf.hard_filter_pass.population_filter.bialleic.LDpruned.vcf.gz --recode --recode-INFO-all --stdout --keep ${i}_Central.downsample >${i}_Central.downsample.vcf
done

python easySFS.py -i AB_Central.downsample.vcf -p ./AB_Central.downsample -a -f  -o sfs_result/AB --proj 42,48
python easySFS.py -i BM_Central.downsample.vcf -p ./BM_Central.downsample -a -f  -o sfs_result/BM --proj 42,38
python easySFS.py -i HN_Central.downsample.vcf -p ./HN_Central.downsample -a -f  -o sfs_result/HN --proj 42,50
python easySFS.py -i NE_Central.downsample.vcf -p ./NE_Central.downsample -a -f  -o sfs_result/NE --proj 42,36
python easySFS.py -i QH_Central.downsample.vcf -p ./QH_Central.downsample -a -f  -o sfs_result/QH --proj 42,40
python easySFS.py -i TW_Central.downsample.vcf -p ./TW_Central.downsample -a -f  -o sfs_result/TW --proj 42,20
cd ..
###########################################02_model#########################################################
mkdir 02_models/ 
cd 02_models/

for pop in {AB,BM,QH,NE,HN,TW}
do
	cd ${pop}
	sed -i "s/pop1/CT/g" 01-scripts/01-run_model_iteration_v2.sh 
	sed -i "s/pop2/${pop}/g" 01-scripts/01-run_model_iteration_v2.sh 
	
	cp ../../01_sfs/sfs_result/${pop}/dadi/CT-${pop}.sfs  03-data 
	for Models in {IM,IM2m,IM2N,IM2N2m,IM2N2mG,AM,AM2m,AM2N,AM2N2m,AM2N2mG,SI,SC,SC2m,SC2N,SC2N2m,SC2N2mG}
	do
	mkdir ${Models}_model
	echo 'echo -e job started at `date`' >${pop}.${Models}.sh
	echo "sh 01-scripts/00.run_dadi_parallel_v2.sh 03-data/CT-${pop}.sfs $Models folded 20" >>${pop}.${Models}.sh
	echo 'echo -e job end at `date`' >>${pop}.${Models}.sh
	done
	cd ..
done 



