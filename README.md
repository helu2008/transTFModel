# transTFModel
This scripts are used to replicate the results in the paper "Tissue-specific variations in transcription factors elucidate complex immune system regulation" to identify variations of TFs and their SNPs that related with gene expression.

## tf-expression.py
To run this script, you will need: 
1. Normalized gene expression values from GTEx Version 8 
2. Covariates files
3. Results from Predixcan: \
`python Predict.py --model_db_path /GTEx_v8/eqtl/mashr/mashr_tissue.db --vcf_genotypes GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf --vcf_mode genotyped --pred_output predixcan_output.txt --pred_summary_output predixcan_summary_output.txt --verbosity 9 --throw --model_db_snp_key varID --on_the_fly_mapping METADATA "{}_{}_{}_{}_b38"`

Example of running tf-expression.py on tissue of Muscle Skeletal: \
`python tf-expression.py -g Muscle_Skeletal.v8.normalized_expression.bed -p predixcan_output_Muscle_Skeletal.txt -c Muscle_Skeletal.v8.covariates.txt-o tf-expression_output.txt` \
The result of hit genes, their TFs as well as the weights will be saved in the output file.

## tf-binding.py
To run this script, you will need: 
1. All input files from tf-expression.py and 
2. genotype .bed, .bim and .fam files from .vcf file using plink: \
`plink --file GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf --make-bed --out nssnp_genotype` 
3. Annotate and filter the nsSNPs using SnpEff and SnpSift: \
For SNP data, we use All_20170403.vcf. There is a similar version available at: https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/All_20170710.vcf.gz \
`java -Xmx8g -jar snpEff.jar -v -stats stats.html GRCh38.78 All_20170710.vcf.gz > annotated_snpeff.vcf` \
`java -jar SnpSift.jar annotate annotated_snpeff.vcf > filtered_snpsift.vcf` 
4. Filter deleterious SNPs using SIFT score: \
`java -jar GRCh38.78/SIFT4G_Annotator.jar -c -i filtered_snpsift.vcf -d GRCh38.78/ -r deleterious_after_filtered.vcf` 

Example of running tf-expression.py on tissue of Muscle Skeletal: \
`python tf-binding.py -g Muscle_Skeletal.v8.normalized_expression.bed -p predixcan_output_Muscle_Skeletal.txt -bed nssnp_genotype.bed -bim nssnp_genotype.bim -fam nssnp_genotype.fam -snpSift filtered_snpsift.vcf -s deleterious_after_filtered.vcf -o tf-binding_output.txt` \
The result of hit genes, the nsSNPs of their TFs as well as the weights will be saved in the output file.
