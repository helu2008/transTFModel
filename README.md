# transTFModel
This scripts are used to replicate the results in the paper "Tissue-specific variations in transcription factors elucidate complex immune system regulation" to identify variations of TFs and their SNPs that related with gene expression.

## tf-expression.py
To run this script, you will need: 
1. Normalized gene expression values from GTEx Version 8 
2. Covariates files
3. Results from Predixcan: \
`python Predict.py --model_db_path /GTEx_v8/eqtl/mashr/mashr_tissue.db --vcf_genotypes GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf --vcf_mode genotyped --pred_output predixcan_output.txt --pred_summary_output predixcan_summary_output.txt --verbosity 9 --throw --model_db_snp_key varID --on_the_fly_mapping METADATA "{}_{}_{}_{}_b38"`

Example of running tf-expression.py on tissue of Muscle Skeletal: \
`python tf-expression.py -g Muscle_Skeletal.v8.normalized_expression.bed -p predixcan_output_Muscle_Skeletal.txt -c Muscle_Skeletal.v8.covariates.txt -o tf-expression_output.txt` \
where -g is for the tissue-specific bed file from GTEx, -p is the PrediXcan imputed expression for that tissue, -c is for the covariates file for the same tissue and -o is the name of the output file. The result of hit genes, their TFs as well as the weights will be saved in the output file. The last section of the output file is the result of the robustness check, which should be used to select the hit genes that remain robust after sub-sampling. 

## tf-binding.py
This is a little more complex and requires a preprocessing step to optimize the execution.
To run this script, you will need: 
1. All input files from tf-expression.py and 
2. genotype vcf file GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf 
3. Annotated deleterious nsSNPs using SnpEff and SnpSift: \
For SNP data, we use All_20170403.vcf. There is a similar version available at: https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/All_20170710.vcf.gz \
`java -Xmx8g -jar snpEff.jar -v -stats stats.html GRCh38.78 All_20170710.vcf.gz > annotated_snpeff.vcf` \
`java -jar SnpSift.jar annotate annotated_snpeff.vcf > filtered_snpsift.vcf` 
4. Filter deleterious SNPs using SIFT score: \
`java -jar GRCh38.78/SIFT4G_Annotator.jar -c -i filtered_snpsift.vcf -d GRCh38.78/ -r deleterious_after_filtered.vcf` 

Example of running tf-binding_preprocess.py (run only once, not tissue specific): \
`python tf-binding_preprocess.py -s deleterious_after_filtered.vcf -g GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf -o filtered_snpsift.vcf  -otf Tf2Pos.txt` \
The result are two files, a vcf of GTEx genotype, including only deleterious SNPs and a mapping between a TF and the location in the vcf file of the associated SNPs.

The next step is to run tf-binding.py using the preprocessed files:
Example of running tf-binding.py on tissue of Muscle Skeletal: \
`python tf-binding.py -e Muscle_Skeletal.v8.normalized_expression.bed -p predixcan_output_Muscle_Skeletal.txt -t Tf2Pos.txt -g filtered_snpsift.vcf -c Muscle_Skeletal.v8.covariates.txt -o tf-binding_output.txt` \

where -e is for the tissue-specific bed file from GTEx (expression), -p is the PrediXcan imputed expression for that tissue, -c is for the covariates file for the same tissue, -t and -g are the outputs of the preprocessing step and -o is the name of the output file. The result of hit genes, the nsSNPs of their TFs as well as the weights will be saved in the output file. The last section of the output file is the result of the robustness check, which should be used to select the hit genes that remain robust after sub-sampling. 
