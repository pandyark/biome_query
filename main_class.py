_author_ = "Riddhika Pandya"
_version_ = "1.0"
'''
    This script takes a .txt file containing a list of the genes needs to be pulled from input vcf file. The program takes 
    the vcf file and checks if it is .gz or not. It will gzip using bgzip and generates an index file using tabix. Then, bcftools
    removes homozygous calls from the vcf files and generates the file with heterozygous calls only. The program then calculates the HWE 
    test statistics using vcftools and generates HEW p-values for the variants in a separate file with .hew extension. GATK varianttotable
    tool then parse the vcf file into a text file. With the help of Pandas, the program then concatenates HWE statistics and parsed vcf files 
    by the variant positions and carries other calculations like Odd ratio and ORS.  
'''
import os
import sys
import subprocess
from os import path
import pandas as pd
import numpy as np



input_gene_file = sys.argv[1]
input_vcf_file = sys.argv[2]

# This block is to run on the cluster
vcftools=""
gatk=""
bcftools=""
tabix=""
bgzip=""
script=""


# Please remove this block comments to run this on local computer
# This block is to run on local computer
gatk="/Users/bowcock_lab/Desktop/Analysis_Softwares/gatk-4.1.5.0/gatk"
vcftools="/Users/bowcock_lab/Desktop/Analysis_Softwares/vcftools-master/src/cpp/vcftools"
bcftools="/Users/bowcock_lab/Desktop/Analysis_Softwares/bcftools-1.10.2/bcftools"
bgzip="/Users/bowcock_lab/Desktop/Analysis_Softwares/tabix-0.2.6/bgzip"
tabix="/Users/bowcock_lab/Desktop/Analysis_Softwares/tabix-0.2.6/tabix"
script="/Users/bowcock_lab/Desktop/Analysis_Softwares/get_gene_coordinates.R"



def check_index_file(file):
    index_file = file.replace(".vcf.gz", ".vcf.gz.tbi")
    if path.exists(index_file) == True:
        pass
    if path.exists(index_file) == False:
        cmd3 = [tabix, "-p", "vcf", file]
        process3 = subprocess.Popen(cmd3)
        process3.wait()


def file_check(vcf_file):
    if vcf_file.endswith(".vcf.gz"):
        return vcf_file
        check_index_file(vcf_file)

    if vcf_file.endswith(".vcf"):
        comp_vcf = vcf_file.replace(".vcf", ".vcf.gz")
        cmd2 = [bgzip, vcf_file]
        process2 = subprocess.Popen(cmd2)
        process2.wait()
        check_index_file(comp_vcf)
        return comp_vcf




try:
    # output_bed_file = input_gene_file.replace(".csv",".bed")
    # process1 = subprocess.Popen(["Rscript",script,input_gene_file,output_bed_file],bufsize=1)
    # process1.communicate()

    output_file_name = (os.path.splitext(input_vcf_file)[0])
    if input_vcf_file.endswith("vcf.gz"):
        output_file_1 = (os.path.splitext(input_vcf_file)[0])
        output_file_final = (os.path.splitext(output_file_1)[0])
        output_vcf_file = output_file_final + "_gene_pull"
        hwe_test_file = output_vcf_file + "_HWE"

    if input_vcf_file.endswith((".vcf")):
        output_file_final = (os.path.splitext(input_vcf_file)[0])
        output_vcf_file = output_file_name + "_gene_pull"

    print(" ")
    print("===================================================================================================")
    print(" ")
    print("Bed file has been generated. The analysis is moving forward to vcftools.")
    print(" ")
    print("===================================================================================================")
    print(" ")

    zip_check_file = file_check(input_vcf_file)


    cmd1 = [vcftools,"--gzvcf",zip_check_file,"--bed",input_gene_file,"--out",output_vcf_file,"--recode","--keep-INFO-all"]
    process1 = subprocess.Popen(cmd1)
    process1.communicate()

    recoded_file = output_vcf_file + ".recode.vcf"
    rm_hom_gt_file = file_check(recoded_file)
    print(rm_hom_gt_file)


    het_vcf = rm_hom_gt_file.replace(".vcf.gz","_Het_Calls.vcf")

    cmd4 = [bcftools,"view","--genotype","het",rm_hom_gt_file]
    process4 = subprocess.Popen(cmd4,stdout=open(het_vcf,"w"))
    process4.wait()

    print(" ")
    print("===================================================================================================")
    print(" ")
    print("The Homozygous call has been removed by vcftools. A New vcf file has been generated with Heterozygous call only.")
    print(" ")
    print("===================================================================================================")
    print(" ")

    #hwe_test_file = het_vcf.replace(".vcf","_HWE")
    hwe_test_file = (os.path.splitext(het_vcf)[0]) + "_HWE"
    logfile = hwe_test_file.replace(".hwe","_logfile.txt")

    cmd5 = [vcftools,"--vcf",het_vcf,"--hardy","--out",hwe_test_file]
    process5 = subprocess.Popen(cmd5,stdout=open(logfile,'w'))
    process5.communicate()

    print(" ")
    print("===================================================================================================")
    print(" ")
    print("Hardy Weinberg Equilibrium Test has been done, a file with HWE stat has been generated.")
    print(" ")
    print("===================================================================================================")
    print(" ")

    vtt_file = het_vcf.replace(".vcf","_VTT.txt")

    cmd6 = [gatk,"VariantsToTable","-V",het_vcf,"-F","CHROM" ,"-F","POS","-F","REF","-F","ALT","-F","TYPE","-F","Gene.refGene",
              "-F","GeneDetail.refGene","-F","AAChange.refGene","-F" ,"Func.refGene","-F","ExonicFunc.refGene","-F" ,"AC","-F","DP",
              "-F","MQ","-F","QD", "-F","AF","-F","CLNDN","-F","CLNSIG","-F","gnomAD_genome_ALL","-F","gnomAD_genome_AFR","-F","gnomAD_genome_AMR",
              "-F","gnomAD_genome_ASJ","-F","gnomAD_genome_EAS","-F","AF_female","-F","gnomAD_genome_FIN","-F","AF_male","-F","gnomAD_genome_NFE",
              "-F","gnomAD_genome_OTH","-F","AF_popmax","-F","AF_raw","-F","AF_sas","-F","AS_FS","-F","ExAC_AFR","-F","ExAC_ALL","-F","ExAC_AMR",
              "-F","ExAC_EAS","-F","ExAC_FIN","-F","ExAC_NFE","-F","ExAC_OTH","-F","ExAC_SAS","-F","Kaviar_AF","-F","ExcessHet",
              "-F","CADD_phred","-F","DANN_score","-F","FATHMM_pred","-F","FATHMM_score","-F","MetaLR_pred","-F","MetaLR_score",
              "-F","MetaSVM_pred","-F","MetaSVM_score","-F","MutationTaster_pred","-F","MutationTaster_score","-F","PROVEAN_pred",
              "-F","PROVEAN_score","-F","Polyphen2_HDIV_pred","-F","Polyphen2_HDIV_score","-F","Polyphen2_HVAR_pred","-F","Polyphen2_HVAR_score",
              "-F","SIFT_pred","-F","SIFT_score","-F","avsnp150","-F","fathmm-MKL_coding_pred","-F","fathmm-MKL_coding_score","-F","HET",
              "-F","HOM-REF","-F","HOM-VAR","-F","NO-CALL","-F","VAR","-F","NSAMPLE","-F","NCALLED","-O",vtt_file]
    process6 = subprocess.Popen(cmd6)
    process6.communicate()

    print(" ")
    print("===================================================================================================")
    print(" ")
    print("vcf file parsing has been done by GATK.")
    print(" ")
    print("====================================================================================================")
    print(" ")


    hwe_file = hwe_test_file + ".hwe"
    vtt_df1 = pd.read_csv(vtt_file, sep='\t',low_memory=True,keep_default_na=False)
    hwe_df2 = pd.read_csv(hwe_file, sep='\t',low_memory=True,keep_default_na=False)
    hwe_df2.drop(hwe_df2.columns[[0,2,3,5,6,7]], axis=1, inplace=True)
    merge_file = (pd.merge(vtt_df1,hwe_df2,on='POS',how='left'))
    merge_file['AF'] = pd.to_numeric(merge_file.AF, errors='coerce')
    merge_file['gnomAD_genome_NFE'] = pd.to_numeric(merge_file.gnomAD_genome_NFE, errors='coerce')
    merge_file['OR'] = merge_file['AF']/merge_file['gnomAD_genome_NFE']
    merge_file['ORS']= (2 * merge_file['HOM-VAR']+ 1 * merge_file['HET'])/(2 * merge_file['NCALLED'])
    merge_output = vtt_file.replace(".txt","_final_Merge.xlsx")
    merge_file.to_excel(merge_output)

    print(" ")
    print("===================================================================================================")
    print(" ")
    print("Program has been run successfully and final file is generated.")
    print(" ")
    print("===================================================================================================")
    print(" ")


except Exception as e:
    print("An Error has occurred generating the  {}".format(e))