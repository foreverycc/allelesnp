# 02-get genotype by querying vcf files
## ==================================================================================
## Function: get_vcf_info_main
###  Aim: To get the genotypes of input SNPs by vcf files
###  INPUT: input snp file (.csv)
###  OUTPUT: input snp file with information of genotypes of vcf files
## ==================================================================================

# outline ---------------------------------------------------------------------------------------------------------

# # 1. read snp file into data frame
# # snp_list should contain both the df and GRange object
# snp_list = read_snpFile()
# 
# # 2. read peak files into a list of GRange object
# vcf_list = gen_vcf_list(vcf_dir, snp_info_gr)
# 
# # 3. add overlapping info to original snp table
# snp_df_proc = add_vcf_info(snp_info_df, vcf_list)

# 0. load packages ------------------------------------------------------------------------------------------------

load_packages_2.3 = function(){
        packages = c('dplyr', 'VariantAnnotation','GenomicRanges')
        load = lapply(packages, require, character.only = T)
}

# function: write.csv0
source("./src/scripts/T1-toolbox.R")

# toolbox: get short names

get_short_fileName = function(fileName_fullPath) {
        str_vec = strsplit(fileName_fullPath, "/")[[1]]
        str_vec[length(str_vec)]
}

# 0. generate output file name ------------------------------------------------------------------------------------

gen_output_file_vcfInfo = function(snp_info_file, output_dir = "./", sample_name = "") {
        # Aim: to get output file name        
        snp_info_vec = strsplit(snp_info_file, "/")[[1]]
        snp_batch_id = gsub(".csv", "", snp_info_vec[length(snp_info_vec)])
        output_file = paste0(output_dir, "/", snp_batch_id, "_", sample_name, "_genotypeInfo.csv")
        
        return (output_file)
        
}

# 1. read snp file ------------------------------------------------------------------------------------------------

source ("./src/scripts/F2.1-get_alleleDist_info.R")
# same as read_inputSNP_file

# Test #
# snp_info_list = read_inputSNP_file("./data/haploreg_files/LUC_Index+LD_SNPs_20160607.csv")


# 2. read vcf files into a list -----------------------------------------------------------------------------------

gen_vcf_list = function(vcf_dir = NA, vcf_file = NA, snp_info_gr) {
# Aim: to generate vcf list by reading vcf files
# Input: 1. vcf dir; 2. snp_info_gr
# Output: vcf list
        
        # Test #
        # vcf_dir = "./data/samples/DDBJ_A549/vcf_files/"
        
        if (!is.na(vcf_dir)) {
                vcf_files = list.files(vcf_dir, ".vcf$", full.names=T)
                vcf_files_short = list.files(vcf_dir, ".vcf$", full.names=F)
        }
        if (!is.na(vcf_file)) {
                vcf_files = vcf_file
                vcf_files_short = get_short_fileName(vcf_file)
        }
        
        vcf_list = replicate(length(vcf_files), list())
        snp_param = ScanVcfParam(which= snp_info_gr)
        
        # read vcf files
        for (i in 1:length(vcf_files)) {
                cat("reading", vcf_files[i], "\n")
                # read vcf files
                vcf_obj_i = read_vcfFile(vcf_files[i], snp_param)
                vcf_list[[i]] = vcf_obj_i
        }
        
        names(vcf_list) = vcf_files_short
        
        return (vcf_list)
}

# helper function
read_vcfFile = function(vcf_file, param){
        compressVcf <- bgzip(vcf_file, tempfile())
        idx <- indexTabix(compressVcf, "vcf")
        tab <- TabixFile(compressVcf, idx)
        vcf <- readVcf(tab, "hg19", param)
}

# Test #
# vcf_list = gen_vcf_list(snp_info_gr = snp_info_list$snp_info_gr)

# 3. process vcf list data -------------------------------------------------------------------------------------
process_vcf_list = function(vcf_list) {
# To process vcf list
# Output: to get a processed df for each vcf elemnt in the list
# - row: SNPs
# - col: seqnames start end width strand paramRangeID REF ALT QUAL FILTER ref_count alt_count
# helper functions: get_vcf_software, sel_hetSnps_vcf, process_vcf_data
        
        # Test #
        # vcf_file = "./data/samples/DDBJ_A549/vcf_files/A549_snv.GATK.formatcor.vcf"
        
        vcf_df_list = replicate(length(vcf_list), list())
        names(vcf_df_list) = names(vcf_list)
        
        for (i in 1:length(vcf_list)) {
                vcf_data_i = vcf_list[[i]]
                # get the software that generates the vcf file
                vcf_software_i = get_vcf_software (names(vcf_list)[[i]])
                # subset vcf with het genotypes
                vcf_data_i_het = sel_hetSnps_vcf(vcf_data_i) 
                # process vcf data to data.frame format
                vcf_data_i_df = process_vcf_data(vcf_data_i_het, vcf_software_i)
                # assign vcf_df_list with processed data.frame
                vcf_df_list[[i]] = vcf_data_i_df
        }
        
        # return to our df
        return (vcf_df_list)
}

 
# 3-1 helper function 
get_vcf_software = function(vcf_fileName) {
# Aim: to get the software that was used to generate the vcf file
        if (grepl("[Gg][Aa][Tt][Kk]", vcf_fileName)) {
                return ("GATK")
        } 
        if (grepl("[Ss][Aa][Mm][Tt][Oo][Oo][Ll][Ss]", vcf_fileName)) {
                return ("Samtools")
        }
        if (grepl("[Bb][Ii][Ss][Ss][Nn][Pp]", vcf_fileName)) {
                return ("BisSNP")
        } 
        else {
                stop("Please add the software name you used in generating vcf file to your vcf file name. 
                     for example: A549.GATK.vcf")
        }
}

# 3-2 helper function 
sel_hetSnps_vcf = function(vcf_data){
# Aim: to select heterozygous SNPs from vcf data
        
        genotype = geno(vcf_data)$GT
        vcf_data_het = vcf_data[genotype[,1] == "0/1"]
        
        return(vcf_data_het)
}       

# 3-3 helper function 
process_vcf_data = function(vcf_data, vcf_software = "GATK") {
# Aim: to transform vcf data to data.frame format
# Input: vcf_data, (vcf annotation format)
# Output: vcf data frame, with genotype, counts of ref and alt alleles added
# helper function: extract_allelic_dist
        
        # calculate ref and alt allele counts
        vcf_allelic_dist = extract_allelic_dist(vcf_data, vcf_software)
        
        ref_count_vec = sapply(vcf_allelic_dist, function(x) {
                as.numeric(strsplit(x, ",")[[1]][1])
        })
        alt_count_vec = sapply(vcf_allelic_dist, function(x) {
                as.numeric(strsplit(x, ",")[[1]][2])
        })
        
        # add ref and alt allele counts information
        vcf_data_gr = rowRanges(vcf_data)
        names(vcf_data_gr) = 1:length(vcf_data_gr)
        vcf_data_df = as.data.frame(vcf_data_gr)
        vcf_data_df = mutate(vcf_data_df, ref_count = ref_count_vec, alt_count = alt_count_vec)
        
        # modify the data format (from DNAString to character)
        vcf_data_df$ALT = sapply(vcf_data_df$ALT, function(x) as.character(unlist(x)))
        
        return(vcf_data_df)
}

# 3-3-1 helper function
extract_allelic_dist = function(vcf_data, vcf_software = "GATK"){
# Aim: to extract allelic distribution in vcf file generated by various softwares
        
        if (vcf_software == "GATK") {
                AD = geno(vcf_data)$AD
                AD_1 = gsub("c\\(","",paste(AD))
                AD_2 = gsub("\\)","", AD_1)
                AD_3 = gsub(":",", ",AD_2)
                AD_4 = gsub(", ", ",", AD_3)
                return(AD_4)
        }
        
        if (vcf_software == "Samtools") {
                DP4 = info(vcf_data)$DP4
                AD = sapply(DP4, function(x) paste0(x[1] + x[2], ",", x[3] + x[4]))
                return(AD)
        }
        
        if (vcf_software == "BisSNP") {
                DP4 = geno(vcf_data)$DP4[ , , c(1,2,3,4)]
                AD = apply(DP4, 1, function(x){paste0(x[1]+x[2],',',x[3]+x[4])})
                return(AD)
        }
}

# Test #
# vcf_df_list = process_vcf_list(vcf_list)

# 4. to add genotype information to snp info data frame -----------------------------------------------------------

add_vcf_info = function(snp_info_df, vcf_df_list) {
# Aim: to get genotype information added
        
        # add snp genotype information
        for (i in 1:length(vcf_df_list)) {
                vcf_df_i = vcf_df_list[[i]]
                vcf_file_i = names(vcf_df_list)[i]
                # select SNPs with het genotypes
                sel_rows = as.character(snp_info_df$rsID) %in% as.character(vcf_df_i$paramRangeID)
                
                # add vcf information
                snp_info_df[, vcf_file_i] = sel_rows
                snp_info_df[sel_rows, paste0(vcf_file_i, "_ref_count")] = vcf_df_i$ref_count
                snp_info_df[sel_rows, paste0(vcf_file_i, "_alt_count")] = vcf_df_i$alt_count
        }
        
        return(snp_info_df)
}



# main function ---------------------------------------------------------------------------------------------------

get_vcf_info_main = function(snp_info_file, vcf_dir = NA, vcf_file = NA, output_dir = "./", sample_name = "", output_file = NA) {
        # 0. load packages
        load_packages_2.3()
        cat("get vcf information for SNPs ... \n")
        
        # 0. generate output file name
        if (is.na(output_file)) {
                output_file = gen_output_file_vcfInfo(snp_info_file = snp_info_file, 
                                                      output_dir = output_dir, 
                                                      sample_name = sample_name)
        }
        cat("    output file name:", output_file, '\n')
        # 1. read snp file into data frame
        snp_info_list = read_inputSNP_file(snp_info_file)
        
        # 2. read vcf file into a list
        vcf_list = gen_vcf_list(vcf_dir = vcf_dir, 
                                vcf_file = vcf_file, 
                                snp_info_gr = snp_info_list$snp_info_gr)

        # 3. process vcf data
        vcf_df_list = process_vcf_list(vcf_list)
        
        # 4. add overlapping info to original snp table
        snp_info_addVcf_df = add_vcf_info (snp_info_df = snp_info_list$snp_info_df, 
                                           vcf_df_list = vcf_df_list)
        
        if (output_file  != F) { # if output_file == F, do not write down file
                write.csv0(snp_info_addVcf_df, output_file)
        }
        
        cat("vcf information added ... \n")
        return(list(snp_info_addVcf_df = snp_info_addVcf_df, output_file = output_file))
}

# Test #
# snp_file_loc = "LUC_Index_SNPs_20160607_riskPop_0.5.csv"
# get_vcf_info_main(snp_info_file = snp_file_loc, vcf_dir = "./data/samples/DDBJ_A549/vcf_files/")



