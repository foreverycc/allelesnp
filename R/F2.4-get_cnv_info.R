# F2.4-get cnv info
## ==================================================================================
## Function: get_cnv_info_main
###  Aim: To get the cnv of input SNPs by vcf/encode cnv data files
###  INPUT: input snp file (.csv), vcf/encode cnv data
###  OUTPUT: input snp file with information of cnv
## ==================================================================================

# outline ---------------------------------------------------------------------------------------------------------

# # 1. read snp file into data frame
# # snp_list should contain both the df and GRange object
# snp_list = read_snpFile()
# 
# # 2. read vcf file 
# vcf_info_list = read_vcfFile(vcf_file_for_cnv, snp_info_gr)
# 
# # 3. add overlapping info to original snp table
# snp_df_proc = add_genotype_info(snp_info_df, vcf_list)

# 0-1. load packages ----------------------------------------------------------------------------------------------

load_packages_2.4 = function(){
        packages = c('dplyr', 'VariantAnnotation')
        load = lapply(packages, require, character.only = T)
}

# function: write.csv0
source("./src/scripts/T1-toolbox.R")

# 1. some helper functions ----------------------------------------------------------------------------------------
# generate output file name
gen_output_file_cnvInfo = function(snp_info_file, output_dir, sample_name) {
        # Aim: to get output file name        
        snp_info_vec = strsplit(snp_info_file, "/")[[1]]
        snp_batch_id = gsub(".csv", "", snp_info_vec[length(snp_info_vec)])
        output_file = paste0(output_dir, "/", snp_batch_id, "_", sample_name, "_cnvInfo.csv")
        
        return (output_file)
        
}

# get other ethic groups
get_other_ethic = function(ethic) {
# Aim: to get teh correct order of ethic groups
        if (ethic == "EUR") return (c("ASN", "AFR"))
        if (ethic == "ASN") return (c("EUR", "AFR"))
        if (ethic == "AFR") return (c("EUR", "ASN"))
}

# merge het snp summary df together in order
merge_het_snp_summary_df = function (het_snp_summary_list) {
# To merge het_snp_summary_df files based on various r2 cutoff
        
        for (i in 1:length(het_snp_summary_list)) {
                if (i == 1) {
                        het_snp_summary_df = het_snp_summary_list[[1]]
                } else {
                        to_add = filter(het_snp_summary_list[[i]], !(query_snp %in% as.character(het_snp_summary_df$query_snp)))
                        het_snp_summary_df = rbind(het_snp_summary_df, to_add)
                }
        }
        
        return (het_snp_summary_df)
}


# 2. get loci cnv  ------------------------------------------------------------------------------------------------

# to source the function: get_loci_cnv_byEthic
source("./src/scripts/F2.4.1-get_loci_cnv_byEthic.R")

# 3. add cnv information on risk SNPs --------------------------------------------------------------------------------

add_cnv_info_riskSnp = function(snp_info_addVcf_df, het_snp_summary_df) {
# Aim: to add cnv information to risk snp table
        
        # Test #
        # snp_info_addVcf_df = ldsnp_info_vcf_list$snp_info_addVcf_df
        
        # modify data: make r2, D., population into single string
        r2 = as.numeric(sapply(as.character(snp_info_addVcf_df$r2), function(x) strsplit(x, ",")[[1]][1]))
        D. = as.numeric(sapply(as.character(snp_info_addVcf_df$D.), function(x) strsplit(x, ",")[[1]][1]))
        query_snp = sapply(as.character(snp_info_addVcf_df$query_snp), function(x) strsplit(x, ",")[[1]][1])
        population = sapply(as.character(snp_info_addVcf_df$population), function(x) strsplit(x, ",")[[1]][1])
        
        snp_info_addVcf_df$r2 = r2
        snp_info_addVcf_df$D. = D.
        snp_info_addVcf_df$query_snp = query_snp
        snp_info_addVcf_df$population = population
        
        snp_info_addVcf_df = mod_vcf_info(snp_info_addVcf_df)
        head(snp_info_addVcf_df)        
        
        # add bootstrap infered cnv 
        snp_info_addVcf_df$ref_cnv = NA
        snp_info_addVcf_df$alt_cnv = NA
        
        for (query_snp in as.character(unique(snp_info_addVcf_df$query_snp))) {
                if (query_snp %in% as.character(het_snp_summary_df$query_snp)) {
                        snp_info_addVcf_df[snp_info_addVcf_df$query_snp == query_snp & snp_info_addVcf_df$D. > 0, 
                                           c("ref_cnv", "alt_cnv")] = 
                                het_snp_summary_df[het_snp_summary_df$query_snp == query_snp, 
                                                   c("allele_1_cnv_count_bootstrap", "allele_2_cnv_count_bootstrap")]
                        # be aware of the order here!!!
                        snp_info_addVcf_df[snp_info_addVcf_df$query_snp == query_snp & snp_info_addVcf_df$D. < 0, 
                                           c("ref_cnv", "alt_cnv")] = 
                                het_snp_summary_df[het_snp_summary_df$query_snp == query_snp,
                                                   c("allele_2_cnv_count_bootstrap", "allele_1_cnv_count_bootstrap")]
                }
        }
        
        return(snp_info_addVcf_df)
        
}

# add cnv info main -----------------------------------------------------------------------------------------------

get_cnv_info_main = function(index_snp_file = "./data/input_snps/LUC_Index_SNPs_20160607.csv", snp_info_file = NA,
                             sample_name = "DDBJ_A549", sample_ethic = "EUR", output_dir = "./", output_file = NA,
                             vcf_file_for_cnv = "data/samples/DDBJ_A549/vcf_files/A549_snv.GATK.formatcor.vcf",
                             r2_cutoff = 0.5, min_ldsnp_num = 3, read_count_cutoff = 500) {
        # Test #
        # index_snp_file = "./data/input_snps/LUC_Index_SNPs_20160607.csv"
        # sample_name = "DDBJ_A549"
        # sample_ethic = "EUR"
        # output_dir = "./"
        # output_file = NA
        # vcf_file_for_cnv = "data/samples/DDBJ_A549/vcf_files/A549_snv.GATK.formatcor.vcf"
        # r2_cutoff = 0.5
        # min_ldsnp_num = 3
        # read_count_cutoff = 500
        
        # 0. preparation
        set.seed(1234)
        load_packages_2.4()
        cat("get cnv information for SNPs ... \n")
        cnv_param_list = list(vcf_file_for_cnv = vcf_file_for_cnv,
                              r2_cutoff = r2_cutoff, 
                              min_ldsnp_num = min_ldsnp_num, 
                              read_count_cutoff = read_count_cutoff)
        
        # 1. calculate cnv by loci and ethic group, and then merge them together
        ethics = c(sample_ethic, get_other_ethic(sample_ethic))
        het_snp_summary_list = replicate(length(ethics), list())
        for (i in 1:length(ethics)) {
                het_snp_summary_list[[i]] = get_loci_cnv_byEthic(index_snp_file = index_snp_file, 
                                                                 sample_name = sample_name, 
                                                                 sample_ethic = ethics[i], 
                                                                 output_dir = output_dir,
                                                                 cnv_param_list = cnv_param_list)
        }
        
        het_snp_summary_df_final = merge_het_snp_summary_df(het_snp_summary_list)
        
        # 2. get and add cnv info to risk SNP info table
        if (is.na(snp_info_file)) {
                ldsnp_info_list = get_ldsnp_info_main(index_snp_file = index_snp_file, 
                                                      output_dir = output_dir)
                snp_info_file = ldsnp_info_list$output_file
        }
        
        ldsnp_info_vcf_list = get_vcf_info_main(snp_info_file = snp_info_file, 
                                                vcf_file = vcf_file_for_cnv, 
                                                output_file = F)
        snp_info_addCnv_df = add_cnv_info_riskSnp(snp_info_addVcf_df = ldsnp_info_vcf_list$snp_info_addVcf_df, 
                                                  het_snp_summary_df = het_snp_summary_df_final)
        
        # 3. write data
        if (is.na(output_file)) {
                output_file = gen_output_file_cnvInfo(snp_info_file = snp_info_file, 
                                                      output_dir = output_dir,
                                                      sample_name = sample_name)
        }
        cat("    output file name:", output_file, '\n')
        
        if (output_file != F) {
                write.csv0(snp_info_addCnv_df, output_file)
        }
        
        cat("cnv information added ... \n")
        
        return(list(snp_info_addCnv_df = snp_info_addCnv_df, output_file = output_file))
        
}


# Test #
# get_cnv_info_main(output_dir = output_dir)

