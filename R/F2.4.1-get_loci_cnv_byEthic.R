
# 1. prepare SNP files --------------------------------------------------------------------------------------------

# 1-1 prepare SNP files for cnv inference
prepare_files_for_cnvInfer = function(index_snp_file, sample_ethic = "EUR", output_dir = "./", sample_name = "",
                                      vcf_file = "data/samples/DDBJ_A549/vcf_files/A549_snv.GATK.formatcor.vcf") {
        # Aim: to generate wgs info added snp file
        # 1. generate ld snp info df from index snp file
        # 2. add wgs info
        # Input: 1. index snp file ; 2. the ethnicity of the sample (not the SNPs)
        
        # Test #
        # index_snp_file = "./data/input_snps/LUC_Index_SNPs_20160607.csv"
        # sample_ethic = "EUR"
        
        source("./src/scripts/F1.1-get_ldsnp_info.R")
        source("./src/scripts/F2.3-get_vcf_info.R")
        
        ldsnp_info_list = get_ldsnp_info_main(index_snp_file = index_snp_file, 
                                              population = sample_ethic, 
                                              output_dir = output_dir,
                                              r2_cutoff = 0.2, 
                                              for_cnv_call = T)
        ldsnp_info_addVcf_list = get_vcf_info_main(snp_info_file = ldsnp_info_list$output_file, 
                                                   vcf_file = vcf_file, 
                                                   output_dir = output_dir, 
                                                   sample_name = sample_name)
        
        return (ldsnp_info_addVcf_list$snp_info_addVcf_df)
        
}


# 1-2. modify vcf information for further analysis
mod_vcf_info = function(snp_info_addVcf_df) {
        # Aim: to modify vcf df to make it more suitable for further analysis
        
        # Test #
        # snp_info_addVcf_df = read.csv("./LUC_Index_SNPs_20160607_cnv_EUR_0.2_genotypeInfo.csv")
        
        snp_info_addVcf_df$allele_1_count = NA
        snp_info_addVcf_df$allele_2_count = NA
        
        # 1. change column names
        names(snp_info_addVcf_df)[grepl("ref_count", names(snp_info_addVcf_df))] = "ref_count"
        names(snp_info_addVcf_df)[grepl("alt_count", names(snp_info_addVcf_df))] = "alt_count"
        # head(snp_info_addVcf_df)
        
        # 2. sort df by query snp
        snp_info_addVcf_df = arrange(snp_info_addVcf_df, query_snp)
        # head(snp_info_addVcf_df)
        
        # 3. modify vcf information by D'
        D_prime_sign = as.numeric(as.character(snp_info_addVcf_df$D.)) > 0
        snp_info_addVcf_df$allele_1_count[D_prime_sign] = snp_info_addVcf_df$ref_count[D_prime_sign]
        snp_info_addVcf_df$allele_2_count[D_prime_sign] = snp_info_addVcf_df$alt_count[D_prime_sign]
        
        snp_info_addVcf_df$allele_1_count[!D_prime_sign] = snp_info_addVcf_df$alt_count[!D_prime_sign]
        snp_info_addVcf_df$allele_2_count[!D_prime_sign] = snp_info_addVcf_df$ref_count[!D_prime_sign]
        
        return(snp_info_addVcf_df)
}


# 2. get loci cnv raw ---------------------------------------------------------------------------------------------
# to source the "get_loci_cnv_raw" function
source("./src/scripts/F2.4.1.1-get_loci_cnv_raw.R")

# get CNV info ---------------------------------------------------------------------------------------------------

get_loci_cnv_byEthic = function(index_snp_file, cnv_param_list, 
                                sample_name = "", sample_ethic = "EUR", output_dir = "./") {
# Aim: to get loci cnv information of a particular ethic group
# Input: index snp file, vcf file, ethic
# Output: het_snp_summary_df_ethic (summarized info of loci cnv based on a particular ethic group)
        
        # Test #
        # index_snp_file = "./data/input_snps/LUC_Index_SNPs_20160607.csv"
        
        # 1. process snp data
        ## 1) get ld snp table; 2) add vcf information
        snp_info_addVcf_df_raw = prepare_files_for_cnvInfer(index_snp_file = index_snp_file, 
                                                            sample_name = sample_name,
                                                            sample_ethic = sample_ethic, 
                                                            output_dir = output_dir,
                                                            vcf_file = cnv_param_list$vcf_file)
        ## modify vcf information for further analysis
        snp_info_addVcf_df = mod_vcf_info(snp_info_addVcf_df_raw)
        
        # 2. get loci cnv by calculating r2(0.5) and 0.2 SNP distribution, then merge the information
        het_snp_summary_df_0.2 = get_loci_cnv_raw(snp_info_addVcf_df = snp_info_addVcf_df, 
                                                  r2_cutoff = 0.2, 
                                                  sample_ethic = sample_ethic, 
                                                  min_ldsnp_num = cnv_param_list$min_ldsnp_num, 
                                                  read_count_cutoff = cnv_param_list$read_count_cutoff)
        # head(as.data.frame(het_snp_summary_df_0.2))
        het_snp_summary_df_r2 = get_loci_cnv_raw(snp_info_addVcf_df = snp_info_addVcf_df, 
                                                 sample_ethic = sample_ethic,
                                                 r2_cutoff = cnv_param_list$r2_cutoff,
                                                 min_ldsnp_num = cnv_param_list$min_ldsnp_num, 
                                                 read_count_cutoff = cnv_param_list$read_count_cutoff)
        
        het_snp_summary_df_ethic = merge_het_snp_summary_df(list(het_snp_summary_df_r2, het_snp_summary_df_0.2))
        
        return(het_snp_summary_df_ethic)
        
}
