# M1-get_assnp
## ==================================================================================
## Function: get_ase_snp
###  Aim: To get allele-specific SNP 
###  INPUT 1: index SNP list (a .csv file), with the 1st col the rsID and the 2nd col Population 
###  INPUT 2: snp information table (a .csv file), contains column names with rsID and population 
###  OUTPUT: potential allele-specific effect (ase) snps
## ==================================================================================


# source scripts for use ------------------------------------------------------------------------------------------

for (script in list.files("./src/scripts/", pattern = "^[FT].*.R", full.names = T)) {
        cat ("source script: ", script, "\n")
        source(script)
}


# build main function ---------------------------------------------------------------------------------------------

get_assnp = function(index_snp_file = NA, 
                     snp_info_file = NA, 
                     sample_name = "DDBJ_A549", 
                     bam_dir = NA,
                     merge_replicates = T,
                     peak_dir = NA,
                     vcf_dir = NA,
                     vcf_file_for_cnv = NA, 
                     sample_ethic = "EUR",
                     het_threshold = 0.1,
                     depth_threshold = 20,
                     base_qual_threshold = 20,
                     mapq_threshold = 20,
                     genotype_by_sample = T, 
                     r2_cutoff = 0.5, 
                     min_ldsnp_num = 3, 
                     read_count_cutoff = 500,
                     use_encode_cnv = F, 
                     ...) {
        # Test #
        # index_snp_file = "./data/input_snps/LUC_Index_SNPs_20160607.csv"
        # snp_info_file = NA
        # sample_name = "DDBJ_A549"
        # sample_ethic = "EUR"
        # bam_dir = "./data/samples/DDBJ_A549/bam_files/"
        # peak_dir = NA
        # vcf_dir = NA
        # vcf_file_for_cnv = NA
        # het_threshold = 0.1
        # depth_threshold = 20
        # base_qual_threshold = 20
        # mapq_threshold = 20
        
        # 0. make a directory for data storation
        index_snp_file_vec = strsplit(index_snp_file, "/")[[1]]
        snp_batch_id = gsub(".csv", "", index_snp_file_vec[length(index_snp_file_vec)])
        
        dir.create(output_dir <- paste0("./", snp_batch_id, "_", sample_name, "_assnp"))
        
        # 1. get ld snp file
        if (!is.na(index_snp_file)) {
                snp_info_list = get_ldsnp_info_main(index_snp_file = index_snp_file, 
                                                    output_dir = output_dir)
                snp_info_file = snp_info_list$output_file
        } 
        
        # 2. get allele distribution
        snp_info_alleleDist_list = get_alleleDist_info_main(snp_info_file = snp_info_file, 
                                                            bam_dir = bam_dir, 
                                                            output_dir = output_dir, 
                                                            sample_name = sample_name,
                                                            base_qual_threshold = base_qual_threshold,
                                                            mapq_threshold = mapq_threshold, 
                                                            merge_replicates = merge_replicates, ...)
        
        # 3. get peak annotation
        if (!is.na(peak_dir)) {
                snp_info_addPeak_list = get_peak_info_main(snp_info_file = snp_info_file, 
                                                           peak_dir = peak_dir, 
                                                           output_dir = output_dir,
                                                           sample_name = sample_name)
        } else {
                snp_info_addPeak_list = list()
                snp_info_addPeak_list$output_file = NA
        }
        
        # 4. get vcf annotation
        if (!is.na(vcf_dir)) {
                snp_info_addVcf_list = get_vcf_info_main(snp_info_file = snp_info_file, 
                                                         vcf_dir = vcf_dir, 
                                                         output_dir = output_dir, 
                                                         sample_name = sample_name)
        } else {
                snp_info_addVcf_list = list()
                snp_info_addVcf_list$output_file = NA
        }
        
        # 5. get cnv annotation
        if (use_encode_cnv) {
                ## use encode cnv
                snp_info_addCnv_list = get_encodeCnv_info_main(snp_info_file = snp_info_file,
                                                                output_dir = output_dir,
                                                                sample_name = sample_name)
        } else if (!is.na(index_snp_file) & !is.na(vcf_file_for_cnv)) {
                ## if have index snp and vcf_file_for_cnv (wgs vcf data file)
                snp_info_addCnv_list = get_cnv_info_main(index_snp_file = index_snp_file, 
                                                         snp_info_file = snp_info_file, 
                                                         output_dir = output_dir, 
                                                         vcf_file_for_cnv = vcf_file_for_cnv,
                                                         sample_name = sample_name, 
                                                         sample_ethic = sample_ethic,
                                                         r2_cutoff = r2_cutoff, 
                                                         min_ldsnp_num = min_ldsnp_num, 
                                                         read_count_cutoff = read_count_cutoff)
        } else {
                snp_info_addCnv_list = list()
                snp_info_addCnv_list$output_file = NA 
        }
        
        # 6. calculate ase-snps
        get_ase_snp_main(snp_info_alleleDist_file = snp_info_alleleDist_list$output_file, 
                         snp_info_peak_file = snp_info_addPeak_list$output_file,
                         snp_info_vcf_file = snp_info_addVcf_list$output_file, 
                         snp_info_cnv_file = snp_info_addCnv_list$output_file,
                         het_threshold = het_threshold,
                         depth_threshold = depth_threshold,
                         output_dir = output_dir, 
                         genotype_by_sample = genotype_by_sample,
                         use_encode_cnv = use_encode_cnv)
        
}

