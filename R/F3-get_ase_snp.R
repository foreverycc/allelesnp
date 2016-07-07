# F3 get_ase_snp
## ==================================================================================
## Function: get ase snps
###  Aim: To get ase snps
###  INPUT: allele-dist snp file, peak annotated snp file, vcf annotated snp file, cnv annotated snp file
###  OUTPUT: ase snp file with peak, vcf and cnv information, and test p-value info added
## ==================================================================================


# 0. toolbox ------------------------------------------------------------------------------------------------------

# function: write.csv0
source("./src/scripts/T1-toolbox.R")

# output_file ----------------------------------------------------------------------------------------------------

gen_output_file_ase = function(snp_info_file, output_dir = "./") {
        # Aim: to get output file name        
        snp_info_vec = strsplit(snp_info_file, "/")[[1]]
        snp_batch_id = gsub(".csv", "", snp_info_vec[length(snp_info_vec)])
        output_file = paste0(output_dir, "/", snp_batch_id, "_ase_snp.csv")
        
        return (output_file)
        
}

# 1. infer genotype from bam files only ---------------------------------------------------------------------------

infer_genotype_from_reads = function(snp_info_alleleDist_df, het_threshold, genotype_by_sample = T) {
# Aim: infer genotypes from bam reads
# 2 types of inference: 1) single file 2) all files in a same sample together
        
        # get genotype by single file reads
        snp_info_alleleDist_df_sel = dplyr::select(snp_info_alleleDist_df, ref, alt)
        genotype_singleBam = apply(snp_info_alleleDist_df_sel, 1, function(x) {
                genotype_infer(as.numeric(x["ref"]), as.numeric(x["alt"]), het_threshold = het_threshold)
        })
        snp_info_alleleDist_df$genotype_singleBam = genotype_singleBam
        
        # get genotype by sample
        snp_info_alleleDist_df_sel2 = dplyr::select(snp_info_alleleDist_df, rsID, ref, alt)
        snp_info_alleleDist_df_summ = 
                snp_info_alleleDist_df_sel2 %>% group_by(rsID) %>% summarise(ref = sum(ref), alt = sum(alt))
        snp_info_alleleDist_df_summ_sel = dplyr::select(snp_info_alleleDist_df_summ, ref, alt)
        ## infer genotype by sample
        if (genotype_by_sample) {
                genotype_sample = apply(snp_info_alleleDist_df_summ[, c("ref", "alt")], 1, function(x) {
                        genotype_infer(as.numeric(x["ref"]), as.numeric(x["alt"]), het_threshold = het_threshold)
                })
                het_snps = as.character(snp_info_alleleDist_df_summ$rsID[genotype_sample])
                
                snp_info_alleleDist_df$genotype_sample = F
                snp_info_alleleDist_df$genotype_sample [snp_info_alleleDist_df$rsID %in% het_snps] = T
        } else {
                snp_info_alleleDist_df$genotype_sample = NA
        }
        # head(snp_info_alleleDist_df)
        return (snp_info_alleleDist_df)
}

# 1-1. helper function: infer genotype based on ref vs alt reads
genotype_infer = function(ref, alt, het_threshold = 0.1) {
# Aim: to infer genotype based on reads counts that contain ref/alt
# Input: ref_v_alt_vec, het_threshold (0 ~ 1, default = 0.1 )
# Output: genotype ("T" (het) or "F" (homozygous) or NA)
        
        ref_v_alt_vec = c(ref, alt)
        
        if (any (is.na(ref_v_alt_vec))) {
                return (NA)
        } else if (het_threshold < 1 & het_threshold > 0) {
                if (all (max (ceiling(ref_v_alt_vec * het_threshold)) < ref_v_alt_vec)) {
                        return (T)
                        break() 
                } else {
                        return (F)
                        break()
                }
        } else return (NA)
}


# 2. add vcf info --------------------------------------------------------------------------------------------------

add_vcf_res = function(snp_info_alleleDist_df, snp_info_vcf_file) {
# Aim: to add vcf information 
        
        if (is.na(snp_info_vcf_file)) {
                # add column "genotype vcf"
                snp_info_alleleDist_df$genotype_vcf = NA
                # get final genotype
                sel_genotype_cols = grepl("^genotype", names(snp_info_alleleDist_df))
                genotype_final = apply(snp_info_alleleDist_df[, sel_genotype_cols], 1, function(x) any(x, na.rm = T))
                snp_info_alleleDist_df$genotype_final = genotype_final
                
                return(snp_info_alleleDist_df)
        }
        
        snp_info_vcf_df = read.csv(snp_info_vcf_file)
        rownames(snp_info_vcf_df) = as.character(snp_info_vcf_df$rsID)
        # head(snp_info_vcf_df)
        
        # add vcf info
        snp_info_vcf_df_sel = snp_info_vcf_df[, grepl(".vcf$", names(snp_info_vcf_df))]
        snp_info_vcf_df_sel$genotype_vcf = apply(snp_info_vcf_df_sel, 1, function(x) any(x, na.rm = T))
        
        sel_vcf_cols = snp_info_vcf_df_sel[as.character(snp_info_alleleDist_df$rsID), "genotype_vcf"]
        snp_info_alleleDist_df = cbind(snp_info_alleleDist_df, "genotype_vcf" = sel_vcf_cols)
        # head(snp_info_alleleDist_df)
        
        # calculate final genotype
        sel_genotype_cols = grepl("^genotype", names(snp_info_alleleDist_df))
        genotype_final = apply(snp_info_alleleDist_df[, sel_genotype_cols], 1, any)
        snp_info_alleleDist_df$genotype_final = genotype_final
        
        rownames(snp_info_alleleDist_df) = 1:nrow(snp_info_alleleDist_df)
        return (snp_info_alleleDist_df)
}


# 3. add peak info ------------------------------------------------------------------------------------------------
add_peak_res = function(snp_info_alleleDist_df, snp_info_peak_file) {
# Aim: to add peak information 
        if (is.na(snp_info_peak_file)) {
                snp_info_alleleDist_df$biofeature_overlap = snp_info_alleleDist_df$biofeature_overlap_num = 
                        snp_info_alleleDist_df$biofeature_overlap_names = NA
                return(snp_info_alleleDist_df)
        }
        
        snp_info_peak_df = read.csv(snp_info_peak_file)
        rownames(snp_info_peak_df) = as.character(snp_info_peak_df$rsID)
        # head(snp_info_peak_df)
        snp_info_peak_df_sel = snp_info_peak_df[, grepl("biofeature", names(snp_info_peak_df))]
        # head(snp_info_peak_df_sel)
        
        # add peak annotation results
        sel_peak_cols = snp_info_peak_df_sel[as.character(snp_info_alleleDist_df$rsID), ]
        snp_info_alleleDist_df = cbind(snp_info_alleleDist_df, sel_peak_cols)        
        head(snp_info_alleleDist_df)
        rownames(snp_info_alleleDist_df) = 1:nrow(snp_info_alleleDist_df)
        
        return (snp_info_alleleDist_df)
}

# 4. add cnv info ------------------------------------------------------------------------------------------------
# 4.1. for local data
add_cnv_res = function(snp_info_alleleDist_df, snp_info_cnv_file) {
# Aim: to add vcf information 
        
        if (is.na(snp_info_cnv_file)) {
                snp_info_alleleDist_df$alt_count = snp_info_alleleDist_df$ref_count = NA
                snp_info_alleleDist_df$alt_cnv = snp_info_alleleDist_df$ref_cnv = 1
                return (snp_info_alleleDist_df)
        }
        
        snp_info_cnv_df = read.csv(snp_info_cnv_file)
        
        rownames(snp_info_cnv_df) = as.character(snp_info_cnv_df$rsID)
        head(snp_info_cnv_df)
        snp_info_cnv_df_sel = snp_info_cnv_df[, c("ref_count", "alt_count", "ref_cnv", "alt_cnv")]
        head(snp_info_cnv_df_sel)
        
        # add cnv results
        sel_cnv_cols = snp_info_cnv_df_sel[as.character(snp_info_alleleDist_df$rsID), ]
        snp_info_alleleDist_df = cbind(snp_info_alleleDist_df, sel_cnv_cols)        
        head(snp_info_alleleDist_df)
        rownames(snp_info_alleleDist_df) = 1:nrow(snp_info_alleleDist_df)
        return (snp_info_alleleDist_df)
}

# 4.2. for encode database only 
add_cnv_res_encode = function(snp_info_alleleDist_df, snp_info_cnv_file) {
# Aim: to add encode-cnv information 

        # Test #
        # snp_info_cnv_file = "LUC_Index+LD_SNPs_20160607_ENCODE_cnvInfo.csv"

        if (is.na(snp_info_cnv_file)) {
                snp_info_alleleDist_df$alt_count = snp_info_alleleDist_df$ref_count = NA
                snp_info_alleleDist_df$alt_cnv = snp_info_alleleDist_df$ref_cnv = 1
                return (snp_info_alleleDist_df)
        }
        
        snp_info_cnv_df = read.csv(snp_info_cnv_file)
        rownames(snp_info_cnv_df) = as.character(snp_info_cnv_df$rsID)
        head(snp_info_cnv_df)
        
        samples = names(snp_info_cnv_df)[which(names(snp_info_cnv_df) == "A549_Rep1"):ncol(snp_info_cnv_df)]
        sample_info_df = as.data.frame(sapply(samples, function(x) unlist(strsplit(x, "_"))))
        
        snp_info_alleleDist_df$encode_cnv_sample = NA
        snp_info_alleleDist_df$encode_cnv_info = NA
        
        for (i in 1:nrow(snp_info_alleleDist_df)) {
                snp_i = as.character(snp_info_alleleDist_df$rsID[i])
                biofeature_i =as.character(snp_info_alleleDist_df$biofeature[i])
                sel = unlist(apply(sample_info_df, 2, function(x) grep (x[1] , biofeature_i)))
                
                if (length(sel) == 0) {
                        snp_cnv_char = NA
                        snp_cnv_sample_char = NA
                        next
                }
                
                snp_cnv_sample_char = paste(names(sel), collapse = ",")
                snp_cnv_info_char = paste(snp_info_cnv_df[snp_i, names(sel)], collapse = ",")

                snp_info_alleleDist_df$encode_cnv_sample[i] = snp_cnv_sample_char
                snp_info_alleleDist_df$encode_cnv_info[i] = snp_cnv_info_char
        }
        
        return (snp_info_alleleDist_df)        
}


# 5. calculate allele-specific effects ----------------------------------------------------------------------------
# 5.1 for local data
calculate_ase = function(snp_info_alleleDist_df, depth_threshold){
        
        # filter snps
        snp_info_alleleDist_df = filter(snp_info_alleleDist_df, genotype_final, ref + alt > depth_threshold)
        if (nrow(snp_info_alleleDist_df) == 0) {
                cat("No allele-specific effects are identified.\n")
                return ()
        }
        
        # raw p-value (binomial p = 0.5)
        p.val.raw = binomial_test(snp_info_alleleDist_df$ref, snp_info_alleleDist_df$alt)
        snp_info_alleleDist_df$p.val.raw = p.val.raw
        
        # p-val (wgs corrected)
        if (all(is.na(snp_info_alleleDist_df$ref_count)) & all(is.na(snp_info_alleleDist_df$alt_count))) {
                snp_info_alleleDist_df$p.val.wgs = NA
        } else {
                binom_p_wgs = snp_info_alleleDist_df$ref_count / 
                        (snp_info_alleleDist_df$ref_count + snp_info_alleleDist_df$alt_count)
                p.val.wgs = binomial_test(snp_info_alleleDist_df$ref, snp_info_alleleDist_df$alt, binom_p_wgs)
                snp_info_alleleDist_df$p.val.wgs = p.val.wgs
        }
        
        # p-val (cnv corrected) 
        binom_p_cnv = snp_info_alleleDist_df$ref_cnv / 
                (snp_info_alleleDist_df$ref_cnv + snp_info_alleleDist_df$alt_cnv)
        p.val.cnv = binomial_test(snp_info_alleleDist_df$ref, snp_info_alleleDist_df$alt, binom_p_cnv)
        snp_info_alleleDist_df$p.val.cnv = p.val.cnv
        
        
        p.val.cnv.bh = p.adjust(p.val.cnv, method = "BH")        
        p.val.cnv.bonf = p.adjust(p.val.cnv, method = "bonferroni")        
        snp_info_alleleDist_df$p.val.cnv.bh = p.val.cnv.bh
        snp_info_alleleDist_df$p.val.cnv.bonf = p.val.cnv.bonf
        
        snp_info_alleleDist_df = arrange(snp_info_alleleDist_df, p.val.cnv)
        
        return (snp_info_alleleDist_df)
}

# 5.2 for encode database only
calculate_ase_encode = function(snp_info_alleleDist_df, depth_threshold) {
        snp_info_alleleDist_df = filter(snp_info_alleleDist_df, genotype_final, ref + alt > depth_threshold)
        if (nrow(snp_info_alleleDist_df) == 0) {
                cat("No allele-specific effects are identified.\n")
                return ()
        }
        
        p.val.raw = binomial_test(snp_info_alleleDist_df$ref, snp_info_alleleDist_df$alt)
        snp_info_alleleDist_df$p.val.raw = p.val.raw
        
        p.val.bh = p.adjust(p.val.raw, method = "BH")        
        p.val.bonf = p.adjust(p.val.raw, method = "bonferroni")        
        snp_info_alleleDist_df$p.val.bh = p.val.bh
        snp_info_alleleDist_df$p.val.bonf = p.val.bonf
        
        return(snp_info_alleleDist_df)
}

# Function: vectorized binom test
binomial_test = function(x, y, p=0.5) {
        xyp = cbind(x, y, p)
        apply(xyp, 1, function(x){
                if (any(is.na(x))) {
                        return (NA)
                }
                if (x[1] == 0 & x[2] == 0) {
                        return (NA) 
                } else {
                        return(binom.test(x[1], x[2]+ x[1], x[3])$p.value)
                }
        })
}

# main function ---------------------------------------------------------------------------------------------------

get_ase_snp_main = function(snp_info_alleleDist_file, 
                            snp_info_peak_file = NA, 
                            snp_info_vcf_file = NA, 
                            snp_info_cnv_file = NA,
                            het_threshold = 0.1, 
                            depth_threshold = 20,
                            genotype_by_sample = T, 
                            output_dir = "./", 
                            output_file = NA,
                            use_encode_cnv = F
                            ) {
        
        # Test #
        # snp_info_alleleDist_file = "./tests/test6/LUC_Index_SNPs_20160607_short_ENCODE_DNase_assnp/LUC_Index_SNPs_20160607_short_riskPop_0.5_ENCODE_DNase_alleleDist.csv"
        # snp_info_peak_file = "./LUC_Index_SNPs_20160607_DDBJ_A549_assnp/LUC_Index_SNPs_20160607_riskPop_0.5_DDBJ_A549_peakAnnotation.csv"
        # snp_info_vcf_file = "./LUC_Index_SNPs_20160607_DDBJ_A549_assnp/LUC_Index_SNPs_20160607_riskPop_0.5_DDBJ_A549_genotypeInfo.csv"
        # snp_info_cnv_file = "./LUC_Index_SNPs_20160607_DDBJ_A549_assnp/LUC_Index_SNPs_20160607_cnv_riskPop_0.5_DDBJ_A549_cnvInfo.csv"
        # het_threshold = 0.1
        # depth_threshold = 20
        # output_dir = "./LUC_Index_SNPs_20160607_DDBJ_A549_assnp/"
        # output_file = NA

        require("dplyr")
        cat("calculate allele-specific effects for SNPs ... \n")
        
        if (is.na(output_file)) {
                output_file = gen_output_file_ase(snp_info_file = snp_info_alleleDist_file, 
                                                  output_dir = output_dir)
        }
        cat("    output file name:", output_file, '\n')
        
        snp_info_alleleDist_df = read.csv(snp_info_alleleDist_file)
        head(snp_info_alleleDist_df)
        # 1. get genotype from allelic distribution
        snp_info_alleleDist_df = infer_genotype_from_reads(snp_info_alleleDist_df = snp_info_alleleDist_df, 
                                                           het_threshold = het_threshold, 
                                                           genotype_by_sample = genotype_by_sample)
        
        # 2. add vcf information to the table
        snp_info_alleleDist_df = add_vcf_res(snp_info_alleleDist_df = snp_info_alleleDist_df, 
                                             snp_info_vcf_file = snp_info_vcf_file)
        
        # 3. add peak information to the table
        snp_info_alleleDist_df = add_peak_res(snp_info_alleleDist_df = snp_info_alleleDist_df, 
                                              snp_info_peak_file = snp_info_peak_file)
                        
        # 4. add cnv information to the table
        if (use_encode_cnv) {
                snp_info_alleleDist_df = add_cnv_res_encode(snp_info_alleleDist_df = snp_info_alleleDist_df, 
                                                            snp_info_cnv_file = snp_info_cnv_file)
        } else {
                snp_info_alleleDist_df = add_cnv_res(snp_info_alleleDist_df = snp_info_alleleDist_df, 
                                                     snp_info_cnv_file = snp_info_cnv_file)
        }
        # head(snp_info_alleleDist_df)
        
        # 5. calculate allele-specificity
        if (use_encode_cnv) {
                snp_info_alleleDist_df = calculate_ase_encode(snp_info_alleleDist_df = snp_info_alleleDist_df, 
                                                              depth_threshold = depth_threshold) 
        } else {
                snp_info_alleleDist_df = calculate_ase(snp_info_alleleDist_df = snp_info_alleleDist_df, 
                                                       depth_threshold = depth_threshold) 
        }
        
        if (output_file != F) { # if output_file == F, do not write down file
                write.csv0(snp_info_alleleDist_df, output_file)
        }
        cat("allele-specific effects information added ... \n")
        
        return (snp_info_alleleDist_df)
}

# Test #
# get_ase_snp_main(snp_info_alleleDist_file, snp_info_peak_file, snp_info_vcf_file, snp_info_cnv_file)
