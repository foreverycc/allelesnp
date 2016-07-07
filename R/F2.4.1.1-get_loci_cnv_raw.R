
# 1. infer loci CNV -----------------------------------------------------------------------------------------------
calculate_loci_cnv = function(snp_info_addVcf_df, sample_ethic = "EUR", r2_cutoff = 0.5, min_ldsnp_num = 3, read_count_cutoff = 500) {
# To calculate CNV number for each locus
# step 1: collect SNP distribution for each locus
# step 2: use the binomial test based algorithm to infer CNV number
# Input: min_ldsnp_num - the min number of LD-SNPs that have wgs count info for one index SNP
# read_count_cutoff - the min reads of the sum of wgs counts for one index SNP
        
        # ----- Test ----- #
        # snp_info_addVcf_df = snp_info_addVcf_df
        # read_count_cutoff = 200
        # ---------------- #
        
        # calculate het-SNP count by LD
        snp_info_addVcf_df_narm = snp_info_addVcf_df[complete.cases(snp_info_addVcf_df) & 
                                                             snp_info_addVcf_df$r2 >= r2_cutoff, ]
        het_snp_summary_df = snp_info_addVcf_df_narm %>% group_by(query_snp) %>% 
                summarise(sum(allele_1_count), sum(allele_2_count))
        colnames(het_snp_summary_df) = c("query_snp", "sum_allele_1_count", "sum_allele_2_count")
        # filter low count het SNPs
        sel_high_count_snps = table(snp_info_addVcf_df_narm$query_snp) > min_ldsnp_num # an input_SNP must have >= 3 het SNPs in LD
        sel_query_snps = names(table(snp_info_addVcf_df_narm$query_snp)[sel_high_count_snps])
        het_snp_summary_df = filter(het_snp_summary_df, 
                                    sum_allele_1_count + sum_allele_2_count > read_count_cutoff &
                                            query_snp %in% sel_query_snps)
        
        # infer CNV
        het_snp_summary_df$ethic = sample_ethic
        het_snp_summary_df$r2 = r2_cutoff
        het_snp_summary_df$allele_1_cnv_count_bootstrap = NA
        het_snp_summary_df$allele_2_cnv_count_bootstrap = NA
        
        ratio_mat = gen_ratio_mat()
        for (i in 1:nrow(het_snp_summary_df)) {
                # ----- debug ----- #
                # print (i)
                # ----------------- #
                snp_info_addVcf_df_narm_selsnp = filter(snp_info_addVcf_df_narm, query_snp == as.character(het_snp_summary_df$query_snp)[i]) 
                allele_count_vec = infer_cnv(snp_info_addVcf_df_narm_selsnp$allele_1_count, snp_info_addVcf_df_narm_selsnp$allele_2_count, ratio_mat)
                
                het_snp_summary_df$allele_1_cnv_count_bootstrap[i] = allele_count_vec[1]
                het_snp_summary_df$allele_2_cnv_count_bootstrap[i] = allele_count_vec[2]
        }
        
        return (het_snp_summary_df)
}

# 1-1. helper function: get ratio matrix
gen_ratio_mat = function(allele_1_num_max = 10, allele_2_num_max = 10) {
# Aim: to generate (copy number) ratio (ref vs alt) matrix 
        ratio_mat = matrix(nrow = allele_1_num_max, ncol = allele_2_num_max)
        rownames(ratio_mat) = 1:allele_1_num_max
        colnames(ratio_mat) = 1:allele_2_num_max
        
        for (allele_1_num in 1:10) {
                for (allele_2_num in 1:10) {
                        ratio_mat[allele_1_num, allele_2_num] = allele_1_num/(allele_1_num + allele_2_num)
                }
        }    
        return (ratio_mat)
}

# 1-2. helper function: to use bootstrap distribution to model the distribution of (sum #ref alleles)/(sum #all alleles)
infer_cnv = function(allele_1_count, allele_2_count, ratio_mat) {
# algorithm to infer CNV
# algorithm is based on statistical model confidence interval
# models include: bootstrap, nb, pois, binom
        
        # ----- Test ----- #
        # allele_1_count = rpois(10, 50)
        # allele_2_count = rpois(10, 30)
        # allele_1_count = rnbinom(n = 10, mu = 50, size = 10)
        # allele_2_count = rnbinom(n = 10, mu = 30, size = 10)
        # ---------------- #
        
        cnv_num_vec = rep(0, 2)
        names(cnv_num_vec) = paste0(rep("cnv", each=2), c("_ref", "_alt"))
        # get confidence interval using bootstrap method
        interval = conf_int_bootstrap(allele_1_count, allele_2_count)
        
        # get cnv number by checking the 1st ratio_matrix's row-col pair that is within the interval
        if (! any (ratio_mat > interval[1] & ratio_mat < interval[2])) {
                cnv_num = c(NA, NA)
        } else {
                cnv_num = which(ratio_mat > interval[1] & ratio_mat < interval[2], arr.ind=T)[1, ]
        }
        (cnv_num_vec[c(1,2)] = cnv_num)
        
        return (cnv_num_vec)
}

# 1-2-1. helper function
conf_int_bootstrap = function (allele_1_count, allele_2_count, alpha = 0.05) {
# Aim: to get confidence interval using bootstrap method
        sim_mean_vec = rep(0,1000)
        for (i in 1:1000 ) {
                sim_mean_vec[i] = sim_bootstrap(allele_1_count, allele_2_count)
        }
        quantile(sim_mean_vec, probs=c(alpha, 1-alpha))
}

# 1-2-1-1. helper function
sim_bootstrap = function(allele_1_count, allele_2_count) {
# Aim: to get count sum ratio using bootstrap method
        n = length(allele_1_count)
        idx_1 = sample(n, n, replace = T)
        idx_2 = sample(n, n, replace = T)
        count_sum_ref = sum (allele_1_count[idx_1])
        count_sum_alt = sum (allele_2_count[idx_2])
        count_sum_all = sum(count_sum_ref, count_sum_alt)
        return (count_sum_ref/count_sum_all)
}


# 2. add cnv info to snp info df ----------------------------------------------------------------------------------
add_cnv_info = function(snp_info_addVcf_df, het_snp_summary_df) {
# To add cnv information into snp_info_addVcf_df table
        
        # ----- Test ----- #
        # snp_info_addVcf_df = snp_info_df_temp
        # het_snp_summary_df = het_snp_summary_df
        # ---------------- #
        
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


# 3. improve cnv results ------------------------------------------------------------------------------------------
backward_cor_cnv = function(snp_info_addCnv_df, OR_cutoff= 0.05) {
        # To correct cnv information based on calculated CNV
        
        # calculate odds ratio using calculated CNV
        or_vec = apply(snp_info_addCnv_df, 1, function(x) {
                if (any(!complete.cases(x))) {
                        return (NA)
                }
                p = as.numeric(x["ref_cnv"])/(as.numeric(x["ref_cnv"]) + as.numeric(x["alt_cnv"]))
                if (as.numeric(x["D."]) < 0 ) {
                        p = as.numeric(x["alt_cnv"])/(as.numeric(x["ref_cnv"]) + as.numeric(x["alt_cnv"]))
                }
                a = as.numeric(x["allele_1_count"])
                b = as.numeric(x["allele_2_count"])
                OR(p, a, b)
        })
        
        snp_info_addCnv_df$OR = or_vec
        temp = snp_info_addCnv_df
        sel_OR = snp_info_addCnv_df$OR < OR_cutoff & !is.na(snp_info_addCnv_df$OR)
        snp_info_addCnv_df$allele_1_count[sel_OR] = temp$allele_2_count[sel_OR]
        snp_info_addCnv_df$allele_2_count[sel_OR] = temp$allele_1_count[sel_OR]
        
        return (list(or_vec = or_vec, snp_info_addCnv_df = snp_info_addCnv_df))
}

# 3-1. helper function: to calculate odds ratio
OR = function(p, a, b) {
        # Function to calculate odds ratio
        q = 1-p
        p^a * q^b / (p^b * q^a)
}

# get loci cnv raw value ------------------------------------------------------------------------------------------

get_loci_cnv_raw = function(snp_info_addVcf_df, sample_ethic = "EUR", 
                            r2_cutoff = 0.5, min_ldsnp_num = 3, read_count_cutoff = 500) {
# Aim: To calculate CNV at loci level
# Input: snp table with vcf (wgs) info added
# Output: het_snp_summary_df (the cnv info summarized by loci)
        
        # calculate cnv based on summing counts of heterozygous SNPs in each loci
        while(T) {
                # 1. calculate loci cnv
                het_snp_summary_df = calculate_loci_cnv(snp_info_addVcf_df, sample_ethic, r2_cutoff = r2_cutoff, 
                                                        min_ldsnp_num = min_ldsnp_num, read_count_cutoff = read_count_cutoff)
                # 2. add cnv information to the snp info table
                snp_info_addCnv_df = add_cnv_info(snp_info_addVcf_df, het_snp_summary_df)
                head(snp_info_addCnv_df)
                # 3. improve results: 
                # if the allele distribution does not match the calculated CNV, switch counts for each allele
                cnv_cor = backward_cor_cnv(snp_info_addCnv_df)
                snp_info_addVcf_df = cnv_cor$snp_info_addCnv_df
                # if all heterozygous SNPs allele distribution match the calculated CNV, break the loop
                print (sum(cnv_cor$or_vec[!is.na(cnv_cor$or_vec)] < 0.05))
                if (sum(cnv_cor$or_vec[!is.na(cnv_cor$or_vec)] < 0.05) == 0) break
        }
        
        return (het_snp_summary_df)
}
