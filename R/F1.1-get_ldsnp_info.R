# F1.1-get ld snp information
## ==================================================================================
## Function: get_ldsnp_info_main
###  Aim: To get HaploReg files from an input SNP file which contains a few SNPs and its corresponding risky population 
###  INPUT: index SNP list (a .csv file), with the 1st col the rsID and the 2nd col Population 
###  OUTPUT: information table of index + LD SNP from HaploReg
## ==================================================================================

### In this version, I use python to webscrape information from HaploReg websites rather than manually get it ###


# outline ---------------------------------------------------------------------------------------------------------

# # read index SNP table
# indexSNP_df = read_indexSNP()
# populations = unique(indexSNP_df$population)
# pop_n = length(populations)
# 
# # find LD SNPs by population
# haploreg_list = replicate(pop_n, list())
# i = 1
# for (pop in populatinos) {
#     indexSNP_df_pop = filter(indexSNP_df, population == pop)
#     haploreg_df_pop = get_haploreg_res(indexSNP_df_pop$rsID)
#     haploreg_list[[i]] = haploreg_df_pop    
#     i = i + 1
# }
# 
# arrange_haploreg_list(haploreg_list)


# toolboxes -------------------------------------------------------------------------------------------------------

# Function to load necessary packages
load_packages_01 = function(){
        packages = c('rPython', 'dplyr', 'rtracklayer','GenomicRanges','FDb.UCSC.snp137common.hg19', "BBmisc")
        load = lapply(packages, require, character.only = T)
}

# function: write.csv0
source("./src/scripts/T1-toolbox.R")

# Function to read data from haploreg output
read.haploreg3 = function(file_loc) {
        return(read.table(file_loc,skip=1, header=T))
}

# Function to read index SNPs to data frame
read_indexSNP = function(file_loc) {
        df = read.csv(file_loc, header = F)
        names(df) = c("rsID", "population")
        return (df)
}


gen_output_file_ldsnp = function(index_snp_file, population, r2, for_cnv_call = F, output_dir = "./") {
        # Aim: to get output file name        
        index_snp_file_vec = strsplit(index_snp_file, "/")[[1]]
        snp_batch_id = gsub(".csv", "", index_snp_file_vec[length(index_snp_file_vec)])
        if (!is.na(population)) {
                if (for_cnv_call) {
                        output_file = paste0(output_dir, "/", snp_batch_id, "_cnv_", population, "_", r2, ".csv")
                } else {
                        output_file = paste0(output_dir, "/", snp_batch_id, "_", population, "_", r2, ".csv")
                }
        }
        if (is.na(population)) {
                output_file = paste0(output_dir, "/", snp_batch_id, "_riskPop_", r2, ".csv")
        }
        return (output_file)
        
}
# Function 1: get index+ld snp information -----------------------------------------------------------------------

get_ldsnp_info = function(index_snp_df_pop, r2_cutoff) {
# Aim: get haploreg results of index + LD snps
# Input: index_snp_df_pop (index SNP df of one population), r2_cutoff
# Output: index + ld snp haploreg results
# helper functions:
#         1. get_haploreg_info
#         2. add_query_snp_info

        # get index SNP and population information
        index_snps = as.character(index_snp_df_pop$rsID)
        population = as.character(index_snp_df_pop$population)[1]
        
        # get the information of the index SNPs from haploreg website
        index_snp_haploreg_df = get_haploreg_info(index_snps, population, r2_cutoff = 1.1)
        
        # sort index SNPs and submit to haploreg to get LD SNPs
        index_snp_df_sorted = arrange(index_snp_haploreg_df, desc(chr), desc(pos))
        index_snps_sorted = as.character(index_snp_df_sorted$rsID)
        ldsnp_haploreg_df = get_haploreg_info(index_snps_sorted, population, r2_cutoff)
        
        # add query snp information
        ldsnp_haploreg_df = add_query_snp_info(ldsnp_haploreg_df, index_snps_sorted, population)
}

# helper function 1-1: to get haploreg file
get_haploreg_info = function(input_snps, population, r2_cutoff) {
# Aim: to get haploreg file through log on to HaploReg website using python
# Input: input snps, population, r2_cutoff
# Output: haploreg data frame
        
        # python code to webscrape haploreg
        python_code = paste0("import urllib, urllib2; ", 
                             "data = urllib.urlencode({'query':'", paste(input_snps, collapse = ","),"', ",
                             "'ldThresh':'", r2_cutoff, "', ", "'ldPop':'", population, "', 'output':'text'}); ", 
                             "info = urllib2.urlopen('http://www.broadinstitute.org/mammals/haploreg/haploreg_v3.php', data); ", 
                             "content = info.read(); ",
                             "open('assnp.tmp.txt', 'w').write(content)")
        
        python.exec(python_code)
        
        # read the haploreg file downloaded by python webscraper
        haploreg_df = read.haploreg3("assnp.tmp.txt")
        
        # delete the temporary file
        system("rm assnp.tmp.txt")
        
        return (haploreg_df)
}

# helper function 1-2: to get haploreg file
add_query_snp_info = function(ldsnp_haploreg_df, index_snps_sorted, population) {
# Aim: to add query snp information to the haploreg data frame
# Input: ldsnp-haploreg-df, index_snps_sorted, population
# Output: ldsnp_haploreg_df with query snp added for each ld block
        
        # step 1: add population info of query SNPs
        ldsnp_haploreg_df$population = population
        
        # step 2: add query SNPs
        query_snp_vec = vector(length = nrow(ldsnp_haploreg_df))
        
        j = 1 # j is the index SNP's idx
        for (i in 1:nrow(ldsnp_haploreg_df)) {
                if (i == 1) {
                        query_snp_vec[i] = index_snps_sorted[j]
                } else if (ldsnp_haploreg_df$chr[i] != ldsnp_haploreg_df$chr[i-1]) {
                        j = j + 1
                        query_snp_vec[i] = index_snps_sorted[j]
                } else if (ldsnp_haploreg_df$chr[i] == ldsnp_haploreg_df$chr[i-1] & 
                           ldsnp_haploreg_df$pos[i] < ldsnp_haploreg_df$pos[i-1]) {
                        j = j + 1
                        query_snp_vec[i] = index_snps_sorted[j]
                } else {
                        query_snp_vec[i] = index_snps_sorted[j]
                }
        }
        
        ldsnp_haploreg_df$query_snp = query_snp_vec
        
        return(ldsnp_haploreg_df)
}


# function 2: arrange haploreg list ------------------------------------------------------------------------------

arrange_haploreg_list = function(haploreg_list, for_cnv_call = F) {
# Aim: to process haploreg_list
# Input: haploreg list
# Output: processed haploreg df
# helper function:
#         1. process_haploreg_df
#         2. merge_haploreg_df
#         3. sel_haploreg_df
        
        # merge df from haploreg list to get the raw df
        haploreg_df_raw = merge_haploreg_df(haploreg_list)
        if (for_cnv_call) return (sel_haploreg_df(haploreg_df_raw))
        
        # compress the raw df by aggregating information of the same SNPs
        haploreg_df_proc = process_haploreg_df(haploreg_df_raw)
        
        # select certain columns and eliminate indel or SNPs with multiple alternative alleles
        haploreg_df_sel = sel_haploreg_df(haploreg_df_proc)
        
        return (haploreg_df_sel)
}

# helper function 2-1: merge haploreg df from list
merge_haploreg_df = function(haploreg_list) {
        
        haploreg_df = haploreg_list[[1]]
        
        if (length(haploreg_list) > 1) {
                for (i in 2:length(haploreg_list)) {
                        haploreg_df = rbind(haploreg_df, haploreg_list[[i]])
                }
        }
        
        return (haploreg_df)    
}

# helper function 2-2: process haploreg df
process_haploreg_df = function(haploreg_df) {
# Aim: to reshape the haploreg df by aggregating information of the same SNPs
# Input: haploreg_df
# Output: processed haploreg_df
        
        # prepare data processing
        haploreg_df$querysnp_info = paste0(haploreg_df$query_snp,'_', haploreg_df$population)
        unique_snps = unique(as.character(haploreg_df$rsID))
        haploreg_df_proc = data.frame() # normal file
        
        # process data    
        for (snp in unique_snps){
                snp_pos = which (haploreg_df$rsID %in% snp)
                
                if (length(snp_pos) == 1) {
                        haploreg_df_proc = rbind(haploreg_df_proc, haploreg_df[snp_pos, ])
                } else if (length(snp_pos) > 1) {
                        # concatnate information if a SNP occurs multiple times
                        base = haploreg_df[snp_pos[1], ]
                        r2 = base$r2
                        D. = base$D.
                        is_query_snp = base$is_query_snp
                        query_snp = base$query_snp
                        population = base$population
                        querysnp_info = base$querysnp_info
                        
                        for (i in 2:length(snp_pos)){
                                temp = haploreg_df[snp_pos[i], ]
                                r2 = c(r2, temp$r2)
                                D. = c(D., temp$D.)
                                is_query_snp = c(is_query_snp, temp$is_query_snp)
                                query_snp = c(query_snp, temp$query_snp)
                                population = c(population, temp$population)
                                querysnp_info = c(querysnp_info, temp$querysnp_info)
                        }
                        
                        ord = order(r2, decreasing = T)
                        base$r2 = paste(r2[ord], collapse = ",")
                        base$D. = paste(D.[ord], collapse = ",")
                        base$is_query_snp = paste(is_query_snp[ord], collapse = ",")
                        base$query_snp = paste(query_snp[ord], collapse = ",")
                        base$population = paste(population[ord], collapse = ",")
                        base$querysnp_info = paste(querysnp_info[ord], collapse = ",")
                        
                        haploreg_df_proc = rbind(haploreg_df_proc, base)
                        
                } else stop("There's a bug in SNP names!\n")
        }
        
        return (haploreg_df_proc)
}

# helper function 2-3: select columns and rows 
sel_haploreg_df = function(haploreg_df_proc) {
# Aim: to get a cleaner table, and get rid of the indels and SNPs with multiple alternative alleles
        
        # select columns
        haploreg_df_sel = dplyr::select(haploreg_df_proc, 
                                        rsID, chr, pos, ref, alt, r2, D., AFR, AMR, ASN, EUR, 
                                        query_snp, population, 
                                        RefSeq_id, RefSeq_name, RefSeq_direction, RefSeq_distance,
                                        dbSNP_functional_annotation)
        # select rows
        haploreg_df_sel = filter(haploreg_df_sel, nchar(as.character(ref)) == 1 & nchar(as.character(alt)) == 1)
        
        # sort df
        haploreg_df_sel = arrange(haploreg_df_sel, chr, pos)
        
        return(haploreg_df_sel)
}

# main function ---------------------------------------------------------------------------------------------------

get_ldsnp_info_main = function(index_snp_file = "./data/input_snps/LUC_Index_SNPs_20160607.csv", 
                          population = NA, r2_cutoff = 0.5, for_cnv_call = F, output_dir = "./", output_file = NA) {
        
        # step 0: load packages and input SNPs, population
        load_packages_01()
        cat("get high-LD SNPs information from index SNPs ... \n")
        # step 0: generate output file name
        if (is.na(output_file)) {
                output_file  = gen_output_file_ldsnp(index_snp_file = index_snp_file, 
                                                     population = population, 
                                                     r2 = r2_cutoff, 
                                                     for_cnv_call = for_cnv_call, 
                                                     output_dir = output_dir)
        }
        cat("output file name:", output_file, '\n')
        
        indexSNP_df = read_indexSNP(index_snp_file)
        if (!is.na(population)) {
                indexSNP_df$population = population # if there's assigned population, add assigned population
        }
        populations = unique(as.character(indexSNP_df$population))
        pop_n = length(populations)
        
        # step 1: find LD SNPs by population
        haploreg_list = replicate(pop_n, list())
        for (i in 1:pop_n) {
                indexSNP_df_pop = filter(indexSNP_df, population == populations[i])
                haploreg_df_pop = get_ldsnp_info(indexSNP_df_pop, r2_cutoff = r2_cutoff)
                haploreg_list[[i]] = haploreg_df_pop
        }
        
        # step 2: get clean data
        haploreg_df = arrange_haploreg_list(haploreg_list, for_cnv_call)
        
        # write file
        if (output_file  != F) { # if output_file == F, do not write down file
                write.csv0(haploreg_df, output_file)
        }
        
        cat("high-LD SNPs added ... \n")
        return (list(ldsnp_info_df = haploreg_df, output_file = output_file))
}

# Test #
# get_ldsnp_info_main(for_cnv_call = F)
# get_ldsnp_info_main(population = "EUR", output_file = "./LUC_Index+ldsnps_SNPs_20160607_EUR.csv")
# get_ldsnp_info_main(population = "ASN", output_file = "./LUC_Index+ldsnps_SNPs_20160607_ASN.csv")
