# F2.1-get allele distribution info
## ==================================================================================
## Function: get_alleleDist_info_main
###  Aim: to get allelic-distribution of input SNPs in functional sequencing files (.bam) ')
###  INPUT: input snp file (.csv), vcf/encode cnv data
###  OUTPUT: input snp file with information of allele distribution in various bam data
## ==================================================================================

# toolbox ---------------------------------------------------------------------------------------------------------

# function: write.csv0
source("./src/scripts/T1-toolbox.R")

# toolbox: convert ascii to quality score
ascii_to_quality_score = function(x) { 
        strtoi(charToRaw(x),16L) 
}


# functions -------------------------------------------------------------------------------------------------------

# step 0-1: load packages -----------------------------------------------------------------------------------------

load_packages_2.1 = function(){
        packages = c('Rsamtools', "BSgenome.Hsapiens.UCSC.hg19", "dplyr")
        load = lapply(packages, require, character.only = T)
}


# step 0-2: function: generate parameter list ---------------------------------------------------------------------

gen_param_list = function(snp_info_file,
                          cell_sel, 
                          output_dir,
                          sample_name,
                          output_file,
                          base_qual_threshold,
                          mapq_threshold, 
                          server,
                          bam_dir,
                          rmdup_file) {
        # Test #
        # snp_info_file = "./data/input_snps/Try_query_SNPs.csv"
        # cell_sel = ""
        # output_dir = "./"
        # sample_name = ""
        # base_qual_threshold = 20
        # mapq_threshold = 20
        # server = "http"
        # bam_dir = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDgf/"
        # rmdup_file = F
        
        # get the analysis id
        snp_info_vec = strsplit(snp_info_file, "/")[[1]]
        snp_batch_id = gsub(".csv", "", snp_info_vec[length(snp_info_vec)])
        
        if (is.na(output_file)) {
                output_file = paste0(output_dir, "/", snp_batch_id, "_", sample_name, "_alleleDist.csv")
        }
                
        param_list = 
                list(snp_info_file = snp_info_file,
                     snp_batch_id = snp_batch_id,
                     cell_sel = cell_sel, 
                     output_file = output_file,
                     base_qual_threshold = base_qual_threshold,
                     mapq_threshold = mapq_threshold, 
                     server = server,
                     bam_dir = bam_dir,
                     rmdup_file = rmdup_file)
        
        return (param_list)
        
}

# step 1: process snp information ---------------------------------------------------------------------------------

read_inputSNP_file = function(snp_info_file) {
        # Aim: to parse the snp information for further use
        # Input: snp_info file, should contain at least the following fields:
        # 1. rsID
        # 2. chr
        # 3. pos
        # 4. ref
        # 5. alt
        # Output: a list that contain basic information of those SNPs
        
        # read and parse snp information file
        snp_info_df = read.csv(snp_info_file, header= T)
        head(snp_info_df)
        snp_info_df = arrange(snp_info_df, as.numeric(gsub("chr", "", chr)), pos)
        snp_id = as.character(snp_info_df$rsID)
        snp_ref = as.character(snp_info_df$ref)
        snp_alt = as.character(snp_info_df$alt)
        chr = as.character(snp_info_df$chr)
        snp_chr = ifelse (grepl("chr", chr), chr, paste0("chr", chr))
        snp_pos = as.numeric(as.character(snp_info_df$pos))
        which = GRanges(seqnames = snp_chr, ranges = IRanges(start = snp_pos,end = snp_pos))
        names(which) = snp_id
        snp_info_gr = which
        
        # return a list of snp information
        return(list(snp_info_gr = snp_info_gr, snp_info_df = snp_info_df,
                    snp_id = snp_id, snp_ref = snp_ref, snp_alt = snp_alt, 
                    snp_chr = snp_chr, snp_pos = snp_pos))
}


# step 2: read in bam files ---------------------------------------------------------------------------------------

read_bam_files = function(snp_info_gr, param_list) {
        # Aim: to read bam reads information on selected snps from:
        # 1. http website
        # 2. ftp website
        # 3. local computer
        # Input: snp_info_gr (grange object of SNP)
        # Output: a list contain bam information
        
        # assign variable values
        server = param_list$server
        bam_dir = param_list$bam_dir
        cell_sel = param_list$cell_sel
        
        # read bam files on http, ftp, and local computer
        if (server == "http") {
                dir_bam = bam_dir
                file_table = read.delim(paste0(dir_bam, "files.txt"), header = F)
                if (all(!grepl(paste0(cell_sel, ".*.bam$"), file_table[, 1]))) {
                        stop("No bam files with the selected cell type is found. 
                             Please modify the remote bam directory or cell types.\n")
                } else {
                        bam_files = as.character(file_table[grep (paste0(cell_sel, ".*.bam$"), file_table[,1]), 1])
                }
        } else if (server == "ftp") {
                dir_bam = gsub("http", "ftp", bam_dir)
                file_table = read.delim(paste0(dir_bam, "files.txt"), header = F)
                if (all(!grepl(paste0(cell_sel, ".*.bam$"), file_table[, 1]))) {
                        stop("No bam files with the selected cell type is found. 
                             Please modify the remote bam directory or cell types.\n")
                } else {
                        bam_files = as.character(file_table[grep (paste0(cell_sel, ".*.bam$"), file_table[,1]), 1])
                }
        } else if (server == "local") {
                dir_bam = bam_dir
                bam_files = as.character(list.files(dir_bam, paste0(cell_sel, ".*.bam$")))
                if (length(bam_files) == 0) {
                        stop("No bam files with the selected cell type is found. 
                             Please modify the local bam directory or cell types.\n")
                }
                } else stop("Please choose the right server: 'http' or 'ftp' or 'local.")
        
        ## get the param
        what = c("qname", "flag","rname","pos", "mapq", "cigar","seq", "qual","strand", "qwidth")
        tag = c("NM", "RG","MD")
        param = ScanBamParam(which = snp_info_gr, what = what, tag = tag) # this function will reorder the snps comparing to haploreg files
        
        ## scan bam files
        bam_super_list = replicate(length(bam_files), list())
        names(bam_super_list) = bam_files
        for (i in 1:length(bam_files)) {
                file = bam_files[i]
                cat("reading", file , '\n')
                bam_file_link = paste0(dir_bam, "/", file)
                bam_super_list[[i]] = scanBam(bam_file_link, param = param)
                gc()
        }
        
        return (list(bam_super_list = bam_super_list, bam_files = bam_files))
}



# step 3: process bam data ----------------------------------------------------------------------------------------

# step 3-1 function: add bam information --------------------------------------------------------------------------

add_bam_sample_info = function(bam_list, bam_lib){
        # Aim: to add some information to bam_list
        # Input: bam_list, bam_lib (library of the bam, normally is just the bam file without .bam)
        # Output: bam_list with some added information
        
        # --------------- test ----------------- #
        #     bam_list = bam_list
        #     bam_lib = bam_lib
        # -------------------------------------- #
        
        # add base quality info
        for (i in 1: length(bam_list)){
                # add mapq info if there's no
                if (!any(grepl("mapq", names(bam_list[[i]])))) {
                        bam_list[[i]]$mapq = NA
                }
                # add tag sample and library information
                bam_list[[i]]$tag$sm = bam_lib
                bam_list[[i]]$tag$lb = rep(bam_lib, length(bam_list[[i]]$seq))
        }
        
        return(bam_list)
}


# step 3-2 function: add bam base quality information -------------------------------------------------------------

add_bam_base_info = function(bam_list, snp_info_list) {
        # Aim: to add base-level information (ref/alt allele, base quality) to bam_list
        # Input: bam_list (the reads information of all SNPs in a bam file), snp_info_list
        # Output: bam_list with base-level information added
        # helper function: parse_base_info_snp
        
        # --------------- test ----------------- #
        #           bam_list = bam_list
        #           snp_chr = snp_info_list$snp_chr
        #           snp_pos = snp_info_list$snp_pos
        #           snp_id = snp_info_list$snp_id
        #           snp_ref = snp_info_list$snp_ref
        #           snp_alt = snp_info_list$snp_alt
        # -------------------------------------- #
        
        # add some snp information
        for (i in 1: length(bam_list)){
                bam_list[[i]]$snp$id = snp_info_list$snp_id[i]
                bam_list[[i]]$snp$chr = snp_info_list$snp_chr[i]
                bam_list[[i]]$snp$pos = snp_info_list$snp_pos[i]
                bam_list[[i]]$snp$ref = snp_info_list$snp_ref[i]
                bam_list[[i]]$snp$alt = snp_info_list$snp_alt[i]
        }
        
        # add some base and base-quality information
        for (i in 1: length(bam_list)){
                # --------------- test ----------------- #
                # print (paste("ith snp:", i))
                # print (bam_list[[i]]$snp$id)
                # -------------------------------------- #
                
                bam_list_i= bam_list[[i]]
                seq_list = bam_list_i$seq
                
                # if there's no reads on that snp in chip-seq, just assign NA
                if (length(seq_list) == 0) {
                        bam_list[[i]]$snp$base = NA
                        bam_list[[i]]$snp$base_qual = NA
                        next
                } else {
                        bam_list[[i]] = parse_base_info_snp(bam_list_i)
                }
        }
        
        return(bam_list)
}


# helper function 3-2-1: parse base information

parse_base_info_snp = function(bam_list_i){
        # Aim: parse base information (base, base_qual) for reads that contain a particular SNP in a bam file
        # Input: bam_list_i (reads info of a SNP in a bam file)
        # Output: bam_list_i with added base-level information
        # helper function: get_relative_loc
        
        seq_list = bam_list_i$seq
        seq_qual_list = bam_list_i$qual
        snp_rel_loc_vec = get_relative_loc(bam_list_i)
        
        base_vec = vector()
        base_qual_vec = vector()
        
        for (j in 1:length(seq_list)){
                # ----- test ----- #
                # print (paste("j:", j))
                # ---------------- #
                seq_j = seq_list[[j]]
                # ----- test ----- #
                # print (seq_j)
                # ---------------- #
                seq_qual_j = seq_qual_list[[j]]
                # ----- test ----- #
                # print (seq_qual_j)
                # ---------------- #
                snp_rel_loc_j = snp_rel_loc_vec[j]
                # ----- test ----- #
                # print (snp_rel_loc_j)
                # ---------------- #
                
                if (snp_rel_loc_j > length(seq_j) | snp_rel_loc_j < 1 | is.na(snp_rel_loc_j)) { 
                        # for some index out-of-boundary problems
                        base_vec = c(base_vec, NA)
                        base_qual_vec = c(base_qual_vec, NA)
                } else {
                        base_vec = c(base_vec, as.character(seq_j[snp_rel_loc_j]))
                        base_qual_vec = c(base_qual_vec, ascii_to_quality_score(as.character(seq_qual_j[snp_rel_loc_j])))
                }
        }
        
        # base vector
        bam_list_i$snp$base = base_vec
        
        # modify quality score
        if ( any(base_qual_vec > (41 + 33), na.rm= T)) {
                bam_list_i$snp$base_qual = base_qual_vec - 64 # for phred score of illumina 1.8-
        } else {
                bam_list_i$snp$base_qual = base_qual_vec - 33 # for phred score of illumina 1.8+
        }
        
        return(bam_list_i)
}

# helper function 3-2-2: to calculate relative location based on cigar information

get_relative_loc = function(bam_list_i) {
        # Aim: to calculate the relative location of a snp in the bam list
        # Input: bam_list_i
        # Output: the relative location of the particular SNP
        
        if (length(bam_list_i$seq) == 0) {
                return (NA)
        } 
        
        snp_rel_loc_vec = vector(length = length(bam_list_i$seq))
        
        # calculate relative location based on cigar
        for (j in 1:length(bam_list_i$seq)) {
                cigar_j = bam_list_i$cigar[j]
                if (cigar_j == paste0(bam_list_i$qwidth[j], "M")) {
                        snp_rel_loc_vec[j] = bam_list_i$snp$pos - bam_list_i$pos[j] + 1 # perfect match / mismatch
                } else if (grepl("^[[:digit:]]+M[[:digit:]]+S$", cigar_j)) { # 3' soft clipping
                        snp_rel_loc_vec[j] = bam_list_i$snp$pos - bam_list_i$pos[j] + 1
                        seq_len = as.numeric(gsub("M[[:digit:]]+S", "", cigar_j))
                        if (snp_rel_loc_vec[j] > seq_len) snp_rel_loc_vec[j] = NA
                } else if (grepl("^[[:digit:]]+S[[:digit:]]+M$", cigar_j)) { # 5' soft clipping
                        soft_clipping_len = as.numeric(gsub("S[[:digit:]]+M", "", cigar_j))
                        #       seq_len = bam_list_i$qwidth[j] - soft_clipping_len
                        snp_rel_loc_vec[j] = bam_list_i$snp$pos - bam_list_i$pos[j] + 1 + soft_clipping_len
                        if (snp_rel_loc_vec[j] > bam_list_i$qwidth[j]) snp_rel_loc_vec[j] = NA
                } else {
                        snp_rel_loc_vec[j] = NA # if there's indel, ignore the read
                }
        }
        
        return (snp_rel_loc_vec)
}


# step 3-3 function: filter bam list  -----------------------------------------------------------------------------

filter_bam_list = function(bam_list, param_list) {
        # Aim: filter bam_list by base_qual, mapq and reads of duplicates
        # Input: bam_list
        # Output: filtered bam_list
        # helper function: filter_bam_reads
        
        # assign variables
        base_qual_threshold = param_list$base_qual_threshold
        mapq_threshold = param_list$mapq_threshold
        rmdup_file = param_list$rmdup_file
        
        for (i in 1: length(bam_list)){
                # ------ debug ------- #
                # print (i) 
                # -------------------- #
                bam_list_i = bam_list[[i]]
                snp_alleles = c(bam_list_i$snp$ref, bam_list_i$snp$alt)
                if (is.null (bam_list_i$index_del_vec)) bam_list_i$index_del_vec = vector()
                # if there's no reads on that SNP in chip-seq, just go next round
                if (length(bam_list_i$seq) == 0){
                        next
                } else {
                        # filter reads of duplicates
                        if (rmdup_file) {
                                # filter reads of duplicates
                                if (any(bam_list_i$flag >= 1024)) {
                                        bam_list_i$index_del_vec = c(bam_list_i$index_del_vec, which(bam_list_i$flag >= 1024))
                                }
                        }
                        
                        # filter reads that has more than 1 mismatch
                        if (any(bam_list_i$tag$NM > 2)) {
                                bam_list_i$index_del_vec = c(bam_list_i$index_del_vec, which(bam_list_i$tag$NM > 2))
                        }
                        # filter reads that does not contain in known snp alleles
                        if (!any(bam_list_i$snp$base %in% snp_alleles)){
                                bam_list_i$index_del_vec = c(bam_list_i$index_del_vec, which(!(bam_list_i$snp$base %in% snp_alleles)))
                        }
                        # filter reads that has reference allele, but has >1 mismatch
                        any_bad_base = any((bam_list_i$snp$base == bam_list_i$snp$ref) & (bam_list_i$tag$NM > 1))
                        if (is.na(any_bad_base) | any_bad_base) {
                                bam_list_i$index_del_vec = c(bam_list_i$index_del_vec,
                                                             which((bam_list_i$snp$base == bam_list_i$snp$ref) & (bam_list_i$tag$NM > 1)))
                                bam_list_i$index_del_vec = c(bam_list_i$index_del_vec, which(is.na(bam_list_i$snp$base)))
                        }
                        # filter reads that has low base quality
                        any_bad_base_qual = any((bam_list_i$snp$base_qual < base_qual_threshold) | (is.na(bam_list_i$snp$base_qual)))
                        if (is.na(any_bad_base_qual) | any_bad_base_qual) {
                                bam_list_i$index_del_vec = c(bam_list_i$index_del_vec, which(bam_list_i$snp$base_qual < base_qual_threshold))
                                bam_list_i$index_del_vec = c(bam_list_i$index_del_vec, which(is.na(bam_list_i$snp$base_qual)))
                        }
                        # filter reads that has low mapping quality
                        any_bad_mapq = any(is.na(bam_list_i$mapq) | bam_list_i$mapq < mapq_threshold)
                        if (is.na(any_bad_mapq) | any_bad_mapq) {
                                bam_list_i$index_del_vec = c(bam_list_i$index_del_vec, which(bam_list_i$mapq < mapq_threshold))
                                bam_list_i$index_del_vec = c(bam_list_i$index_del_vec, which(is.na(bam_list_i$mapq)))
                        }
                        ## summarize filtered reads
                        if (length(bam_list_i$index_del_vec) == 0) {
                                next
                        } else {
                                bam_list_i$index_del_vec = unique(bam_list_i$index_del_vec)
                                bam_list[[i]] = filter_bam_reads(bam_list_i, bam_list_i$index_del_vec)
                        }
                }
        }
        
        return (bam_list)
}

# helper function 3-3-1: filter bam reads 
filter_bam_reads = function(bam_list_i, index_del) {
        # Aim: filter bam reads by index to delete
        # Input: bam_list_i, index_del (index to delete)
        # Output: bam_list_i with deleted indexes
        
        index_remain = setdiff(1:length(bam_list_i$rname), index_del)
        bam_list_i$qname = bam_list_i$qname[index_remain]
        bam_list_i$flag = bam_list_i$flag[index_remain]
        bam_list_i$rname = bam_list_i$rname[index_remain]
        bam_list_i$strand = bam_list_i$strand[index_remain]
        bam_list_i$pos = bam_list_i$pos[index_remain]
        bam_list_i$qwidth = bam_list_i$qwidth[index_remain]
        bam_list_i$mapq = bam_list_i$mapq[index_remain]
        bam_list_i$cigar = bam_list_i$cigar[index_remain]
        bam_list_i$seq = bam_list_i$seq[index_remain]
        bam_list_i$qual = bam_list_i$qual[index_remain]
        bam_list_i$tag$NM = bam_list_i$tag$NM[index_remain]
        bam_list_i$tag$RG = bam_list_i$tag$RG[index_remain]
        bam_list_i$tag$MD = bam_list_i$tag$MD[index_remain]
        bam_list_i$tag$LB = bam_list_i$tag$LB[index_remain]
        bam_list_i$snp$base = bam_list_i$snp$base[index_remain]
        bam_list_i$snp$base_qual = bam_list_i$snp$base_qual[index_remain]
        
        return(bam_list_i)
}


# step 3-4 function: calculate reads stats ------------------------------------------------------------------------

add_bam_reads_stats = function(bam_list){
        # Aim: to add ref-v.s.-alt statistics to bam files
        # Input: bam_list
        # Output: bam_list with reads statistics added
        
        for (i in 1:length(bam_list)){
                # --------------- debug ---------------- #
                #       print (i)
                # -------------------------------------- #
                bam_list_i = bam_list[[i]]
                snp_ref_i = bam_list_i$snp$ref
                snp_alt_i = bam_list_i$snp$alt
                
                # add snp ref vs alt info
                bam_list_i$snp$ref_v_alt = c(sum(bam_list_i$snp$base == snp_ref_i, na.rm = T),
                                             sum(bam_list_i$snp$base == snp_alt_i, na.rm = T))
                # add snp ref vs alt info with duplicates removed
                uniq_pos_ref = unique(bam_list_i$pos[bam_list_i$snp$base == snp_ref_i])
                ref_rmdup = ifelse (is.na(uniq_pos_ref), 0, length(uniq_pos_ref))
                ref_rmdup = ifelse (!length(ref_rmdup), 0, ref_rmdup)
                uniq_pos_alt = unique(bam_list_i$pos[bam_list_i$snp$base == snp_alt_i])
                alt_rmdup = ifelse (is.na(uniq_pos_alt), 0, length(uniq_pos_alt))
                alt_rmdup = ifelse (!length(alt_rmdup), 0, alt_rmdup)
                bam_list_i$snp$ref_v_alt_rmdup = c(ref_rmdup, alt_rmdup)
                
                names(bam_list_i$snp$ref_v_alt) = c("ref", "alt")
                
                # add ref vs alt info for each library
                bam_sm_i = as.character(bam_list_i$tag$sm)
                bam_lib_rva_df = data.frame()
                
                # if there's no reads on that snp in chip-seq, just assign a 1-d data frame
                bam_lib_rva_df = data.frame(t(c(bam_list_i$snp$ref_v_alt, bam_list_i$snp$ref_v_alt_rmdup)))
                colnames(bam_lib_rva_df) = c('ref', 'alt', 'ref_rmdup', 'alt_rmdup')
                bam_lib_rva_df$rsID = bam_list_i$snp$id
                bam_lib_rva_df$biofeature = bam_sm_i
                # reorder
                bam_lib_rva_df = bam_lib_rva_df[, c('rsID','biofeature', 'ref','alt', 'ref_rmdup', 'alt_rmdup')]
                bam_list_i$stats = bam_lib_rva_df
                bam_list[[i]] = bam_list_i
        }
        
        return(bam_list = bam_list)
}


# 3-5 function: generate allele distribution table ----------------------------------------------------------------

gen_allele_distribution_table = function(bam_list){
        # Aim: to generate allele distribution table by gathering together the SNPs' distributions
        # Input: bam_list
        # Output: allele distribution table of the bam file
        
        bam_stats_list = lapply(bam_list, function(x) x$stats)
        stats_df = data.frame()
        # merge stats tables of each snp
        for (i in 1: length(bam_stats_list)) {
                stats_df = rbind(stats_df, bam_stats_list[[i]])
        }
        
        return(stats_df)
}


# 3-6 function: merge replicates ----------------------------------------------------------------------------------

merge_replicates_table = function(ad_table) {
# Aim: to merge replicates
        ad_table$biofeature = sapply(ad_table$biofeature, function(x) {
                gsub("[Rr][Ee][Pp][1-9]", replacement = "", x)
        })
        ad_table_summ = ad_table %>% group_by(rsID, biofeature) %>% 
                summarise(ref = sum(ref), alt = sum(alt), 
                          ref_rmdup = sum(ref_rmdup), alt_rmdup = sum(alt_rmdup))
        
        return(ad_table_summ)
}

### --------------------------------------------------------------------------------- ###
### ------- main function: to check allele distribution of snps in batch mode ------- ###
### --------------------------------------------------------------------------------- ###

get_alleleDist_info_main = function(
        snp_info_file,
        cell_sel = "", 
        output_dir = "./",
        sample_name = "",
        output_file = NA, 
        base_qual_threshold = 20,
        mapq_threshold = 20, 
        server = "local",
        bam_dir = "./",
        rmdup_file = F,
        merge_replicates = F
) {
        # ---------------------- test --------------------- #
        # snp_info_file = "./data/input_snps/LUC_Index_SNPs_20160607.csv"
        # cell_sel = ""
        # output_dir = "./"
        # sample_name = "DDBJ_A549"
        # base_qual_threshold = 20
        # mapq_threshold = 20
        # server = "local"
        # bam_dir = "./data/samples/DDBJ_A549/bam_files/"
        # rmdup_file = F
        # ------------------------------------------------- #
        
        ### -------------------- step 0-1: load the packages -------------------- ###
        load_packages_2.1()
        cat("get bam information for SNPs ... \n")
        
        ### ------------------ step 0-2: generate parmater list ------------------- ###
        param_list = gen_param_list(snp_info_file = snp_info_file,
                                    cell_sel = cell_sel, 
                                    output_dir = output_dir,
                                    sample_name = sample_name,
                                    output_file = output_file,
                                    base_qual_threshold = base_qual_threshold,
                                    mapq_threshold = mapq_threshold, 
                                    server = server,
                                    bam_dir = bam_dir,
                                    rmdup_file = rmdup_file)
        
        cat("output file name:", param_list$output_file, '\n')
        
        ### ----------------------- step 1: get the snp ----------------------- ###
        snp_info_list = read_inputSNP_file(snp_info_file= snp_info_file)
        
        ### ------- step 2: read bam files that contain the snp region -------- ###
        if (server != "local") {
                dir.create('./bai_files') # create a directory and store the bai files
                setwd('./bai_files')
                cat(".bai files will be downloaded in ", getwd(), '\n')
                bam_data_list = read_bam_files(snp_info_gr = snp_info_list$snp_info_gr, param_list)
                setwd("../")
        } else {
                bam_data_list = read_bam_files(snp_info_gr = snp_info_list$snp_info_gr, param_list)
        }
        
        ### --------------------- step 3: process bam data --------------------- ###
        # extract bam data
        bam_files = bam_data_list$bam_files
        bam_super_list = bam_data_list$bam_super_list
        
        # initialize a data frame for allele-distribution
        ad_table = data.frame()
        
        # for each bam file, process data...
        for (i in 1: length(bam_files)) {
                # ----- test ----- #
                # print(i)
                # ---------------- #
                
                cat("processing reads ...", bam_files[[i]], "\n")
                
                # step 3-1: add bam information
                bam_lib = gsub(".bam","", bam_files[[i]])
                bam_list = bam_super_list[[i]]
                bam_list_add_bam_sample_info = add_bam_sample_info(bam_list= bam_list, 
                                                                   bam_lib= bam_lib)
                
                # step 3-2: add base quality information
                bam_list_add_bam_base_q = add_bam_base_info(bam_list= bam_list_add_bam_sample_info, 
                                                            snp_info_list = snp_info_list)
                
                # step 3-3: filter low base/mapping quality reads --- after this step, we can get relatively clean data
                bam_list_filter = filter_bam_list(bam_list= bam_list_add_bam_base_q, param_list)
                
                # step 3-4: calculate allele-distribution statistics
                bam_list_add_stats = add_bam_reads_stats(bam_list= bam_list_filter)
                
                # step 3-5: generate allele-distribution table
                ad_table_i = gen_allele_distribution_table(bam_list= bam_list_add_stats) 
                
                # add allele-distribution table together
                ad_table = rbind(ad_table, ad_table_i)
        }
        
        if (merge_replicates & nrow(ad_table) > 0) {
                ad_table = merge_replicates_table(ad_table)
        }
        
        ad_table = filter(ad_table, !(ref == 0 & alt == 0))
        # write down allele-distribution table
        if (param_list$output_file != F) { # if output_file == F, do not write down file
                write.csv0(ad_table, param_list$output_file)
        }
        
        cat("allele-distribution for SNPs written ... \n")
        
        return(list(snp_info_alleleDist_df = ad_table, output_file = param_list$output_file))
}


### ------------ main function for checking uw dgf data ------------ ###
get_alleleDist_info_encodeDGF = 
        function(snp_info_file,
                 cell_sel = "", 
                 output_dir = "./", 
                 sample_name = "ENCODE_DGF",
                 server = "http",
                 bam_dir = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDgf/", 
                 ...) {
                get_alleleDist_info_main(
                        snp_info_file = snp_info_file,
                        server = server, 
                        bam_dir = bam_dir,
                        cell_sel = cell_sel, 
                        output_dir = output_dir, 
                        sample_name = sample_name, ...)
        } 

### ----------- main function for checking uw dnase data ---------- ###
get_alleleDist_info_encodeDNase = 
        function(snp_info_file,
                 cell_sel = "", 
                 output_dir = "./", 
                 sample_name = "ENCODE_DNase",
                 server = "http",
                 bam_dir = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/", 
                 ...) {
                get_alleleDist_info_main(
                        snp_info_file = snp_info_file,
                        server = server, 
                        bam_dir = bam_dir,
                        cell_sel = cell_sel, 
                        output_dir = output_dir, 
                        sample_name = sample_name, ...)
        } 

### ---------- main function for checking haib tfbs data ---------- ###
get_alleleDist_info_encodeHaibTFBS = 
        function(snp_info_file,
                 cell_sel = "", 
                 output_dir = "./", 
                 sample_name = "ENCODE_HaibTFBS",
                 server = "http",
                 bam_dir = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/", 
                 ...) {
                get_alleleDist_info_main(
                        snp_info_file = snp_info_file,
                        server = server, 
                        bam_dir = bam_dir,
                        cell_sel = cell_sel, 
                        output_dir = output_dir, 
                        sample_name = sample_name, ...)
        }

### ---------- main function for checking sydh tfbs data ---------- ###
get_alleleDist_info_encodeSydhTFBS = 
        function(snp_info_file,
                 cell_sel = "", 
                 output_dir = "./", 
                 sample_name = "ENCODE_SydhTFBS",
                 server = "http",
                 bam_dir = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/", 
                 ...) {
                get_alleleDist_info_main(
                        snp_info_file = snp_info_file,
                        server = server, 
                        bam_dir = bam_dir,
                        cell_sel = cell_sel, 
                        output_dir = output_dir, 
                        sample_name = sample_name, ...)
        }

### -------- main function for checking uchicago tfbs data -------- ###
get_alleleDist_info_encodeUchicagoTFBS = 
        function(snp_info_file,
                 cell_sel = "", 
                 output_dir = "./", 
                 sample_name = "ENCODE_UChicagoTFBS",
                 server = "http",
                 bam_dir = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUchicagoTfbs/", 
                 ...) {
                get_alleleDist_info_main(
                        snp_info_file = snp_info_file,
                        server = server, 
                        bam_dir = bam_dir,
                        cell_sel = cell_sel, 
                        output_dir = output_dir, 
                        sample_name = sample_name, ...)
        }

### ------- main function for checking broad histone data -------- ###
get_alleleDist_info_encodeBroadHistone = 
        function(snp_info_file,
                 cell_sel = "", 
                 output_dir = "./", 
                 sample_name = "ENCODE_BroadHistone",
                 server = "http",
                 bam_dir = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/", 
                 ...) {
                get_alleleDist_info_main(
                        snp_info_file = snp_info_file,
                        server = server, 
                        bam_dir = bam_dir,
                        cell_sel = cell_sel, 
                        output_dir = output_dir, 
                        sample_name = sample_name, ...)
        }

### --------- main function for checking uw histone data --------- ###
get_alleleDist_info_encodeUwHistone = 
        function(snp_info_file,
                 cell_sel = "", 
                 output_dir = "./", 
                 sample_name = "ENCODE_UWHistone",
                 server = "http",
                 bam_dir = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwHistone/", 
                 ...) {
                get_alleleDist_info_main(
                        snp_info_file = snp_info_file,
                        server = server, 
                        bam_dir = bam_dir,
                        cell_sel = cell_sel, 
                        output_dir = output_dir, 
                        sample_name = sample_name, ...)
        }

### --------------------------------------------------------------------------------- ###
### --------- identify allele-specific effects in luc risk-associated snps ---------- ###
### --------------------------------------------------------------------------------- ###

### ------------- test ------------- ###
# snp_file_loc = "./data/haploreg_files/LUC_Index+LD_SNPs_20160607.csv"
# snp_file_loc = "./LUC_Index_SNPs_20160607_short_ENCODE_DNase_assnp/LUC_Index_SNPs_20160607_short_riskPop_0.5.csv"
# get_alleleDist_info_local(snp_info_file = snp_file_loc, bam_dir = "./data/samples/DDBJ_A549/bam_files/")
# get_alleleDist_info_encodeDNase(snp_info_file= "./data/input_snps/Try_query_SNPs.csv", cell_sel = "A549", merge_replicates= T)
# get_alleleDist_info_encodeDNase(snp_info_file= snp_file_loc, cell_sel = "A549", merge_replicates= T)
# get_alleleDist_info_encodeHaibTFBS(snp_info_file = snp_file_loc, cell_sel="A549")
# get_alleleDist_info_encodeDGF(snp_info_file = "./data/input_snps/Try_query_SNPs.csv")
### ------------- test ------------- ###

### run analysis 
# setwd("/volumes/macintosh_hd_2/research/00-projects/00-assnp_encode")
# snp_file_loc = "./data/luc_ldsnp_haploreg_arrange_summ_table.csv"
# get_alleleDist_info_encodeDGF(snp_info_file = snp_file_loc)
# get_alleleDist_info_encodeDNase(snp_info_file = snp_file_loc)
# get_alleleDist_info_encodeHaibTFBS(snp_info_file = snp_file_loc)
# get_alleleDist_info_encodeSydhTFBS(snp_info_file = snp_file_loc)
# get_alleleDist_info_encodeUchicagoTFBS(snp_info_file = snp_file_loc)
# get_alleleDist_info_encodeBroadHistone (snp_info_file = snp_file_loc)
# get_alleleDist_info_encodeUwHistone(snp_info_file = snp_file_loc)
