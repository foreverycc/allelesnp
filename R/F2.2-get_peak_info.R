# F2.2-get snps in biofeatures
## ==================================================================================
## Function: get_peak_info_main
###  Aim: To check whether the input SNPs are within certain peaks (biofeatures) or not
###  INPUT: input snp file (.csv)
###  OUTPUT: input snp file with information of peaks added
## ==================================================================================

# outline ---------------------------------------------------------------------------------------------------------

# # 1. read snp file into data frame
# # snp_list should contain both the df and GRange object 
# snp_list = read_snpFile()
# 
# # 2. read peak files into a list of GRange object
# peak_gr_list = read_peakFiles(peak_dir)
# 
# # 3. Find whether the SNPs overlap with each peak
# overlap_mat = get_overlap_mat(snp_gr, peak_gr_list)
# 
# # 4. add overlapping info to original snp table 
# snp_df_proc = add_overlap_info(snp_df, overlap_mat)

# 0. load packages ------------------------------------------------------------------------------------------------

load_packages_2.2 = function(){
        packages = c('dplyr', 'rtracklayer','GenomicRanges')
        load = lapply(packages, require, character.only = T)
}

# function: write.csv0
source("./src/scripts/T1-toolbox.R")

# 0. generate output file name ------------------------------------------------------------------------------------

gen_output_file_peakInfo = function(snp_info_file, output_dir = "./", sample_name = "") {
# Aim: to get output file name        
        snp_info_vec = strsplit(snp_info_file, "/")[[1]]
        snp_batch_id = gsub(".csv", "", snp_info_vec[length(snp_info_vec)])
        output_file = paste0(output_dir, "/", snp_batch_id, "_", sample_name, "_peakAnnotation.csv")
        
        return (output_file)
        
}

# 1. read snp file ------------------------------------------------------------------------------------------------

source ("./src/scripts/F2.1-get_alleleDist_info.R")
# same as read_inputSNP_file

# Test #
# snp_info_list = read_inputSNP_file("./data/haploreg_files/LUC_Index+LD_SNPs_20160607.csv")


# 2. read peak files ----------------------------------------------------------------------------------------------

read_peak_files = function(peak_dir){
# Aim: to read peak files (.bed, .narrowPeak, .broadPeak) into GRange object
# Input: peak dir
# Output: peak GRange object
        
        files = list.files(path=peak_dir, pattern="[.P][be][ea][dk]$")
        n = length(files)
        peak_gr_list = replicate(n,list(vector()))
        
        # read peak files
        for (i in 1:n) {
                # -------- Debug -------- #
                cat ("reading", i, "th/", n, 'peak file\n')
                # ----------------------- #
                if (length(grep ("bed", files[i])) == 1) {
                        file_i = paste0(peak_dir, files[i])
                        peak_gr_list[[i]] = import(file_i)
                } else if (length(grep ("narrowPeak", files[i])) == 1){
                        file_i = paste0(peak_dir, files[i])
                        peak_gr_list[[i]] = import(file_i, format = "bedGraph")
                } else if (length(grep ("broadPeak", files[i])) == 1){
                        file_i = paste0(peak_dir, files[i])
                        peak_gr_list[[i]] = import(file_i, format = "bedGraph")
                }
        }
        
        # assign names
        chip_library = gsub("broadPeak", "", gsub(".narrowPeak", "", gsub(".bed", "", files)))
        names(peak_gr_list) = chip_library

        return (peak_gr_list)
}

# Test #
# peak_gr_list = read_peak_files("./data/samples/DDBJ_A549/peak_files/")


# 3. get overlap table --------------------------------------------------------------------------------------------

get_overlap_mat = function(snp_info_gr, peak_gr_list){
# Aim: get overlap matrix between snps and peaks
# Input: snp GRange object, peak GRange object
# Output: overlap matrix
        
        # construct the matrix
        snp_num = length(snp_info_gr)
        biofeature_num = length(peak_gr_list)
        overlap_mat = as.data.frame(matrix(nrow = snp_num, ncol = biofeature_num + 1))
        overlap_mat[,1] = names(snp_info_gr)
        colnames(overlap_mat) = c("SNP", names(peak_gr_list))
        
        # fill in the matrix
        for (i in 1:biofeature_num){
                # -------- Debug -------- #
                # print(i)
                # ----------------------- #
                overlap_mat[, i+1] = overlapsAny(snp_info_gr, peak_gr_list[[i]])
        }
        
        return(overlap_mat)
}

# Test #
# head(get_overlap_mat(snp_info_list$snp_info_gr, peak_gr_list))


# 4. add overlapping info -----------------------------------------------------------------------------------------

add_overlap_info = function(snp_info_df, overlap_mat) {
# Aim: to add overlapping information to snp table
# Input: snp-info-df, overlap-mat
# Output: snp-info-df with peak overlapping information added
        
        # find overlap
        if (ncol(overlap_mat) == 2) {
                overlaps = overlap_mat[2]
        } else overlaps = overlap_mat[, 2:ncol(overlap_mat)]
        
        biofeature_overlap = apply(overlaps, 1, any)
        biofeature_overlap_num = apply(overlaps, 1, sum)
        
        # list out the biofeatures overlap with that SNP
        biofeature_overlap_names =
                apply(overlaps, 1, function(x){
                        pos = which(x == T)
                        temp = names(x)[pos[1]] 
                        if (length(pos) > 1){
                                for (i in 2:length(pos)){
                                        temp = paste0(temp,',',names(x)[pos[i]])
                                }
                        }
                        return (temp)})
        
        # combine data
        snp_info_addPeak_df = cbind(snp_info_df, 
                                    biofeature_overlap, 
                                    biofeature_overlap_num, 
                                    biofeature_overlap_names)

        return (snp_info_addPeak_df)
        
}


# main function ---------------------------------------------------------------------------------------------------

get_peak_info_main = function(snp_info_file, peak_dir, output_dir = "./", sample_name = "", output_file = NA) {
        # 0. load packages
        load_packages_2.2()
        cat("get peak information for SNPs ... \n")
        
        # 0. generate output file name
        if (is.na(output_file)) {
                output_file = gen_output_file_peakInfo(snp_info_file = snp_info_file, 
                                                       output_dir = output_dir, 
                                                       sample_name = sample_name)
        }
        cat("    output file name:", output_file, '\n')
        # 1. read snp file into data frame
        # snp_list should contain both the df and GRange object
        snp_info_list = read_inputSNP_file(snp_info_file)

        # 2. read peak files into a list of GRange object
        peak_gr_list = read_peak_files(peak_dir)

        # 3. Find whether the SNPs overlap with each peak
        overlap_mat = get_overlap_mat(snp_info_gr = snp_info_list$snp_info_gr, 
                                      peak_gr_list = peak_gr_list)

        # 4. add overlapping info to original snp table
        snp_info_addPeak_df = add_overlap_info(snp_info_df = snp_info_list$snp_info_df, 
                                               overlap_mat = overlap_mat)
        
        if (output_file  != F) { # if output_file == F, do not write down file
                write.csv0 (snp_info_addPeak_df, output_file)
        }
        cat("peak information added ... \n")
        return (list(snp_info_addPeak_df = snp_info_addPeak_df, output_file = output_file))
        
}

# Test #
# snp_file_loc = "LUC_Index_SNPs_20160607_riskPop_0.5.csv"
# get_peak_info_main(snp_info_file = snp_file_loc, peak_dir = "./data/samples/DDBJ_A549/peak_files/")





