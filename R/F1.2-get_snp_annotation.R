# F1.2-get snp annotation and distribution
## ==================================================================================
## Function: get_snp_annotation
###  Aim: To get snp annotation (coding, intron, intergenic, UTR, promoter) distribution; and provide ANNOVAR analysis file
###  INPUT: haploreg snp information files
###  OUTPUT1: distribution table of index SNP, LD SNP, and all SNP;
###  OUTPUT2: information of key SNP (SYM, NSM, NSN, UTR, Promoter)
###  OUTPUT3: files for annovar analysis
## ==================================================================================

gen_output_file_list = function(snp_info_file, output_dir = "./") {
        snp_info_file_vec = strsplit(snp_info_file, "/")[[1]]
        snp_batch_id = gsub(".csv", "", snp_info_file_vec[length(snp_info_file_vec)])
        
        # 1. snp distribution files
        snp_dist_file = paste0(output_dir, "/", snp_batch_id, "_snp_distribution.csv")
        indexsnp_dist_file = paste0(output_dir, "/", snp_batch_id, "_indexsnp_distribution.csv")
        ldsnp_dist_file = paste0(output_dir, "/", snp_batch_id, "_ldsnp_distribution.csv")
        
        # 2. snp information files
        indexsnp_info_file = paste0(output_dir, "/", snp_batch_id, "_indexsnp_info.csv")
        sym_snp_info_file = paste0(output_dir, "/", snp_batch_id, "_sym_snp_info.csv")
        nsm_snp_info_file = paste0(output_dir, "/", snp_batch_id, "_nsm_snp_info.csv")
        nsn_snp_info_file = paste0(output_dir, "/", snp_batch_id, "_nsn_snp_info.csv")
        utr_snp_info_file = paste0(output_dir, "/", snp_batch_id, "_utr_snp_info.csv")
        promoter_snp_info_file = paste0(output_dir, "/", snp_batch_id, "_promoter_snp_info.csv")
        
        # 3. annovar input files
        indexsnp_annovar_file = paste0(output_dir, "/", snp_batch_id, "_indexsnp_annovar.txt")
        ldsnp_annovar_file = paste0(output_dir, "/", snp_batch_id, "_ld_snp_annovar.txt")
        exon_snp_annovar_file = paste0(output_dir, "/", snp_batch_id, "_exon_snp_annovar.txt")
        utr_snp_annovar_file = paste0(output_dir, "/", snp_batch_id, "_utr_snp_annovar.txt")
        
        return(list(snp_dist_file = snp_dist_file, 
                    indexsnp_dist_file = indexsnp_dist_file,
                    ldsnp_dist_file = ldsnp_dist_file,
                    indexsnp_info_file = indexsnp_info_file,
                    sym_snp_info_file = sym_snp_info_file,
                    nsm_snp_info_file = nsm_snp_info_file,
                    nsn_snp_info_file = nsn_snp_info_file,
                    utr_snp_info_file = utr_snp_info_file,
                    promoter_snp_info_file = promoter_snp_info_file,
                    indexsnp_annovar_file = indexsnp_annovar_file,
                    ldsnp_annovar_file = ldsnp_annovar_file,
                    exon_snp_annovar_file = exon_snp_annovar_file,
                    utr_snp_annovar_file = utr_snp_annovar_file))
}

# 1. get annotation -----------------------------------------------------------------------------------------------

add_annotation = function(snp_info_df) {
        snp_info_df$annotation_edit = sapply(snp_info_df$dbSNP_functional_annotation, get_anno)
        sel_promoter = snp_info_df$RefSeq_direction == 5 & snp_info_df$RefSeq_distance < 1000 & 
                snp_info_df$dbSNP_functional_annotation == "."
        snp_info_df$annotation_edit[sel_promoter] = "Promoter"
        
        return (snp_info_df)
}

get_anno = function(anno) {
        if (length(grep("NSN",anno)) == 1) {
                return ("Exon_Nonsense")
        } else if (length(grep ("NSM",anno)) == 1) {
                return("Exon_Nonsynonymous") 
        } else if (length(grep ("SYN", anno)) == 1) {
                return ("Exon_Synonymous")
        } else if (length(grep ("U3", anno)) == 1) {
                return ("3-UTR")
        } else if (length(grep ("U5", anno)) == 1){
                return ("5-UTR")
        } else if (length(grep ("INT", anno)) == 1) {
                return ("Intron")
        } else if (length(grep (".", anno)) == 1) {
                return ("Intergenic")
        }
}


# write distribution info -----------------------------------------------------------------------------------------

write_dist_info = function(snp_info_df, output_file_list) {
        # add a feature of whether the SNP is query snp
        snp_info_df$is_query_snp = apply(snp_info_df, 1, function(x) grepl(x["rsID"], x["query_snp"]))
        
        # write down index snp information table
        write.csv(filter(snp_info_df, is_query_snp), output_file_list$indexsnp_info_file)
        
        # write all distribution
        write.csv(table(snp_info_df$annotation_edit), output_file_list$snp_dist_file)
        # write tag snp distribution
        write.csv(table(snp_info_df %>% filter(is_query_snp) %>% dplyr::select(annotation_edit)),
                  output_file_list$indexsnp_dist_file)
        # write ld snp distribution
        write.csv(table(snp_info_df %>% filter(!is_query_snp) %>% dplyr::select(annotation_edit)),
                  output_file_list$ldsnp_dist_file)
}



# write coding snp and promoter snp table -------------------------------------------------------------------------

write_coding_snp_table = function(snp_info_df, output_list) {
        
        # write down doding snp information
        ## synonymous 
        write.csv(filter(snp_info_df, annotation_edit == "Exon_Synonymous"), output_file_list$sym_snp_info_file)
        ## non-synonymous
        write.csv(filter(snp_info_df, annotation_edit == "Exon_Nonsynonymous"), output_file_list$nsm_snp_info_file)
        ## nonsense
        write.csv(filter(snp_info_df, annotation_edit == "Exon_Nonsense"), output_file_list$nsn_snp_info_file)
        ## UTR
        write.csv(filter(snp_info_df, grepl("UTR", annotation_edit)), output_file_list$utr_snp_info_file)
        ## promoter
        write.csv(filter(snp_info_df, annotation_edit == "Promoter"), output_file_list$promoter_snp_info_file)
        
}



# write file for annovar analysis ---------------------------------------------------------------------------------

write_annovar_files = function(snp_info_df, output_file_list) {
        
        snp_info_df$is_query_snp = apply(snp_info_df, 1, function(x) grepl(x["rsID"], x["query_snp"]))
        
        snp_info_df_indexsnp = filter(snp_info_df, is_query_snp)
        sink (output_file_list$indexsnp_annovar_file)
        convert_annovar(snp_info_df_indexsnp)
        sink()
        
        snp_info_df_ldsnp = filter(snp_info_df, !is_query_snp)
        sink (output_file_list$ldsnp_annovar_file)
        convert_annovar(snp_info_df_ldsnp)
        sink()
        
        snp_info_df_exon = filter(snp_info_df, grepl("Exon", annotation_edit))
        sink (output_file_list$exon_snp_annovar_file)
        convert_annovar(snp_info_df_exon)
        sink()    
        
        snp_info_df_utr = filter(snp_info_df, grepl("UTR", annotation_edit))
        sink (output_file_list$utr_snp_annovar_file)
        convert_annovar(snp_info_df_utr)
        sink()
}


convert_annovar = function(snp_info_df){
        for (i in 1:dim(snp_info_df)[1]) {
                cat (paste(gsub(pattern="chr","", snp_info_df$chr[i]), # chromosome
                           snp_info_df$pos[i], # start site
                           (snp_info_df$pos[i] + nchar(as.character(snp_info_df$ref[i])) -1), # end site
                           snp_info_df$ref[i],
                           snp_info_df$alt[i],
                           paste0(as.character(snp_info_df$rsID[i]), "_",as.character(snp_info_df$query_snp[i]), "-",
                                  as.character(snp_info_df$population[i]), "-",as.character(snp_info_df$r2[i])),
                           sep ="\t"), '\n')
        }
}


# main function ---------------------------------------------------------------------------------------------------

get_snp_annotation_main = function(snp_info_file = "./data/haploreg_files/LUC_Index+LD_SNPs_20160607.csv",
                                   output_dir = "./", output_file = NA) {
        
        require("dplyr")
        
        # 0. generate output file list
        output_file_list = gen_output_file_list(snp_info_file = snp_info_file, output_dir = output_dir)
        
        # 1. get annotation from the data annotation
        snp_info_df = read.csv(snp_info_file)
        snp_info_df = add_annotation(snp_info_df)
        head(snp_info_df)
        
        # 2. write down distribution table
        write_dist_info(snp_info_df, output_file_list)        
        
        # 3. write down coding/promoter snp table
        write_coding_snp_table(snp_info_df, output_file_list)      
        
        # 4. write down snp table for annovar annotation
        write_annovar_files(snp_info_df, output_file_list)
        
}
        
# Test #
# get_snp_annotation_main()

