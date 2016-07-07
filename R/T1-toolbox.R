# T1 toolbox

# toolbox: write.csv with no row names
write.csv0 = function(...) {
        return(write.csv(..., row.names = F))
}

# toolbox: get # of lines of code
get_code_lines_counts = function() {
        line_count = 0
        for (file in list.files("./src/scripts/", ".R", full.names = T)) {
                print (file)
                bash_output = system(paste('wc -l', file), intern = T)
                bash_output_vec = strsplit(bash_output, split = " ")[[1]]
                line_count = line_count + as.numeric(bash_output_vec[length(bash_output_vec) - 1])
        }
        cat ("Altogether", line_count, "lines of code.")
}

# get_code_lines_counts()
