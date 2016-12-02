fname <- 'SCLC.muTect.vcf'
con <- file(fname, open='r')

# skip initial comment lines
skip_ct = 0
l <- readLines(con, n=1)
while (substring(l, 2, 2) == '#') {
    skip_ct = skip_ct + 1
    l <- readLines(con, n=1)
}

# get header names
if (substring(l, 1, 1) == '#') {
    header <- unlist(strsplit(substring(l, 2), split='\t'))
} else {
    header <- unlist(strsplit(l, split='\t'))
}
skip_ct = skip_ct + 1
print(skip_ct)

vcf_df <- read.table(fname, header=FALSE, sep='\t', skip=skip_ct)
colnames(vcf_df) <- header

# TODO: tranformation

close(con)
