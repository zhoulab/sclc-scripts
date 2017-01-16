printf <- function(...) cat(sprintf(...),"\n")
library("optparse")
library("data.table")
option_list = list(
    make_option(c("-i", "--gene_pairs_num_file"), type="character",
                help="gene pairs filepath"),
    make_option("--out_all", type="character",
                help="output filepath"),
    make_option(c("-n", "--num_samples"), type="integer",
                help="number of samples")
);
opt_parser = OptionParser(option_list=option_list);
args = parse_args(opt_parser);

printf("Loading dataset from %s...", args$gene_pairs_num_file)
data <- fread(args$gene_pairs_num_file, header=T, sep='\t')
printf("Starting analysis...")
get_fisher <- function(df) {
    n4 <- as.numeric(df["Common"])
    n2 <- as.numeric(df["Gene2Samples"]) - n4
    n3 <- as.numeric(df["Gene1Samples"]) - n4
    n1 <- args$num_samples - (n2 + n3 + n4)
    mat <- matrix(c(n1,n2,n3,n4), ncol=2)
    f <- fisher.test(mat)
    return(f$p.value)
}
data$P_value <- apply(data, 1, get_fisher)
printf("Adjusting p-values...")
data$Adjusted_p <- p.adjust(data$P_value, method="BH")
data <- data[order(data$P_value), ]
printf("Writing output file: %s", args$out_all)
write.table(format(data, scientific=T, digits=3),
            file=args$out_all, sep="\t", quote=F, row.names=F)
