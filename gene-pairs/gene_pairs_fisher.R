printf <- function(...) cat(sprintf(...))
library("optparse")
library("data.table")
option_list = list(
    make_option(c("-f", "--file"), type="character", default=NULL,
                help="dataset file name", metavar="character"),
    make_option("--out_all", type="character", default="out_all.txt",
                help="output file name [default= %default]", metavar="character"),
    make_option("--out_filtered", type="character", default="out_filtered.txt",
                help="output file name [default= %default]", metavar="character"),
    make_option(c("-a", "--alpha"), type="double", default=NULL,
                help="alpha level", metavar="character"),
    make_option(c("-n", "--num_samples"), type="integer", default=NULL,
                help="number of samples", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
args = parse_args(opt_parser);

printf("Reading lines from %s...\n", args$file)
printf("Output file destination: %s...\n", args$out_all)
printf("Filtered file destination: %s...\n", args$out_filtered)

ptm = proc.time()
print("Loading dataset...")
data <- fread(args$file, header=T, sep='\t')
print("Starting analysis...")
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
data$Adjusted_p <- p.adjust(data$P_value, method="BH")
data <- data[order(data$P_value), ]
write.table(format(data, scientific=T, digits=3),
            file=args$out_all, sep="\t", quote=F, row.names=F)
write.table(format(data[data$P_value < args$alpha, ], scientific=T, digits=3),
            file=args$out_filtered, sep="\t", quote=F, row.names=F)

printf("\n")
print (proc.time() - ptm)
