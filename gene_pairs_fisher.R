printf <- function(...) cat(sprintf(...))
library("optparse")
option_list = list(
    make_option(c("-f", "--file"), type="character", default=NULL,
                help="dataset file name", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="out.txt",
                help="output file name [default= %default]", metavar="character"),
    make_option(c("-a", "--alpha"), type="double", default=NULL,
                help="alpha level", metavar="character")
);
 
opt_parser = OptionParser(option_list=option_list);
args = parse_args(opt_parser);

printf("Reading lines from %s...\n", args$file)
printf("Output file destination: %s...\n", args$out)

print("Loading dataset...")
data <- read.table(args$file, header=T, sep='\t')
print("Starting analysis...")
get_fisher <- function(df) {
    n4 <- as.numeric(df["Common"])
    n2 <- as.numeric(df["Gene2Samples"]) - n4
    n3 <- as.numeric(df["Gene1Samples"]) - n4
    n1 <- 228 - (n2 + n3 + n4)
    mat <- matrix(c(n1,n2,n3,n4), ncol=2)
    f <- fisher.test(mat)
    return(f$p.value)
}
data$P_value <- apply(data, 1, get_fisher)

data$Adjusted_p <- p.adjust(data$P_value, method="BH")
data <- data[order(data$P_value),]
# TODO: alpha level
write.table(format(data, scientific=T, digits=3),
            file=args$out, sep="\t", quote=F, row.names=F)
