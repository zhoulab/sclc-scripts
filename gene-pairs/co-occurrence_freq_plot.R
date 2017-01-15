library("optparse")
library("ggplot2")
option_list <- list(
    make_option(c("-i", "--gene_pairs_file"), type="character",
                help="gene pairs filepath"),
    make_option(c("-o", "--out_file"), type="character",
                help="output PDF plot filepath"),
    make_option("--freq_table_out", type="character",
                help="frequency table output filepath")
);
opt_parser <- OptionParser(option_list=option_list);
args <- parse_args(opt_parser);

pretty_breaks_log <- function(n) {
    function(x) {
        axisTicks(log10(range(x, na.rm=TRUE)), log=TRUE, n=n)
    }
}

data <- read.table(args$gene_pairs_file, header=T, sep="\t")
pdf(args$out_file)
p <- ggplot(data, aes(Co_Occurrence)) +
     geom_histogram(binwidth=1) +
     scale_y_log10(breaks=pretty_breaks_log(n=5)) +
     xlab("Co-occurrence odds ratio") +
     ylab("Frequency") +
     ggtitle("Frequency of Co-Occurrences") +
     theme_bw() +
     theme(plot.title = element_text(hjust = 0.5),
           panel.grid.major=element_blank(),
           panel.grid.minor=element_blank())
p
dev.off()

num_breaks <- as.integer(max(data$Co_Occurrence)) + 1
hist_data <- hist(data$Co_Occurrence, breaks=num_breaks, plot=F)
breaks <- 0:(num_breaks - 1)
freq_table <- data.frame(breaks=breaks, counts=hist_data$counts)
write.table(freq_table, args$freq_table_out, quote=F, row.names=F, sep="\t")
