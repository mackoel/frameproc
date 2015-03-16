#!/usr/bin/Rscript --silent

library(optparse)

option_list <- list(
	make_option(c("-a", "--lowgreen"), type="double", default=0,
	help="threshold [default %default]",
	metavar="number"),
	make_option(c("-b", "--highgreen"), type="double", default=255,
	help="threshold [default %default]",
	metavar="number"),
	make_option(c("-c", "--lowred"), type="double", default=0,
	help="threshold [default %default]",
	metavar="number"),
	make_option(c("-d", "--highred"), type="double", default=255,
	help="threshold [default %default]",
	metavar="number"),
	make_option(c("-k", "--breaks"), type="integer", default=100,
	help="breaks [default %default]",
	metavar="number"),
	make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
	help="Print extra output [default=false]")
)
opt_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, description = "Post-proc object table", epilogue = "Send feedback to mackoel@gmail.com")
flaggs <- parse_args(opt_parser, args = commandArgs(trailingOnly = TRUE), print_help_and_exit = TRUE, positional_arguments = TRUE)
opts <- flaggs$options
args <- flaggs$args
infile <- args[1]
outfiles <- unlist(strsplit(args[2], ",", fixed = TRUE))
outfile <- outfiles[1]
outgraph <- outfiles[2]
a <- opts$lowgreen
b <- opts$highgreen
c <- opts$lowred
d <- opts$highred
q <- read.csv(infile)
colnames(q) <- sub("eea1", "green", colnames(q))
colnames(q) <- sub("egfr", "red", colnames(q))
colnames(q) <- sub("cort", "red", colnames(q))
pdf(file = outgraph)
hist(q$red_mean, breaks=seq(from=0,to=255,length=opts$breaks), main="red", xlim=c(0, opts$breaks))
hist(q$green_mean, breaks=seq(from=0,to=255,length=opts$breaks), main="green", xlim=c(0, opts$breaks))
q <- q[(q$green_mean >= a) & (q$green_mean < b) & (q$red_mean >= c) & (q$red_mean < d),]
write.csv(q, file = outfile)
hist(q$red_mean, breaks=seq(from=0,to=255,length=opts$breaks), main="red", xlim=c(0, opts$breaks))
hist(q$green_mean, breaks=seq(from=0,to=255,length=opts$breaks), main="green", xlim=c(0, opts$breaks))
res <- dev.off()

