#!/usr/bin/Rscript --silent

library(optparse)

run_tests_once <- function(string_tmpl, time, ID, a, b, c, d, FILTER=FALSE) {
	fds <- list.files(path = infile, pattern = string_tmpl, all.files = FALSE, full.names = TRUE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE)
	if (FILTER) {
		res <- lapply(fds, FUN=function(fd) system(paste0("../frame_proc_filt.R -a ", a, " -b ", b, " -c ", c, " -d ", d, " ", basename(fd), " ", paste0(ID, "-", basename(fd)), " ", paste0(ID, "-", basename(fd), "-hist.pdf"))))
		res <- lapply(fds, FUN=function(fd) system(paste0("../frame_proc_obj.R -m ", time, " ", paste0(ID, "-", basename(fd)), " ", paste0(ID, "-obj-", basename(fd)))))
	} else {
		res <- lapply(fds, FUN=function(fd) system(paste0("../frame_proc_obj.R -m ", time, " ", basename(fd), " ", paste0(ID, "-obj-", basename(fd)))))
	}
	cat("OBJ run:", max(unlist(res)), "\n")
	obj_sum = sum(res == 0)
	res <- lapply(fds, FUN=function(fd) if (file.exists(paste0(ID, "-obj-", basename(fd)))) { system(paste0("../frame_proc_cells.R -m ", time, " ", paste0(ID, "-obj-", basename(fd)), " ", paste0(ID, "-cl-", basename(fd)))) } )
	cat("CL run:", max(unlist(res)), "\n")
	cl_sum = sum(res == 0)
	cat("Time=", time, "obj_sum=", obj_sum, "cl_sum=", cl_sum, "\n")
	if (cl_sum < 1) {
		return(NULL)
	}
	datatab <- do.call("rbind", lapply(fds, function(fd) if (file.exists(paste0(ID, "-cl-", basename(fd)))) { read.csv(paste0(ID, "-cl-", basename(fd))) } ))
	expres <- data.frame(t(colMeans(datatab)))
	return(expres)
}

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
	make_option(c("-f", "--filter"), action="store_true", default=FALSE,
	help="filter [default %default]",
	metavar="number"),
	make_option(c("-i", "--id"), type="character", default="id",
	help="ID [default %default]",
	metavar="character"),
	make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
	help="Print extra output [default=false]")
)
opt_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, description = "Post-proc object table", epilogue = "Send feedback to mackoel@gmail.com")
flaggs <- parse_args(opt_parser, args = commandArgs(trailingOnly = TRUE), print_help_and_exit = TRUE, positional_arguments = TRUE)
opts <- flaggs$options
args <- flaggs$args
infile <- args[1]
outfile <- args[2]
outgraph <- args[3]
ID <- opts$id
a <- opts$lowgreen
b <- opts$highgreen
c <- opts$lowred
d <- opts$highred
FILTER <- opts$filter
expres <- NULL
string_tmpl <- paste("^C.*-", "CC_merged_counts.csv$", sep="")
das <- run_tests_once(string_tmpl, "CC", ID, a, b, c, d, FILTER)
if (!is.null(das)) {
	expres <- das
}

string_tmpl <- paste("^1.*-", "15_merged_counts.csv$", sep="")
das <- run_tests_once(string_tmpl, 15, ID, a, b, c, d, FILTER)
if (!is.null(das)) {
	expres <- rbind(expres, das)
}

string_tmpl <- paste("^3.*-", "30_merged_counts.csv$", sep="")
das <- run_tests_once(string_tmpl, 30, ID, a, b, c, d, FILTER)
if (!is.null(das)) {
	expres <- rbind(expres, das)
}

string_tmpl <- paste("^6.*-", "60_merged_counts.csv$", sep="")
das <- run_tests_once(string_tmpl, 60, ID, a, b, c, d, FILTER)
if (!is.null(das)) {
	expres <- rbind(expres, das)
}

if (is.null(expres)) {
	cat("NULL TABLE!\n")
	q()
}

write.csv(file=outfile, expres, row.names=F)

pdf(file = outgraph)

plot(expres$time,expres$dco,main="Mean intensity green vs red",xlab="Time",ylab="Mean correlation per cell", pch=15, type = 'o')
minman = min(c(expres$manders1, expres$manders2))
maxman = max(c(expres$manders1, expres$manders2))
plot(expres$time,expres$manders1,main="Mean Manders coeffs",xlab="Time",ylab="Mean Coeff", pch=15, type = 'o', ylim = c(minman, maxman))
lines(expres$time,expres$manders2,main="Mean Manders coeffs",xlab="Time",ylab="Mean Coeff", pch=17, type = 'o')
legend("topright", c("green(red>0)", "red(green>0)"), pch = c(15, 17), border = NA)
plot(expres$time,expres$mco,main="Number of pixels",xlab="Time",ylab="Ratio overlap/mask per cell", pch=15, type = 'o')
plot(expres$time,expres$nco,main="Median intensity green vs red",xlab="Time",ylab="Mean correlation per cell", pch=15, type = 'o')
plot(expres$time,expres$cco,main="Number of pixels in mask green vs red",xlab="Time",ylab="Mean correlation per cell", pch=15, type = 'o')

#counts <- table(expres$vol_counts1, expres$vol_counts2, expres$vol_counts3, expres$vol_counts4, expres$vol_counts5, expres$vol_counts6, expres$vol_counts7, expres$vol_counts8, expres$vol_counts9, expres$vol_counts10)
#expres
#counts <- table(expres$vol_counts1, expres$vol_counts2)
#counts
#barplot(counts, main="Share of overlap", xlab="Time", col=c("darkblue","red"), legend = rownames(counts), beside=TRUE)
cat("EGFR\n");
kounts <- cbind(expres$red_counts1, expres$red_counts2, expres$red_counts3, expres$red_counts4, expres$red_counts5, expres$red_counts6, expres$red_counts7, expres$red_counts8, expres$red_counts9, expres$red_counts10)
rwm <- rowSums(kounts)
kounts <- kounts/rwm
barplot(main="Ratio overlap/red", names.arg=c("CC", "15 min", "30 min", "60 min"), t(kounts),beside=T)
expres[,c("red_counts1", "red_counts2", "red_counts3", "red_counts4", "red_counts5", "red_counts6", "red_counts7", "red_counts8", "red_counts9", "red_counts10")] <- kounts
cat("EEA1\n");
kounts <- cbind(expres$green_counts1, expres$green_counts2, expres$green_counts3, expres$green_counts4, expres$green_counts5, expres$green_counts6, expres$green_counts7, expres$green_counts8, expres$green_counts9, expres$green_counts10)
rwm <- rowSums(kounts)
kounts <- kounts/rwm
barplot(main="Ratio overlap/green", names.arg=c("CC", "15 min", "30 min", "60 min"), t(kounts),beside=T)
expres[,c("green_counts1", "green_counts2", "green_counts3", "green_counts4", "green_counts5", "green_counts6", "green_counts7", "green_counts8", "green_counts9", "green_counts10")] <- kounts
cat("VOL\n");
kounts <- cbind(expres$vol_counts1, expres$vol_counts2, expres$vol_counts3, expres$vol_counts4, expres$vol_counts5, expres$vol_counts6, expres$vol_counts7, expres$vol_counts8, expres$vol_counts9, expres$vol_counts10)
rwm <- rowSums(kounts)
kounts <- kounts/rwm
barplot(main="Ratio overlap/obj", names.arg=c("CC", "15 min", "30 min", "60 min"), t(kounts),beside=T)
expres[,c("vol_counts1", "vol_counts2", "vol_counts3", "vol_counts4", "vol_counts5", "vol_counts6", "vol_counts7", "vol_counts8", "vol_counts9", "vol_counts10")] <- kounts
cat("CVOL\n");
kounts <- cbind(expres$cvol_counts1, expres$cvol_counts2, expres$cvol_counts3, expres$cvol_counts4, expres$cvol_counts5, expres$cvol_counts6, expres$cvol_counts7, expres$cvol_counts8, expres$cvol_counts9, expres$cvol_counts10)
rwm <- rowSums(kounts)
kounts <- kounts/rwm
barplot(main="Ratio red/obj", names.arg=c("CC", "15 min", "30 min", "60 min"), t(kounts),beside=T)
expres[,c("cvol_counts1", "cvol_counts2", "cvol_counts3", "cvol_counts4", "cvol_counts5", "cvol_counts6", "cvol_counts7", "cvol_counts8", "cvol_counts9", "cvol_counts10")] <- kounts
cat("EVOL\n");
kounts <- cbind(expres$evol_counts1, expres$evol_counts2, expres$evol_counts3, expres$evol_counts4, expres$evol_counts5, expres$evol_counts6, expres$evol_counts7, expres$evol_counts8, expres$evol_counts9, expres$evol_counts10)
rwm <- rowSums(kounts)
kounts <- kounts/rwm
barplot(main="Ratio green/obj", names.arg=c("CC", "15 min", "30 min", "60 min"), t(kounts),beside=T)
expres[,c("evol_counts1", "evol_counts2", "evol_counts3", "evol_counts4", "evol_counts5", "evol_counts6", "evol_counts7", "evol_counts8", "evol_counts9", "evol_counts10")] <- kounts

barplot(main="Ratio overlap/red", names.arg=c("CC", "15 min", "30 min", "60 min"), t(cbind(expres$red_counts1, expres$red_counts2+expres$red_counts3+expres$red_counts4+ expres$red_counts5, expres$red_counts6+expres$red_counts7+expres$red_counts8+expres$red_counts9, expres$red_counts10)),beside=T)
barplot(main="Ratio overlap/green", names.arg=c("CC", "15 min", "30 min", "60 min"), t(cbind(expres$green_counts1, expres$green_counts2+expres$green_counts3+expres$green_counts4+expres$green_counts5, expres$green_counts6+expres$green_counts7+expres$green_counts8+expres$green_counts9, expres$green_counts10)),beside=T)
barplot(main="Ratio overlap/obj", names.arg=c("CC", "15 min", "30 min", "60 min"), t(cbind(expres$vol_counts1, expres$vol_counts2+expres$vol_counts3+expres$vol_counts4+expres$vol_counts5, expres$vol_counts6+expres$vol_counts7+expres$vol_counts8+expres$vol_counts9, expres$vol_counts10)),beside=T)
barplot(main="Ratio red/obj", names.arg=c("CC", "15 min", "30 min", "60 min"), t(cbind(expres$cvol_counts1, expres$cvol_counts2+expres$cvol_counts3+expres$cvol_counts4+expres$cvol_counts5, expres$cvol_counts6+expres$cvol_counts7+expres$cvol_counts8+expres$cvol_counts9, expres$cvol_counts10)),beside=T)
barplot(main="Ratio green/obj", names.arg=c("CC", "15 min", "30 min", "60 min"), t(cbind(expres$evol_counts1, expres$evol_counts2+expres$evol_counts3+expres$evol_counts4+expres$evol_counts5, expres$evol_counts6+expres$evol_counts7+expres$evol_counts8+expres$evol_counts9, expres$evol_counts10)),beside=T)

counts <- data.frame(t(expres))
colnames(counts) <- expres$time

barplot(main="CC", names.arg=c("overlap/red", "overlap/green", "overlap/obj", "red/obj", "green/obj"), cbind(counts[7:16,]$"-1", counts[17:26,]$"-1", counts[27:36,]$"-1", counts[37:46,]$"-1", counts[47:56,]$"-1"),beside=T)
barplot(main="15 min", names.arg=c("overlap/red", "overlap/green", "overlap/obj", "red/obj", "green/obj"), cbind(counts[7:16,]$"15", counts[17:26,]$"15", counts[27:36,]$"15", counts[37:46,]$"15", counts[47:56,]$"15"),beside=T)
barplot(main="30 min", names.arg=c("overlap/red", "overlap/green", "overlap/obj", "red/obj", "green/obj"), cbind(counts[7:16,]$"30", counts[17:26,]$"30", counts[27:36,]$"30", counts[37:46,]$"30", counts[47:56,]$"30"),beside=T)
barplot(main="60 min", names.arg=c("overlap/red", "overlap/green", "overlap/obj", "red/obj", "green/obj"), cbind(counts[7:16,]$"60", counts[17:26,]$"60", counts[27:36,]$"60", counts[37:46,]$"60", counts[47:56,]$"60"),beside=T)

counts2 <- counts
counts2[8,] <- counts[8,] + counts[9,] + counts[10,] + counts[11,]
counts2[9,] <- counts[12,] + counts[13,] + counts[14,] + counts[15,]
counts2[10,] <- counts[16,]
counts2[11,] <- counts[17,]
counts2[12,] <- counts[18,] + counts[19,] + counts[20,] + counts[21,]
counts2[13,] <- counts[22,] + counts[23,] + counts[24,] + counts[25,]
counts2[14,] <- counts[26,]
counts2[15,] <- counts[27,]
counts2[16,] <- counts[28,] + counts[29,] + counts[30,] + counts[31,]
counts2[17,] <- counts[32,] + counts[33,] + counts[34,] + counts[35,]
counts2[18,] <- counts[36,]
counts2[19,] <- counts[37,]
counts2[20,] <- counts[38,] + counts[39,] + counts[40,] + counts[41,]
counts2[21,] <- counts[42,] + counts[43,] + counts[44,] + counts[45,]
counts2[22,] <- counts[46,]
counts2[23,] <- counts[47,]
counts2[24,] <- counts[48,] + counts[49,] + counts[50,] + counts[51,]
counts2[25,] <- counts[52,] + counts[53,] + counts[54,] + counts[55,]
counts2[26,] <- counts[56,]

counts2 <- counts2[1:26,]


barplot(main="CC", names.arg=c("overlap/red", "overlap/green", "overlap/obj", "red/obj", "green/obj"), cbind(counts2[7:10,]$"-1", counts2[11:14,]$"-1", counts2[15:18,]$"-1", counts2[19:22,]$"-1", counts2[23:26,]$"-1"),beside=T)
barplot(main="15 min", names.arg=c("overlap/red", "overlap/green", "overlap/obj", "red/obj", "green/obj"), cbind(counts2[7:10,]$"15", counts2[11:14,]$"15", counts2[15:18,]$"15", counts2[19:22,]$"15", counts2[23:26,]$"15"),beside=T)
barplot(main="30 min", names.arg=c("overlap/red", "overlap/green", "overlap/obj", "red/obj", "green/obj"), cbind(counts2[7:10,]$"30", counts2[11:14,]$"30", counts2[15:18,]$"30", counts2[19:22,]$"30", counts2[23:26,]$"30"),beside=T)
barplot(main="60 min", names.arg=c("overlap/red", "overlap/green", "overlap/obj", "red/obj", "green/obj"), cbind(counts2[7:10,]$"60", counts2[11:14,]$"60", counts2[15:18,]$"60", counts2[19:22,]$"60", counts2[23:26,]$"60"),beside=T)
cat("SUCCESS\n");
res <- dev.off()

