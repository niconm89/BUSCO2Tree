### Program to filter neutral sites from aligned CDS.
### Place fasta input files (one per gene) inside an input directory and, if codon usage bias estimates are wanted, also place a table file with the Most Frequent Codons info
### Run this script from Terminal e.g. with:
### 'Rscript --vanilla neutral_filter.R -i Input -l no -g no -c no -m MFC_dros_nuc.txt -o Output


### LIBRARIES
library("optparse")
library("stringr")
library(seqinr)
library(Biostrings)


### ARGUMENTS
option_list <- list(
  make_option(c("--input", "-i"), default=NULL,
              help = "input directory name (with a fasta file per gen)"),
  make_option(c("--length_check", "-l"), default='no',
              help = "want to check seqs length? if 'yes' files will be tested for multiple of 3 and equal length, and main programs won't be executed"),
  make_option(c("--non_gact", "-g"), default='no',
              help = "non-GACT sites may be there? if 'yes' unresolved codons will be filter out as output and main programs won't be executed"),
  make_option(c("--cub", "-c"), default='yes',
              help = "want codon usage bias estimates? 'yes' or 'no'"),
  make_option(c("--mfc", "-m"), default=NULL,
              help = "most frequent codons file name"),
  make_option(c("--output", "-o"), default="Output",
              help = "output directory name")
)
opt <- parse_args(OptionParser(option_list=option_list))
# Validating arguments
if (is.null(opt$input) || !file.exists(opt$input)) {
  stop(paste0("Input file not provided or does not exist. --input ", opt$input))
}
if (opt$length_check != "yes" && opt$length_check != "no") {
	stop(paste0("Must choose 'yes' or 'no' for --length_check ", opt$length_check))
}
if (opt$length_check == "yes") {
	print("As seqs length need to be checked -c and -m parameters are not being evaluated")
}
if (opt$non_gact != "yes" && opt$non_gact != "no") {
	stop(paste0("Must choose 'yes' or 'no' for --non_gact ", opt$non_gact))
}
if (opt$non_gact == "yes") {
	print("As non-GACT sites need to be removed -c and -m parameters are not being evaluated")
}
if (opt$length_check == "no" && opt$non_gact == "no") {
 	if (opt$cub != "yes" && opt$cub != "no") {
			stop(paste0("Must choose 'yes' or 'no' for --cub ", opt$cub))
	}
	if (opt$cub == "yes") {
		if (is.null(opt$mfc) || !file.exists(opt$mfc)) {
			stop(paste0("Most frequent codons file not provided or does not exist. --mfc", opt$mfc))
		}
	}
}


### CONSTANTS
# Read Files
input_dir <- opt$input
input_files <- list.files(input_dir)
mfc_file <- opt$mfc
# Write Files
output_dir <- opt$output


### MAIN PROGRAM: 4-FOLD DEGENERATED SITES FILTERING
dir.create(output_dir)
N <- length(input_files)
FFDC <- c ("CT", "GT", "TC", "CC", "AC", "GC", "CG", "AG", "GG")

#for (i in 1:N) {
#	myseqs <- readDNAStringSet(file=paste(input_dir, input_files[i], sep="/"), format="fasta")
#	sps <- names(myseqs)
#	M <- length(sps)
#	pos <- vector(mode="list", length=M)
#	for (j in 1:M) {
#		myNchar <- s2c(toString(myseqs[j]))
#		n <- width(myseqs)[1]/3
#		myC2char <- character(n)
#		for (k in 1:n) {	
#			myC2char[k] <- paste(myNchar[(k*3-2):(k*3-1)], collapse="")
#		}
#		pos[[j]] <- which(myC2char %in% FFDC)
#		rm(myNchar, myC2char)
#	}
#	rm(myseqs)
#	pos <- Reduce(intersect, pos)
#	myseqs <- read.fasta(file=paste(input_dir, input_files[i], sep="/"), as.string=TRUE)
#	for (j in 1:M) {
#		mystring <- s2c(myseqs[[j]])
#		myseqs[[j]] <- c2s(mystring[pos*3])
#	}
#	rm(mystring)
#	write.fasta(myseqs, names=names(myseqs), as.string=TRUE, file.out=paste(output_dir, input_files[i], sep="/"))
#	rm(myseqs)
#	print(paste0("Filtering Sites progress ", round(i*100/N, 4), "%"))
#}



### CODON USAGE BIAS PROGRAM
if (opt$cub == "yes") {
	MFC_table <- read.table(mfc_file, header=T)
	average <- numeric(N)
	for (i in 1:N) {
		myseqs <- readDNAStringSet(file=paste(input_dir, input_files[i], sep="/"), format="fasta")
		sps <- names(myseqs)
		M <- length(sps)
		print(paste(input_dir, input_files[i], sep="/"))
		myprots <- translate(myseqs)
		resultado <- numeric(M)
		for (j in 1:M) {
			myNchar <- s2c(toString(myseqs[j]))
			myPchar <- s2c(toString(myprots[j]))
			pos <- which(myPchar == "M" | myPchar == "W")
			if (length(pos) > 0) {
				myPchar <- myPchar[-pos]
				myNchar <- myNchar[-c(pos*3-2, pos*3-1, pos*3)]
			}
			n <- length(myPchar)
			myCchar <- character(n)
			for (k in 1:n) {	
				myCchar[k] <- paste(myNchar[(k*3-2):(k*3)], collapse="")
			}
			rm(myNchar)
			res <- logical(n)
			for (k in 1:n) {
				for (l in 1:21) {
					if (myPchar[k] == as.character(MFC_table[l,1]) && myCchar[k] == as.character (MFC_table[l,2])) {
						res[k] <- T
					}
				}			
			}
			resultado[j] <- sum(res)/n
		}
		rm(myseqs, myprots)
		average[i] <- mean(resultado)
		print(paste0("Codon Usage Bias progress ", round(i*100/N, 4), "%"))
	}	
	print("Creating final report")
	tabla <- data.frame(id=input_files, bias=average, low_bias_0.1quantil=logical(N), low_bias_0.2quantil=logical(N))
	for (i in 1:N) {
		if (tabla$bias[i] <= quantile(average, probs=0.1)) {
			tabla$low_bias_0.1quantil[i] <- T
		}
		if (tabla$bias[i] <= quantile(average, probs=0.2)) {
			tabla$low_bias_0.2quantil[i] <- T
		}
	}
	write.table(tabla, file=paste(output_dir, "Mean_Codon_Usage_Bias.txt", sep="/"), row.names=FALSE, quote=FALSE, sep="\t")
}
