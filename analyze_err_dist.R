#!/usr/bin/env Rscript
#
# visualizes the mismatch distribution table generated by merge-illumina-pairs script
#

suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))


# command line options
option_list <- list(
		make_option(c("-o", "--output_file_prefix"), default="unknown",
				help = "Output file prefix for visualization files [default \"%default\"]"),
		make_option("--title", default="No Title",
				help="Title for the output figure [default '%default']")
)

parser <- OptionParser(usage = "analyze_results.R [options] input_file", option_list=option_list,
		description="An R program to visualize error distribution along the reads")

arguments <- parse_args(parser, positional_arguments = TRUE)
options <- arguments$options

# check if the positional argument is set
if(length(arguments$args) != 1) {
	cat("Incorrect number of required positional arguments\n\n")
	print_help(parser)
	stop()
} else {
	input_file <- arguments$args
}

# check if the input file is accessible
if(file.access(input_file) == -1){
	stop(sprintf("Specified file '%s' does not exist", input_file))
}

# load data frame.
input_file <- '/Users/meren/Desktop/MBL/PrimerShift/test'
df <- as.data.frame(read.csv(input_file, header=TRUE, sep="\t"))
row.names <- df$oligo
col.names <- colnames(df)

#df <- transform(df, oligo=reorder(oligo, -diff) ) 
#df$oligo

k <- cor(df[df$bin == 'observed', ]$frequency, df[df$bin == 'null', ]$frequency, method='kendall')

P <- function(){
	p = ggplot(df, aes(x = pos, y = frequency, group = source, color = source))
	p <- p + geom_line(size=1, alpha=.75)
	p <- p + theme(axis.text.x = element_text(size = 16, angle=90, vjust=0.5))
	p <- p + theme(axis.text.y = element_text(size = 12))
	p <- p + theme(legend.position = 'bottom', legend.text=element_text(size=14))
	p <- p + labs(x='', y='Frequency')
	p <- p + scale_x_continuous(breaks=c(df$pos))
	p <- p + ggtitle(options$title)
	p <- p + coord_cartesian(ylim=c(0, max(df$frequency) * 1.1))
	p <- p + theme(plot.title = element_text(hjust=0, vjust=1))
	#p <- p + facet_grid(source ~ .)
	
	print(p)
}

P()

# gen PDF
pdf_output <- paste(options$output_file_prefix,".pdf",sep="")
pdf(pdf_output, width = 16, height = 6)
P()
sprintf("Lines PDF: '%s'", pdf_output)
