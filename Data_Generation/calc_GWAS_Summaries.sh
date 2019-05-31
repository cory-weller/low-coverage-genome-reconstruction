#!/usr/bin/env bash

module load R/3.3.0

while read population_type population_stem N_founders Population_Replicate Replicate N_variable_sites PVE MAF MAF_exact focal_snp_chr focal_snp_pos focal_snp_freq N_case_individuals N_control_individuals; do
    Rscript - ${population_type} ${population_stem} ${N_founders} ${Population_Replicate} ${Replicate} ${N_variable_sites} ${PVE} ${MAF} ${MAF_exact} ${focal_snp_chr} ${focal_snp_pos} ${focal_snp_freq} ${N_case_individuals} ${N_control_individuals} <<EOF

#!/usr/bin/env Rscript
library(data.table)
library(foreach)

args <- commandArgs(trailingOnly=TRUE)

# args <- c("outbred_dgrp","outbred","128","1","1","3620137","0.1","0.5","0.495082653274744","2L","19967459","0.504917346725256","2251","2327")

population_type <- args[[1]]
population_stem <- args[[2]]
N_founders <- as.numeric(args[[3]])
Population_Replicate <- as.numeric(args[[4]])
Replicate <- as.numeric(args[[5]])
N_variable_sites <- as.numeric(args[[6]])
PVE <- as.numeric(args[[7]])
MAF <- as.numeric(args[[8]])
MAF_exact <- as.numeric(args[[9]])
focal_snp_chr <- args[[10]]
focal_snp_pos <- as.numeric(args[[11]])
focal_snp_freq <- as.numeric(args[[12]])
N_case_individuals <- as.numeric(args[[13]])
N_control_individuals <- as.numeric(args[[14]])


read.freqs <- function(population_type, N_founders, MAF, PVE, Population_Replicate, Replicate) {
    #zipFileName <- paste(population_type, "_", N_founders, "_", MAF, "MAF", "_", PVE, "PVE", ".zip", sep="")
    #fread(paste("unzip -p ", zipFileName, " ", Population_Replicate, ".", Replicate, ".freqs", sep=""), colClasses=c("character","integer","double","double","double","double"), key=c("CHROM","POS"), header=TRUE)
    fread(paste(Population_Replicate, ".", Replicate, ".freqs", sep=""), colClasses=c("character","integer","double","double","double","double"), key=c("CHROM","POS"), header=TRUE)
}

get.chisq <- function(DT) {
    DT[,chisq := (((case.Ref)+(control.Ref)+((case.Alt))+((control.Alt))) * (((control.Ref)*(case.Alt)) - ((case.Ref)*(control.Alt)))**2) / (((case.Ref)+(control.Ref)) * (((case.Ref)+(case.Alt))) * ((control.Ref)+(control.Alt)) * (((case.Alt))+(control.Alt)))]
}

get.P <- function(DT) {
    DT[,P := pchisq(chisq, 1, lower.tail=FALSE)]
}

get.lambdaGC1000 <- function(chisq.array, N_case_individuals, N_control_individuals) {
	expected.median.chisq <- qchisq(0.5, df=1)
	unadjusted.lambda <- median(chisq.array) / expected.median.chisq
    lambda.1000 <- 1 + (unadjusted.lambda-1) * (1/N_case_individuals + 1/N_control_individuals)/(1/1000 + 1/1000)
	return(lambda.1000)
}

get.max_chisq <- function(DT) {
    return(max(DT[,chisq]))
}

get.top_snps <- function(DT, max_chisq) {
    return(DT[chisq==max_chisq])
}

get.distance_top_to_focal <- function(focal_snp_chr, focal_snp_pos, top_chr, top_pos) {
    if(focal_snp_chr == top_chr) {
        return(abs(top_pos - focal_snp_pos))
    } else {
        # If on different chromosome, snps are unlinked
        return(Inf)
    }
}

get.snp.clusters <- function(DT, focal_snp_chr, focal_snp_pos, threshold_vector=c(0, 1e3, 10e3, 100e3)) {
    # Iterate through distance thresholds. Do not run in parallel with %dopar%.
    o <- foreach(peak_width=threshold_vector, .combine="rbind") %do% {
    	# Initialize output data.table
    	# DT.results <- copy(DT[0, c("CHROM","POS","chisq","P")])
        DT.results <- copy(DT[0, c("CHROM","POS","chisq")])
        #DT.tmp <- copy(DT[chisq > 19.51142, c("CHROM","POS","chisq","P")])
        #DT.tmp <- copy(DT[chisq > 19.51142, c("CHROM","POS","chisq")])
        DT.tmp <- copy(DT[chisq>2, c("CHROM","POS","chisq")])

        if(nrow(DT.tmp) == 0) {
            DT.results <- data.table("CHROM"=NA, "POS"=NA, "chisq"=NA, "peak_width"=peak_width, "rank"=NA)
            return(DT.results)
        } else {

    	# Iteratively extract windows centered on most significant SNP until all clusters have been accounted for
        i <- 1
    	while(nrow(DT.tmp) > 0) {
    		# Determine most significant row; if ties, randomly select one

            if(i > 25) {break}

    		most_signif_chisq <- max(DT.tmp[,chisq])
    		most_signif_row <- copy(DT.tmp[chisq == most_signif_chisq][sample(.N, size=1)])
            DT.results <- rbindlist(list(DT.results, most_signif_row))

            # Remove SNPs +/- peak_width from this most significant SNP
            # SNPs on different chromosomes + SNPs at least <threshold> bp away
            DT.tmp <- DT.tmp[ CHROM != most_signif_row[,CHROM] | (CHROM == most_signif_row[,CHROM] & abs(POS-most_signif_row[,POS]) > peak_width)]

            i <- i + 1
        }
        DT.results[, "peak_width" := peak_width]
        DT.results[, rank := 1:.N]
        return(DT.results)
        }
    }

    o[, containsFocal := ifelse(CHROM==focal_snp_chr & abs(POS - focal_snp_pos) <= peak_width, TRUE, FALSE)][]
    return(o)
}

get.summary <- function(population_type, N_founders, MAF, PVE, Population_Replicate, Replicate, focal_snp_chr, focal_snp_pos, N_case_individuals, N_control_individuals) {
    dat.freqs <- read.freqs(population_type, N_founders, MAF, PVE, Population_Replicate, Replicate)
    get.chisq(dat.freqs)
    # get.P(dat.freqs)

    # If the population is inbred (homozygous), correct chisq statistic
    if(population_type=="inbred" ) {
        dat.freqs[, chisq := (chisq/2)]
    }

    # calculate chi-square rank
    dat.freqs[, chisq_rank := frank(-chisq, ties.method="random")]

    lambdaGC1000 <- get.lambdaGC1000(dat.freqs[,chisq], N_case_individuals, N_control_individuals)

    max_chisq <- get.max_chisq(dat.freqs)
    focal_chisq <- dat.freqs[.(focal_snp_chr, focal_snp_pos), chisq]
    focal_chisq_rank <- dat.freqs[.(focal_snp_chr, focal_snp_pos), chisq_rank]


    equally_top_snps <- get.top_snps(dat.freqs, max_chisq)

    top_snp <- equally_top_snps[sample(.N, size=1)]
    top_snp_chr <- top_snp[,CHROM]
    top_snp_pos <- top_snp[,POS]

    distance_from_top_to_focal <- get.distance_top_to_focal(focal_snp_chr, focal_snp_pos, top_snp_chr, top_snp_pos)

    clusters <- get.snp.clusters(dat.freqs, focal_snp_chr, focal_snp_pos)

output_summary <- data.table(
        population_type, N_founders, MAF, PVE, Population_Replicate, Replicate, focal_snp_chr, focal_snp_pos, N_case_individuals, N_control_individuals,
        "lambdaGC1000" = lambdaGC1000,
        "max_chisq" = max_chisq,
        "focal_chisq" = focal_chisq,
        "focal_rank" = focal_chisq_rank,
        "N_equally_top_snps" = nrow(equally_top_snps),
        "top_snp_chr" = top_snp[,CHROM],
        "top_snp_pos" = top_snp[,POS],
        "distance_from_top_to_focal_snp" = distance_from_top_to_focal
        )

write.table(cbind(population_type, N_founders, MAF, PVE, Population_Replicate, Replicate, clusters), file="clusters.tab", append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

write.table(output_summary, file="summaries.tab", append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

}

get.summary(population_type, N_founders, MAF, PVE, Population_Replicate, Replicate, focal_snp_chr, focal_snp_pos, N_case_individuals, N_control_individuals)

EOF

done < <(head -n 11 results.tab | tail -n 10)

# If zipped, change to <(unzip -p $zipfilename results.tab)
