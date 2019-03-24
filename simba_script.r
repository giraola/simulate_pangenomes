library(simba)
library(seqinr)
library(vegan)

# Perform simulations

set.seed(3241)

sim0 <- simpg(ref='pan_genome_reference.fa',norg=10,ngenes=100,dir_out='s0',ne=1e+06)
sim1 <- simpg(ref='pan_genome_reference.fa',norg=10,ngenes=100,dir_out='s1',ne=2e+06)
sim2 <- simpg(ref='pan_genome_reference.fa',norg=10,ngenes=100,dir_out='s2',ne=5e+06)
sim3 <- simpg(ref='pan_genome_reference.fa',norg=10,ngenes=100,dir_out='s3',ne=1e+07)
sim4 <- simpg(ref='pan_genome_reference.fa',norg=10,ngenes=100,dir_out='s4',ne=2e+07)
sim5 <- simpg(ref='pan_genome_reference.fa',norg=10,ngenes=100,dir_out='s5',ne=5e+07)
sim6 <- simpg(ref='pan_genome_reference.fa',norg=10,ngenes=100,dir_out='s6',ne=1e+08)
sim7 <- simpg(ref='pan_genome_reference.fa',norg=10,ngenes=100,dir_out='s7',ne=2e+08)
sim8 <- simpg(ref='pan_genome_reference.fa',norg=10,ngenes=100,dir_out='s8',ne=5e+08)


# Accessory genome

s0 <- as.vector(as.matrix(vegdist(sim0$panmatrix,method='jaccard')))
s1 <- as.vector(as.matrix(vegdist(sim1$panmatrix,method='jaccard')))
s2 <- as.vector(as.matrix(vegdist(sim2$panmatrix,method='jaccard')))
s3 <- as.vector(as.matrix(vegdist(sim3$panmatrix,method='jaccard')))
s4 <- as.vector(as.matrix(vegdist(sim4$panmatrix,method='jaccard')))
s5 <- as.vector(as.matrix(vegdist(sim5$panmatrix,method='jaccard')))
s6 <- as.vector(as.matrix(vegdist(sim6$panmatrix,method='jaccard')))
s7 <- as.vector(as.matrix(vegdist(sim7$panmatrix,method='jaccard')))
s8 <- as.vector(as.matrix(vegdist(sim8$panmatrix,method='jaccard')))


# Jaccard Boxplot (Fig. 1A)

pdf('boxplot_accessory.pdf')

boxplot(s0,s1,s2,s3,s4,s5,s6,s7,s8)

dev.off()

# Mutations

simus <- paste0('sim',c(0,1,2,3,4,5,6,7,8))

allres <- NULL

for (s in simus) {
	
	cb       <- combn(paste0('genome', 1:10), 2)
	sim      <- get(s)
	tdep     <- coalescent.intervals(sim$coalescent)$total.depth
	observed <- expected <- rep(NA, dim(cb)[2])

	for (j in 1:dim(cb)[2]) {
  
 		gg          <- cb[, j]
 		dd          <- cophenetic(sim$coalescent)[gg[1], gg[2]]
 		expected[j] <- 1 - exp( -attr(sim, 'mu') * dd * (attr(sim, 'ne') / tdep) )
  
  		# Which ogs have both genomes?
  
  		ogs <- which(colSums(sim$panmatrix[c(gg[1], gg[2]), ])==2)
  		rr  <- rl <- rep(NA, length(ogs))

  		for (i in 1:length(ogs)) {
    
    		fi    <- paste0(attr(sim, 'dir_out'),'/', names(ogs[i]), '.fasta')
    		ge    <- read.fasta(fi, seqtype = 'DNA')[paste0(gg, '_', names(ogs[i]))]
    		m     <- do.call(rbind, ge)
    		dm    <- dim(m)
    		subs  <- which(apply(m, 2, function(X) {X[1] != X[2]}))
    		rr[i] <- length(subs)
    		rl[i] <- dm[2]
    	}

    	observed[j] <- sum(rr) / sum(rl)
    }

	allres <- rbind(allres,data.frame(genomeA = cb[1, ],genomeB = cb[2, ],observed = observed,expected = expected))
}


fit <- lm(observed ~ expected, data = df)


# Figure 1B

pdf('lm.pdf')

plot(observed ~ expected, data=df, pch=19, col='grey')

abline(fit, lty=2)

dev.off()
