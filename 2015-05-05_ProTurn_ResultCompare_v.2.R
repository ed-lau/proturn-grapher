############# INSTRUCTIONS #############

# Note: v.2.x differs from v.1.x in that we now consider peptide-level comparisons and Welch's t-test first, before protein-level comparison.
# Protein-level comparison now includes the maximal Welch's t test p value as an estimate of protein-level significance.
# Confidence intervals of comparisons are now also included.
# Note: Reads in the ProTurn output files AFTER analysis by ProTurn Grapher (without reoptimization).
# ProTurn Grapher is necessary because it corrects the dk value from ProTurn output.
########################################


############## USER INPUT ##############
#R2_threshold <- 0.81             # R2 threshold to further filter out some peptides from Grapher output
#SE_threshold <- 0.1              # Standard error threshold to further filter out some peptides from Grapher output

#home_directory <- "~/Documents/Ping Lab/R Projects/proturn-grapher/"     # Working directory. Need back slach "/" at the end
#home_directory <- "~/Documents/Ping Lab/Project Files/2015 Isoform Turnover/" # Working directory. Need back slach "/" at the end
home_directory <- "~/Documents/Ping Lab/Project Files/2015 Paraquat Turnover/" # Working directory. Need back slach "/" at the end


#dataset1.directory <- "Data/ctrl/ctrl mouse heart cyto/"                              # The Proturn grapher output file for dataset 1, e.g., control hearts
#dataset2.directory <- "Data/ctrl/ctrl mouse heart mito/"                              # The Proturn grapher output file for dataset 2, e.g., iso hearts
dataset1.directory <- "Data/3. MitoTEMPO Mito/"                              # The Proturn grapher output file for dataset 1, e.g., control hearts
dataset2.directory <- "Data/4. Tempol Mito/"                              # The Proturn grapher output file for dataset 2, e.g., iso hearts



local_fasta_location <- "~/Documents/Ping Lab/Project Files/2015 Isoform Turnover/Fasta/Swissprot_Mouse_16689entries_20141219.fasta"

protein_comparison_output <- paste("ProTurn_Protein_Compare_", format(Sys.time(),"%Y-%m-%d %H.%M"), ".txt",sep="")                              # Name of the result file (comparing two proteins)
peptide_comparison_output <- paste("ProTurn_Peptide_Compare_", format(Sys.time(),"%Y-%m-%d %H.%M"), ".txt",sep="")              # Name of the peptide result file (comparing each common peptide for all proteins)
# Load annotation files
annotation.location <- "~/Documents/Ping Lab/Project Files/2015 Paraquat Turnover/Annotation file/Uniprot_musmusculus_annotations.csv"      # Annotation file

#################################################

############ LOAD LIBRARIES ###############
require("BSDA")
require("httr")
require("data.table")
require("dplyr")

############ LOAD INDIVIDUAL PEPTIDE OUTPUT FILES ###############
setwd(home_directory)
annot <- read.csv(annotation.location)

# Get hl-data.out
hl.data.1 <- read.table(paste(dataset1.directory,"hl-data.out",sep=""), header=TRUE, fill=TRUE)
hl.data.2 <- read.table(paste(dataset2.directory,"hl-data.out",sep=""), header=TRUE, fill=TRUE)

# Get Proturn Grapher output
peptide.output.1 <- read.table(paste(dataset1.directory,"Proturn_output.txt",sep=""), header=TRUE, fill=TRUE)
peptide.output.2 <- read.table(paste(dataset2.directory,"Proturn_output.txt",sep=""), header=TRUE, fill=TRUE)


proteins.in.output.1 <-  as.matrix(peptide.output.1$Uniprot)
proteins.in.output.2 <-  as.matrix(peptide.output.2$Uniprot)

protein.list <- unique(
                rbind(proteins.in.output.1,proteins.in.output.2)
                )

################################################

############## FUNCTIONS ##############
# This function calcualtes dk from the solution for partial derivatives dA/dk ; dA is SS.

# This is the normal KL function, which we will use to visualize the Proturn output k.
Rescaled_Model <-function(x){
        z <- 0
        k <- 2^(subset$k[j])
        N <- subset$N[j]
        kp <- subset$kp[j]
        pss <- subset$pss[j]
        a <- subset$a[j]
        for (n in 0:N) {
                b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
                bp <- (k/(k-n*kp))*b
                y<- (bp*exp(-n*kp*x)+exp(-k*x)*(1/(N+1)-bp))
                z <- y + z}
        final <- (z-1)/(((1-pss)^N)-1)
        return(final)
}

# Calculate Fractional Synthesis from the hl-data.out subset.
Calculate_FS <- function(x){
        A0. <- subset$a[j]
        Ainf. <- subset$a[j]*(1-subset$pss[j])^subset$N[j]
        FS. <- (x-A0.)/(Ainf.-A0.)
        return(FS.)
}

## This function gets individual peptide sequences through Uniprot (slow)
Get_sequence <- function(Uniprot,peptide){

        peptide <- gsub( " *\\(.*?\\) *", "", peptide)        
        fasta_address <- paste("http://www.uniprot.org/uniprot/", Uniprot, ".fasta",sep="")
        fasta <- readLines(fasta_address)
        fasta <- paste(
                fasta[2:length(fasta)]
                ,collapse="")
        protein_length <- nchar(fasta)

        position <- regexpr(peptide,fasta,fixed=TRUE)
        starting_pos <- position[1]
        result <- list(protein_length,starting_pos)
        }

# This function reads the local fasta file
Read_local_fasta <- function(file) {
        # Read the file line by line
        fasta<-readLines(file)
        # Identify header lines
        ind<-grep(">", fasta)
        # Identify the sequence lines
        s<-data.frame(ind=ind, from=ind+1, to=c((ind-1)[-1], length(fasta)))
        # Process sequence lines
        seqs<-rep(NA, length(ind))
        for(i in 1:length(ind)) {
                seqs[i]<-paste(fasta[s$from[i]:s$to[i]], collapse="")
                }
        #Create a data frame 
                DF<-data.frame(name=gsub(">", "", fasta[ind]), sequence=seqs)
        # Return the data frame as a result object from the function
        return(DF)
        }

## This function gets individual peptide sequences from the local fasta file (fast)
Get_local_sequence <- function(Uniprot,peptide) {
        
        peptide <- gsub( " *\\(.*?\\) *", "", peptide)        
        
        Local_sequence <- Local_fasta[grep(Uniprot,Local_fasta$name),2]
        
        protein_length <- nchar(as.character(Local_sequence))
        
        position <- regexpr(peptide,Local_sequence,fixed=TRUE)
        starting_pos <- position[1]
        result <- list(protein_length,starting_pos)
}


###################################################

############# Read the Local FASTA file ############
# Read the FASTA file
Local_fasta <- Read_local_fasta(local_fasta_location)
#####################################################

################ MAIN LOOP: PEPTIDE AND PROTEIN COMPARISON OUTPUT ####################
peptide.output <- data.frame()
protein.output <- data.frame()

# Graphical output parameters
pdf(file="Proturn_Peptide_Compare.pdf")
par(mfrow=c(5,3), mar=c(2.2,2.2,2,2))

#for(i in 1:20){
for (i in 1:nrow(protein.list)) {
        
        print(paste("Now running Protein # ", i, " of ", nrow(protein.list), ". ", round(i/nrow(protein.list)*100,2), "% done."))
        subset.output.1 <- peptide.output.1[ which(peptide.output.1$Uniprot == protein.list[i]), ]
        subset.output.2 <- peptide.output.2[ which(peptide.output.2$Uniprot == protein.list[i]), ]
        
        # Annotate gene name and protien name
        GN <- annot[ which(annot$Entry == protein.list[i]),"Gene.names"]
        PN <- annot[ which(annot$Entry == protein.list[i]),"Protein.names"]
        
        # Add a new column which concatenates peptide and charge to become a unique identifier for peptide 
        
        # First for Sample 1
        if (nrow(subset.output.1) > 0) {
                subset.output.1 <- mutate(subset.output.1, concat = paste(Peptide,z,sep=""))
                subset.output.1$protein_length <- NA
                subset.output.1$peptide_pos <- NA
                for (j in 1:nrow(subset.output.1)){
                        
                        # If you are using combined organelle files, also paste the first character of the modified IDs
                        #subset.output.1$concat[j] <- paste(subset.output.1$Peptide[j],subset.output.1$z[j],substr(subset.output.1$ID[j],1,1), sep="")
                        
                        # Getting protein length and sequence length for each peptide
                        sequence <- Get_local_sequence(subset.output.1$Uniprot[j],subset.output.1$Peptide[j])
                        subset.output.1$protein_length[j] <- unlist(sequence[1])
                        subset.output.1$peptide_pos[j] <- unlist(sequence[2])
                }
        }
        
        # Then for Sample 2
        if (nrow(subset.output.2) > 0) {
                subset.output.2 <- mutate(subset.output.2, concat = paste(Peptide,z,sep=""))
                subset.output.2$protein_length <- NA
                subset.output.2$peptide_pos <- NA
                for (j in 1:nrow(subset.output.2)){
                       
                        # If you are using combined organelle files, also paste the first character of the modified IDs
                        #subset.output.2$concat[j] <- paste(subset.output.2$Peptide[j],subset.output.2$z[j],substr(subset.output.2$ID[j],1,1), sep="")
                        
                        # Getting protein length and sequence length for each peptide
                        sequence <- Get_local_sequence(subset.output.2$Uniprot[j],subset.output.2$Peptide[j])
                        subset.output.2$protein_length[j] <- unlist(sequence[1])
                        subset.output.2$peptide_pos[j] <- unlist(sequence[2])
                }
        }
        
        # If there are exact peptide-charge matches between sample 1 and sample 2, subset out the common peptides.
        if (nrow(subset.output.1) > 0 & nrow(subset.output.2) >0 ) {
                # Subset out the list of common peptides from subset.output.1
                common.peptides.1 <- subset.output.1[ which(subset.output.1$concat %in% subset.output.2$concat), ]
                common.peptides.1 <- common.peptides.1[order(common.peptides.1$concat),]
                colnames(common.peptides.1) <- paste(colnames(common.peptides.1),".1",sep="")
                # Subset out the list of common peptides from subset.output.2
                common.peptides.2 <- subset.output.2[ which(subset.output.2$concat %in% subset.output.1$concat), ]
                common.peptides.2 <- common.peptides.2[order(common.peptides.2$concat),]
                colnames(common.peptides.2) <- paste(colnames(common.peptides.2),".2",sep="")
                
                # Combine the two tables into one!
                common.peptides <- cbind(common.peptides.1,common.peptides.2)
                
                # Calculate ratios of all common peptides between Sample 1 and Sample 2
                common.peptides$ratio <- common.peptides$k.2-common.peptides$k.1
                
                # If there are common peptides between the two Samples for this protein, plot graphs and calculate statistics.
                if  (nrow(common.peptides) > 0){
                        
                        # Preparing the first plot
                        plot(-10, -10, xlab="Pos", ylab="log2 Ratio", ylim=c(-2,2), xlim=c(0,max(common.peptides$protein_length.1)), 
                             cex=0.8, cex.axis=0.6, cex.lab=0.6, pch=16, ps=28, lwd=2, lty=1, mgp=c(1.55,0.48,0), las=1)
                        mtext(paste(common.peptides$Uniprot.1[1],
                                    GN)
                              ,side = 3, line = 0.5, cex=0.4);
                        
                        # Loop through every common peptides and calculate the Welch's t test ratio for each peptide pair
                        for (j in 1:nrow(common.peptides)){                                
                                # R's if-else syntax. Matching colors to modifications.
                                if(grepl("(42.0106)",common.peptides$Peptide.1[j]) == TRUE) {
                                        color <- "red"
                                } else if(grepl("(114.042927)",common.peptides$Peptide.1[j]) == TRUE){
                                        color <- "purple"
                                } else if (grepl("(79.9663)",common.peptides$Peptide.1[j]) == TRUE){
                                        color <- "green"
                                } else if (grepl("\\(",common.peptides$Peptide.1[j]) == TRUE){
                                        color <- "grey"
                                } else {
                                        color <- "black"
                                }
                                
                                # Plot out the segment ratio
                                segments(x0=common.peptides$peptide_pos.1[j],
                                         y0=common.peptides$ratio[j],
                                         x1=common.peptides$peptide_pos.1[j] + nchar(as.character(gsub( " *\\(.*?\\) *", "", common.peptides$Peptide.1[j]))),
                                         y1=common.peptides$ratio[j],
                                         lwd=4, col=color)
                                
                                # This is a function from the BSDA package
                                welch.t <- tsum.test(mean.x=common.peptides$k.1[j],   s.x=common.peptides$dk.1[j], n.x=common.peptides$DP.1[j],
                                                     mean.y=common.peptides$k.2[j], s.y=common.peptides$dk.2[j], n.y=common.peptides$DP.2[j])
                                
                                common.peptides$welch.t[j] <- welch.t$p.value
                                
                                # Write out individual peptide ratios and p value
                                
                                peptide_summary_table <- cbind(protein.list[i], as.character(GN), as.character(PN), common.peptides$protein_length.1[j],
                                                           as.character(common.peptides$Peptide.1[j]),common.peptides$z.1[j], common.peptides$peptide_pos.1[j],
                                                           common.peptides$k.1[j],common.peptides$k.2[j],
                                                           common.peptides$ratio[j], common.peptides$welch.t[j])
                                
                                colnames(peptide_summary_table) <- c("Uniprot","GN","PN","protein.length","sequence","z","pos","k1","k2","ratio","p")
                                
                                peptide.output <- rbind(peptide.output, peptide_summary_table)
                                }
                        
                        if (nrow(common.peptides) >2 ) {
                                T_shared <- tryCatch({t.test(common.peptides$k.1, common.peptides$k.2, paired=TRUE)},error=function(cond){return(NA)})
                                T_shared_results <- tryCatch({T_shared$p.value},error=function(cond){return(NA)})
                                
                                W_shared <- wilcox.test(common.peptides$k.1, common.peptides$k.2, paired=TRUE)
                                W_shared_results <- W_shared$p.value
                                }
                        else {
                                T_shared_results <- NA
                                W_shared_results <- NA
                        }
                        
                        # Write out the protein ratios
                protein_summary_table <- cbind(protein.list[i], as.character(GN), as.character(PN), nrow(common.peptides), median(common.peptides$k.1),
                                                       median(common.peptides$k.2), median(common.peptides$k.2)-median(common.peptides$k.1),max(as.numeric(common.peptides$welch.t)),
                                                       T_shared_results,W_shared_results)
                
                colnames(protein_summary_table) <- c("Uniprot","GN","PN","peptides","k1","k2","ratio","max.welch.p","t.p","w.p")
                protein.output <- rbind(protein.output, protein_summary_table)
                }      
        }
}


# Write the output file.
write.table(peptide.output, file=peptide_comparison_output, quote=FALSE, sep="\t", row.names=FALSE)
write.table(protein.output, file=protein_comparison_output, quote=FALSE, sep="\t", row.names=FALSE)
dev.off()

################## SECOND GRAPH (PROTEIN FIT) ####################

# Graphical output parameters
pdf(file="Proturn_fit.pdf")
par(mfrow=c(5,4), mar=c(2.2,2.2,2,2))

#Loop through all Uniprot IDs in the protein list, and extract all the peptides from either result files that belong to the protein.
#for (i in 6){
for (i in 1:nrow(protein.list)) {
        
        
        print(paste("Now graphing the turnover curves for protein # ", i, " of ", nrow(protein.list), ". ", round(i/nrow(protein.list)*100,2), "% done."))
        subset.output.1 <- peptide.output.1[ which(peptide.output.1$Uniprot == protein.list[i]), ]
        subset.output.1 <- mutate(subset.output.1, concat = paste(Peptide,z,sep=""))
        subset.output.2 <- peptide.output.2[ which(peptide.output.2$Uniprot == protein.list[i]), ]
        subset.output.2 <- mutate(subset.output.2, concat = paste(Peptide,z,sep=""))
        
        if (nrow(subset.output.1) > 0 & nrow(subset.output.2) >0 ) {
                # Subset out the list of common peptides from subset.output.1
                common.peptides.1 <- subset.output.1[ which(subset.output.1$concat %in% subset.output.2$concat), ]
                common.peptides.1 <- common.peptides.1[order(common.peptides.1$concat),]
                # Subset out the list of common peptides from subset.output.2
                common.peptides.2 <- subset.output.2[ which(subset.output.2$concat %in% subset.output.2$concat), ]
                common.peptides.2 <- common.peptides.2[order(common.peptides.2$concat),]
        
        
        # Create an empty plot
        plot(-2, -2, xlab="Time", ylab="FS", ylim=c(-0.1,1.1), xlim=c(0,15), cex=0.8, cex.axis=0.6, cex.lab=0.6, pch=16, ps=28, lwd=2, lty=1, mgp=c(1.55,0.48,0), las=1)
        mtext(protein.list[i],side = 3, line = 0.5, cex=0.4);
        
        # Plot the turnover graphs (kinetic curve and data points for file 1 and file 2) for comparison
        
        # For file 1:
        if (nrow(common.peptides.1) > 0) {
                subset <- common.peptides.1
                for(j in 1:nrow(common.peptides.1)) {
                        # Plot out each fitted kinetic curve
                        curve(Rescaled_Model(x), from=0, to=15, col="red", add=TRUE)
                        # Find all the data points linked to the peptide IDs in the protein being considered
                        # NOTE, changed from == to %in% because == didn't seem to work when the IDs were no longer numbers?
                        hl.data.peptide.subset <- hl.data.1[ which(hl.data.1$ID %in% subset$ID[j]),1:3] 
                        # Calculate their fractional synthesis
                        FS <- sapply(hl.data.peptide.subset$A0, Calculate_FS)
                        hl.data.peptide.subset <- cbind(hl.data.peptide.subset,FS)
                        # Plot out all data points
                        points(hl.data.peptide.subset$t,hl.data.peptide.subset$FS, col="red", ylim=c(-0.1,1.1), xlim=c(0,15), cex=0.8)      
                }
        }
        # For file 2:
        if (nrow(common.peptides.2) > 0) {
                subset <- common.peptides.2
                for(j in 1:nrow(common.peptides.2)){
                        # Plot out each fitted kinetic curve
                        curve(Rescaled_Model(x), from=0, to=15, col="blue", add=TRUE)
                        # Find all the data points linked to the peptide IDs in the protein being considered
                        hl.data.peptide.subset <- hl.data.2[ which(hl.data.2$ID %in% subset$ID[j]),1:3]
                        # Calculate their fractional synthesis
                        FS <- sapply(hl.data.peptide.subset$A0, Calculate_FS)
                        hl.data.peptide.subset <- cbind(hl.data.peptide.subset,FS)
                        # Plot out all data points
                        points(hl.data.peptide.subset$t,hl.data.peptide.subset$FS, col="blue", ylim=c(-0.1,1.1), xlim=c(0,15), cex=0.8)
                }
        }
        }
        
}
# Close graphical stream
dev.off()

