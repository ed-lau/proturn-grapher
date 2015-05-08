############# INSTRUCTIONS #############
#Note: Reads in the ProTurn output files AFTER analysis by ProTurn Grapher, because it corrects the dk value from ProTurn output.
#If using old version of ProTurn (before 2.0.0.5) hl.txt and hl-data.txt in Excel, remove the non-peptides in hl.txt, fill all the holes, find gene names, and label the
#following columns: ID, Uniprot, GN, Peptide, DP, Charge, Isotopomer, SS (Sums of squares), a, pss, kp, N, k, and R2; ID, Time and A0.
# kp, pss, N and a. N is from Commerford and must be rounded.
########################################


############## USER INPUT ##############
R2_threshold <- 0.81                                                    # R2 threshold to further filter out some peptides from Grapher output
SS_threshold <- 0.1                                                     # Standard error threshold to further filter out some peptides from Grapher output

home_directory <- "~/Documents/Ping Lab/Heavy Water/ProTurn Output/"     # Working directory
dataset1.directory <- "test1/"                              # The Proturn grapher output file for dataset 1, e.g., control hearts
dataset2.directory <- "test2/"                              # The Proturn grapher output file for dataset 2, e.g., iso hearts

comparison_output <- "ProTurn_Compare.txt"                              # Name of the result file.
# Load annotation files
annotation.location <- "Annotation file/10090_annotations.csv"      # Annotation file

#################################################

############ LOAD INDIVIDUAL PEPTIDE OUTPUT FILES ###############
setwd(home_directory)
annot <- read.csv(file=annotation.location)

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

summary <- paste(
        "Uniprot", "GeneName","ProteinName",
        "peptides.1","k.median.1","k.mad.1","k.cv.1", 
        "peptides.2","k.median.2","k.mad.2","k.cv.2",
        "ratio", "Wilcoxon.test", "T test",
        sep = "\t")
write(summary, file=comparison_output, append=T)
################################################

############## FUNCTIONS ##############
# This function calcualtes dk from the solution for partial derivatives dA/dk ; dA is SS.
Calculate_dk <-function(x){
        z <- 0
        for (n in 0:dt$N[c]) {
                b <- factorial(dt$N[c])/(factorial(n)*factorial(dt$N[c]-n))*(1-dt$pss[c])^(dt$N[c]-n)*(dt$pss[c])^n
                bp <- (dt$k[c]/(dt$k[c]-n*dt$kp[c]))*b
                y<- dt$a[c]*((n*dt$kp[c])/(dt$k[c]*(dt$k[c]-n*dt$kp[c]))*bp*(exp(-dt$k[c]*x)-exp(-n*dt$kp[c]*x))-x*(1/(dt$N[c]+1)-bp)*exp(-dt$k[c]*x))                                                              
                z <- y + z}
        return(dt$SS[c]/z)}

Calculate_FS <- function(x){
        A0. <- dt$a[c]
        Ainf. <- dt$a[c]*(1-dt$pss[c])^dt$N[c]
        FS. <- (x-A0.)/(Ainf.-A0.)
        return(FS.)
}

# This is the normal KL function, which we will use to visualize the Proturn output k.

# This is the normal KL function, which we will use to visualize the Proturn output k.
Rescaled_Model <-function(x){
        z <- 0
        k <- subset$k[j]
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

########################################

######### GRAPHICAL OUTPUT PARAMETERS ##########
pdf(file="Proturn_fit.pdf")
par(mfrow=c(5,4), mar=c(2.2,2.2,2,2))
################################################
#Loop through all Uniprot IDs in the protein list, and extract all the peptides from either result files that belong to the protein.
for (i in 1:nrow(protein.list)) {
        
        print(paste("Now running Protein # ", i, " of ", nrow(protein.list), ". ", round(i/nrow(protein.list)*100,2), "% done."))
        subset.output.1 <- peptide.output.1[ which(peptide.output.1$Uniprot == protein.list[i]), ]
        subset.output.2 <- peptide.output.2[ which(peptide.output.2$Uniprot == protein.list[i]), ]
        
        protein.median.1 <- median(subset.output.1$k)
        protein.median.2 <- median(subset.output.2$k)

        # Calculate median absolute deviation
        if (length(subset.output.1$k) > 1) {
                protein.mad.1 <- mad(subset.output.1$k)
                }
        else {protein.mad.1 <- NA}

        if (length(subset.output.2$k) > 1) {
                protein.mad.2 <- mad(subset.output.2$k)
        }
        else { protein.mad.2 <- NA}
        
        # Calculate protein CV
        protein.cv.1 <- protein.mad.1/protein.median.1
        protein.cv.2 <- protein.mad.2/protein.median.2
        
        # Calculate protein ratio
        
        protein.ratio <- protein.median.2/protein.median.1

        GN <- annot[ which(annot$Uniprot == protein.list[i]),5]
        PN <- annot[ which(annot$Uniprot == protein.list[i]),3]
        
        # Perform Mann Whitney U test
        if (nrow(subset.output.1) >= 3 && nrow(subset.output.2) >= 3) {
                Utest <- wilcox.test(subset.output.1$k, subset.output.2$k)
                Ttest <- t.test(subset.output.1$k, subset.output.2$k)
        }
        else{
                Utest <- wilcox.test(1:3,1:3)
                Ttest <- t.test(1:3,1:3)
        }
        
        # Create an empty plot
        plot(-2, -2, xlab="Time", ylab="FS", ylim=c(-0.1,1.1), xlim=c(0,15), cex=0.8, cex.axis=0.6, cex.lab=0.6, pch=16, ps=28, lwd=2, lty=1, mgp=c(1.55,0.48,0), las=1, add=TRUE)
        mtext(protein.list[i],side = 3, line = 0.5, cex=0.4);
        # Plot the turnover graphs for comparison (no data points yet as that require the hl-data.out file)
        for(j in 1:nrow(subset.output.1)){
                subset <- subset.output.1
                if (nrow(subset) > 0){
                        curve(Rescaled_Model(x), from=0, to=15, col="red", add=TRUE)
                        }
        }
        
        for(j in 1:nrow(subset.output.2)){
                subset <- subset.output.2
                if (nrow(subset) > 0) {
                curve(Rescaled_Model(x), from=0, to=15, col="blue", add=TRUE)
                }
        }
        
        # Write out summary for output
        summary <- paste(protein.list[i], GN, PN, 
                         nrow(subset.output.1), protein.median.1, protein.mad.1, protein.cv.1, 
                         nrow(subset.output.2), protein.median.2, protein.mad.2, protein.cv.2, 
                         Utest$p.value, Ttest$p.value,
                         sep = "\t")
        write(summary, file=comparison_output, append=T)
        }

# Close graphical stream
dev.off()########################################