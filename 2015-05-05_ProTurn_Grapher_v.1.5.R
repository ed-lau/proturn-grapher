############# INSTRUCTIONS #############
#Note: Read in hl-out and hl-data.out directly. 
#If using old version of ProTurn (before 2.0.0.5) hl.txt and hl-data.txt in Excel, remove the non-peptides in hl.txt, fill all the holes, find gene names, and label the
#following columns: ID, Uniprot, GN, Peptide, DP, Charge, Isotopomer, SS (Sums of squares), a, pss, kp, N, k, and R2; ID, Time and A0.
# kp, pss, N and a. N is from Commerford and must be rounded.
########################################

############## USER INPUT ##############
R2_threshold <- 0.9 # 0.9    # R2 Threshold
SE_threshold <- 0.01 # 0.01    # Standard Error of Estimate Threshold

home_directory <- "~/Documents/Ping Lab/Project Files/2014 Drosophila Turnover/Data/yw mito/"
#home_directory <- "~/Documents/Ping Lab/Project Files/2015 Isoform Turnover/Data/hmdp/c57 iso"
#home_directory <- "~/Documents/Ping lab/Project Files/2015 Paraquat Turnover/Data/4. Tempol Mito"


hl.out_location <- "hl.out"
hl.data.out_location <- "hl-data.out"
output_file <- "ProTurn_Output.txt"

Refit <- FALSE           # Should the Grapher attempt to gather the data points for each peptide and refit the kinetic curve in R.
Is_fly <- TRUE           # Is this a fly (free fitting) experiment; if so, will require the plateaued time point.
plateau_point <- 42      # Plateau time point

#local_fasta_location <- "~/Documents/Ping Lab/Project Files/2015 Isoform Turnover/Fasta/Swissprot_Mouse_16689entries_20141219.fasta"
#local_fasta_location <- "~/Documents/Ping Lab/Project Files/2015 Isoform Turnover/Fasta/Swissprot_Human_20187entries_20141228.fasta"
local_fasta_location <- "~/Documents/Ping Lab/Project Files/2014 Drosophila Turnover/Fasta/Uniprot_Drosophila_42360entries_20150717.fasta"

annotation.location <- "~/Documents/Ping Lab/Project Files/2015 Paraquat Turnover/Annotation file/Uniprot_musmusculus_annotations.csv"
#annotation.location <- "~/Documents/Ping Lab/Heavy Water/ProTurn Output/Annotation file/9606_annotations.csv"

########################################

############## LOADING FILES ##############
setwd(home_directory)
annot <- read.csv(file=annotation.location)

oput <- paste("ID","Uniprot","Peptide","z","R2","SS","SE","k","dk","DP","a","pss","kp","N", sep = "\t")
write(oput, file=output_file, append=T)

# Preparing the graph
pdf(file="Proturn_fit.pdf")
par(mfrow=c(5,4), mar=c(2.2,2.2,2,2))

# Reading hl.out and hl-data.out
hl.out <- read.table(hl.out_location, header=TRUE, fill=TRUE)
dp <- read.table(hl.data.out_location, header=TRUE, fill=TRUE)
dp <- dp[,1:3]

# Calculate the Standard Error of Estimate from SS (Sums of Squares) based on SE = (SS/DP-1)^0.5
hl.out$SE <- (hl.out$SS/(hl.out$DP-1))^0.5

# Subsetting the hl.out table to house only those peptides that pass the threshold
dt <- hl.out[ which(hl.out$R2 >= R2_threshold | hl.out$SE <= SE_threshold),] # Subsetting by R2 and SS

########################################

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
Model <- function(t){
        z <- 0
        for (n in 0:dt$N[c]) {
                k <- dt$k[c]
                N <- dt$N[c]
                kp <- dt$kp[c]
                pss <- dt$pss[c]
                a <- dt$a[c]
                
                b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
                bp <- (k/(k-n*kp))*b
                y<- a*(bp*exp(-n*kp*t)+exp(-k*t)*(1/(N+1)-bp))                                                                  
                z <- y + z}
        return(z)
        }

# This is the normal KL function, which we will use to visualize the upper bound.
UpperModel <-function(x){
        z <- 0
        for (n in 0:dt$N[c]) {
                k <- dt$k[c]
                N <- dt$N[c]
                kp <- dt$kp[c]
                pss <- dt$pss[c]
                a <- dt$a[c]
                
                b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
                bp <- (Upper/(Upper-n*kp))*b
                y<- a*(bp*exp(-n*kp*x)+exp(-Upper*x)*(1/(dt$N[c]+1)-bp))                                                                  
                z <- y + z}
        return(z)}

# This is the normal KL function, which we will use to visualize the lower bound.
LowerModel <-function(x){
        z <- 0
        for (n in 0:dt$N[c]) {
                b <- factorial(dt$N[c])/(factorial(n)*factorial(dt$N[c]-n))*(1-dt$pss[c])^(dt$N[c]-n)*(dt$pss[c])^n
                bp <- (Lower/(Lower-n*dt$kp[c]))*b
                y<- dt$a[c]*(bp*exp(-n*dt$kp[c]*x)+exp(-Lower*x)*(1/(dt$N[c]+1)-bp))                                                                  
                z <- y + z}
        return(z)}

# This is the normal KL function, which we will use to visualize the Proturn output k.
Rescaled_Model <-function(x){
        z <- 0
        for (n in 0:dt$N[c]) {
                b <- factorial(dt$N[c])/(factorial(n)*factorial(dt$N[c]-n))*(1-dt$pss[c])^(dt$N[c]-n)*(dt$pss[c])^n
                bp <- (dt$k[c]/(dt$k[c]-n*dt$kp[c]))*b
                y<- (bp*exp(-n*dt$kp[c]*x)+exp(-dt$k[c]*x)*(1/(dt$N[c]+1)-bp))
                z <- y + z}
        zz <- (z-1)/(((1-dt$pss[c])^dt$N[c])-1)
        return(zz)}

Rescaled_UpperModel <-function(x){
        z <- 0
        for (n in 0:dt$N[c]) {
                b <- factorial(dt$N[c])/(factorial(n)*factorial(dt$N[c]-n))*(1-dt$pss[c])^(dt$N[c]-n)*(dt$pss[c])^n
                bp <- (Upper/(Upper-n*dt$kp[c]))*b
                y<- (bp*exp(-n*dt$kp[c]*x)+exp(-Upper*x)*(1/(dt$N[c]+1)-bp))
                z <- y + z}
        zz <- (z-1)/(((1-dt$pss[c])^dt$N[c])-1)
        return(zz)}

Rescaled_LowerModel <-function(x){
        z <- 0
        for (n in 0:dt$N[c]) {
                b <- factorial(dt$N[c])/(factorial(n)*factorial(dt$N[c]-n))*(1-dt$pss[c])^(dt$N[c]-n)*(dt$pss[c])^n
                bp <- (Lower/(Lower-n*dt$kp[c]))*b
                y<- (bp*exp(-n*dt$kp[c]*x)+exp(-Lower*x)*(1/(dt$N[c]+1)-bp))
                z <- y + z}
        zz <- (z-1)/(((1-dt$pss[c])^dt$N[c])-1)
        return(zz)}

# This is the Refitting Function - given a k value, plug in all time points and A0 value into the KL function, and return the 1-Refitting_R2 (Residual)
# Note, here we are using the rescaled (FS) values for fitting, so that we can combine 
Refitting_Function <- function(ki){
        current_k <<- ki                # May have to set the current k from within the function as a global parameter for use in the Core function.
        Refitting_Predicted <- sapply(ds$t,Refitting_Function_Core)
        Refitting_R2 <- 1- (sum((ds$FS-Refitting_Predicted)^2))/(sum((ds$FS-mean(ds$FS))^2))
        return(1-Refitting_R2)
}

## This is the core function for reoptimization with KL - given a k and and t, plot out the predicted A0
Refitting_Function_Core <-function(x){
        z <- 0
        N <- dt$N[c]
        kp <- dt$kp[c]
        pss <- dt$pss[c]
        a <- dt$a[c]
        for (n in 0:N) {
      
                b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
                bp <- (current_k/(current_k-n*kp))*b
                #y<- a*(bp*exp(-n*kp*x)+exp(-current_k*x)*(1/(N+1)-bp))                                                                  
                y<- (bp*exp(-n*kp*x)+exp(-current_k*x)*(1/(N+1)-bp))
                z <- y + z}
                final <- (z-1)/((1-pss)^N-1)
                return(final)
                }
        

# This function takes in time information and returns the KL equation y value using the now optimized k
Refitted_Model <-function(x){
        z <- 0
        for (n in 0:dt$N[c]) {
                b <- factorial(dt$N[c])/(factorial(n)*factorial(dt$N[c]-n))*(1-dt$pss[c])^(dt$N[c]-n)*(dt$pss[c])^n
                bp <- (Optimize$par/(Optimize$par-n*dt$kp[c]))*b
                y<- dt$a[c]*(bp*exp(-n*dt$kp[c]*x)+exp(-Optimize$par*x)*(1/(dt$N[c]+1)-bp))                                                                  
                z <- y + z}
        return(z)
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
########################################

############## MAIN #####################
# Looping through the hl.out table, find the currect ID number, then consider the ID number in hl.data.out

# This is the regular graphing block, which uses the ProTurn output k to graph out the kinetic curve.
if (Refit == FALSE){
for (c in 1:nrow(dt)) { 
        print(paste("Now running peptide ", c, " of ", nrow(dt), ". ", round(c/nrow(dt)*100,2), "% done."))
        id <- dt$ID[c] 
        ds <- dp[ which(dp$ID == id),]  # Subsetting that particular ID
 
        dkTable <- sapply(ds$t,Calculate_dk)
        dk <- min(abs(dkTable))
        
        ### Calculate Fractional Synthesis
        FS <- sapply(ds$A0, Calculate_FS)
        ds <- cbind(ds,FS)
        
        if (Is_fly == TRUE & max(ds$t) != plateau_point){
                print("Skipping peptide c because plateau time point is absent")
                next
                }
        
        Upper <- dt$k[c]+dk
        Lower <- dt$k[c]^2/(dt$k[c]+dk)
        Model_Predicted <- sapply(ds$t,Model)
        Model_R2 <- 1- (sum((ds$FS-Model_Predicted)^2))/(sum((ds$FS-mean(ds$FS))^2))
                   
        # This block plots out the regular curve (a, from A.0 to A.inf)
        plot(ds$t, ds$A0, xlab="", ylab="m0", ylim=c(min(ds$A0)-0.1,max(ds$A0)+0.1), 
        xlim=c(-1,max(dp$t)+1), cex=0.8, cex.axis=0.6, cex.lab=0.6, pch=16, ps=28, lwd=2, lty=1, mgp=c(1.55,0.48,0), las=1);
 
        mtext(paste(
        "a:", round(dt$a[c],3), 
        "pss:", round(dt$pss[c],3), 
        "R2:", round(dt$R2[c],3), 
        "SE:", round(dt$SE[c],3)),
        side = 3, line = 0, cex=0.4);

        mtext(paste(
        "k:", round(dt$k[c],3),
        "(", round(Lower,3),"-", round(Upper,3), ")",
        "N:", round(dt$N[c],3), "kp:", round(dt$kp[c], 3)), 
        side = 3, line = 0.54, cex=0.4);

         mtext(paste(
                dt$Uniprot[c],dt$Peptide[c], dt$z[c], "+"),side = 3, line = 1.08, cex=0.4);
         mtext(paste(
                "Time(days)"), side=1, line=1, cex=0.4);

         curve(Model(x), from=0, to=max(dp$t), col="red", add=TRUE)
         curve(UpperModel(x), from=0, to=max(dp$t), col=rgb(100,0,0,50,maxColorValue=255), add=TRUE)
         curve(LowerModel(x), from=0, to=max(dp$t), col=rgb(100,0,0,50,maxColorValue=255), add=TRUE)

        # This block plots out the rescaled curve (Fractional Synthesis from 0 to 1)
#plot(ds$t, ds$FS, xlab="",
#     ylab="m0", ylim=c(min(ds$FS)-0.1,max(ds$FS)+0.1), xlim=c(-1,max(dp$t)+1), cex=0.8, cex.axis=0.6, cex.lab=0.6, pch=16, ps=28, lwd=2, lty=1, mgp=c(1.55,0.48,0), las=1);
#mtext(paste("a:", round(dt$a[c],3), "pss:", round(dt$pss[c],3), "r:", round((dt$R2[c])^0.5,3), "SE:", round((dt$SS[c]/(dt$DP[c]-3))^0.5,3)), side = 3, line = 0, cex=0.4);
#mtext(paste("k:", round(dt$k[c],3), "(", round(Lower,3),"-", round(Upper,3), ")", "N:", round(dt$N[c],3), "kp:", round(dt$kp[c], 3)), side = 3, line = 0.54, cex=0.4);
#mtext(paste(dt$Uniprot[c],dt$Peptide[c], dt$z[c], "+"),side = 3, line = 1.08, cex=0.4);
#mtext(paste("Time(days)"), side=1, line=1, cex=0.4);
#curve(Rescaled_Model(x), from=0, to=max(dp$t), col="red", add=TRUE)
#curve(Rescaled_UpperModel(x), from=0, to=max(dp$t), col=rgb(100,0,0,50,maxColorValue=255), add=TRUE)
#curve(Rescaled_LowerModel(x), from=0, to=max(dp$t), col=rgb(100,0,0,50,maxColorValue=255), add=TRUE)

         oput <- paste(dt$ID[c],dt$Uniprot[c], dt$Peptide[c], dt$z[c], round(dt$R2[c],3), round(dt$SS[c],4), round(dt$SE[c],4), round(log2(dt$k[c]),3), round(log2(dt$k[c]+dk)-log2(dt$k[c]),3), dt$DP[c], round(dt$a[c],3), round(dt$pss[c],4), dt$kp[c], dt$N[c], sep = "\t")

        write(oput, file=output_file, append=T)
}
}

# This is the refitting graphing block
if (Refit == TRUE){
        for (c in 1:nrow(dt)) { 
                print(paste("Now refitting peptide ", c, " of ", nrow(dt), ". ", round(c/nrow(dt)*100,2), "% done."))
                id <- dt$ID[c] 
                ds <- dp[ which(dp$ID == id),]  # Subsetting that particular ID
                
                dkTable <- sapply(ds$t,Calculate_dk)
                dk <- min(abs(dkTable))
                
                ### Calculate Fractional Synthesis
                FS <- sapply(ds$A0, Calculate_FS)
                ds <- cbind(ds,FS)
                
                Upper <- dt$k[c]+dk
                Lower <- dt$k[c]^2/(dt$k[c]+dk)
                Model_Predicted <- sapply(ds$t,Model)
                Model_R2 <- 1- (sum((ds$FS-Model_Predicted)^2))/(sum((ds$FS-mean(ds$FS))^2))
                
              
                # This is the optimization function that tries out different k values and plug into the Refitting function. Optimzation$par will give you the now optimized K. 
                Optimize <- optim(0.3,Refitting_Function) #optim(start value, fxn) # Use optim() for Nelder-Mead
                if(Optimize$par > 10){
                     Optimize$par <- 10}
                 
                # This is the KL function that takes in different values of k and return 1 minus R2 value, which the optim function will try to minimize
                

                
                Refitted_Predicted <- sapply(ds$t,Refitted_Model)
                Refitted_R2 <- 1- (sum((ds$A0-Refitted_Predicted)^2))/(sum((ds$A0-mean(ds$A0))^2))
                
                print(dt$k[c])
                print(Optimize$par)

                # This block plots out the regular curve (a, from A.0 to A.inf) using the reoptimized k
                plot(ds$t, ds$A0, xlab="", ylab="m0", ylim=c(min(ds$A0)-0.1,max(ds$A0)+0.1), 
                     xlim=c(-1,max(dp$t)+1), cex=0.8, cex.axis=0.6, cex.lab=0.6, pch=16, ps=28, lwd=2, lty=1, mgp=c(1.55,0.48,0), las=1);
                
                mtext(paste(
                        "a:", round(dt$a[c],3), 
                        "pss:", round(dt$pss[c],3), 
                        "R2:", round(dt$R2[c],3), 
                        "SE:", round(dt$SE[c],3)),
                      side = 3, line = 0, cex=0.4);
                
                mtext(paste(
                        "k:", round(dt$k[c],3),
                        "(", round(Lower,3),"-", round(Upper,3), ")",
                        "N:", round(dt$N[c],3), "kp:", round(dt$kp[c], 3)), 
                      side = 3, line = 0.54, cex=0.4);
                
                mtext(paste(
                        dt$Uniprot[c],dt$Peptide[c], dt$z[c], "+"),side = 3, line = 1.08, cex=0.4);
                mtext(paste(
                        "Time(days)"), side=1, line=1, cex=0.4);
                
                curve(Model(x), from=0, to=max(dp$t), col="red", add=TRUE)
                curve(UpperModel(x), from=0, to=max(dp$t), col=rgb(100,0,0,50,maxColorValue=255), add=TRUE)
                curve(LowerModel(x), from=0, to=max(dp$t), col=rgb(100,0,0,50,maxColorValue=255), add=TRUE)
                
                # This block plots out the rescaled curve (Fractional Synthesis from 0 to 1)
                #plot(ds$t, ds$FS, xlab="",
                #     ylab="m0", ylim=c(min(ds$FS)-0.1,max(ds$FS)+0.1), xlim=c(-1,max(dp$t)+1), cex=0.8, cex.axis=0.6, cex.lab=0.6, pch=16, ps=28, lwd=2, lty=1, mgp=c(1.55,0.48,0), las=1);
                #mtext(paste("a:", round(dt$a[c],3), "pss:", round(dt$pss[c],3), "r:", round((dt$R2[c])^0.5,3), "SE:", round((dt$SS[c]/(dt$DP[c]-3))^0.5,3)), side = 3, line = 0, cex=0.4);
                #mtext(paste("k:", round(dt$k[c],3), "(", round(Lower,3),"-", round(Upper,3), ")", "N:", round(dt$N[c],3), "kp:", round(dt$kp[c], 3)), side = 3, line = 0.54, cex=0.4);
                #mtext(paste(dt$Uniprot[c],dt$Peptide[c], dt$z[c], "+"),side = 3, line = 1.08, cex=0.4);
                #mtext(paste("Time(days)"), side=1, line=1, cex=0.4);
                #curve(Rescaled_Model(x), from=0, to=max(dp$t), col="red", add=TRUE)
                #curve(Rescaled_UpperModel(x), from=0, to=max(dp$t), col=rgb(100,0,0,50,maxColorValue=255), add=TRUE)
                #curve(Rescaled_LowerModel(x), from=0, to=max(dp$t), col=rgb(100,0,0,50,maxColorValue=255), add=TRUE)
                
                oput <- paste(dt$ID[c],dt$Uniprot[c], dt$Peptide[c], dt$z[c], round(dt$R2[c],3), round(dt$SS[c],4), round(dt$SE[c],4), round(log2(dt$k[c]),3), round(log2(dt$k[c]+dk)-log2(dt$k[c]),3), dt$DP[c], round(dt$a[c],3), round(dt$pss[c],4), dt$kp[c], dt$N[c], sep = "\t")
                
                write(oput, file=output_file, append=T)
        }     
}        
dev.off()

############ Code for protein summary graphs ################################################
library("beeswarm")
op <- read.table(output_file, header=TRUE, fill=TRUE)
attach(op)
temp_order <- aggregate(op$k, list(Uniprot), median, na.rm = TRUE) 
ordered_Uniprot <- factor(Uniprot, levels = temp_order[order(temp_order$x), 1])
protein_count <- unique(unlist(op$Uniprot, use.names=FALSE))
pdf(file="Proturn_summary.pdf", width=length(protein_count)*0.2) # Automatically adjusts width based on number of proteins from output
beeswarm(op$k~ordered_Uniprot, method=c("swarm"),cex=0.8, ylim=c(-3,0.5), las=3)#pwcol=Organ
bxplot(op$k~ordered_Uniprot, data=op, add=1)
dev.off()
################################################################################################

############ OUTPUT SUMMARY OF PROTEIN MEDIANS AND DEVIATIONS, PLOT OUT POSITIONAL K ###############

# Reading the ProTurn Output file from the Grapher
peptide.output <- read.table(output_file, header=TRUE, fill=TRUE)
protein.list <- unique(peptide.output$Uniprot)


# Read the FASTA file
Local_fasta <- Read_local_fasta(local_fasta_location)

# Preparing the Graphical output file
pdf(file="Proturn_Peptide_Positions.pdf")
par(mfrow=c(5,3), mar=c(2.2,2.2,2,2))

# Preparing the Summary file (headers)
summary <- paste("Uniprot", "GeneName","ProteinName","No. Peptides",
                 "k.mean", "k.sd", "k.mean.cv",
                 "k.median","k.mad","k.median.cv",
                 "normality",
                 sep = "\t")
write(summary, file="Proturn_protein_summary.txt", append=T)


#Looping over unique Uniprot ID to get the median of k values from all peptides, etc.
for (i in 1:length(protein.list)) {
        print(paste("Now summarizing Protein # ", i, " of ", length(protein.list), ". ", round(i/length(protein.list)*100,2), "% done."))
        
        # Get the annotations
        GN <- annot[ which(annot$Entry == as.character(protein.list[i])),"Gene.names"]
        PN <- annot[ which(annot$Entry == as.character(protein.list[i])),"Protein.names"]
        
        # Subset the current protein being considered
        subset <- peptide.output[ which(peptide.output$Uniprot == protein.list[i]), ]
        
        # Calculate median and mean of the k values of all peptides for this protein
        protein.median <- median(subset$k)
        protein.mean <- mean(subset$k)
        
        
        # If there are more than one peptide in this protein, calculate median absolute deviation and standard deviation
        if (nrow(subset) > 1) {
                protein.mad <- mad(subset$k)
                protein.sd <- sd(subset$k)
                }
        else {
                protein.mad <- NA 
                protein.sd <- NA
        }
        
 
        
        # If there are three or more peptides in this protein, test whether they are normally distributed
        if (nrow(subset) >= 3 & sd(subset$k) != 0) {
                normality <- shapiro.test(subset$k)
        }
        else { normality <- shapiro.test(c(1,2,1)) }    # Just forcing a p value of 0 (assume non-normal distributed)
        
        # Calculate CV
        protein.median.cv <- protein.mad/protein.median
        protein.mean.cv <- protein.sd/protein.mean

        # Write out the Summary
        summary <- paste(subset$Uniprot[1], GN, PN, 
                         nrow(subset), protein.mean, protein.sd, 
                         protein.mean.cv, protein.median, protein.mad, 
                         protein.median.cv, normality$p.value, 
                         sep = "\t")
        write(summary, file="Proturn_protein_summary.txt", append=T)
        
        ######## Plot the position vs. k graph
        
        # Filling in the protein length and peptide sequence length for each peptide in the protein
        subset$protein_length <- NA
        subset$peptide_pos <- NA
        
        for (j in 1:nrow(subset)) {
                # Getting protein length and sequence length for each peptide
                sequence <- Get_local_sequence(subset$Uniprot[j],subset$Peptide[j])
                subset$protein_length[j] <- unlist(sequence[1])
                subset$peptide_pos[j] <- unlist(sequence[2])
                }
        
        # Preparing the plot
        plot(-10, -10, xlab="Pos", ylab="log2 k", ylim=c(-8,1), xlim=c(0,max(subset$protein_length)), 
             cex=0.8, cex.axis=0.6, cex.lab=0.6, pch=16, ps=28, lwd=2, lty=1, mgp=c(1.55,0.48,0), las=1)
        mtext(paste(subset$Uniprot[1],GN),
              side = 3, line = 0.5, cex=0.4)
        
        for (j in 1:nrow(subset)){
                # Conditional coloring for modified peptides
                
                
                # R's if-else syntax. Matching colors to modifications.
                if(grepl("(42.0106)",subset$Peptide[j]) == TRUE) {
                        color <- "red"
                        } else if(grepl("(114.042927)",subset$Peptide[j]) == TRUE){
                                color <- "purple"
                                } else if (grepl("(79.9663)",subset$Peptide[j]) == TRUE){
                                        color <- "green"
                                        } else if (grepl("\\(",subset$Peptide[j]) == TRUE){
                                                color <- "grey"
                                        } else {
                                                color <- "black"
                                        }
                        
                
                # Plot out the segment ratio
                segments(x0=subset$peptide_pos[j],
                         y0=subset$k[j],
                         x1=subset$peptide_pos[j] + nchar(as.character(gsub( " *\\(.*?\\) *", "", subset$Peptide[j]))),
                         y1=subset$k[j],
                         lwd=5, col=color)
                }
        }
dev.off()



