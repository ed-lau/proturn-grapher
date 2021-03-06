############# INSTRUCTIONS #############
#Note: Read in hl-out and hl-data.out directly. 
#If using old version of ProTurn (before 2.0.0.5) hl.txt and hl-data.txt in Excel, remove the non-peptides in hl.txt, fill all the holes, find gene names, and label the
#following columns: ID, Uniprot, GN, Peptide, DP, Charge, Isotopomer, SS (Sums of squares), a, pss, kp, N, k, and R2; ID, Time and A0.
# kp, pss, N and a. N is from Commerford and must be rounded.
########################################

############## USER INPUT ##############
R2_threshold <- 0.9    # R2 Threshold
SE_threshold <- 0.01     # Standard Error of Estimate Threshold

#home_directory <- "~/Documents/Ping Lab/R Projects/proturn-grapher/test_refit"
home_directory <- "~/Documents/Ping Lab/Project Files/2015 Paraquat Turnover/6. Paraquat Cyto"

hl.out_location <- "hl.out"
hl.data.out_location <- "hl-data.out"
output_file <- "ProTurn_Refitted_Output.txt"

annotation.location <- "~/Documents/Ping Lab/Heavy Water/ProTurn Output/Annotation file/10090_annotations.csv"
########################################

############## LOADING FILES ##############
setwd(home_directory)
annot <- read.csv(file=annotation.location)

oput <- paste("Uniprot","GN","PN","Peptides", "DP", "Refitted_k", "Refitted_dk", "Refitted_R2", sep = "\t")
write(oput, file=output_file, append=T)

pdf(file="Proturn_Combined_fit.pdf")
par(mfrow=c(5,4), mar=c(2.2,2.2,2,2))
hl.out <- read.table(hl.out_location, header=TRUE, fill=TRUE)
dp <- read.table(hl.data.out_location, header=TRUE, fill=TRUE)
dp <- dp[,1:3]

# Calculate the Standard Error of Estimate from SS (Sums of Squares) based on SE = (SS/DP-1)^0.5
SE <- (hl.out$SS/(hl.out$DP-1))^0.5
hl.out <- cbind(hl.out,SE)
rm(SE)
# Subsetting the hl.out table to house only those peptides that pass the threshold
dt <- hl.out[ which(hl.out$R2 >= R2_threshold | hl.out$SE <= SE_threshold),] # Subsetting by R2 and SS
#FS <- 1:nrow(dt)
#dt <- cbind(dt,FS)

########################################

############## FUNCTIONS ##############
# This function calcualtes dk from the solution for partial derivatives dA/dk ; dA is SS.


Refit_Calculate_FS <- function(x){
        
        a <- dt[which(dt$ID == ds$ID[i]),"a"]
        pss <- dt[which(dt$ID == ds$ID[i]),"pss"]
        N <- dt[which(dt$ID == ds$ID[i]),"N"]
        
        Ainf. <- a*(1-pss)^N
        FS. <- (x-a)/(Ainf.-a)
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
        N <- mean(all_peptides_of_protein$N)
        kp <- mean(all_peptides_of_protein$kp)
        pss <- mean(all_peptides_of_protein$pss)

        for (n in 0:N) {
      
                b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
                bp <- (current_k/(current_k-n*kp))*b
                                                             
                y<- (bp*exp(-n*kp*x)+exp(-current_k*x)*(1/(N+1)-bp))
                z <- y + z}
                final <- (z-1)/((1-pss)^N-1)
                return(final)
                }
        


Calculate_Refitted_dk <-function(x){
        z <- 0
        N <- mean(all_peptides_of_protein$N)
        kp <- mean(all_peptides_of_protein$kp)
        pss <- mean(all_peptides_of_protein$pss)
        
        for (n in 0:N) {
                b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
                bp <- (Optimize$minimum/(Optimize$minimum-n*kp))*b
                y<- ((n*kp)/(Optimize$minimum*(Optimize$minimum-n*kp))*bp*(exp(-Optimize$minimum*x)-exp(-n*kp*x))-x*(1/(N+1)-bp)*exp(-Optimize$minimum*x))                                                              
                z <- y + z}
        dA.over.dk <- (z-1)/((1-pss)^N-1)
        dk <- Refitted_SS/dA.over.dk
        return(dk)}

# This function takes in time information and returns the KL equation y value using the now optimized k
Refitted_Model <-function(x){
        z <- 0
        N <- mean(all_peptides_of_protein$N)
        kp <- mean(all_peptides_of_protein$kp)
        pss <- mean(all_peptides_of_protein$pss)
        
        for (n in 0:N) {
                b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
                bp <- (Optimize$minimum/(Optimize$minimum-n*kp))*b
                y<- (bp*exp(-n*kp*x)+exp(-Optimize$minimum*x)*(1/(N+1)-bp))                                                                  
                z <- y + z}
        zz <- (z-1)/((1-pss)^N-1)
        return(zz)
        }

Refitted_UpperModel <-function(x){
        z <- 0
        N <- mean(all_peptides_of_protein$N)
        kp <- mean(all_peptides_of_protein$kp)
        pss <- mean(all_peptides_of_protein$pss)
        
        for (n in 0:N) {
                b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
                bp <- (Refitted_Upper/(Refitted_Upper-n*kp))*b
                y<- (bp*exp(-n*kp*x)+exp(-Refitted_Upper*x)*(1/(N+1)-bp))                                                                  
                z <- y + z}
        zz <- (z-1)/((1-pss)^N-1)
        return(zz)
}

Refitted_LowerModel <-function(x){
        z <- 0
        N <- mean(all_peptides_of_protein$N)
        kp <- mean(all_peptides_of_protein$kp)
        pss <- mean(all_peptides_of_protein$pss)
        
        for (n in 0:N) {
                b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
                bp <- (Refitted_Lower/(Refitted_Lower-n*kp))*b
                y<- (bp*exp(-n*kp*x)+exp(-Refitted_Lower*x)*(1/(N+1)-bp))                                                                  
                z <- y + z}
        zz <- (z-1)/((1-pss)^N-1)
        return(zz)
}

########################################

############## MAIN #####################
# Looping through the hl.out table, find the currect protein, then consider all the ID numbers in hl.data.out
protein_list <- unique(dt$Uniprot)
# This is the combined refitting graphing block
#for (c in 1:35){        
for (c in 1:length(protein_list)) { 
                print(paste("Now refitting protein ", c, " of ", length(protein_list), ". ", round(c/length(protein_list)*100,2), "% done."))
                
                # Get the annotations for the protein being considered
                GN <- annot[ which(annot$Uniprot == as.character(protein_list[c])),5]
                PN <- annot[ which(annot$Uniprot == as.character(protein_list[c])),3]
                
                # Subset all peptides for the protein being currently considered
                all_peptides_of_protein <- dt[ which(dt$Uniprot == protein_list[c]),]  # Subsetting that particular ID
                
                # Subset all data points for all peptides for the protein being currently considered
                tally <- dp$ID  %in% all_peptides_of_protein$ID
                ds <- dp[tally,]
                
                        #dkTable <- sapply(ds$t,Calculate_dk)
                        #dk <- min(abs(dkTable))
                        #Upper <- dt$k[c]+dk
                        #Lower <- dt$k[c]^2/(dt$k[c]+dk)
                        #Model_Predicted <- sapply(ds$t,Model)
                        #Model_R2 <- 1- (sum((ds$FS-Model_Predicted)^2))/(sum((ds$FS-mean(ds$FS))^2))
                
                # Calculate the Fractional Synthesis for each data point
                for (i in 1:nrow(ds)) {
                        ds$FS[i] <- Refit_Calculate_FS(ds$A0[i])
                }
                
                # This is the optimization function that tries out different k values and plug into the Refitting function. Optimzation$par will give you the now optimized K. 
                Optimize <- optimize(Refitting_Function,c(0,5),tol=0.001)
                # This is the KL function that takes in different values of k and return 1 minus R2 value, which the optim function will try to minimize
                
                # Calculate the R2 of the refitting, the sum of squares of the refitting, and the standard deviation
                Refitted_Predicted <- sapply(ds$t,Refitted_Model)
                Refitted_R2 <- 1- (sum((ds$FS-Refitted_Predicted)^2))/(sum((ds$FS-mean(ds$FS))^2))
                Refitted_SS <- sum((ds$FS-Refitted_Predicted)^2)
                Refitted_SD <- (mean((ds$FS-Refitted_Predicted)^2))^0.5
                
                # Calculate dk based on the equation for dA/dk (Using SS as dA)
                dk <- min(abs(sapply(ds$t,Calculate_Refitted_dk)))
                
                # Calculate the Upper bound and Lower bound of the fitted k using dk
                Refitted_Upper <- Optimize$minimum + dk
                Refitted_Lower <- Optimize$minimum^2/(Optimize$minimum+dk)
                
                print(protein_list[c])
                print(Optimize$minimum)

                # This block plots out the regular curve (a, from A.0 to A.inf) using the reoptimized k
                plot(ds$t, ds$FS, xlab="", ylab="FS", ylim=c(-0.1,1.1), 
                     xlim=c(-1,15), cex=0.8, cex.axis=0.6, cex.lab=0.6, pch=16, ps=28, lwd=2, lty=1, mgp=c(1.55,0.48,0), las=1);
                
                
                mtext(paste(
                        "Protein", protein_list[c],
                        "Refitted k:", round(Optimize$minimum,3), 
                        "Refitted R2:", round(Refitted_R2,3)),
                      side = 3, line = 0.54, cex=0.4);
                
                mtext(paste(
                        "Time(days)"), side=1, line=1, cex=0.4);
                
                curve(Refitted_Model(x), from=0, to=max(dp$t), col="red", add=TRUE)
                curve(Refitted_UpperModel(x), from=0, to=max(dp$t), col=rgb(100,0,0,50,maxColorValue=255), add=TRUE)
                curve(Refitted_LowerModel(x), from=0, to=max(dp$t), col=rgb(100,0,0,50,maxColorValue=255), add=TRUE)
                
                oput <- paste(protein_list[c],GN,PN,nrow(all_peptides_of_protein),nrow(ds), Optimize$minimum, dk, Refitted_R2, sep = "\t")
                
                write(oput, file=output_file, append=T)
        }     
       
dev.off()
