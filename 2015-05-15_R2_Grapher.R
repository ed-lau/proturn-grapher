
############## USER INPUT ##############

# Plotting steps
steps <- 100

# Plotting range
R2_start <- 1.0
R2_end <- 0.5
SE_start <- 0.0
SE_end <- 0.25
 
#home_directory <- "~/Documents/Ping Lab/R Projects/proturn-grapher/mitotempo cyto"
home_directory <- "~/Documents/Ping Lab/Project Files/2015 Isoform Turnover/Data/iso/iso mouse heart mito"
#home_directory <- "~/Documents/Ping lab/Project Files/2015 Paraquat Turnover/6. Paraquat Cyto New"

output_file <- "ProTurn_Output.txt"

########################################

############## LOAD PACKAGES ###########

require("dplyr")
require("lattice")
########################################

peptide.output <- read.table(output_file, header=TRUE, fill=TRUE)

protein_matrix <- matrix(,nrow=steps,ncol=steps)
peptide_matrix <- matrix(,nrow=steps,ncol=steps)
cv_matrix <- matrix(,nrow=steps,ncol=steps)

for (i in 0:steps) {
        print(paste("Now running step ", i, ", ", i/steps*100, "% done."))
        
        for (j in 0:steps){
        R2_threshold <- R2_start + (i/steps)*(R2_end-R2_start)
        SE_threshold <- SE_start + (j/steps)*(SE_end-SE_start)

        # Subsetting the output table to house only those peptides that pass the threshold
        peptide.subset <- peptide.output [ which(peptide.output$R2 >= R2_threshold | peptide.output $SE <= SE_threshold),] # Subsetting by R2 and SS
        
        # Summarize mean CV each uniprot.
        summary <- summarize(group_by(peptide.subset,Uniprot), mean.k = mean (k), sd.k = sd (k)) 
        summary <- mutate(summary, cv.k = sd.k/mean.k)
        mean.cv <- mean(summary$cv.k,na.rm=TRUE)
        # How many proteinsa are left with these R2 and SE criteria ...
        protein_matrix[i,j] <- length(unique(peptide.subset$Uniprot))
        peptide_matrix[i,j] <- nrow(peptide.subset)
        cv_matrix[i,j] <- mean.cv
        }
}


# begin generating my 3D shape
R2 <- seq(R2_start,R2_end,length=steps)
SE <- seq(SE_start,SE_end,length=steps)
plotting_matrix <- expand.grid(R2=R2,SE=SE)
plotting_matrix$protein_count <- c(protein_matrix)  # Flatten out the matrix
plotting_matrix$protein_count[plotting_matrix$protein_count < -1] <- -1

plotting_matrix$peptide_count <- c(peptide_matrix)  # Flatten out the matrix
plotting_matrix$peptide_count[plotting_matrix$peptide_count < -1] <- -1

plotting_matrix$cv <- c(cv_matrix)  # Flatten out the matrix
plotting_matrix$cv[plotting_matrix$cv < -1] <- -1

# Combining both 
plotting_matrix$combined <- (1-(plotting_matrix$cv))^3*(plotting_matrix$peptide_count)
# end generating 3D shape

# Plot out 3d shape with lattice
wireframe(combined ~ R2 * SE, plotting_matrix, shade = TRUE, aspect = c(1, 0.5),
          light.source = c(10,10,10), main = "Protein Count",
          scales = list(z.ticks=5,arrows=FALSE, col="black", font=10, tck=0.5),
          screen = list(z = 45, x = -60, y = 0))

plotting_matrix[which(plotting_matrix$combined == max(plotting_matrix$combined)),]
