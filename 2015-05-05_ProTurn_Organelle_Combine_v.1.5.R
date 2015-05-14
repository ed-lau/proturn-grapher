# THIS FILE SIMPLY COMBINES THE RESULTS FROM CYTO, MITO, AND DEBRIS INTO ONE.


############## USER INPUT ##############

home_directory <- "~/Documents/Ping Lab/Project Files/2015 Isoform Turnover"

hl.out_location_1 <- "ctrl/ctrl mouse heart cyto/hl.out"
hl.out_location_2 <- "ctrl/ctrl mouse heart debris/hl.out"
hl.out_location_3 <- "ctrl/ctrl mouse heart mito/hl.out"
hl.data.out_location_1 <- "ctrl/ctrl mouse heart cyto/hl-data.out"
hl.data.out_location_2 <- "ctrl/ctrl mouse heart debris/hl-data.out"
hl.data.out_location_3 <- "ctrl/ctrl mouse heart mito/hl-data.out"

combined_hl.out_1 <- "ctrl/hl.out"                      # Have to make the directory manually
combined_hl.data.out_1 <- "ctrl/hl-data.out"


hl.out_location_4 <- "iso/iso mouse heart cyto/hl.out"
hl.out_location_5 <- "iso/iso mouse heart debris/hl.out"
hl.out_location_6 <- "iso/iso mouse heart mito/hl.out"
hl.data.out_location_4 <- "iso/iso mouse heart cyto/hl-data.out"
hl.data.out_location_5 <- "iso/iso mouse heart debris/hl-data.out"
hl.data.out_location_6 <- "iso/iso mouse heart mito/hl-data.out"

combined_hl.out_2 <- "iso/hl.out"                      # Have to make the directory manually
combined_hl.data.out_2 <- "iso/hl-data.out"

annotation.location <- "~/Documents/Ping Lab/Heavy Water/ProTurn Output/Annotation file/10090_annotations.csv"


########################################
setwd(home_directory)

hl.out.1 <- read.table(hl.out_location_1, header=TRUE, fill=TRUE)
hl.out.1$ID <- paste("c",hl.out.1$ID,sep="")
hl.out.2 <- read.table(hl.out_location_2, header=TRUE, fill=TRUE)
hl.out.2$ID <- paste("d",hl.out.2$ID,sep="")
hl.out.3 <- read.table(hl.out_location_3, header=TRUE, fill=TRUE)
hl.out.3$ID <- paste("m",hl.out.3$ID,sep="")

hl.out <- rbind(hl.out.1,hl.out.2,hl.out.3)
write.table(hl.out, file=combined_hl.out_1, sep= "\t", row.names=FALSE, col.names=TRUE, append=T)
                       
hl.data.out.1 <- read.table(hl.data.out_location_1, header=TRUE, fill=TRUE)
hl.data.out.1$ID <- paste("c",hl.data.out.1$ID,sep="")
hl.data.out.1 <- hl.data.out.1[,1:3]
hl.data.out.2 <- read.table(hl.data.out_location_2, header=TRUE, fill=TRUE)
hl.data.out.2$ID <- paste("d",hl.data.out.2$ID,sep="")
hl.data.out.2 <- hl.data.out.2[,1:3]
hl.data.out.3 <- read.table(hl.data.out_location_3, header=TRUE, fill=TRUE)
hl.data.out.3$ID <- paste("m",hl.data.out.3$ID,sep="")
hl.data.out.3 <- hl.data.out.3[,1:3]

hl.data.out <- rbind(hl.data.out.1,hl.data.out.2,hl.data.out.3)
write.table(hl.data.out[1:3], file=combined_hl.data.out_1, sep= "\t", row.names=FALSE, col.names=TRUE, append=T)


########################################
hl.out.4 <- read.table(hl.out_location_4, header=TRUE, fill=TRUE)
hl.out.4$ID <- paste("c",hl.out.4$ID,sep="")
hl.out.5 <- read.table(hl.out_location_5, header=TRUE, fill=TRUE)
hl.out.5$ID <- paste("d",hl.out.5$ID,sep="")
hl.out.6 <- read.table(hl.out_location_6, header=TRUE, fill=TRUE)
hl.out.6$ID <- paste("m",hl.out.6$ID,sep="")

hl.out <- rbind(hl.out.4,hl.out.5,hl.out.6)
write.table(hl.out, file=combined_hl.out_2, sep= "\t", row.names=FALSE, col.names=TRUE, append=T)

hl.data.out.4 <- read.table(hl.data.out_location_4, header=TRUE, fill=TRUE)
hl.data.out.4$ID <- paste("c",hl.data.out.4$ID,sep="")
hl.data.out.4 <- hl.data.out.4[,1:3]
hl.data.out.5 <- read.table(hl.data.out_location_5, header=TRUE, fill=TRUE)
hl.data.out.5$ID <- paste("d",hl.data.out.5$ID,sep="")
hl.data.out.5 <- hl.data.out.5[,1:3]
hl.data.out.6 <- read.table(hl.data.out_location_6, header=TRUE, fill=TRUE)
hl.data.out.6$ID <- paste("m",hl.data.out.6$ID,sep="")
hl.data.out.6 <- hl.data.out.6[,1:3]

hl.data.out <- rbind(hl.data.out.4,hl.data.out.5,hl.data.out.6)
write.table(hl.data.out[1:3], file=combined_hl.data.out_2, sep= "\t", row.names=FALSE, col.names=TRUE, append=T)