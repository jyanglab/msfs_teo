# Jinliang Yang 11.5.2010
# last update:
# Purpose: GWAS version3$ fine genetic map

p2g <- function (predictdf=test, train="gam_training_nam.map.txt") {
  #===========================================================================================#
  # Run GAM function to convert Physical position to Genetic position
  # Estimate Gentetic positions from Physical positions by using Integrated map (Liu et al. PLoS Genetics) 
  # Usage: p2g(input_file_name)
  # Input file:  Three columns: 1.marker_name 2.chromosome (Integer) 3. Physical position (Colname must be Physical);
  # Author: Sanzhen Liu
  # Update date: 8/25/09
  #===========================================================================================#
  library(mgcv)
  # read marker file:
  marker <- predictdf
  # read P-G data:
  d <- train
  #d <- read.delim(train, header=T) # subject to change
  # Plot and output the predicted data:
  marker_Pos<-NULL
  for (chr in sort(unique(marker[,2]))) {
    # extract marker data for certain chromosome:
    d.curr<-d[d[,2]==chr,]
    # GAM fitting
    out = gam(Genetic~s(Physical,k=50),data=d.curr)
    # Prediction:
    markername <- data.frame(marker[marker[,2]==chr,1])
    phy <- data.frame(marker[marker[,2]==chr,3])
    colnames(phy)="Physical"
    gen <- data.frame(predict.gam(out,phy))
    gen2 <- round(gen, 5)
    Pchr <- cbind(markername,chr,phy,gen2)
    marker_Pos <- rbind(marker_Pos,Pchr)
  }
  # change column name:
  colnames(marker_Pos)=c("marker","chr","physical","genetic")
  # output
  return(marker_Pos)
}


