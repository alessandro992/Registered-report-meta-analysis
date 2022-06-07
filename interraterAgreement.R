#+ setup, include=FALSE

library(tidyverse)
library(irr)





firstCoder <- read_delim("RR-first-support.csv", delim = ",")
secondCoder <- read_delim("RR-second-support.csv", delim = ",")

# Check of the first 10% (remove it from the comment to calculate the interrater agreement of the first 10%)
# ratersData <- list(firstCoder[1:11,], secondCoder[1:11,])
# 20% overall
ratersData <- list(firstCoder[1:3,], secondCoder[1:3,])

variableNames <- names(secondCoder)[c(1,4,10:19,21:26,28:30,32:54,61)]


# kappas <- list()
 #for(i in variableNames){
  # kappas[[i]] <- kappa2(cbind(select(.data = ratersData[[1]], contains(variableNames[i])), select(.data = ratersData[[2]], contains(variableNames[i]))))
 #}
# names(kappas) <- variableNames

agreement <- list()
for(i in variableNames){
  agreement[[i]] <- agree(cbind(ratersData[[1]][i], ratersData[[2]][i]))
}
names(agreement) <- variableNames


#+ include=TRUE
#'# Inter-rater agreement being in nature
agreement

rm(list = ls())
#+ include=FALSE
firstCoder <- read_delim("RR-first-nature.csv", delim = ",")
secondCoder <- read_delim("RR-second-nature.csv", delim = ",")

# Check of the first 10% (remove it from the comment to calculate the interrater agreement of the first 10%)
# ratersData <- list(firstCoder[1:11,], secondCoder[1:11,])
# 20% overall
ratersData <- list(firstCoder[1:3,], secondCoder[1:3,])

variableNames <- names(secondCoder)[c(1,4,10:19,21:26,28:30,32:54,61)]


# kappas <- list()
#for(i in variableNames){
# kappas[[i]] <- kappa2(cbind(select(.data = ratersData[[1]], contains(variableNames[i])), select(.data = ratersData[[2]], contains(variableNames[i]))))
#}
#names(kappas) <- variableNames

agreement <- list()
for(i in variableNames){
  agreement[[i]] <- agree(cbind(ratersData[[1]][i], ratersData[[2]][i]))
}
names(agreement) <- variableNames


#+ include=TRUE
#'# Inter-rater agreement emotional social support
agreement