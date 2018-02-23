### The master file
### Change variables and so on here.
#setwd("/home/mark/diss_may/")

nLGM <- 100
phi <- 0.3
burnIn <- 2000
nIter <- 12000
nChain <- 3

#Sets up the proper node - this code was necessary on the supercomputer stampede
#args <- commandArgs(trailingOnly = TRUE)
#name <- ifelse(length(args)>0, args[1], "nodelist.txt")
#assign("hostnames", scan(name, what="", sep="\n"), envir = .GlobalEnv)

# begin workflow
source("functions.R")
parameters <- read.table("item.iparms.txt")
names(parameters) <- c("Item ID", "a", "b", "c", "content")
source("data_generation.R")
source("inits_mirt.R")
source("chainy_ffbs_gamma_TMat_fixed.R")
source("estimation.R")

# Save the file properly
files=list.files()
num=length(grep("^m_diss",files))
if(num==0) {
  print("not exits")
}else {
  print(num)
}
fileName=paste("m_diss_",nLGM,"_",phi*10,"_",num+1,".Rdata",sep="")
save.image(file=fileName)
