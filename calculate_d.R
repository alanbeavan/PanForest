library(caper)
library(phytools)
library(getopt)
library(future)
library(flock)
library(stringr)

#Get call input
spec <- matrix(c('path', 'a', 1, "character",
                 'phylogeny', 't', 1, "character",
                 'gene_pa', 'g', 1, "character",
                 'cores', 'c', 1, "integer",
                 'output', 'o', 1, "character"), byrow=TRUE, ncol=4)
opt <- getopt(spec)
#hardcoded
#opt = list("~/OneDrive - The University of Nottingham/Post_doc_notts/Software/pangenome_rf/test_results/", "strain_phylogeny_mod.phy", "collapsed_matrix.csv", 1, "d_test")
#names(opt) = c("path", "phylogeny", "gene_pa", "cores", "output")
#setwd(opt$path)

#Read in
outstr <- paste(opt$output, "_nodes_in.csv", sep="")
genes  <- read.csv(outstr, check.names=TRUE) #"coincident_nodes_in.csv")

#Read in tree
tree <- read.tree(opt$phylogeny)
tree$tip.label <- make.names(tree$tip.label) #ensure tree tip names will match annot rownames
#Ensure no zero branch lengths
if (!is.na(match(0, tree$edge.length))) {
  print("Phylogeny contains pairs of tips on zero branch lengths, cannot currently simulate")
  quit()
}
#Root phylogeny for input into comparative.data
treeRt <- midpoint.root(tree)
treeRt <- makeLabel(treeRt)

#Gene_pa read in
#genepa <- read.csv(opt$gene_pa, header=T, row.names=1)
#load in line by line; if rowname of line %in% names(genes), put into annot
con = file(opt$gene_pa, "r") #GOOD
header=readLines(con, n=1) #read in line-GOOD
header=gsub("[.]", "_", header)
header.sp = strsplit(header,",") #split line into list-GOOD

header.sp[[1]] <- make.names(header.sp[[1]]) #puts "X" into the first 3 columns - not sure this is necessary
annot <- matrix(ncol=length(header.sp[[1]])) #creates an empty list of length = nrow matrix
colnames(annot) <- header.sp[[1]] #annotates them
flag <- 1
print("Read in gene_pa file..")
while(TRUE) {
  line=readLines(con, n=1) #read in line
  if (length(line)==0) {
    break #file done
  }
  line.sp = strsplit(line,",") #split line into list
  if (line.sp[[1]][1] %in% names(genes)) {
    print(line.sp[[1]][1]) #its a gene of interest, keep around
    if (flag == 1) {
      annot[1,] <- as.character(line.sp[[1]])
      flag <- 0
    } else {
      annot <- rbind(annot, as.character(line.sp[[1]]))#, stringsAsFactors=FALSE)
    }
  }
}
close(con)
annot <- as.data.frame(annot)
rownames(annot) <- annot[,1]
annot[,1:3] <- NULL
#Make annot table
#genepa[,1:14] <- NULL
#rownames(genepa) <- gsub(" ", "", rownames(genepa))
#annot <- genepa[(rownames(genepa) %in% names(genes)),]
#genepa <- NULL
annot <- t(annot)
annot <- as.data.frame(annot)
annot$Id <- rownames(annot)
#Double check that each column has 2 states
remove <- vector()
for(a in 1:length(colnames(annot))) {
  if (length(table(annot[a])) == 1) {
    print("Encountered the following gene that only exists in one state in the data.")
    print(colnames(annot)[a])
    print("I'm removing it, but you should be cautious as to what this means...")
    remove <- append(remove, a)
  }
}
#Remove any columns with only one state
if (length(remove)>0) {
  for(a in 1:length(remove)) {
    annot[remove[a]] <- NULL
  }
}


####IAMHERE!!!!!!!!####


#Make comparative data object
dataset <- comparative.data(phy = treeRt,
                            data = annot,
                            names.col = Id,
                            vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
#Estimate D (Fritz & Purvis 2010) for all columns in annot
#Detect available cores and decrease the number of parallel runs to that if its < opt$cores
availcores <- availableCores()
cores <- 1
if (availcores < opt$cores) {
  cores <- as.integer(availcores)
} else {
  cores <- opt$cores
}
print("Cores is set to:")
print(cores)
outstr <- paste(opt$output, "_nodes.tsv", sep="")
write(paste("ID","Result",sep="\t"),file=outstr,append=FALSE)
parallelCluster <- parallel::makeCluster(cores, type="FORK")
mkWorker <- function(dataset) {
  library(flock)
  #Make sure each value is passed
  #force(dataset)
  #Define function
  calcD <- function(dataset, binvarry) {
    if(binvarry != "Id") {
      result <- eval(parse(text=paste("caper::phylo.d(data=dataset, binvar=",binvarry,", permut=1000)", sep="")))
      line <- paste(binvarry,result$DEstimate, sep="\t")
      locked_towrite <- flock::lock("./.lock")
      write(line,file=outstr,append=TRUE)
      flock::unlock(locked_towrite)
      rm(result)
    }
  }
  #Define & return worker function
  worker <- function(binvarry) {
    calcD(dataset, binvarry)
  }
  return(worker)
}
#results <- parallel::parLapply(parallelCluster, colnames(annot), mkWorker(dataset))
parallel::parLapply(parallelCluster, colnames(annot), mkWorker(dataset))
# Shutdown cluster neatly
if(!is.null(parallelCluster)) {
  parallel::stopCluster(parallelCluster)
  parallelCluster <- c()
}