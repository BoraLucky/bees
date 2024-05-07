library(SNPRelate)
args <- commandArgs(TRUE)

vcf.infile <- args[1]
metadata.infile<- args[2]

transformNegativeFSTtozero<- TRUE
keepLowerTriangularMatrixOnly<- TRUE

##################################### Parsing VCF file#####################################################
snpgdsVCF2GDS(vcf.infile, "allpopFst.gds",  method="biallelic.only")
genofile <- snpgdsOpen("allpopFst.gds")
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))


metadata=read.table(metadata.infile,header = T,sep = "\t",stringsAsFactors = F)
metadata=metadata[order(factor(metadata$samples, levels = sample.id)),]
pop_code=metadata$pop
poplevels=unique(pop_code)


res= outer(X= poplevels , Y= poplevels, 
           FUN = function(h,k){
             paste(h,k,sep = "/")
           }
)
colnames(res)=poplevels
rownames(res)=poplevels
as.data.frame(res)


for(i in poplevels) {
  for(j in poplevels) {
    if(grepl("/", res[i,j])) { 
      popelem= unlist(strsplit(res[i,j],"/"))

      flag<- pop_code %in% c(popelem[1],popelem[2])
      samp.sel<- sample.id[flag]
      pop.sel<- pop_code[flag]
      
      if (popelem[1]==popelem[2]){result="0"}else{
        results = snpgdsFst(genofile, sample.id=samp.sel, population=as.factor(pop.sel), 
                            autosome.only=FALSE, method="W&C84",verbose=TRUE)
        result = results$Fst
      }
      res[i,j]=as.character(result)
      res[j,i]=as.character(result) #copy the result
    }
  }
}


res=data.frame(apply(res, 2, function(x) as.numeric(as.character(x))))
rownames(res)=poplevels


if (transformNegativeFSTtozero==TRUE){
  for(i in colnames(res)) {
    for(j in colnames(res)) {
      if (res[i,j]<0){res[i,j]<-0} 
    }
  }
}


if(keepLowerTriangularMatrixOnly == TRUE){
  res[upper.tri(res)]<-0
}


outfile <- paste("pairFstMatrix",".txt", sep="")
write.table(x = res, file = outfile, sep = "\t", dec = ".", 
            quote = F, row.names = T,col.names = NA)

numericMatrix=res

mirror.matrix <- function(x) {
  xx <- as.data.frame(x);
  xx <- rev(xx);
  xx <- as.matrix(xx);
  xx;
}

#----Rotate matrix 270 clockworks----
rotate270.matrix <- function(x) {
  mirror.matrix(t(x))
}

Matrix <- rotate270.matrix(numericMatrix)
ColorRamp <- colorRampPalette(c("white", "steelblue1", "blue3"))
outfileGraphic <- paste("ALLpop_pairFstMatrix",".png", sep="")
outfileGraphic_pdf <- paste("ALLpop_pairFstMatrix",".pdf", sep="")

#save graphic--------------------------------
pdf(outfileGraphic_pdf, width = 10, height = 10)  
smallplot <- c(0.874, 0.9, 0.18, 0.83)
bigplot <- c(0.13, 0.85, 0.14, 0.87)
old.par <- par(no.readonly = TRUE)
# draw legend --------------------------------
par(plt = smallplot)
# get legend values --------------------------------
Min <- min(Matrix, na.rm=TRUE)
Max <- max(Matrix, na.rm=TRUE)
binwidth <- (Max - Min) / 64
y <- seq(Min + binwidth/2, Max - binwidth/2, by = binwidth)
z <- matrix(y, nrow = 1, ncol = length(y))

image(1, y, z, col = ColorRamp(64),xlab="", ylab="", axes=FALSE)

# adjust axis if only one value exists
if(Min == Max){
  axis(side=4, las = 2, cex.axis=0.8, at=Min, labels=round(Min, 2))
} else {
  axis(side=4, las = 2, cex.axis=0.8)
}

box()
mtext(text=expression(bold(F[ST])), side=4, line=2.5, cex=1.1)

#draw main graphic ---------------------------
a <- ncol(numericMatrix)
b <- nrow(numericMatrix)

x <- c(1:a)
y <- c(1:b)

par(new = TRUE, plt = bigplot)

image(x,y,as.matrix(Matrix), col=ColorRamp(64),
      main=expression(bold(Matrix~of~pairwise~F[ST])), xlab="",ylab="", axes=FALSE)
box()

#add labels
Labels=poplevels
if(is.null(Labels)){
  axis(1, at = c(1:a))
  axis(2, at = c(1:b), labels=c(b:1))
  mtext(side = 1, at =(a/2), line = 2.5, text = "Population", cex=1,font=2)
  mtext(side = 2, at =(b/2), line = 2.7, text = "Population", cex=1,font=2)
} else{
  axis(1, at = c(1:a), labels=Labels[1:length(Labels)], cex.axis=1.5,las=2)
  axis(2, at = c(1:b), labels=Labels[length(Labels):1], cex.axis=1.5,las=2)
}

par(old.par)  #reset graphic parameters
dev.off()
############## END FSTPAIRWISE TABLE PLOT ##############
snpgdsClose(genofile)


########ref: https://github.com/rgiannico/RpairwiseFST/################