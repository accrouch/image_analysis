#Make Dictionary form AT&T database
#Author: Matt Wheeler (modified code from Matt Shotwell)
#Date: 9/2/17

#Function to read pgm files

pgm <- function(file) {
  if(!is.character(file) || !file.exists(file))
    stop("invalid file argument")
  con <- file(file, "rb")
  # Check for magic number (ascii "P5")
  if(scan(con, character(0), 1, quiet=TRUE)!="P5")
    stop("file not PGM")
  # Read width, height, maxval
  width  <- scan(con, integer(0), 1, quiet=TRUE)
  height <- scan(con, integer(0), 1, quiet=TRUE)
  maxval <- scan(con, integer(0), 1, quiet=TRUE)
  # Read data
  size <- ifelse(maxval > 255, 2, 1)
  data <- readBin(con, "integer", height*width, size, FALSE)
  close(con)
  return(matrix(data, height, width, byrow=TRUE))
}

## Read in the face PGM files
pix <- vector()
grp <- vector()
cnt <- 1
Nsub <- 40
Nimg <- 10
N <- Nsub * Nimg
for( sub in 1:Nsub ) {
  for( img in 1:Nimg ) {
    fname = paste("faces/s",sub,img,".pgm",sep="")
    cat(paste("loading", fname, "\n"))
    dat <- pgm(fname)
    
    tmp <- t(matrix(dat,nrow=112,ncol=92))[,112:1] #orient the images correctly to view
    pix <- c(pix, c(tmp))
    grp <- c(grp, rep(cnt, prod(dim(dat))))
    cnt <- cnt + 1
  }
}
rm(fname, dat, cnt, sub, img)



facepixmat <- unstack(pix, pix ~ grp)


## Read in the face PGM files
pix <- vector()
grp <- vector()
cnt <- 1

N <- Nsub * Nimg

for( img in 1:20 ) {
    fname = paste("cats/c",img,".pgm",sep="")
    cat(paste("loading", fname, "\n"))
    dat <- pgm(fname)
    
    tmp <- t(matrix(dat,nrow=112,ncol=92)[112:1,]) #orient the images correctly to view
    pix <- c(pix, c(tmp))
    grp <- c(grp, rep(cnt, prod(dim(dat))))
    cnt <- cnt + 1
    image( tmp, col=gray(seq(0,1,length.out=1024)),axes=F)
    readline("Press any key to continue:")
}
catpixmat <- unstack(pix, pix ~ grp)
image( matrix(catpixmat[,20],nrow=92,ncol=112), col=gray(seq(0,1,length.out=1024)),axes=F)

face.pca.fnc <- prcomp(t(facepixmat))
cats.pca.fnc <- prcomp(t(catpixmat))

par(mfrow=c(2,1),mai=c(0,0,0,0))
image( matrix(face.pca.fnc$rotation[,2],nrow=92,ncol=112),col=gray(seq(0,1,length.out=1024)),axes=F)
image( matrix(cats.pca.fnc$rotation[,2],nrow=92,ncol=112),col=gray(seq(0,1,length.out=1024)),axes=F)
########################################################################################################
#TRAIN the MODEL
########################################################################################################
library(MASS)
h.matrix.coef <- matrix(0,nrow=200,ncol=13)

for (ii in 1:200){
  h.matrix.coef[ii,] = coef(lm(facepixmat[,ii] ~ face.pca.fnc$rotation[,1:5] + cats.pca.fnc$rotation[,1:7]))
}

c.matrix.coef <- matrix(0,nrow=20,ncol=13)
for (ii in 1:20){
  c.matrix.coef[ii,] = coef(lm(catpixmat[,ii] ~ face.pca.fnc$rotation[,1:5] + cats.pca.fnc$rotation[,1:7]))
}


library(e1071)
classification <- as.factor(c(rep('COMPUTER-SCIENTIST',200),rep('CAT',20)))



a = rbind(h.matrix.coef,c.matrix.coef)



a1 = a[,1]
a2 = a[,2]
a3 = a[,3]
a4 = a[,4]
a5 = a[,5]
a6 = a[,6]
a7 = a[,7]
a8 = a[,8]
a9 = a[,9]
a10 = a[,10]
a11 = a[,11]
a12 = a[,12]
a13 = a[,13]
my.svm <-svm(classification~a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12+a13,probability=TRUE)



#now we test the SVM on new data 
temp <- matrix(0,nrow=5,ncol=13)
####################################

par(mfrow=c(2,2),mar=c(1,1,0,0))

dat <- pgm("cats/kiper.pgm")
kiper <- t(matrix(dat,nrow=112,ncol=92)[112:1,])
image( kiper,col=gray(seq(0,1,length.out=1024)),axes=F)
temp[1,] <- coef(lm(c(kiper)~face.pca.fnc$rotation[,1:5] + cats.pca.fnc$rotation[,1:7]))
dat <- pgm("cats/zmuda.pgm")
zmuda<- t(matrix(dat,nrow=112,ncol=92)[112:1,])
temp[2,] <- coef(lm(c(zmuda)~face.pca.fnc$rotation[,1:5] + cats.pca.fnc$rotation[,1:7]))
image( zmuda,col=gray(seq(0,1,length.out=1024)),axes=F)
dat <- pgm("cats/gcat.pgm")
gcat<- t(matrix(dat,nrow=112,ncol=92)[112:1,])
temp[3,] <- coef(lm(c(gcat)~face.pca.fnc$rotation[,1:5] + cats.pca.fnc$rotation[,1:7]))
image( gcat,col=gray(seq(0,1,length.out=1024)),axes=F)
dat <- pgm("cats/mc2.pgm")
gcat<- t(matrix(dat,nrow=112,ncol=92)[112:1,])
temp[4,] <- coef(lm(c(gcat)~face.pca.fnc$rotation[,1:5] + cats.pca.fnc$rotation[,1:7]))
image( gcat,col=gray(seq(0,1,length.out=1024)),axes=F)
dat <- pgm("cats/steve.pgm")
gcat<- t(matrix(dat,nrow=112,ncol=92)[112:1,])
temp[5,] <- coef(lm(c(gcat)~face.pca.fnc$rotation[,1:5] + cats.pca.fnc$rotation[,1:7]))
image( gcat,col=gray(seq(0,1,length.out=1024)),axes=F)


#now we test the SVM on new data 
newdata <- data.frame(a1=temp[,1],a2=temp[,2],a3=temp[,3],a4=temp[,4],a5=temp[,5],a6=temp[,6],a7=temp[,7],a8=temp[,8],a9=temp[,9],a10=temp[,10],a11=temp[,11],a12=temp[,12],a13=temp[,13])

predict(my.svm,newdata,probability=TRUE,decision.values = TRUE)


par(mfrow=c(2,2),mar=c(1,1,0,0))

dat <- pgm("cats/kiper.pgm")
kiper <- t(matrix(dat,nrow=112,ncol=92)[112:1,])
image( kiper,col=gray(seq(0,1,length.out=1024)),axes=F)
text(0.5,0.5,"Not a Cat",col=2,cex=5)
temp[1,] <- coef(lm(c(kiper)~face.pca.fnc$rotation[,1:5] + cats.pca.fnc$rotation[,1:7]))
dat <- pgm("cats/zmuda.pgm")
zmuda<- t(matrix(dat,nrow=112,ncol=92)[112:1,])
temp[2,] <- coef(lm(c(zmuda)~face.pca.fnc$rotation[,1:5] + cats.pca.fnc$rotation[,1:7]))
image( zmuda,col=gray(seq(0,1,length.out=1024)),axes=F)
text(0.5,0.5,"Not a Cat",col=2,cex=5)
dat <- pgm("cats/gcat.pgm")
gcat<- t(matrix(dat,nrow=112,ncol=92)[112:1,])
temp[3,] <- coef(lm(c(gcat)~face.pca.fnc$rotation[,1:5] + cats.pca.fnc$rotation[,1:7]))
image( gcat,col=gray(seq(0,1,length.out=1024)),axes=F)
text(0.5,0.5,"A cat",col=2,cex=5)
dat <- pgm("cats/mc2.pgm")
gcat<- t(matrix(dat,nrow=112,ncol=92)[112:1,])
temp[4,] <- coef(lm(c(gcat)~face.pca.fnc$rotation[,1:5] + cats.pca.fnc$rotation[,1:7]))
image( gcat,col=gray(seq(0,1,length.out=1024)),axes=F)
text(0.5,0.5,"A cat",col=2,cex=5)
