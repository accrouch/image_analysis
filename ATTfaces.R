#Make Dictionary form AT&T database
#Author: Matt Wheeler (modified code from Matt Shotwell)
#Date: 9/2/17
library(glmnet)
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

dat<- pgm("faces/bailer_spray.pgm")
bailer3 <- matrix(dat + 0*rnorm(92*112,0,10),92,112)[92:1,112:1]
image(bailer3,col=gray(seq(0,1,length.out=1024)),axes=F)

## Read in the PGM files
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


par(mfrow=c(2,2),mar=c(1,1,0,0))

pixmat <- unstack(pix, pix ~ grp)
w <- 112; h <- 92  #images are 92x112 pixles 
#look at one of the images to see that we are ok
image( matrix(pixmat[,83],nrow=92,ncol=112), col=gray(seq(0,1,length.out=1024)),axes=F)
image( matrix(pixmat[,183],nrow=92,ncol=112), col=gray(seq(0,1,length.out=1024)),axes=F)
image( matrix(pixmat[,50],nrow=92,ncol=112), col=gray(seq(0,1,length.out=1024)),axes=F)
image( matrix(pixmat[,330],nrow=92,ncol=112), col=gray(seq(0,1,length.out=1024)),axes=F)

#################################################################
#Break up the images into patches, which will then be used to form 
#the basis for our dictionary, which is created using a principal components
#analysis. 
#################################################################
pix <- vector()
grp <- vector()
cnt <- 1

h.slice <- c(0,seq(28,112,28)) #cut up the image horizontally
v.slice <- c(0,seq(23,92,23))  #cut the image vertically
for (ii in 1:400){ #For each image
  # break it up into squares
    for (jj in 1:(length(h.slice)-1)){
      for (kk in 1:(length(v.slice)-1)){
        h.start <- h.slice[jj]+1
        h.end   <- h.slice[jj+1]
        
        v.start <- v.slice[kk] +1
        v.end   <- v.slice[kk+1]
        
        tmp <- matrix(pixmat[,ii],nrow=92,ncol=112)
        tmp <- tmp[v.start:v.end,h.start:h.end]
        pix <- c(pix,c(tmp))
        grp <- c(grp, rep(cnt, prod(dim(tmp))))
        cnt <- cnt + 1
        # For testing making sure I got what I asked for
        # which are small 23x28 slices of the origional picture
        #print( dim(tmp))
        #image(tmp,col=gray(seq(0,1,length.out=1024)),axes=F)
        #invisible(readline(prompt="Press [enter] to continue"))
      }
    }
}

smallpixmat <- unstack(pix, pix ~ grp)

############################################################
#make the principal components and the functions
pca.fnc <- prcomp(t(smallpixmat))
pca.fnc.big <- prcomp(t(pixmat))


par(mfrow=c(4,5), mai = c(0,0,0,0))

# plot the first 20 principal components
for(ii in 1:20){
  image(matrix(pca.fnc$rotation[,ii],23,28),col=gray(seq(0,1,length.out=1024)),axes=F)
}
# plot the next 20 principal components
par(mfrow=c(4,5), mai = c(0,0,0,0))
for(ii in 21:40){
  image(matrix(pca.fnc$rotation[,ii],23,28),col=gray(seq(0,1,length.out=1024)),axes=F)
}




############################################################################################
# We now have all of our "Dictionary Elements"
# Woo hooo!!!
############################################################################################
#break up the bailer image
#
#
############################################################################################
# break it up into squares
#################################################################
pix <- vector()
grp <- vector()
cnt <-1
for (jj in 1:(length(h.slice)-1)){
  for (kk in 1:(length(v.slice)-1)){
    h.start <- h.slice[jj]+1
    h.end   <- h.slice[jj+1]
    
    v.start <- v.slice[kk] +1
    v.end   <- v.slice[kk+1]
    
    tmp <- matrix(bailer3,nrow=92,ncol=112)
    tmp <- tmp[v.start:v.end,h.start:h.end]
    pix <- c(pix,c(tmp))
    grp <- c(grp, rep(cnt, prod(dim(tmp))))
    cnt <- cnt+1
    # For testing making sure I got what I asked for
    # which are small 23x28 slices of the origional picture
    #print( dim(tmp))
    #image(tmp,col=gray(seq(0,1,length.out=1024)),axes=F)
    #invisible(readline(prompt="Press [enter] to continue"))
  }
}



bailer.chunks <- unstack(pix, pix ~ grp)
############################################
#attempt to predict new 23x28 chunks of data from the noisy image
############################################
cnt <- 1
h.slice <- c(0,seq(28,112,28)) #cut up the image horizontally
v.slice <- c(0,seq(23,92,23))  #cut the image vertically
t101 <- pgm("faces/s8673409.pgm")
bailer.new1 <- t(matrix(t101,112,91))[,112:1]
bailer.new <- matrix(0,92,112)
cnt <-1
for (jj in 1:(length(h.slice)-1)){
    for (kk in 1:(length(v.slice)-1)){
      h.start <- h.slice[jj]+1
      h.end   <- h.slice[jj+1]
      
      v.start <- v.slice[kk] +1
      v.end   <- v.slice[kk+1]

      new.est <- glmnet(cbind(diag(644),pca.fnc$rotation[,1:322]),bailer.chunks[,cnt],penalty.factor=c(rep(1.15,644),rep(0.8,50),rep(1.1,50),rep(1.15,222)))
    
      tmp <- predict(new.est,cbind(diag(644)*0,pca.fnc$rotation[,1:322]))
      bailer.new[v.start:v.end,h.start:h.end] <- tmp[,floor(ncol(tmp))]
      cnt <- cnt + 1
      # For testing making sure I got what I asked for
      # which are small 23x28 slices of the origional picture
      #print( dim(tmp))
      #image(tmp,col=gray(seq(0,1,length.out=1024)),axes=F)
      #invisible(readline(prompt="Press [enter] to continue"))
    }
}

##image(bailer.new1,axes=F,col=gray(seq(0,1,length.out=1024)))
par(mfrow=c(2,1),mai=c(0,0,0,0))


image(bailer.new,axes=F,col=gray(seq(0,1,length.out=1024)))
image(bailer3,axes=F,col=gray(seq(0,1,length.out=1024)))
cnt <- cnt + 1
#check your residuals !!!!
image(bailer3-bailer.new,axes=F,col=gray(seq(0,1,length.out=1024)))

quantile(bailer3-bailer.new,.01)
#Use the residual error to develop a new template
#This would normally done automatically using advanced methods that
#We are not going to discuss here. 
max(bailer3-bailer.new)
temp <- bailer3-bailer.new
temp[temp >-60] = 0
temp[temp != 0] = 255
image(temp,axes=F,col=gray(seq(0,1,length.out=1024))) #this is what we want to extract!!


### Final Analysis
cnt<-1
bailer.final <- matrix(0,nrow = 92,112)
for (jj in 1:(length(h.slice)-1)){
  for (kk in 1:(length(v.slice)-1)){
    h.start <- h.slice[jj]+1
    h.end   <- h.slice[jj+1]
    
    v.start <- v.slice[kk] +1
    v.end   <- v.slice[kk+1]
    
    gbob <-  temp[v.start:v.end,h.start:h.end]
    
   
    new.est <- glmnet(cbind(c(gbob),pca.fnc$rotation),bailer.chunks[,cnt],penalty.factor=c(0,rep(1,644)))
    tmp <- predict(new.est,cbind(c(gbob)*0,pca.fnc$rotation))
    bailer.final[v.start:v.end,h.start:h.end] <- tmp[,floor(ncol(tmp))]
    cnt <- cnt + 1
    # For testing making sure I got what I asked for
    # which are small 23x28 slices of the origional picture
    #print( dim(tmp))
    #image(tmp,col=gray(seq(0,1,length.out=1024)),axes=F)
    #invisible(readline(prompt="Press [enter] to continue"))
  }
}
par(mfrow=c(2,1),mai=c(0,0,0,0))
image(bailer.final,axes=F,col=gray(seq(0,1,length.out=1024)))
image(bailer3,axes=F,col=gray(seq(0,1,length.out=1024)))


