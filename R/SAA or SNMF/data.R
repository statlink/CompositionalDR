library(Rfast)
library(robCompositions)
library(ToolsForCoDa)
library(xtable)
library(coda.base)
library(CompositionalDR)
source("saakl.eg2.R")


R <- 20
nam <- c("FQP", "FEG", "KLD", "KLD2") 
ti <- matrix(nrow = R, ncol = 4)
colnames(ti) <- nam
obj <- times <- list()
for ( i in 1:12 )  obj[[ i ]] <- times[[ i ]] <- ti
names(obj) <- names(times) <- c( "economics", "pottery", "chorizonDL", 
"gemas", "gjovik", "honey", "nutrients", "payments", "bmi_activity", 
"waste", "bronze", "FADN" )


### economics
economics <- read.table("economics.txt", header = TRUE)
x <- as.matrix(economics[,2:6])
x <- x/rowsums(x)
group <- economics[,7]
mod <- saa.qp(x, 3)
ternary(mod$W, colour = group)

for ( i in 1:R ) {
  mod1 <- saa.qp(x, 3, maxiter = 50000)
  mod2 <- saa.eg(x, 3, maxiter = 50000)
  mod3 <- saakl.eg(x, 3, maxiter = 50000)
  mod3a <- saakl.eg2(x, 3, maxiter = 1000)
  times[[ 1 ]][i, ] <- c( mod1$runtime[3], mod2$runtime[3], mod3$runtime[3], mod3a$runtime[3])
  obj[[ 1 ]][i, ] <- c( mod1$obj, mod2$obj, mod3$obj, mod3a$obj) 
  times[[ 1 ]][i, ] <-  c( mod1$runtime[3], mod2$runtime[3], mod3$runtime[3], mod3a$runtime[3]) 
}
save.image("experiments.RData")

## pottery
pottery <- read.table("pottery.txt", header = TRUE)
x <- as.matrix(pottery[, -1])
x <- x/rowsums(x)
group <- as.numeric( as.factor(pottery[, 1]) )
mod <- saa.qp(x,3)
ternary(mod$W, colour = group)

for ( i in 1:R ) {
  mod1 <- saa.qp(x, 3, maxiter = 50000)
  mod2 <- saa.eg(x, 3, maxiter = 50000)
  mod3 <- saakl.eg(x, 3, maxiter = 50000)
  mod3a <- saakl.eg2(x, 3, maxiter = 1000)
  times[[ 2 ]][i, ] <- c( mod1$runtime[3], mod2$runtime[3], mod3$runtime[3], mod3a$runtime[3])
  obj[[ 2 ]][i, ] <- c( mod1$obj, mod2$obj, mod3$obj, mod3a$obj) 
  times[[ 2 ]][i, ] <-  c( mod1$runtime[3], mod2$runtime[3], mod3$runtime[3], mod3a$runtime[3]) 
}
save.image("experiments.RData")

### Package robCompositions
## chorizonDL
x <- chorizonDL 
x <- as.matrix(x[, -c(1:3)])
x <- x[, 1:53]
x <- x[, -c(3, 10, 20, 23, 29, 31, 33, 37, 43, 48, 53)]
x <- x/rowsums(x)
mod <- saa.qp(x, 3)
ternary(mod$W)

for ( i in 1:R ) {
  mod1 <- saa.qp(x, 3, maxiter = 50000)
  mod2 <- saa.eg(x, 3, maxiter = 50000)
  mod3 <- saakl.eg(x, 3, maxiter = 50000, clip_exp = 20)
  #mod3a <- saakl.eg2(x, 3, maxiter = 50000, clip_exp = 20)
  times[[ 3 ]][i, ] <- c( mod1$runtime[3], mod2$runtime[3], mod3$runtime[3], mod3a$runtime[3])
  obj[[ 3 ]][i, ] <- c( mod1$obj, mod2$obj, mod3$obj, mod3a$obj) 
  times[[ 3 ]][i, ] <-  c( mod1$runtime[3], mod2$runtime[3], mod3$runtime[3], mod3a$runtime[3]) 
  save.image("experiments.RData")
}
save.image("experiments.RData")

### gemas
x <- gemas[, -c(1:8)]
x <- as.matrix(x)
x <- x/ rowsums(x)
id <- which( is.na( rowsums(x) ) )
x <- x[-id, ]
group <- as.numeric( as.factor(gemas[-id, 8]) )
mod <- saa.qp(x, 3)
ternary(mod$W, colour = group)

for ( i in 1:R ) {
  mod1 <- saa.qp(x, 3, maxiter = 50000)
  mod2 <- saa.eg(x, 3, maxiter = 50000)
  mod3 <- saakl.eg(x, 3, maxiter = 50000)
  mod3a <- saakl.eg2(x, 3, maxiter = 1000)
  times[[ 4 ]][i, ] <- c( mod1$runtime[3], mod2$runtime[3], mod3$runtime[3], mod3a$runtime[3])
  obj[[ 4 ]][i, ] <- c( mod1$obj, mod2$obj, mod3$obj, mod3a$obj) 
  times[[ 4 ]][i, ] <-  c( mod1$runtime[3], mod2$runtime[3], mod3$runtime[3], mod3a$runtime[3]) 
}
save.image("experiments.RData")

## gjovik
x <- gjovik
x <- as.matrix(x[, -c(1:10)])
x <- x/ rowsums(x)
id <- which( is.na( rowsums(x) ) )
x <- x[-id, ]
group <- as.numeric( as.factor(gjovik[-id, 10]) )
mod <- saa.eg(x, 3)
ternary(mod$W, colour = group)

for ( i in 1:R ) {
  mod1 <- saa.qp(x, 3, maxiter = 50000)
  mod2 <- saa.eg(x, 3, maxiter = 50000)
  mod3 <- saakl.eg(x, 3, maxiter = 50000)
  mod3a <- saakl.eg2(x, 3, maxiter = 1000)
  times[[ 5 ]][i, ] <- c( mod1$runtime[3], mod2$runtime[3], mod3$runtime[3], mod3a$runtime[3])
  obj[[ 5 ]][i, ] <- c( mod1$obj, mod2$obj, mod3$obj, mod3a$obj) 
  times[[ 5 ]][i, ] <-  c( mod1$runtime[3], mod2$runtime[3], mod3$runtime[3], mod3a$runtime[3]) 
}
save.image("experiments.RData")

## honey
x <- honey
group <- as.numeric(x[, 1])
x <- as.matrix( x[, -c(1:5)] )
x[which( is.na(x)) ] <- 0
x <- x / rowsums(x)
mod <- saa.eg(x, 3)
ternary(mod$W, colour = group)

for ( i in 1:R ) {
  mod1 <- saa.qp(x, 3, maxiter = 50000)
  mod2 <- saa.eg(x, 3, maxiter = 50000)
  mod3 <- saakl.eg(x, 3, maxiter = 50000)
  mod3a <- saakl.eg2(x, 3, maxiter = 1000)
  times[[ 6 ]][i, ] <- c( mod1$runtime[3], mod2$runtime[3], mod3$runtime[3], mod3a$runtime[3])
  obj[[ 6 ]][i, ] <- c( mod1$obj, mod2$obj, mod3$obj, mod3a$obj) 
  times[[ 6 ]][i, ] <-  c( mod1$runtime[3], mod2$runtime[3], mod3$runtime[3], mod3a$runtime[3]) 
}
save.image("experiments.RData")

## nutrients
x <- as.data.frame(nutrients)
x <- x[, -c(1:12)]
x <- as.matrix( x[, c(1:30, 36)] )
x <- na.omit(x)
x <- x / rowsums(x)
mod <- saa.eg(x, 3)
ternary(mod$W)

for ( i in 1:R ) {
  mod1 <- saa.qp(x, 3, maxiter = 50000)
  mod2 <- saa.eg(x, 3, maxiter = 50000)
  mod3 <- saakl.eg(x, 3, maxiter = 50000)
  mod3a <- saakl.eg2(x, 3, maxiter = 1000)
  times[[ 7 ]][i, ] <- c( mod1$runtime[3], mod2$runtime[3], mod3$runtime[3], mod3a$runtime[3])
  obj[[ 7 ]][i, ] <- c( mod1$obj, mod2$obj, mod3$obj, mod3a$obj) 
  times[[ 7 ]][i, ] <-  c( mod1$runtime[3], mod2$runtime[3], mod3$runtime[3], mod3a$runtime[3]) 
}
save.image("experiments.RData")

## payments
x <- as.matrix( payments[, -c(1:6)] )
x <- x / rowsums(x)
group <- as.numeric(payments[, 2]) - 1
mod <- saa.qp(x, 3)
ternary(mod$W, colour = group)

for ( i in 1:R ) {
  mod1 <- saa.qp(x, 3, maxiter = 50000)
  mod2 <- saa.eg(x, 3, maxiter = 50000)
  mod3 <- saakl.eg(x, 3, maxiter = 50000)
  mod3a <- saakl.eg2(x, 3, maxiter = 1000)
  times[[ 8 ]][i, ] <- c( mod1$runtime[3], mod2$runtime[3], mod3$runtime[3], mod3a$runtime[3])
  obj[[ 8 ]][i, ] <- c( mod1$obj, mod2$obj, mod3$obj, mod3a$obj) 
  times[[ 8 ]][i, ] <-  c( mod1$runtime[3], mod2$runtime[3], mod3$runtime[3], mod3a$runtime[3]) 
}
save.image("experiments.RData")

## coda.base
x <- bmi_activity
x <- as.matrix(x[, 3:7])
x <- x/rowsums(x)
group <- as.numeric( as.factor(bmi_activity[, 8]) )
mod <- saa.qp(x, 3)
ternary(mod$W, colour = group)

for ( i in 1:R ) {
  mod1 <- saa.qp(x, 3, maxiter = 50000)
  mod2 <- saa.eg(x, 3, maxiter = 50000)
  mod3 <- saakl.eg(x, 3, maxiter = 50000)
  mod3a <- saakl.eg2(x, 3, maxiter = 1000)
  times[[ 9 ]][i, ] <- c( mod1$runtime[3], mod2$runtime[3], mod3$runtime[3], mod3a$runtime[3])
  obj[[ 9 ]][i, ] <- c( mod1$obj, mod2$obj, mod3$obj, mod3a$obj) 
  times[[ 9 ]][i, ] <-  c( mod1$runtime[3], mod2$runtime[3], mod3$runtime[3], mod3a$runtime[3]) 
}
save.image("experiments.RData")

x <- waste
x <- as.matrix(x[, 5:9])
x <- x/rowsums(x)
group <- as.numeric( as.factor(waste[, 4]) )
mod <- saa.qp(x, 3)
ternary(mod$W, colour = group)

for ( i in 1:R ) {
  mod1 <- saa.qp(x, 3, maxiter = 50000)
  mod2 <- saa.eg(x, 3, maxiter = 50000)
  mod3 <- saakl.eg(x, 3, maxiter = 50000)
  mod3a <- saakl.eg2(x, 3, maxiter = 1000)
  times[[ 10 ]][i, ] <- c( mod1$runtime[3], mod2$runtime[3], mod3$runtime[3], mod3a$runtime[3])
  obj[[ 10 ]][i, ] <- c( mod1$obj, mod2$obj, mod3$obj, mod3a$obj) 
  times[[ 10 ]][i, ] <-  c( mod1$runtime[3], mod2$runtime[3], mod3$runtime[3], mod3a$runtime[3]) 
}
save.image("experiments.RData")

## From here:  https://github.com/michaelgreenacre/CODAinPractice/blob/master/Shang_West_East.csv
bronze <- read.csv("bronze.csv")
x <- as.matrix(bronze[, -c(1:4)])
x <- x / rowsums(x)
group <- as.numeric( as.factor( bronze[, 4] ) )
mod <- saa.qp(x, 3)
ternary(mod$W, colour = group)

for ( i in 1:R ) {
  mod1 <- saa.qp(x, 3, maxiter = 50000)
  mod2 <- saa.eg(x, 3, maxiter = 50000)
  mod3 <- saakl.eg(x, 3, maxiter = 50000)
  mod3a <- saakl.eg2(x, 3, maxiter = 1000)
  times[[ 11 ]][i, ] <- c( mod1$runtime[3], mod2$runtime[3], mod3$runtime[3], mod3a$runtime[3])
  obj[[ 11 ]][i, ] <- c( mod1$obj, mod2$obj, mod3$obj, mod3a$obj) 
  times[[ 11 ]][i, ] <-  c( mod1$runtime[3], mod2$runtime[3], mod3$runtime[3], mod3a$runtime[3]) 
}
save.image("experiments.RData")

## FADN
x <- as.matrix( read.csv("y.csv") )
x <- x / rowsums(x)
mod <- saa.qp(x, 3)
ternary(mod$W)

for ( i in 1:R ) {
  mod1 <- saa.qp(x, 3, maxiter = 50000)
  mod2 <- saa.eg(x, 3, maxiter = 50000)
  mod3 <- saakl.eg(x, 3, maxiter = 50000)
  mod3a <- saakl.eg2(x, 3, maxiter = 1000)
  times[[ 12 ]][i, ] <- c( mod1$runtime[3], mod2$runtime[3], mod3$runtime[3], mod3a$runtime[3])
  obj[[ 12 ]][i, ] <- c( mod1$obj, mod2$obj, mod3$obj, mod3a$obj) 
  times[[ 12 ]][i, ] <-  c( mod1$runtime[3], mod2$runtime[3], mod3$runtime[3], mod3a$runtime[3]) 
}
save.image("experiments.RData")

for ( i in 1:12 )  {
  ind <- which( obj[[ i ]] == 0 ) 
  if ( length(ind) > 0 ) obj[[ i ]][ind] <- NA
}
obj[[ 3 ]][, 4] <- NA
times[[ 3 ]][, 4] <- NA

mat <- matrix( nrow = 36, ncol = 4 )
colnames(mat) <- colnames(ti)
id <- matrix(1:36, ncol = 12)
rownames(mat) <- rep( c("Mean", "Min", "Max"), 12)
for ( j in 1:12 ) {
  a1 <- apply( obj[[ j ]], 2, mean, na.rm = TRUE )
  a2 <- apply( obj[[ j ]], 2, min, na.rm = TRUE )
  a3 <- apply( obj[[ j ]], 2, max, na.rm = TRUE )
  mat[id[, j], ] <- rbind(a1, a2, a3)  
}
xtable(mat, digits = 6)


mat <- matrix( nrow = 36, ncol = 4 )
colnames(mat) <- colnames(ti)
id <- matrix(1:36, ncol = 12)
rownames(mat) <- rep( c("Mean", "Min", "Max"), 12)
for ( j in 1:12 ) {
  mat[id[, j], ] <- rbind( colmeans(times[[ j ]]), colMinsMaxs(times[[ j ]]) )  
}
xtable(mat, digits = 6)





