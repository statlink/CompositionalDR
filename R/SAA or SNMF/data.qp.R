library(Rfast)
library(robCompositions)
library(ToolsForCoDa)
library(xtable)
library(coda.base)
library(CompositionalDR)
source("saakl.eg2.R")

### economics
economics <- read.table("economics.txt", header = TRUE)
x <- as.matrix(economics[,2:6])
x <- x/rowsums(x)
group <- economics[,7]
D <- dim(x)[2]
sse <- numeric(D - 2)
for ( k in 3:D )  sse[k - 2] <- saa.qp(x, k, maxiter = 50000)$obj 
png(filename = "data1.png", width = 3000, height = 3000, res = 600)
plot(sse, type = "b", pch = 16, col = 4, lwd = 2, xlab = "Projected dimensionality (k)", 
ylab = "Frobenius norm", xaxt = "n", cex.lab = 1.3, cex.axis = 1.3)
axis(1, at = 1:(D-2), lab = 3:D)
dev.off()

## pottery
pottery <- read.table("pottery.txt", header = TRUE)
x <- as.matrix(pottery[, -1])
x <- x/rowsums(x)
D <- dim(x)[2]
sse <- numeric(D - 2)
for ( k in 3:D )  sse[k - 2] <- saa.qp(x, k, maxiter = 50000)$obj 
png(filename = "data2.png", width = 3000, height = 3000, res = 600)
plot(sse, type = "b", pch = 16, col = 4, lwd = 2, xlab = "Projected dimensionality (k)", 
ylab = "Frobenius norm", xaxt = "n", cex.lab = 1.3, cex.axis = 1.3)
axis(1, at = 1:(D-2), lab = 3:D)
dev.off()

### Package robCompositions
## chorizonDL
x <- chorizonDL 
x <- as.matrix(x[, -c(1:3)])
x <- x[, 1:53]
x <- x[, -c(3, 10, 20, 23, 29, 31, 33, 37, 43, 48, 53)]
x <- x/rowsums(x)
D <- dim(x)[2]
sse <- numeric(18)
for ( k in 3:20 )  sse[k - 2] <- saa.qp(x, k, maxiter = 50000)$obj 
png(filename = "data3.png", width = 3000, height = 3000, res = 600)
plot(sse[1:18], type = "b", pch = 16, col = 4, lwd = 2, xlab = "Projected dimensionality (k)", 
ylab = "Frobenius norm", xaxt = "n", cex.lab = 1.3, cex.axis = 1.3)
axis(1, at = 1:18, lab = 3:20)
dev.off()

### gemas
x <- gemas[, -c(1:8)]
x <- as.matrix(x)
x <- x/ rowsums(x)
id <- which( is.na( rowsums(x) ) )
x <- x[-id, ]
group <- as.numeric( as.factor(gemas[-id, 8]) )
D <- dim(x)[2]
sse <- numeric(D - 2)
for ( k in 3:D )  sse[k - 2] <- saa.qp(x, k, maxiter = 50000)$obj 
png(filename = "data4.png", width = 3000, height = 3000, res = 600)
plot(sse, type = "b", pch = 16, col = 4, lwd = 2, xlab = "Projected dimensionality (k)", 
ylab = "Frobenius norm", xaxt = "n", cex.lab = 1.3, cex.axis = 1.3)
axis(1, at = 1:(D-2), lab = 3:D)
dev.off()

## gjovik
x <- gjovik
x <- as.matrix(x[, -c(1:10)])
x <- x/ rowsums(x)
id <- which( is.na( rowsums(x) ) )
x <- x[-id, ]
group <- as.numeric( as.factor(gjovik[-id, 10]) )
D <- dim(x)[2]
sse <- numeric(10)
for ( k in 3:12 )  sse[k - 2] <- saa.qp(x, k, maxiter = 50000)$obj 
png(filename = "data5.png", width = 3000, height = 3000, res = 600)
plot(sse, type = "b", pch = 16, col = 4, lwd = 2, xlab = "Projected dimensionality (k)", 
ylab = "Frobenius norm", xaxt = "n", cex.lab = 1.3, cex.axis = 1.3)
axis(1, at = 1:10, lab = 3:12)
dev.off()

## honey
x <- honey
group <- as.numeric(x[, 1])
x <- as.matrix( x[, -c(1:5)] )
x[which( is.na(x)) ] <- 0
x <- x / rowsums(x)
D <- dim(x)[2]
sse <- numeric(D - 2)
for ( k in 3:D )  sse[k - 2] <- saa.qp(x, k, maxiter = 50000)$obj 
png(filename = "data6.png", width = 3000, height = 3000, res = 600)
plot(3:D, sse, type = "b", pch = 16, col = 4, lwd = 2, xlab = "Projected dimensionality (k)", 
ylab = "Frobenius norm", xaxt = "n", cex.lab = 1.3, cex.axis = 1.3)
axis(1, at = 1:(D-2), lab = 3:D)
dev.off()

## nutrients
x <- as.data.frame(nutrients)
x <- x[, -c(1:12)]
x <- as.matrix( x[, c(1:30, 36)] )
x <- na.omit(x)
x <- x / rowsums(x)
D <- dim(x)[2]
sse <- numeric(13)
for ( k in 3:15 )  sse[k - 2] <- saa.qp(x, k, maxiter = 50000)$obj 
png(filename = "data7.png", width = 3000, height = 3000, res = 600)
plot(sse, type = "b", pch = 16, col = 4, lwd = 2, xlab = "Projected dimensionality (k)", 
ylab = "Frobenius norm", xaxt = "n", cex.lab = 1.3, cex.axis = 1.3)
axis(1, at = 1:13, lab = 3:15)
dev.off()

## payments
x <- as.matrix( payments[, -c(1:6)] )
x <- x / rowsums(x)
group <- as.numeric(payments[, 2]) - 1
D <- dim(x)[2]
sse <- numeric(D - 2)
for ( k in 3:D )  sse[k - 2] <- saa.qp(x, k, maxiter = 50000)$obj 
png(filename = "data8.png", width = 3000, height = 3000, res = 600)
plot(3:D, sse, type = "b", pch = 16, col = 4, lwd = 2, xlab = "Projected dimensionality (k)", 
ylab = "Frobenius norm", xaxt = "n", cex.lab = 1.3, cex.axis = 1.3)
axis(1, at = 1:(D-2), lab = 3:D)
dev.off()

## coda.base
x <- bmi_activity
x <- as.matrix(x[, 3:7])
x <- x/rowsums(x)
group <- as.numeric( as.factor(bmi_activity[, 8]) )
D <- dim(x)[2]
sse <- numeric(D - 2)
for ( k in 3:D )  sse[k - 2] <- saa.qp(x, k, maxiter = 50000)$obj 
png(filename = "data9.png", width = 3000, height = 3000, res = 600)
plot(3:D, sse, type = "b", pch = 16, col = 4, lwd = 2, xlab = "Projected dimensionality (k)", 
ylab = "Frobenius norm", xaxt = "n", cex.lab = 1.3, cex.axis = 1.3)
axis(1, at = 1:(D-2), lab = 3:D)
dev.off()

x <- waste
x <- as.matrix(x[, 5:9])
x <- x/rowsums(x)
group <- as.numeric( as.factor(waste[, 4]) )
D <- dim(x)[2]
sse <- numeric(D - 2)
for ( k in 3:D )  sse[k - 2] <- saa.qp(x, k, maxiter = 50000)$obj 
png(filename = "data10.png", width = 3000, height = 3000, res = 600)
plot(3:D, sse, type = "b", pch = 16, col = 4, lwd = 2, xlab = "Projected dimensionality (k)", 
ylab = "Frobenius norm", xaxt = "n", cex.lab = 1.3, cex.axis = 1.3)
axis(1, at = 1:(D-2), lab = 3:D)
dev.off()

## From here:  https://github.com/michaelgreenacre/CODAinPractice/blob/master/Shang_West_East.csv
bronze <- read.csv("bronze.csv")
group <- as.numeric( as.factor( bronze[, 4] ) )
x <- as.matrix(bronze[, -c(1:4)])
x <- x / rowsums(x)
D <- dim(x)[2]
sse <- numeric(D - 2)
for ( k in 3:D )  sse[k - 2] <- saa.qp(x, k, maxiter = 50000)$obj 
png(filename = "data11.png", width = 3000, height = 3000, res = 600)
plot(3:D, sse, type = "b", pch = 16, col = 4, lwd = 2, xlab = "Projected dimensionality (k)", 
ylab = "Frobenius norm", xaxt = "n", cex.lab = 1.3, cex.axis = 1.3)
axis(1, at = 1:(D-2), lab = 3:D)
dev.off()

## FADN
x <- as.matrix( read.csv("y.csv") )
x <- x / rowsums(x)
D <- dim(x)[2]
sse <- numeric(D - 2)
for ( k in 3:D )  sse[k - 2] <- saa.qp(x, k, maxiter = 50000)$obj 
png(filename = "data12.png", width = 3000, height = 3000, res = 600)
plot(3:D, sse, type = "b", pch = 16, col = 4, lwd = 2, xlab = "Projected dimensionality (k)", 
ylab = "Frobenius norm", xaxt = "n", cex.lab = 1.3, cex.axis = 1.3)
axis(1, at = 1:(D-2), lab = 3:D)
dev.off()





