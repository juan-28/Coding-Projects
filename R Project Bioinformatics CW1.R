library(MASS)
library(brnn)
library(nnet)

setwd("~/Desktop/bio cw1")

# data preparation #
x=read.csv('ecoli.two.class.csv', row.names = 1)
inp=x[,which(colnames(x)!='class')]
# data visualization using PCA #
pca=prcomp(inp)
head (pca)
ratio=round(100*sum(pca$sdev[1:2]^2)/sum(pca$sdev^2), digits =2)
print(ratio)
my.col=c('blue','red')[x$class]
plot(pca$x[,1:2], col =my.col, pch = 19, xlab = 'PC1', ylab = 'PC2', main = 'PCA Model')
barplot(pca$sdev, las=2, main= 'PC Variance')
legend("topright", paste("Var prop=0.75"))


# applying LDA model to data set #
ldamodel = lda(class~.,x)

pred.lda=predict(ldamodel)
head(pred.lda$posterior)
prob=pred.lda$posterior[,2]


# cross validation for LDA #
cv=sample(1:5, nrow(x),replace = TRUE)
prob=rep(0,nrow(x))
for (fold in 1:5)
{
  test=which(cv==fold)
  model=lda(class ~.,x[-test,])
  prob[test]=predict(model,
                     newdata=inp[test,])$posterior[,2]
}


# visualizing LDA model #
par(mfrow=c(1,1))
hist(pred.lda$x[,1],nclass = 30, col = "blue", main="LDA",xlab='LDA Projection',cex.axis=1.3,cex.lab=1.5,prob=
       TRUE)
lines(density(pred.lda$x[,1]), col ="red",lwd=2)
conf.lda=table(x$class, 1*(prob>0.5))
head(conf.lda)


# 5 measurements lda #
spe=conf.lda[1,1]/sum(conf.lda[1,])
sen=conf.lda[2,2]/sum(conf.lda[2,])
tot=sum(diag(conf.lda))/sum((conf.lda))
npp=conf.lda[1,1]/sum(conf.lda[,1])
ppp=conf.lda[2,2]/sum(conf.lda[,2])
print(spe)
print(sen)
print(tot)
print(npp)
print(ppp)


# Neural network algorithm #
inp.nn=x[,which(colnames(x)!='class')]
input= as.matrix(inp.nn)
label=1*(x$class == 1)
model.nn=brnn(input,label,neurons = 17)
pred.nn=predict(model.nn)
conf.nn= table(label, 1* (pred.nn>0.5))
print(conf.nn)


# density analysis for neural network #
par(mfrow=c(2,1))
hist(pred.nn[which(x$class==1)], nclass=30,xlab='group1')
hist(pred.nn[which(x$class==2)], nclass=30,xlab='group2')


#5 measurements for neural network#
spe.nn=conf.nn[1,1]/sum(conf.nn[1,])
sen.nn=conf.nn[2,2]/sum(conf.nn[2,])
tot.nn=sum(diag(conf.nn))/sum((conf.nn))
npp.nn=conf.nn[1,1]/sum(conf.nn[,1])
ppp.nn=conf.nn[2,2]/sum(conf.nn[,2])
print(spe.nn)
print(sen.nn)
print(tot.nn)
print(npp.nn)
print(ppp.nn)


# cross validation neural network #
cv.nn=sample(1:5, nrow(x),replace = TRUE)
prob=rep(0,nrow(x))
for (fold in 1:5)
{
  test=which(cv.nn==fold)
  model=brnn(class ~.,x[-test,], neurons=17)
  prob[test]=predict(model,
                     newdata=inp[test,],type='raw')
}


# Random forest algortithm #
library(randomForest)
model.forest=randomForest(class~., data =x, importance=TRUE, ntree=50)
barplot(importance(model.forest), beside = TRUE,legend=TRUE)
y.hat=predict(model.forest,inp)
conf=table(x$class, 1*(y.hat>1.5))
print(conf)


1# density analysis for random forest #
par(mfrow=c(2,1))
hist(y.hat[which(x$class==1)], nclass=30,xlab='group1')
hist(y.hat[which(x$class==2)], nclass=30,xlab='group2')


# 5 measurements random forest #
spe.rf=conf[1,1]/sum(conf[1,])
sen.rf=conf[2,2]/sum(conf[2,])
tot.rf=sum(diag(conf))/sum((conf))
npp.rf=conf[1,1]/sum(conf[,1])
ppp.rf=conf[2,2]/sum(conf[,2])
print(spe.rf)
print(sen.rf)
print(tot.rf)
print(npp.rf)
print(ppp.rf)


# ROC analysis for LDA #
library(ROCR)
roc.obj=prediction(prob,x$class)
curve.obj=performance(roc.obj, 'tpr', 'fpr')
auc.obj=performance(roc.obj, measure = 'auc', fpr.stop=0.1)
auc=round(auc.obj@y.values[[1]], digits = 3)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(curve.obj, lwd=2, xlim=c(0,0.1))
legend('topleft', paste('AUC', auc))


