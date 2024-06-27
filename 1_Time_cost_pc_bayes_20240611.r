####
library(data.table)
gen<-fread(".../PIC/genotypes.txt",sep=",",header=TRUE)
phe <- read.csv(".../PIC/phenotypes.txt",header=TRUE)  

####
dim(phe)
phe[1:10,]

dim(gen)
gen[1:3,1:3]
####

yused=2

#### delete na 
position=which(phe[,yused] != ".")
y=as.numeric(phe[position,yused])
y[1:10]
length(y)
####

g=as.matrix(gen[,-1])
dim(g)
G=g[position,]
dim(G)

#### maf 0.05
freq=apply(G, 2, function(ggg) {
                                mu = mean(ggg)
                                f=mu/2
                                return(f)
})

pos=which(freq >= 0.95 | freq <= 0.05)
G=G[,-pos]
dim(G)



#### 
#### 
#### 
#### 
#### input X_ori y_ori
x_ori=apply(G, 2, function(ggg) {
                            mu = mean(ggg)
                            f=mu/2
                            (ggg - 2*f) 
})
y_ori=y

dim(x_ori)
length(y_ori)

#### set CV						   
## times
repeat_times <-1
## fold
cvnum <- 5

sample_size=length(y_ori)
cv_assign_ori<- cut(1:sample_size, cvnum, labels=FALSE)


#### create rx ry tx ty
			   
i=1
set.seed(i)
cv_assign <- sample(cv_assign_ori, length(cv_assign_ori))						   
						   
j=1
index <- which(cv_assign == j)  
ty=y_ori[index] 

##
rx=x_ori[-index,]
tx=x_ori[index,]
ry=y_ori[-index]
ty=y_ori[index]					   

##
library(BGLR)
##



##
##
##						   
## traditional bayesA
start_1=Sys.time()   ## START

set.seed(1)
eta<-list(list(X=rx,model="BayesA"))
fm<- BGLR(y=ry,ETA=eta,nIter=4500,burnIn=500,verbose=FALSE)
bata=fm$ETA[[1]]$b
py=tx%*%bata
ra=cor(py,ty)
end_1=Sys.time()     ## END

print(ra)            ## print cor
print(start_1)       ## print start time
print(end_1)         ## print end time
print(end_1-start_1) ## time cost

##
##
##						   
## PC bayesA

start_2=Sys.time()   ## START
xxt=tcrossprod(rx)
ea=eigen(xxt)
D=ea$values
U=ea$vectors
#delete the last
U_passed=U[,-nrow(U)]
D_passed=D[-nrow(U)]
#vmatrix
delta=sqrt(D_passed)
s_delta=1/delta
vmatrix_middle=crossprod(rx,U_passed)
vmatrix_trans=apply(vmatrix_middle,1,function(vmatrix_middle_each){  (vmatrix_middle_each*s_delta) })
#pcs
pcs_ori_rx=(apply(U_passed,1,function(u_each){  (u_each*delta) }))
pc_rx=t(pcs_ori_rx)	 
pc_tx=tcrossprod(tx,vmatrix_trans)
#bayes
set.seed(1)
eta<-list(list(X=pc_rx,model="BayesA"))
fm<- BGLR(y=ry,ETA=eta,nIter=4500,burnIn=500,verbose=FALSE)
bata=fm$ETA[[1]]$b
py=pc_tx%*%bata
ra=cor(py,ty)
end_2=Sys.time()   ## END

print(ra)            ## print cor
print(start_2)       ## print start time
print(end_2)         ## print end time
print(end_2-start_2) ## time cost