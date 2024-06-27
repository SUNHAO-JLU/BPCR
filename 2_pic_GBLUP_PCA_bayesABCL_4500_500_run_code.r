
source("Base_rx_tx_ry_ty_return_r_Bayes_abcL.r")
###

library(data.table)
gen<-fread(".../PIC/genotypes.txt",sep=",",header=TRUE)
phe <- read.csv(".../PIC/phenotypes.txt",header=TRUE)   

dim(gen)
dim(phe)
head(phe)


##
## genotype & phenotype

yused=2
position=which(phe[,yused] != ".")
y_ori=as.numeric(phe[position,yused])
y_ori[1:3]
sample_size=length(y_ori)
sample_size


gen[1:3,1:3]
G=as.matrix(gen[,-1])
dim(G)
G=G[position,]
dim(G)
freq=apply(G, 2, function(ggg) {
  mu = mean(ggg)
  f=mu/2
})
pos=which(freq>=0.95 | freq<=0.05)
G=G[,-pos]
print(dim(G))




x_ori=apply(G,2,function(ggg){ mu=mean(ggg)
                               (ggg-mu)
                             })
						   
## kinship for hiblup				   
corg=tcrossprod(x_ori)
sample_size=nrow(corg)
print (sample_size)
kinship=corg/(1/sample_size*sum(diag(corg)))


## write_hiblup_GRM
write.table(kinship,"ww.txt",row.names=FALSE,col.names=FALSE)
write.table(c(1:sample_size) ,"wwid.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
system( ".../hiblup     --trans-xrm   --xrm-txt ww.txt    --xrm-id  wwid.txt --out ww_grm_hiblup ") 
	
#############################
#############################
#SET CV
##times
repeat_times <-10
##fold
cvnum <- 5

cv_assign_ori<- cut(1:sample_size, cvnum, labels=FALSE)

prediction_cor=matrix(-2,repeat_times*cvnum,5)

library(BGLR)



for(i in 1:repeat_times){
    # i=1
	  set.seed(i)
	  cv_assign <- sample(cv_assign_ori, length(cv_assign_ori))
    
	for(j in 1:cvnum){
	# j=1		
		y_NA <- y_ori
		index <- which(cv_assign == j)       ## tx
		y_NA[index] <- NA
		ty=y_ori[index] 
		
		pheno_NA=data.frame(c(1:length(y_NA)),y_NA)
		colnames(pheno_NA)=c("id","phe")
        write.table(pheno_NA,"pheno_NA.txt",quote = FALSE,col.names=TRUE,row.names=FALSE)

## hiblup	
system( " .../hiblup   --single-trait   --pheno pheno_NA.txt   --pheno-pos 2  --xrm  ww_grm_hiblup   --vc-method AI  --out gblup_cv")		
##
pre=read.table("gblup_cv.rand",header=TRUE)	 
head(pre)
id_na=pheno_NA[index,1]  
pre_result=data.frame( index,ty,pre[match(id_na,pre[,1]),])
head(pre_result)
rg= cor(pre_result[,2],pre_result[,4])
print(rg)

prediction_cor[(i - 1) * cvnum+j,1]=rg
system( "rm  gblup_cv.rand " )
rg=-999

##############
##############

rx=x_ori[-index,]
tx=x_ori[index,]
ry=y_ori[-index]
ty=y_ori[index]

## get PCS
rx=x_ori[-index,]
xxt=tcrossprod(rx)
ea=eigen(xxt)
D=ea$values
U=ea$vectors

U_passed=U[,-nrow(U)]
D_passed=D[-nrow(U)]

delta=sqrt(D_passed)
s_delta=1/delta
vmatrix_middle=crossprod(rx,U_passed)
vmatrix_trans=apply(vmatrix_middle,1,function(vmatrix_middle_each){  (vmatrix_middle_each*s_delta) })
#pcs
pcs_ori_rx=(apply(U_passed,1,function(u_each){  (u_each*delta) }))
pc_rx=t(pcs_ori_rx)	 
pc_tx=tcrossprod(tx,vmatrix_trans)



results=give_rx_tx_ry_ty_return_r_bayes_abcL(pc_rx,pc_tx,ry,ty)

prediction_cor[(i - 1) * cvnum+j,2]=results[1]
prediction_cor[(i - 1) * cvnum+j,3]=results[2]
prediction_cor[(i - 1) * cvnum+j,4]=results[3]
prediction_cor[(i - 1) * cvnum+j,5]=results[4]

print ((i - 1) * cvnum+j)
print ((i - 1) * cvnum+j)
print ((i - 1) * cvnum+j)

}
}
prediction_cor
apply(prediction_cor,2,mean)
apply(prediction_cor,2,sd)
