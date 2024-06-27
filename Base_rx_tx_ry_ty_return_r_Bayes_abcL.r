library(BGLR)

give_rx_tx_ry_ty_return_r_bayes_abcL=function(rx,tx,ry,ty){

##bayes A
 set.seed(1)
 eta<-list(list(X=rx,model="BayesA"))
 fm<- BGLR(y=ry,ETA=eta,nIter=4500,burnIn=500,verbose=FALSE)
 bata=fm$ETA[[1]]$b
 py=tx%*%bata
 ra=cor(py,ty)
 print(ra)
 
 py=NULL
 bata=NULL

 
##bayes B 
 set.seed(1)
 eta<-list(list(X=rx,model="BayesB"))
 fm<- BGLR(y=ry,ETA=eta,nIter=4500,burnIn=500,verbose=FALSE)
 bata=fm$ETA[[1]]$b
 py=tx%*%bata
 rb=cor(py,ty)
 print(rb)
 
 py=NULL
 bata=NULL

 
##bayes C 
 set.seed(1)
 eta<-list(list(X=rx,model="BayesC"))
 fm<- BGLR(y=ry,ETA=eta,nIter=4500,burnIn=500,verbose=FALSE)
 bata=fm$ETA[[1]]$b
 py=tx%*%bata
 rc=cor(py,ty)
 print(rc)
 
 py=NULL
 bata=NULL

 
##bayes L 
 set.seed(1)
 eta<-list(list(X=rx,model="BL"))
 fm<- BGLR(y=ry,ETA=eta,nIter=4500,burnIn=500,verbose=FALSE)
 bata=fm$ETA[[1]]$b
 py=tx%*%bata
 rL=cor(py,ty)
 print(rL)
 
 py=NULL
 bata=NULL

result=c(ra,rb,rc,rL)

return(result)

}