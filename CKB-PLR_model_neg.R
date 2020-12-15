#!/usr/bin/Rscript

library("flexsurv")
library("rstpm2")
library("foreign")
load("model_for_negative_traning.Rdata")
load("model_for_mor_traning.Rdata")

load("data_for_negative_validation.Rdata")


stpm2cif_nev <- function(model_neg,model_mor,newdata,X) {
    stpm2if <- function(X) {
        ## setup up newdata
        newdat <- data.frame(as.data.frame(lapply(newdata, rep, 21)))
		newdat$hccpyear=X
		newdat$dupyear=X
		
		## calculate cause-specific survival
        Sk1 <- predict(model_neg, newdata=newdat)
		Sk2 <- predict(model_mor, newdata=newdat)
        S <- apply(cbind(Sk1,Sk2),1,prod) #±©Â¶·çÏÕ=Éú´æ·çÏÕ*Î´·¢²¡·çÏÕ#
        ## calculate cause-specific hazard
        S*predict(model_neg,newdata=newdat,type="hazard") #±©Â¶·çÏÕ*·¢²¡·çÏÕ#
    }
    sapply(X, function(xi)
           ifelse(xi==0,0,
                  integrate(stpm2if,0,xi)$value))
}


results_neg <- as.data.frame(matrix(NA,ncol=2,nrow=1))
colnames(results_neg) <- c("studyid","cumhaz")

results_neg$studyid="ID"
results_neg$cumhaz=stpm2cif_nev(model_neg,model_mor,val_neg[1,],10)	

write.csv(results_neg,"Post_estimation_of_negative_validation.csv",row.names=F)
