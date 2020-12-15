#!/usr/bin/Rscript

library("flexsurv")
library("rstpm2")
library("foreign")

load("model_for_positive_traning.Rdata")
load("model_for_mor_traning.Rdata")

load("data_for_positive_validation.Rdata")

stpm2cif_pos <- function(model_pos,model_mor,newdata,X) {
    stpm2if <- function(X) {
        ## setup up newdata
        newdat <- data.frame(as.data.frame(lapply(newdata, rep, 21)))
		newdat$hccpyear=X
		newdat$dupyear=X
		
		## calculate cause-specific survival
        Sk1 <- predict(model_pos, newdata=newdat)
		Sk2 <- predict(model_mor, newdata=newdat)
        S <- apply(cbind(Sk1,Sk2),1,prod) #±©Â¶·çÏÕ=Éú´æ·çÏÕ*Î´·¢²¡·çÏÕ#
        ## calculate cause-specific hazard
        S*predict(model_pos,newdata=newdat,type="hazard") #±©Â¶·çÏÕ*·¢²¡·çÏÕ#
    }
    sapply(X, function(xi)
           ifelse(xi==0,0,
                  integrate(stpm2if,0,xi)$value))
}

results_pos <- as.data.frame(matrix(NA,ncol=2,nrow=1))
colnames(results_pos) <- c("studyid","cumhaz")

results_pos$studyid="ID"
results_pos$cumhaz=stpm2cif_pos(model_pos,model_mor,val_pos[1,],10)	

write.csv(results_pos,"Post_estimation_of_positive_validation.csv",row.names=F)
