suppressMessages(require(ranger))
suppressMessages(require(PRROC))

OUTLIER_CUT = -785.42  # quantile(a[,87],0.99), all training data
#99% -785.4232

options<-commandArgs(trailingOnly=T)


if(length(options) < 2) stop("Invalid argument number\n\nRscript [].r [model file : data.rb] [indep test set]\n")
if(!is.na(options[1])) model_name=options[1]
if(!is.na(options[2])) test_fn=options[2]

load(model_name)    # data.rb

## test data, no out value 
blind_data <- read.table(test_fn)


blind_seq <- blind_data[,1];blind_data <- blind_data[,-1]
#nd <- ncol(blind_data)
#colnames(blind_data)[nd]="Outvalue"

final_blind_data.pred <- predict(data.rb,blind_data)
tmp_df <- data.frame(seq=blind_seq, prediction=round(final_blind_data.pred$predictions,5))
write.table(tmp_df,quote=F,col.names=F,row.names=F)



