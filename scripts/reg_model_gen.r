require(ranger)
require(PRROC)

options<-commandArgs(trailingOnly=T)
if(length(options) < 1) stop("Invalid argument number\n\nRscript [].r [seq_ddg_values.feat.mat]\n")
if(!is.na(options[1])) data_fn=options[1]
if(!is.na(options[2])) model_name=options[2]

data <- read.table(data_fn) # tScale 85-dim vec 
sc = data[,ncol(data)]; data <- data[,-ncol(data)]
energy = data[,ncol(data)]; data <- data[,-ncol(data)]
ddg = data[,ncol(data)]; data <- data[,-ncol(data)]
data_seq <- data[,1];data <- data[,-1]

tot_n = nrow(data)

combined_data = cbind(data,ddg)
nd = ncol(combined_data)

if(model_name == "ddg"){
	# ddg prediction model, change the col name to energy or sc 
	data.rb <- ranger(ddg ~., combined_data,num.trees = 100)
}else if(model_name == "energy"){
	data.rb <- ranger(energy ~., combined_data,num.trees = 100)
}else if(model_name == "sc"){
	data.rb <- ranger(sc ~., combined_data,num.trees = 100)
}else{
	stop("Model name incorrect!")
}

print(data.rb)

df <- data.frame(pred = data.rb$predictions, real = combined_data[,nd])
write.table(df, file=paste(model_name,".oob_pred.txt",sep=""),quote=F,col.names=F,row.names=F)
name_str=paste(model_name,".pred_model.ranger",sep="")
print(paste(name_str,"model saved as a file"))
save(data.rb, file=name_str)



#[1] "1 44 ROC 0.778755543301915 LeftPosScore 0.5856603937019 Error 0.194811931671842"
#[1] "2 28 ROC 0.885156019496771 Error 0.123839748037717"
#[1] "3 5 ROC 0.851630351520609 Error 0.129200557940455"
#[1] "4 30 ROC 0.828119611092685 Error 0.171449792144822"
#[1] "5 5 ROC 0.840267061042524 Error 0.120274707740581"



