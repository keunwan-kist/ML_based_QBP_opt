# Keunwan Park

options<-commandArgs(trailingOnly=T)
if(length(options) < 1) stop("Invalid argument number\n\nRscript [].r [seq file: sequence header]\n")
if(!is.na(options[1])) seq_fn=options[1]

suppressMessages(require(Peptides))
suppressMessages(require(hash))

AA = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")

tScale_h = hash()
for(aa in AA){
	tScale_h[aa] = unlist(tScales(seq = aa))	# 5-dim vector
}

conv2tScaleVec <- function(s)
{
	sv = unlist(strsplit(s, split= ""))
	out_vec <- c()
	for(s in sv){
		out_vec <- c(out_vec,values(tScale_h[s])[,1])
	}
	out_vec
}


#a<-read.table(seq_fn,header=T)
a<-read.table(seq_fn)

#seq_vec = a[,"sequence"]
seq_vec = a[,1]

all_feat_mat <- matrix(0,nrow(a),5*nchar(as.character(seq_vec[1])))

cnt=1
for(seq in seq_vec){
	all_feat_mat[cnt,] = conv2tScaleVec(seq)		
	cnt = cnt + 1
}

combined_data = data.frame(seq_vec,all_feat_mat)

write.table(combined_data, quote=F,col.names=F,row.names=F)
