i=1
BurnIn = as.numeric(read.table("BurnIn.txt"))
BurnIn = ifelse(BurnIn==0,"F","T")
while (file.exists(paste0("~/scratch/dynastyResults/dynasty_nWF_populationsize_mutationload_prevalence_burnIn",BurnIn, "_",i,".txt"))){ 
 i=i+1;
}
write.table(i,file="index2.txt",row.names = F,col.names = F)
