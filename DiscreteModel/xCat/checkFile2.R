i=1
while (file.exists(paste0("~/scratch/dynastyResults/dynasty_mutationload_prevalence_", i, ".txt"))){
 i=i+1;
}
write.table(i,file="index2.txt",row.names = F,col.names = F)
