i=1
while (file.exists(paste0("~/scratch/randomResults/random_nWF_ills_mutationload_difFitness_", i, ".txt"))){
 i=i+1;
}
write.table(i,file="index.txt",row.names = F,col.names = F)
