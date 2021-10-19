i=1
while (file.exists(paste0("~/scratch/randomResults/random_populationsize_mutationload_prevalence_", i, "Lambda.txt"))){
  i=i+1;
}
write.table(i,file="LambdaCoef.txt",row.names = F,col.names = F)

