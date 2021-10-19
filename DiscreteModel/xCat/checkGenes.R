i=1
while (file.exists(paste0("~/scratch/randomResults/random_populationsize_mutationload_prevalence_", i, "Genes.txt"))){
  i=i+1;
}
write.table(i,file="Num_genes.txt",row.names = F,col.names = F)

