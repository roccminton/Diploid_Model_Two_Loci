i=1
while (file.exists(paste0("~/scratch/dynastyResults/dynasty_grandparents4_",i,".txt"))){
 i=i+1;
}
write.table(i,file="index3.txt",row.names = F,col.names = F)
