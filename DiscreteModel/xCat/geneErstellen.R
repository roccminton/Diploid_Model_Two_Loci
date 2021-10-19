setwd("~/scratch/")
#Chromosomenlaengen
c_length=rep(0,22); c_length[1]=2490; c_length[2]=2422; c_length[3]=1983; c_length[4]=1902;
c_length[5]=1815; c_length[6]=1708; c_length[7]=1593; c_length[8]=1451; c_length[9]=1384;
c_length[10]=1338; c_length[11]=1351; c_length[12]=1333; c_length[13]=1144; c_length[14]=1070;
c_length[15]=1020; c_length[16]=903; c_length[17]=833; c_length[18]=803; c_length[19]=586;
c_length[20]=644; c_length[21]=467; c_length[22]=508; c_length=c_length*10^5
#Enden der Chromosomen  
ends=rep(0,22);
tmp=0;
for (k in 1:22){
  tmp=tmp+c_length[k];
  ends[k]=tmp;
}
#Gene setzten inkl nonCoding regions
#Chromosomen in Abschnitte unterteilen zw 1500 und 30000 bp
base=0; j=1;
i=1
list=matrix(c(0,0,2,0),1,4);
for(i in 1:22){
  while(base+33000<ends[i]){
    list[j,1]=base;
    base=base+as.integer(runif(1,1650,33000));
    list[j,2]=base-1;
    j=j+1;     
    list=rbind(list,matrix(c(0,0,2,0),1,4));
  }
  list[j,1]=base;
  base=ends[i]+1;
  list[j,2]=base-1;
  if(i!=22){
    j=j+1;
    list=rbind(list,matrix(c(0,0,2,0),1,4));
  }
}
#zufaellig 300 Gene zuordnen
sam=sample(seq(1,length(list[,1]),by=1),300);
sam=sort(sam);
for(i in sam){
  list[i,3]=1;
  list[i,4]=1.2e-8;
}
write.table(list, "list2_2complete.txt")

#komprimieren der Liste
i<-1;
while(i<length(list[,1])){
  if(list[i,3]==list[i+1,3]){
    list[i+1,1]<-list[i,1]
    list<-list[-c(i),]
  } else{i<-i+1}  
  print(length(list[,1]))
}
#abspeichern
write.table(list[,1], "start2_2.txt",row.names = F)
write.table(list[,2], "end2_2.txt",row.names = F)
write.table(list[,3], "genomicType2_2.txt",row.names = F)
write.table(list[,4], "mutationRate2_2.txt",row.names = F)

