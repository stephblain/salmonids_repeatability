#Functions




#Calculates observed chi square (lazily) and simulates null distribution for comparison
#Input: Matrix or dataframe with numeric values
#Output: observed chi square, BCa confidence interval, p value based on simulations, histogram simulated chi stats
simulate.chisq<-function(mat.x){
  
  a<-chisq.test(mat.x) #calculate observed stat
  #Randomly generate some tables using Patefield 1981 algorithm
  null.x<-r2dtable(10000,rowSums(mat.x),colSums(mat.x))
  
  z<-chisq.test(null.x[[1]]) #calculate first chisq
  y<-z$statistic #save to object y
  for(x in 2:10000){ #for each simulated matrix
    z<-chisq.test(null.x[[x]]) #calculate chisq
    y<-c(y,z$statistic)} #add chisq to object y
  b<-bca(y,conf.level=0.95) #function from library coxed
  b<-paste(b[1],b[2],sep=",")
  d<-y[y>a$statistic]
  p.value<-1-length(d)/length(y)
  print(paste("Observed stat:",a$statistic)) #report test stat
  print(paste("Chi/df:",a$statistic/a$parameter))
  print(paste("df:",a$parameter))
  #print(paste("BCa Confidence interval:",b)) #report CI
  print(paste("simulated p-value:",p.value))
  hist(y,main="Simulated chi stats")} #histogram of simulated stats

#format data to run k-sample Anderson-Darling test
#morphs.tr is morphs dataframe formatted for trait of interest
#trait_col is the name of the column of trait of interest means
#genus.ID is "Sal" for Salvelinus or "Cor" for Coregonus
#if adTest=T the performs A-D test, otherwise outputs formatted dataframe for plotting
#before running tests, calculates mean values for each ecotype from each lake with multiple samples
phenotype_repeat<-function(morphs.tr,trait_col,n_col,genus.ID,adTest){
  
  # #TEMP TEST DATA
  # morphs.tr<-morphs.fl
  # trait_col<-"Total_length"
  # n_col<-"length_n"
  # genus.ID<-"Sal"
  # adTest<-T
  
  
  
  morphs.tr$ID<-paste(morphs.tr$Paper_ID,morphs.tr$Lake,sep="_") #make dataframe for trait analysis
  morphs.tr<-morphs.tr[substr(morphs.tr$Species,1,3)==genus.ID,] #keep only one genus
  #start IDing lakes and selecting replicate
  morphs.tr.ID<-data.frame(table(morphs.tr$ID)) #get paper IDS
  morphs.tr.ID$Lake<-t(data.frame(strsplit(as.character(morphs.tr.ID$Var1),"_")))[,2] #extract lake name from paper ID
  morphs.tr.ID<-merge(morphs.tr.ID,df,by="Lake")[,c("Var1","Freq","Lake","Lacustrine_morphs")] #get # of morphs per lake
  morphs.tr.ID$Var1<-as.character(morphs.tr.ID$Var1)
  morphs.tr.ID<-morphs.tr.ID[morphs.tr.ID$Freq==morphs.tr.ID$Lacustrine_morphs,]#only keep cases where all morphs have a mean recorded
  lakes.rep<-as.character(data.frame(table(morphs.tr.ID$Lake)[table(morphs.tr.ID$Lake)>1])$Var1) #find lakes with repeated measures
  morphs.tr<-morphs.tr[morphs.tr$ID%in%morphs.tr.ID$Var1,] #only keep selected IDs
  
  morphs.tr<-merge(morphs.tr,df,by="Lake") #add info about lakes
  #only keep data from papers where all morphs have been identified
  temp2<-morphs.tr[substr(morphs.tr$Species,1,3)==genus.ID, #only keep one genus
                   c("ID","Lake",trait_col,"Species","Lacustrine_morphs","Morph_zone","Morph_diet")]
  
  colnames(temp2)[3]<-"trait"
  temp2<-arrange(temp2,Lake,ID,trait)
  # for(i in 1:length(lakes.rep)){ #pick one randomly from lakes with repeated measures
  #   Var1s<-morphs.tr.ID[morphs.tr.ID$Lake==lakes.rep[i],"Var1"]
  #   morphs.tr.ID<-morphs.tr.ID[!morphs.tr.ID$Var1%in%Var1s[Var1s!=sample(Var1s,1)],]
  # }
  if(length(lakes.rep)>0){
    for(i in 1:length(lakes.rep)){ 
      Var1s<-morphs.tr.ID[morphs.tr.ID$Lake==lakes.rep[i],"Var1"]
      Var.dat<-morphs.tr[morphs.tr$ID%in%Var1s,c("ID",trait_col,n_col)]
      colnames(Var.dat)<-c("ID","trait","tr_n")
      Var.dat<-arrange(Var.dat,ID,trait)
      n<-length(Var.dat$ID)/nlevels(as.factor(Var.dat$ID))
      k<-vector()
      for(j in 1:n){ #make vector with means of each ecotype
        k.i<-seq(j,nrow(Var.dat),n)
        k.mean<-sum(Var.dat$trait[k.i]*Var.dat$tr_n[k.i])/sum(Var.dat$tr_n[k.i])
        k<-c(k,k.mean)}
      temp2<-temp2[temp2$Lake!=lakes.rep[1],]
      temp2<-rbind(temp2,data.frame(ID=paste("01E",lakes.rep[i],sep="_"),
                                    Lake=lakes.rep[i],trait=k,
                                    Species=temp2[temp2$Lake==lakes.rep[i],"Species"][1],
                                    Lacustrine_morphs=temp2[temp2$Lake==lakes.rep[i],"Lacustrine_morphs"][1],
                                    Morph_zone=temp2[temp2$Lake==lakes.rep[i],"Morph_zone"][1:n],
                                    Morph_diet=temp2[temp2$Lake==lakes.rep[i],"Morph_diet"][1:n]))}}
  temp3<-temp2[grep("01E",temp2$ID),]
  temp2<-temp2[!temp2$Lake%in%lakes.rep,]
  temp2<-rbind(temp2,temp3); rm(temp3)
  temp2$ID2<-paste(temp2$ID,1:nrow(temp2),sep="_")
  
  for(i in unique(temp2$ID)){
    df.i<-temp2[temp2$ID==i,]
    df.i$trait-mean(df.i$trait)
  }
  
  
  
  temp1<-as.data.frame(pivot_wider(temp2,names_from=Lake,values_from=trait)) #means in columns with each lake as a column name
  temp1
  if(adTest==T){x<-ad.test(temp1[,7:ncol(temp1)],method="simulated")
  print(x)
  print(1-x$ad[1,3])}else{ #run ad test on numeric columns
    temp2 } }


byMinusOne<-function(x){x*(-1)} #multiply by -1

#Test for association between different morphs using presence-absence matrices
#input counts matrix, with lakes as columns and ecotypes as rows
pr.ab.tests<-function(mat.x){
  mat.x<-t(mat.x)
  if(is.numeric(mat.x)){
    mat.x[mat.x%in%c(2,3,4)]<-1 #make matrix presence/absence
    
    #site (row) frequencies preserved, column (species) frequencies proportional to observed
    print("V RATIO OUTPUT")
    print(oecosimu(mat.x,V.ratio,"r1",nsimul=10000,alternative="two.sided"))
    print("NESTEDNESS OUTPUT")
    print(oecosimu(mat.x,nestednodf,"r1",nsimul=10000,alternative="greater"))}else{
      print("matrix not numeric")}}


