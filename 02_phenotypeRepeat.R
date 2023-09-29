setwd("C:/Users/Steph/OneDrive - Texas A&M University/Salmonids/writing/morphs/dryad_upload/")

###Read in packages
pkgs<-c("tidyverse","coxed","gridExtra","matrixcalc","forcats",
        "vegan","rgeos","rnaturalearth", "rnaturalearthdata","sf","bipartite",
        "kSamples","ggpubr","metaDigitise","metafor","ape","rfishbase")

lapply(pkgs,library,character.only=T); rm(pkgs)
source("C:/Users/Steph/OneDrive - Texas A&M University/Salmonids/salmonids/functions_morphs_GEB.R")

###Colours
col.cor<-"#8cb5a7"
col.sal<-"#b58c93"
cols.sal<-c("#241c1d","#5b464a","#917076","#c4a3a9")
cols.cor<-c("#2a3632","#546d64","#709186","#a3c4b9")


###Read in data

df<-read.csv("lakes.csv",stringsAsFactors=F)
m<-read.csv("ecotypes_sum.csv",stringsAsFactors=F)
morphs<-read.csv("ecotypes.csv",stringsAsFactors=F,fileEncoding="latin1")

m$Genus<-sub("_.*","",m$Species)


##################################################
#repeatability of phenotype distributions
##################################################

#prep gill raker dataframe
morphs$GRc_e<-as.numeric(morphs$GRc_e); morphs$GRc_n<-as.numeric(morphs$GRc_n)
morphs$GRc_mean<-as.numeric(morphs$GRc_mean)
morphs.GR<-morphs[!is.na(morphs$GRc_n)&!is.na(morphs$GRc_e)&!is.na(morphs$GRc_mean),]


#phenotype repeatability - gill rakers
phenotype_repeat(morphs.GR,"GRc_mean","GRc_n","Cor",T)
phenotype_repeat(morphs.GR,"GRc_mean","GRc_n","Sal",T)


#convert fork etc lengths to total length with fish base conversion
morphs$length_e<-as.numeric(morphs$length_e); morphs$length_n<-as.numeric(morphs$length_n)
morphs$length_mean<-as.numeric(morphs$length_mean)
morphs.fl<-morphs[!is.na(morphs$length_n)&!is.na(morphs$length_e)&!is.na(morphs$length_mean)&
                    morphs$length_type%in%c("fork","total","standard"),]
morphs.fl$Total_length<-NA

morphs.fl$length_type[morphs.fl$length_type=="fork"]<-"FL"
morphs.fl$length_type[morphs.fl$length_type=="total"]<-"TL"
morphs.fl$length_type[morphs.fl$length_type=="standard"]<-"SL"


#make new vector for variances
morphs.fl$Total_length_var<-morphs.fl$length_e
morphs.fl$Total_length_var<-as.numeric(morphs.fl$Total_length_var)
#convert SD, SE and 95CI to Variance
morphs.fl$Total_length_var[morphs.fl$length_e_type=="SD"]<-
  morphs.fl$Total_length_var[morphs.fl$length_e_type=="SD"]^2

morphs.fl$Total_length_var[morphs.fl$length_e_type=="SE"]<-
  se_to_sd(morphs.fl$Total_length_var[morphs.fl$length_e_type=="SE"],
           morphs.fl$length_n[morphs.fl$length_e_type=="SE"])^2

morphs.fl$Total_length_var[morphs.fl$length_e_type=="95CI"]<-
  CI95_to_sd(morphs.fl$Total_length_var[morphs.fl$length_e_type=="95CI"],
             morphs.fl$length_n[morphs.fl$length_e_type=="95CI"])^2


#fish base function to extract conversion estimates
to_TL<-rfishbase::length_length(sub("_"," ",levels(as.factor(morphs.fl$Species))))
to_TL$a==0 # all a values are zero - need b's only
to_TL$b<-as.numeric(to_TL$b)
to_TL$Length1_2<-paste(to_TL$Length1,to_TL$Length2,sep="_")

sp1<-levels(as.factor(morphs.fl$Species))[1]

for(sp1 in levels(as.factor(morphs.fl$Species))){
  for(lt1 in c("FL","SL")){
    morphs.lt1.sp<-morphs.fl[morphs.fl$Species==sp1&morphs.fl$length_type==lt1,]
    to_TL1<-to_TL[to_TL$Species==sub("_"," ",sp1)&to_TL$Length1_2%in%
                    c(paste(lt1,"TL",sep="_"),paste("TL",lt1,sep="_")),]
    if(nrow(to_TL1)>1){to_TL1<-to_TL1[sample(1:nrow(to_TL1),1),]}
    if(nrow(morphs.lt1.sp)>0&nrow(to_TL1)>0){
      if(to_TL1$Length1_2==paste("TL",lt1,sep="_")){
        #multiply by constant to adjust mean
        morphs.fl[morphs.fl$Species==sp1&morphs.fl$length_type==lt1,]$Total_length<-
          to_TL1$b*morphs.lt1.sp$length_mean
        #multiple by square constant to adjust variance
        morphs.fl[morphs.fl$Species==sp1&morphs.fl$length_type==lt1,]$Total_length_var<-
          to_TL1$b^2*morphs.lt1.sp$Total_length_var}
      if(to_TL1$Length1_2==paste(lt1,"TL",sep="_")){
        #multiply by 1/constant
        morphs.fl[morphs.fl$Species==sp1&morphs.fl$length_type==lt1,]$Total_length<-
          morphs.lt1.sp$length_mean/to_TL1$b
        #multiple by 1/constant squared
        morphs.fl[morphs.fl$Species==sp1&morphs.fl$length_type==lt1,]$Total_length_var<-
          morphs.lt1.sp$Total_length_var*(1/to_TL1$b)^2}}}}
morphs.fl[morphs.fl$length_type=="TL","Total_length"]<-morphs.fl[morphs.fl$length_type=="TL","length_mean"]

#manual typo fixes
morphs.fl[morphs.fl$Lake=="Muddus",c("Total_length","Paper_ID")]
morphs.fl[morphs.fl$Lake=="Muddus"&morphs.fl$Paper_ID%in%c("40B","55B"),
          c("Total_length")]<-
  morphs.fl[morphs.fl$Lake=="Muddus"&morphs.fl$Paper_ID%in%c("40B","55B"),
            c("Total_length")]*10

#evaluate repeatability for length
phenotype_repeat(morphs.fl,"Total_length","length_n","Sal",T)
phenotype_repeat(morphs.fl,"Total_length","length_n","Cor",T)

morphs.fl.pr.S<-phenotype_repeat(morphs.fl,"Total_length","length_n","Sal",F)
morphs.fl.pr.S$Lake<-factor(morphs.fl.pr.S$Lake,
                            levels=rev(unique(arrange(morphs.fl.pr.S,Lacustrine_morphs,trait)$Lake)))
morphs.fl.pr.C<-phenotype_repeat(morphs.fl,"Total_length","length_n","Cor",F)
morphs.fl.pr.C$Lake<-factor(morphs.fl.pr.C$Lake,
                            levels=rev(unique(arrange(morphs.fl.pr.C,Lacustrine_morphs,trait)$Lake)))

#plot by species complex
A<-ggplot(data=morphs.fl.pr.S,
          aes(y=Lake,x=trait,shape=Species))+
  geom_point(size=4,colour=col.sal,stroke=2)+
  scale_shape_manual(values=c(0,1,2,5))+
  theme_classic()+xlab("Total Length")+ylab("Lake")+
  theme(panel.grid.major.y=element_line(size=.1,color="grey80"))

B<-ggplot(data=morphs.fl.pr.C,
          aes(y=Lake,x=trait,shape=Species))+
  geom_point(size=4,color=col.cor,stroke=2)+
  scale_shape_manual(values=c(0,1,2,5))+
  theme_classic()+xlab("Total Length")+ylab("Lake")+
  
  theme(panel.grid.major.y=element_line(size=.1,color="grey80"))

ggarrange(A,B,ncol=2,labels=c("A","B"))

morphs.GR.S.pr<-phenotype_repeat(morphs.GR,"GRc_mean","GRc_n","Sal",F)
morphs.GR.S.pr$Lake<-factor(morphs.GR.S.pr$Lake,
                            levels=rev(unique(arrange(morphs.GR.S.pr,Lacustrine_morphs,trait)$Lake)))
morphs.GR.C.pr<-phenotype_repeat(morphs.GR,"GRc_mean","GRc_n","Cor",F)
morphs.GR.C.pr$Lake<-factor(morphs.GR.C.pr$Lake,
                            levels=rev(unique(arrange(morphs.GR.C.pr,Lacustrine_morphs,trait)$Lake)))


#And gill rakers by species complex
A<-ggplot(data=morphs.GR.S.pr,
          aes(y=Lake,x=trait,colour=Species))+
  geom_point(size=4)+
  theme_classic()+xlab("Gill Raker Count")+ylab("Lake")+
  scale_color_manual(values=cols.sal)+
  theme(panel.grid.major.y=element_line(size=.1,color="grey80"))
B<-ggplot(data=morphs.GR.C.pr,
          aes(y=Lake,x=trait,colour=Species))+
  geom_point(size=4)+
  theme_classic()+xlab("Gill Raker Count")+ylab("Lake")+
  scale_color_manual(values=cols.cor)+
  theme(panel.grid.major.y=element_line(size=.1,color="grey80"))
ggarrange(A,B,ncol=2,labels=c("A","B"))


#get list of lakes used for phenotype repeatability
lakes.list.trait<-unique(c(as.character(morphs.fl.pr.C$Lake),as.character(morphs.fl.pr.S$Lake),
                           as.character(morphs.GR.C.pr$Lake),as.character(morphs.GR.S.pr$Lake)))


##################################################
#Phenotype repeatability - common assemblage types
##################################################

##Start with habitat and diet assemblage histograms

m.zone<-subset(m,!is.na(m$Morph_zone))
m.zone$Morph_zone<-substr(m.zone$Morph_zone,1,3)
#check lakes
temp<-cbind(table(m.zone$Lake),df$Lacustrine_morphs[df$Lake%in%m.zone$Lake],df$Lake[df$Lake%in%m.zone$Lake])
m.zone<-subset(m.zone,m.zone$Lake%in%temp[,3][temp[,1]==temp[,2]])
y<-data.frame(Lake="x",Morphs="x")

for(l in levels(as.factor(m.zone$Lake))){
  x<-data.frame(Lake=l,Morphs=paste(sort(
    as.character(m.zone[m.zone$Lake==l,"Morph_zone"])),collapse="\n"))
  y<-rbind(y,x)
}
y<-y[-1,]
z.z<-y%>%group_by(Morphs)%>%mutate(counts=n())
z.z<-distinct(left_join(z.z,m,by="Lake")[,c("Lake","Morphs","counts","Genus","Species")])

A<-ggplot(z.z,aes(x=reorder(Morphs,-counts),fill=Species))+
  geom_bar(stat="count")+theme_classic()+
  scale_fill_manual(values=c(cols.cor,cols.sal))+
  theme(legend.position="none")+
  xlab("Ecotypes present - Habitat")+ylab("Count")

#diet histogram
m.diet<-subset(m,!is.na(m$Morph_diet))
m.diet$Morph_diet<-substr(m.diet$Morph_diet,1,3)
#check lakes
temp<-cbind(table(m.diet$Lake),df$Lacustrine_morphs[df$Lake%in%m.diet$Lake],
            df$Lake[df$Lake%in%m.diet$Lake])
m.diet<-subset(m.diet,m.diet$Lake%in%temp[,3][temp[,1]==temp[,2]])
y<-data.frame(Lake="x",Morphs="x")

for(l in levels(as.factor(m.diet$Lake))){
  x<-data.frame(Lake=l,Morphs=paste(sort(
    as.character(m.diet[m.diet$Lake==l,"Morph_diet"])),collapse="\n"))
  y<-rbind(y,x)
}
y<-y[-1,]
z.d<-y%>%group_by(Morphs)%>%mutate(counts=n())

z.d<-distinct(left_join(z.d,m,by="Lake")[,c("Lake","Morphs","counts","Genus","Species")])

B<-ggplot(z.d,aes(x=reorder(Morphs,-counts),fill=Species))+
  geom_bar(stat="count")+theme_classic()+
  scale_fill_manual(values=c(cols.cor,cols.sal))+
  theme(legend.position=c(.8,.7))+ylim(0,24)+
  xlab("Ecotypes present - Diet")+ylab("Count")

ggarrange(A,B,labels=c("A","B"))


z.z.fl<-subset(z.z,z.z$Lake%in%c(morphs.fl$Lake))
table(z.z.fl$Morphs,z.z.fl$Genus)


z.bp<-subset(z.d,z.d$Morphs=="ben\npla")
z.bp

GR.pb<-phenotype_repeat(morphs.GR,"GRc_mean","GRc_n","Cor",F)

A<-ggplot(data=GR.pb,
          aes(y=fct_reorder(Lake,Lacustrine_morphs,mean,.desc=T),
              x=trait,colour=Species,shape=Morph_diet))+
  geom_point(size=4)+
  theme_classic()+xlab("Gill Raker Count")+ylab("Lake")+
  scale_color_manual(values=cols.cor)+
  theme(panel.grid.major.y=element_line(size=.1,color="grey80"))#+

#repeatability in subset groups

morphs.fl.lpp<-morphs.fl[morphs.fl$Morph_zone%in%c("littoral/benthic","pelagic","profundal"),]
phenotype_repeat(morphs.fl.lpp,"Total_length","length_n","Cor",T)
phenotype_repeat(morphs.fl.lpp,"Total_length","length_n","Sal",T)

morphs.fl.bpp<-morphs.fl[morphs.fl$Morph_diet%in%c("planktivore","benthivore","piscivore"),]
phenotype_repeat(morphs.fl.bpp,"Total_length","length_n","Cor",T)
phenotype_repeat(morphs.fl.bpp,"Total_length","length_n","Sal",T)



morphs.GR.lpp<-morphs.GR[morphs.GR$Morph_zone%in%c("littoral/benthic","pelagic","profundal"),]
phenotype_repeat(morphs.GR.lpp,"GRc_mean","Cor",T)

morphs.fl.bpp<-morphs.fl[morphs.fl$Morph_diet%in%c("planktivore","benthivore","piscivore"),]
phenotype_repeat(morphs.fl.bpp,"Total_length","length_n","Cor",T)
phenotype_repeat(morphs.fl.bpp,"Total_length","length_n","Sal",T)




##################################################
#estimate repeatability as interaction relative to main effects
##################################################


#add error type column for total length
morphs.fl$Total_length_etype<-"Variance"
colnames(morphs)

#make dataframe for results
res.tr<-data.frame()
#cycle through each trait
for(tr.i in c("gill raker count","total body length","C13","N15","age")){
  
  #for each trait, assign correct column numbers for mean, sample size, error, error type
  if(tr.i=="gill raker count"){
    tr_cols<-6:9
    tr<-morphs[c(1:5,tr_cols)]}
  if(tr.i=="total body length"){ #use total length, based on adjustments
    tr_cols<-c(28,19,29,30)
    tr<-morphs.fl[c(1:5,tr_cols)]}
  if(tr.i=="C13"){
    tr_cols<-10:13
    tr<-morphs[c(1:5,tr_cols)]}
  if(tr.i=="N15"){
    tr_cols<-14:17
    tr<-morphs[c(1:5,tr_cols)]}
  if(tr.i=="age"){
    tr_cols<-24:27
    tr<-morphs[c(1:5,tr_cols)]}
  
  #cycle through ecological categories
  for(eco.id in c("diet","habitat")){
    
    #rename columns, set as numeric
    colnames(tr)[6:9]<-c("tr.mean","tr.n","tr.e","tr.etype")
    tr$tr.mean<-as.numeric(tr$tr.mean)
    tr$tr.n<-as.numeric(tr$tr.n)
    tr$tr.e<-as.numeric(tr$tr.e)
    
    tr<-tr[!is.na(tr$tr.n),] # remove if no sample size
    tr<-tr[!is.na(tr$tr.e),] # remove if no error
    
    #remove if morph sample size less than 5
    tr<-tr[tr$tr.n>=5,]
    tr<-tr[!is.na(tr$tr.etype),]
    
    #convert SD, SE and 95CI to Variance
    tr$tr.e[tr$tr.etype=="SD"]<-tr$tr.e[tr$tr.etype=="SD"]^2
    tr$tr.e[tr$tr.etype=="SE"]<-se_to_sd(tr$tr.e[tr$tr.etype=="SE"],
                                         tr$tr.n[tr$tr.etype=="SE"])^2
    tr$tr.e[tr$tr.etype=="95CI"]<-CI95_to_sd(tr$tr.e[tr$tr.etype=="95CI"],
                                             tr$tr.n[tr$tr.etype=="95CI"])^2
    tr$tr.etype<-"Variance"
    
    if(eco.id=="diet"){colnames(tr)[colnames(tr)=="Morph_diet"]<-"Morph_niche"}
    if(eco.id=="habitat"){colnames(tr)[colnames(tr)=="Morph_zone"]<-"Morph_niche"}
    
    tr<-tr[tr$Morph_niche!="",] #only keep values with associated niche label
    
    #cycle through each lake
    for(i in levels(as.factor(tr$Lake))){
      
      df.i<-subset(tr,tr$Lake==i)
      #only keep papers with all morphs included
      df.i<-subset(df.i,df.i$Paper_ID%in%
                     names(table(df.i$Paper_ID)==df[df$Lake==i,"Lacustrine_morphs"]))
      
      
      
      #calc weighted mean and variance for all records of an ecotype
      #weighted mean if (1) 2 populations of one ecotype in a lake
      #(2) 2 papers with records of an ecotype
      res.w<-data.frame()
      for(j in levels(as.factor(df.i$Morph_niche))){
        df.j<-subset(df.i,df.i$Morph_niche==j) #select all records for an ecotype
        tr.mean.w<-sum(df.j$tr.mean*df.j$tr.n)/sum(df.j$tr.n) #calculate weighted mean
        tr.e.w<-(1/(sum(df.j$tr.n)^2))*sum(df.j$tr.n^2*df.j$tr.e) #calculated weighted variance
        tr.n.w<-sum(df.j$tr.n)
        #add new line to dataframe with weighted mean, variance, etc
        res.w<-rbind(res.w,data.frame(Lake=i,Morph_niche=j,Genus=substr(df.i$Species,1,3),
                                      trait=tr.i,niche.id=eco.id,
                                      tr.mean.w,tr.e.w,tr.n.w))
      }
      res.tr<-rbind(res.tr,res.w)
    }}}


res.tr<-subset(res.tr,res.tr$Lake%in%names(table(res.tr$Lake)[table(res.tr$Lake)>1]))
res.tr<-distinct(res.tr)

#calculate coefficient of variation
res.tr$tr.cv<-sqrt(res.tr$tr.e.w)/res.tr$tr.mean.w

#get common ecotypes
res.tr<-subset(res.tr,res.tr$Morph_niche%in%
                 c("profundal","pelagic","littoral/benthic","planktivore","benthivore","piscivore"))


temp1<-distinct(res.tr[,c("trait","Genus","Lake")])

#get distinct ecotype by lake levels
temp2<-distinct(res.tr[,c("Morph_niche","Lake")])
#get lakes list, subset res.tr by 2+ morphs present
res.tr<-subset(res.tr,res.tr$Lake%in%names(table(temp2$Lake)[table(temp2$Lake)>1]))


#make a table for results output, sorted by genus and trait
out.res.tr<-res.tr%>%tidyr::expand(Genus,trait,niche.id)%>%
  #filter out trait by genus combos with low samples size
  filter(!(Genus=="Cor"&trait=="age")&!(Genus=="Sal"&trait=="gill raker count"))

#prep empty variables and lists
out.res.tr$D.mean<-0;out.res.tr$D.cv<-0
plots.mean.res.tr<-list();plots.cv.res.tr<-list();plots.res.tr<-list()

#fill in results output
for(i in 1:nrow(out.res.tr)){
  df.i<-res.tr%>%filter(Genus==out.res.tr$Genus[i]&trait==out.res.tr$trait[i]&
                          niche.id==out.res.tr$niche.id[i])
  
  #remove lakes with only one ecotype
  temp1<-data.frame(addmargins(table(df.i$Lake,df.i$Morph_niche)))
  df.i<-subset(df.i,df.i$Lake%in%temp1[temp1$Var2=="Sum"&temp1$Freq>1,"Var1"])
  
  #test of repeatability for means
  lm1<-lm(tr.mean.w~Morph_niche*Lake,weights=1/tr.e.w,data=df.i)
  aov1<-anova(lm1)
  #3 = interaction term, 1 = diet category
  D.mean.i<-aov1$`Sum Sq`[3]/(aov1$`Sum Sq`[1]+aov1$`Sum Sq`[3])
  out.res.tr$D.mean[i]<-D.mean.i
  
  #test of repeatability for cv
  lm2<-lm(tr.cv~Morph_niche*Lake,data=df.i)
  aov2<-anova(lm2)
  #3 = interaction term, 1 = diet category
  D.cv.i<-aov2$`Sum Sq`[3]/(aov2$`Sum Sq`[1]+aov2$`Sum Sq`[3])
  out.res.tr$D.cv[i]<-D.cv.i
  
  #get correct colour for the genus
  col.genus<-if_else(df.i$Genus[1]=="Cor",col.cor,col.sal)
  #make plot and add to list
  
  pl.mean<-ggplot(data=df.i,aes(x=Morph_niche,y=tr.mean.w,colour=Genus))+
    geom_point(size=2)+
    geom_line(aes(group=Lake),size=1)+
    ylab(paste("Mean",df.i$trait[1]))+
    xlab(paste("Ecotype -",df.i$niche.id[1]))+
    theme_classic()+scale_colour_manual(values=col.genus)+
    theme(legend.position="none")
  
  plots.mean.res.tr[[i]]<-pl.mean
  
  pl.cv<-ggplot(data=df.i,aes(x=Morph_niche,y=tr.cv,colour=Genus))+
    geom_point(size=2)+
    geom_line(aes(group=Lake),size=1)+
    ylab(paste("CV",df.i$trait[1]))+
    xlab(paste("Ecotype -",df.i$niche.id[1]))+
    theme_classic()+scale_colour_manual(values=col.genus)+
    theme(legend.position="none")
  
  plots.cv.res.tr[[i]]<-pl.cv
  
  plots.res.tr[[i]]<-ggarrange(pl.mean,pl.cv,labels=c("A","B"))
  
  
  
}

out.res.tr


#make supp mat figures

ggarrange(plots.mean.res.tr[[3]],plots.mean.res.tr[[4]],
          plots.mean.res.tr[[7]],plots.mean.res.tr[[8]],
          plots.mean.res.tr[[5]],plots.mean.res.tr[[6]],
          plots.mean.res.tr[[1]],plots.mean.res.tr[[2]],
          labels=LETTERS[1:8],ncol=2,nrow=4)

ggarrange(plots.mean.res.tr[[15]],plots.mean.res.tr[[16]],
          plots.mean.res.tr[[9]],plots.mean.res.tr[[10]],
          plots.mean.res.tr[[11]],plots.mean.res.tr[[12]],
          plots.mean.res.tr[[13]],plots.mean.res.tr[[14]],
          labels=LETTERS[1:8],ncol=2,nrow=4)

ggarrange(plots.cv.res.tr[[3]],plots.cv.res.tr[[4]],
          plots.cv.res.tr[[7]],plots.cv.res.tr[[8]],
          plots.cv.res.tr[[5]],plots.cv.res.tr[[6]],
          plots.cv.res.tr[[1]],plots.cv.res.tr[[2]],
          labels=LETTERS[1:8],ncol=2,nrow=4)

ggarrange(plots.cv.res.tr[[15]],plots.cv.res.tr[[16]],
          plots.cv.res.tr[[9]],plots.cv.res.tr[[10]],
          plots.cv.res.tr[[11]],plots.cv.res.tr[[12]],
          plots.cv.res.tr[[13]],plots.cv.res.tr[[14]],
          labels=LETTERS[1:8],ncol=2,nrow=4)


#fig 5



#select benthivore and planktivore
#only keep lakes that have both ecotypes
res.tl.sal.d.bpl<-subset(res.tr,res.tr$Morph_niche%in%c("benthivore","planktivore")&
                           res.tr$Genus=="Sal"&res.tr$trait=="total body length")
res.tl.sal.d.bpl<-subset(res.tl.sal.d.bpl,res.tl.sal.d.bpl$Lake%in%
                           data.frame(table(res.tl.sal.d.bpl$Lake))
                         [data.frame(table(res.tl.sal.d.bpl$Lake))[,"Freq"]>1,"Var1"])

A1<-ggplot(data=res.tl.sal.d.bpl,
           aes(x=Morph_niche,y=tr.mean.w,colour=Genus))+
  geom_point(size=2)+
  geom_line(aes(group=Lake),size=1)+
  ylab("Total Length")+xlab("Ecotype - Diet")+
  theme_classic()+scale_colour_manual(values=c(col.sal))+
  theme(legend.position="none")

res.tl.cor.d.2<-subset(res.tr,res.tr$Morph_niche%in%c("benthivore","planktivore")&
                         res.tr$Genus=="Cor"&res.tr$trait=="total body length")
res.tl.cor.d.2<-subset(res.tl.cor.d.2,res.tl.cor.d.2$Lake%in%
                         data.frame(table(res.tl.cor.d.2$Lake))
                       [data.frame(table(res.tl.cor.d.2$Lake))[,"Freq"]>1,"Var1"])


A2<-ggplot(data=res.tl.cor.d.2,
           aes(x=Morph_niche,y=tr.mean.w,colour=Genus))+
  geom_point(size=2)+
  geom_line(aes(group=Lake),size=1)+
  ylab("Total length")+xlab("Ecotype - Diet")+
  theme_classic()+scale_colour_manual(values=c(col.cor))+
  theme(legend.position="none")

res.tl.sal.d.bpi<-subset(res.tr,res.tr$Morph_niche%in%c("benthivore","piscivore")&
                           res.tr$Genus=="Sal"&res.tr$trait=="total body length")
res.tl.sal.d.bpi<-subset(res.tl.sal.d.bpi,res.tl.sal.d.bpi$Lake%in%
                           data.frame(table(res.tl.sal.d.bpi$Lake))
                         [data.frame(table(res.tl.sal.d.bpi$Lake))[,"Freq"]>1,"Var1"])

B<-ggplot(data=res.tl.sal.d.bpi,
          aes(x=Morph_niche,y=tr.mean.w,colour=Genus))+
  geom_point(size=2)+
  geom_line(aes(group=Lake),size=1)+
  ylab("Total Length")+xlab("Ecotype - Diet")+
  theme_classic()+scale_colour_manual(values=c(col.sal))+
  theme(legend.position="none")

res.tl.sal.d.pipl<-subset(res.tr,res.tr$Morph_niche%in%c("piscivore","planktivore")&
                            res.tr$Genus=="Sal"&res.tr$trait=="total body length")
res.tl.sal.d.pipl<-subset(res.tl.sal.d.pipl,res.tl.sal.d.pipl$Lake%in%
                            data.frame(table(res.tl.sal.d.pipl$Lake))
                          [data.frame(table(res.tl.sal.d.pipl$Lake))[,"Freq"]>1,"Var1"])

C<-ggplot(data=res.tl.sal.d.pipl,
          aes(x=Morph_niche,y=tr.mean.w,colour=Genus))+
  geom_point(size=2)+
  geom_line(aes(group=Lake),size=1)+
  ylab("Total Length")+xlab("Ecotype - Diet")+
  theme_classic()+scale_colour_manual(values=c(col.sal))+
  theme(legend.position="none")

morphs.fl.pr<-rbind(morphs.fl.pr.S,morphs.fl.pr.C)
morphs.fl.pr$Genus<-t(data.frame(strsplit(morphs.fl.pr$Species,"_")))[,1]
TL.plot<-ggplot(data=morphs.fl.pr,
                aes(y=Lake,x=trait,colour=Genus,shape=Genus))+
  geom_point(size=3,stroke=2)+
  scale_colour_manual(values=c(col.cor,col.sal))+
  scale_shape_manual(values=c(0,2))+
  theme_classic()+xlab("Total Length")+ylab("Lake")+
  theme(panel.grid.major.y=element_line(size=.1,color="grey80"))

ggarrange(TL.plot,
          ggarrange(A1,B,C,A2,
                    nrow=2,ncol=3,
                    labels=c("B","C","D","E")),
          common.legend=T,legend="right",
          ncol=2,widths=c(1,1.2),
          labels=c("A",""))



##Fig 4


res.gr.cor.h.bpe<-subset(res.tr,res.tr$Morph_niche%in%c("littoral/benthic","pelagic")&
                           res.tr$Genus=="Cor"&res.tr$trait=="gill raker count")
res.gr.cor.h.bpe<-subset(res.gr.cor.h.bpe,res.gr.cor.h.bpe$Lake%in%
                           data.frame(table(res.gr.cor.h.bpe$Lake))
                         [data.frame(table(res.gr.cor.h.bpe$Lake))[,"Freq"]>1,"Var1"])

I<-ggplot(data=res.gr.cor.h.bpe,
          aes(x=Morph_niche,y=tr.mean.w,colour=Genus))+
  geom_point(size=2)+
  geom_line(aes(group=Lake),size=1)+
  ylab("Gill raker count")+xlab("Ecotype - Habitat")+
  theme_classic()+scale_colour_manual(values=c(col.cor,col.sal))+
  theme(legend.position="none")


res.gr.cor.h.bpr<-subset(res.tr,res.tr$Morph_niche%in%c("littoral/benthic","profundal")&
                           res.tr$Genus=="Cor"&res.tr$trait=="gill raker count")
res.gr.cor.h.bpr<-subset(res.gr.cor.h.bpr,res.gr.cor.h.bpr$Lake%in%
                           data.frame(table(res.gr.cor.h.bpr$Lake))
                         [data.frame(table(res.gr.cor.h.bpr$Lake))[,"Freq"]>1,"Var1"])

J<-ggplot(data=res.gr.cor.h.bpr,
          aes(x=Morph_niche,y=tr.mean.w,colour=Genus))+
  geom_point(size=2)+
  geom_line(aes(group=Lake),size=1)+
  ylab("Gill raker count")+xlab("Ecotype - Habitat")+
  theme_classic()+scale_colour_manual(values=c(col.cor,col.sal))+
  theme(legend.position="none")


res.gr.cor.h.pepr<-subset(res.tr,res.tr$Morph_niche%in%c("pelagic","profundal")&
                            res.tr$Genus=="Cor"&res.tr$trait=="gill raker count")
res.gr.cor.h.pepr<-subset(res.gr.cor.h.pepr,res.gr.cor.h.pepr$Lake%in%
                            data.frame(table(res.gr.cor.h.pepr$Lake))
                          [data.frame(table(res.gr.cor.h.pepr$Lake))[,"Freq"]>1,"Var1"])

K<-ggplot(data=res.gr.cor.h.pepr,
          aes(x=Morph_niche,y=tr.mean.w,colour=Genus,sp))+
  geom_point(size=2)+
  geom_line(aes(group=Lake),size=1)+
  ylab("Gill raker count")+xlab("Ecotype - Habitat")+
  theme_classic()+scale_colour_manual(values=c(col.cor,col.sal))+
  theme(legend.position="none")



morphs.gr.pr<-rbind(morphs.GR.C.pr,morphs.GR.S.pr)
morphs.gr.pr$Genus<-t(data.frame(strsplit(morphs.gr.pr$Species,"_")))[,1]
GR.plot<-ggplot(data=morphs.gr.pr,
                aes(y=Lake,x=trait,colour=Genus,shape=Genus))+
  geom_point(size=3,stroke=2)+
  scale_colour_manual(values=c(col.cor,col.sal))+
  scale_shape_manual(values=c(0,2))+
  theme_classic()+xlab("Gill Raker Count")+ylab("Lake")+
  theme(panel.grid.major.y=element_line(size=.1,color="grey80"))



ggarrange(GR.plot,
          ggarrange(I,J,K,
                    nrow=2,ncol=2,
                    labels=c("B","C","D")),
          common.legend=T,legend="right",
          ncol=2,widths=c(1,1),
          labels=c("A",""))


