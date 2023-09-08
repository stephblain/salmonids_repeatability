setwd("C:/Users/Steph/OneDrive - Texas A&M University/Salmonids/data/")

###Read in packages
pkgs<-c("tidyverse","ggplot2","coxed","gridExtra","matrixcalc","forcats",
        "vegan","rgeos","rnaturalearth", "rnaturalearthdata","sf","bipartite",
        "kSamples","ggpubr","metaDigitise","metafor","ape","rfishbase")

lapply(pkgs,library,character.only=T); rm(pkgs)
source("../salmonids/functions_phenotypes.R")

###Colours
col.cor<-"#8cb5a7"
col.sal<-"#b58c93"
cols.sal<-c("#241c1d","#5b464a","#917076","#c4a3a9")
cols.cor<-c("#2a3632","#546d64","#709186","#a3c4b9")
ecotype_cols<-c("#5C697B","#69788C","#76879E","#8496B0","#8496B0","#90A0B7","#9CABBF")


###Read in data

df<-read.csv("lakes_summary.csv",stringsAsFactors=F)
m<-read.csv("morphs_summary.csv",stringsAsFactors=F)
morphs<-read.csv("morphs.csv",stringsAsFactors=F,fileEncoding="latin1")
phylo<-read.tree("CreteLafreniere2012_pruned.tre")

#take log of area and depth
df$logArea<-log(df$Area)
df$logDepthMax<-log(df$Depth_max)

#remove data not from Coregonus or Salvelinus
rmLake<-m$Lake[m$Species=="Prosopium_coulterii"][1]
m<-m[m$Lake!=rmLake,];morphs<-morphs[morphs$Lake!=rmLake,];df<-df[df$Lake!=rmLake,]
#remove useless columns
df$X.1<-NULL; m$X<-NULL; morphs$X.1<-NULL


#output wd for figures
setwd("C:/Users/Steph/OneDrive - Texas A&M University/Salmonids/figures/morphs")


df$morphs_level<-df$Lacustrine_morphs
df$morphs_level[df$morphs_level>3]<-"4+"
df$morphs_level<-factor(df$morphs_level,levels=c("2","3","4+"))


##############plot niche categories - supp mat
m$Genus<-sub("_.*","",m$Species)
m$Morph_diet<-factor(m$Morph_diet,levels=c("benthivore","planktivore","piscivore",
                                           "generalist invertivore","omnivore"))
h.d<-ggplot(data=subset(m,!is.na(m$Morph_diet)),aes(x=Morph_diet,fill=Species))+
  geom_bar()+
  scale_fill_manual(values=c(cols.cor,cols.sal))+
  theme_classic()+xlab("Ecotype - Diet")+ylab("Number of Ecotypes")+
  theme(legend.position="none")+
  ylim(0,68)

m$Morph_zone<-factor(m$Morph_zone,levels=c("littoral/benthic","pelagic",
                                           "profundal","generalist","shallow"))
h.z<-ggplot(data=subset(m,!is.na(m$Morph_zone)),aes(x=Morph_zone,fill=Species))+
  geom_bar()+  scale_fill_manual(values=c(cols.cor,cols.sal))+
  theme_classic()+xlab("Ecotype - Habitat")+ylab("Number of Ecotypes")+
  theme(legend.position=c(0.8,0.7))+
  ylim(0,68)

ggarrange(h.d,h.z,ncol=2,labels=c("A","B"))

####Make dataframes for graphing

#Morph ID diet

m.d<-m[m$Morph_diet!="",c(1:3,6)]

t_m.d<-data.frame(addmargins(table(m.d$Lake,m.d$Morph_diet)))
colnames(t_m.d)<-c("Lake","Ecotype","Freq")
temp<-t_m.d[t_m.d$Ecotype=="Sum",]
temp<-merge(temp,df,by="Lake")[,c(1,2,3,16)]
#get vector of lakes where number of morphs with ID is the number in the lake
temp<-temp[temp$Freq==temp$Lacustrine_morphs,]
t_m.d<-t_m.d[t_m.d$Lake%in%temp$Lake,]
t_m.d<-t_m.d[t_m.d$Ecotype!="Sum",]
temp<-temp[order(temp$Lacustrine_morphs),]

#reorder factor levels for Lake by number of morphs
t_m.d$Lake<-factor(t_m.d$Lake,levels=temp$Lake)
t_m.d$Lake<-fct_rev(t_m.d$Lake)
t_m.d<-merge(t_m.d,temp,by="Lake")[,c(1:3,6)]
t_m.d<-merge(t_m.d,df,by="Lake")[,c(1:4,8:9,12,13)]

colnames(t_m.d)[1:4]<-c("Lake","Ecotype","Freq","Lacustrine_morphs")


t_m.d.p<-t_m.d
t_m.d.p$Ecotype<-factor(t_m.d.p$Ecotype,levels=c("benthivore","planktivore","piscivore",
                                                 "generalist.invertivore",
                                                 "generalist invertivore","omnivore"))

temp<-data.frame(table(m$Lake,m$Species))
temp<-temp[temp$Freq!=0,]
colnames(temp)<-c("Lake","Species","Freq2")
t_m.d.p<-merge(t_m.d.p,temp,by="Lake")[,1:9]
t_m.d.p$Ecotype[t_m.d.p$Ecotype=="generalist invertivore"]<-"generalist.invertivore"

#Get factor order for Lakes, not rearranging dataframe
t_m.d.wide<-pivot_wider(t_m.d.p,names_from=Ecotype,values_from=Freq)
t_m.d.wide<-arrange(t_m.d.wide,Lacustrine_morphs,desc(benthivore),desc(planktivore),
                    desc(piscivore),desc(generalist.invertivore),desc(omnivore))

t_m.d.p$Lake<-factor(t_m.d.p$Lake,levels=rev(t_m.d.wide$Lake))


#convert axis labels to two lines
levels(t_m.d.p$Ecotype)<-gsub("st.in","st\nin",levels(t_m.d.p$Ecotype))
t_m.d.p$Freq<-as.character(t_m.d.p$Freq)
t_m.d.p$Freq<-as.character(t_m.d.p$Freq)
t_m.d.p$Freq<-factor(t_m.d.p$Freq,levels=c(0,1,2,3,4))


t_m.d.p$Genus<-substr(t_m.d.p$Species,1,3)
t_m.d.cor<-subset(t_m.d.p,t_m.d.p$Genus%in%c("Cor","Pro"))
t_m.d.sal<-subset(t_m.d.p,t_m.d.p$Genus=="Sal")

t_m.d.cor.2<-t_m.d.cor[t_m.d.cor$Lacustrine_morphs==2,]
t_m.d.sal.2<-t_m.d.sal[t_m.d.sal$Lacustrine_morphs==2,]



#Morph ID zone

m.z<-m[m$Morph_zone!="",c(1:3,5)]


t_m.z<-data.frame(addmargins(table(m.z$Lake,m.z$Morph_zone)))
colnames(t_m.z)<-c("Lake","Ecotype","Freq")
temp<-t_m.z[t_m.z$Ecotype=="Sum",]
temp<-merge(temp,df,by="Lake")[,c(1,2,3,16)]
#get vector of lakes where number of morphs with ID is the number in the lake
temp<-temp[temp$Freq==temp$Lacustrine_morphs,]
t_m.z<-t_m.z[t_m.z$Lake%in%temp$Lake,]
t_m.z<-t_m.z[t_m.z$Ecotype!="Sum",]
temp<-temp[order(temp$Lacustrine_morphs),]

#reorder factor levels for Lake by number of morphs
t_m.z$Lake<-factor(t_m.z$Lake,levels=temp$Lake)
t_m.z$Lake<-fct_rev(t_m.z$Lake)
t_m.z<-merge(t_m.z,temp,by="Lake")[,c(1:3,6)]
t_m.z<-merge(t_m.z,df,by="Lake")[,c(1:4,8:9,12,13)]

colnames(t_m.z)[1:4]<-c("Lake","Ecotype","Freq","Lacustrine_morphs")


t_m.z.p<-t_m.z
t_m.z.p$Ecotype<-factor(t_m.z.p$Ecotype,levels=c("littoral.benthic","littoral/benthic",
                                                 "pelagic","profundal","generalist","shallow"))
t_m.z.p$Ecotype[t_m.z.p$Ecotype=="littoral/benthic"]<-"littoral.benthic"


temp<-data.frame(table(m$Lake,m$Species))
temp<-temp[temp$Freq!=0,]
colnames(temp)<-c("Lake","Species","Freq2")
t_m.z.p<-merge(t_m.z.p,temp,by="Lake")[,1:9]


#Get factor order for Lakes, not rearranging dataframe
t_m.z.wide<-pivot_wider(t_m.z.p,names_from=Ecotype,values_from=Freq)
t_m.z.wide<-arrange(t_m.z.wide,Lacustrine_morphs,desc(littoral.benthic),desc(pelagic),
                    desc(profundal),desc(generalist),desc(shallow))

t_m.z.p$Lake<-factor(t_m.z.p$Lake,levels=rev(t_m.z.wide$Lake))


#convert axis labels to two lines
levels(t_m.z.p$Ecotype)<-gsub("l.b","l/\nb",levels(t_m.z.p$Ecotype))
t_m.d.p$Freq<-as.character(t_m.d.p$Freq)
t_m.z.p$Freq<-as.character(t_m.z.p$Freq)
t_m.z.p$Freq<-factor(t_m.z.p$Freq,levels=c(0,1,2,3,4))


t_m.z.p$Genus<-substr(t_m.z.p$Species,1,3)
t_m.z.cor<-subset(t_m.z.p,t_m.z.p$Genus%in%c("Cor","Pro"))
t_m.z.sal<-subset(t_m.z.p,t_m.z.p$Genus=="Sal")

t_m.z.cor.2<-t_m.z.cor[t_m.z.cor$Lacustrine_morphs==2,]
t_m.z.sal.2<-t_m.z.sal[t_m.z.sal$Lacustrine_morphs==2,]

#lakes list for habitat or diet IDs
lakes_with_ID<-unique(c(levels(t_m.z$Lake),levels(t_m.d$Lake)))

cols.tile<-c("white",col.sal,"#8A6E73","#5E5153","#333333")

a1<-ggplot(data=t_m.d.sal,aes(y=Lake,x=Ecotype,fill=Freq))+
  geom_tile(size=1,colour="white")+theme_classic()+
  scale_fill_manual(values=cols.tile)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Ecotype - Diet")
b1<-ggplot(data=t_m.z.sal,aes(y=Lake,x=Ecotype,fill=Freq))+
  geom_tile(size=1,colour="white")+theme_classic()+
  scale_fill_manual(values=cols.tile,drop=F)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Ecotype - Habitat")


#two ecotypes
c1<-ggplot(data=t_m.d.sal.2,aes(y=Lake,x=Ecotype,fill=Freq))+
  geom_tile(size=0.5,colour="white")+theme_classic()+
  scale_fill_manual(values=cols.tile)+
  theme(legend.position="none")+
  xlab("Ecotype - Diet")
d1<-ggplot(data=t_m.z.sal.2,aes(y=Lake,x=Ecotype,fill=Freq))+
  geom_tile(size=0.5,colour="white")+theme_classic()+
  scale_fill_manual(values=cols.tile)+
  theme(legend.position=c(0.9,0.85))+
  xlab("Ecotype - Habitat")


t_m.d.sal.pa<-t_m.d.sal
t_m.d.sal.pa$Freq[t_m.d.sal.pa$Freq%in%2:4]<-1
e1<-ggplot(data=t_m.d.sal.pa,aes(y=Lake,x=Ecotype,fill=Freq))+
  geom_tile(size=0.5,colour="white")+theme_classic()+
  scale_fill_manual(values=c("white",col.sal))+
  theme(legend.position="none")+
  xlab("Habitat")
t_m.z.sal.pa<-t_m.z.sal
t_m.z.sal.pa$Freq[t_m.z.sal.pa$Freq%in%2:4]<-1
f1<-ggplot(data=t_m.z.sal.pa,aes(y=Lake,x=Ecotype,fill=Freq))+
  geom_tile(size=0.5,colour="white")+theme_classic()+
  scale_fill_manual(values=c("white",col.sal))+
  theme(legend.position="none")+
  xlab("Habitat")



t_m.d.cor$Freq<-as.character(t_m.d.cor$Freq)
t_m.z.cor$Freq<-as.character(t_m.z.cor$Freq)
cols.tile.cor<-c("white",col.cor,"#6E8A80","#515E5A","#333333")

a2<-ggplot(data=t_m.d.cor,aes(y=Lake,x=Ecotype,fill=Freq))+
  geom_tile(size=1,colour="white")+theme_classic()+
  scale_fill_manual(values=cols.tile.cor)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Ecotype - Diet")

b2<-ggplot(data=t_m.z.cor,aes(y=Lake,x=Ecotype,fill=Freq))+
  geom_tile(size=1,colour="white")+theme_classic()+
  scale_fill_manual(values=cols.tile.cor)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Ecotype - Habitat")

#Figure 2 and Figure 3
fig2<-ggarrange(a2,a1,ncol=2,nrow=1,labels=c("A","B"),align="v")
fig3<-ggarrange(b2,b1,ncol=2,nrow=1,labels=c("A","B"),align="v")

# ggsave("figure2.png",plot=fig2,dpi=300,width=22,height=14,units="cm")
# ggsave("figure3.png",plot=fig3,dpi=300,width=22,height=14,units="cm")



t_m.d.cor.pa<-t_m.d.cor
t_m.d.cor.pa$Freq[t_m.d.cor.pa$Freq%in%2:4]<-1
e2<-ggplot(data=t_m.d.cor.pa,aes(y=Lake,x=Ecotype,fill=Freq))+
  geom_tile(size=0.5,colour="white")+theme_classic()+
  scale_fill_manual(values=c("white",col.cor))+
  theme(legend.position="none")+
  xlab("Habitat")
t_m.z.cor.pa<-t_m.z.cor
t_m.z.cor.pa$Freq[t_m.z.cor.pa$Freq%in%2:4]<-1
f2<-ggplot(data=t_m.z.cor.pa,aes(y=Lake,x=Ecotype,fill=Freq))+
  geom_tile(size=0.5,colour="white")+theme_classic()+
  scale_fill_manual(values=c("white",col.cor))+
  theme(legend.position="none")+
  xlab("Habitat")

#presence/absence of ecotypes - supp mat
ggarrange(e2,e1,ncol=2,nrow=1,labels=c("A","B"),align="v")
ggarrange(f2,f1,ncol=2,nrow=1,labels=c("A","B"),align="v")


t_m.d.cor.2$Freq<-as.character(t_m.d.cor.2$Freq)
t_m.z.cor.2$Freq<-as.character(t_m.z.cor.2$Freq)

c2<-ggplot(data=t_m.d.cor.2,aes(y=Lake,x=Ecotype,fill=Freq))+
  geom_tile(size=0.5,colour="white")+theme_classic()+
  scale_fill_manual(values=cols.tile.cor)+
  theme(legend.position="none")+
  xlab("Ecotype - Diet")

d2<-ggplot(data=t_m.z.cor.2,aes(y=Lake,x=Ecotype,fill=Freq))+
  geom_tile(size=0.5,colour="white")+theme_classic()+
  scale_fill_manual(values=cols.tile.cor)+
  theme(legend.position=c(0.9,0.85))+
  xlab("Ecotype - Habitat")

#two ecotypes only - supp mat
ggarrange(c1,d1,c2,d2,ncol=2,nrow=2,labels=c("A","B","C","D"),align="v")

#clean up environment
rm(t_m.d.cor,t_m.d.cor.2,t_m.d.cor.pa,
   t_m.z.cor,t_m.z.cor.2,t_m.z.cor.pa,
   t_m.d.sal,t_m.d.sal.2,t_m.d.sal.pa,
   t_m.z.sal,t_m.z.sal.2,t_m.z.sal.pa,
   t_m.d.p,t_m.z.p,m.d,m.z,
   a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,fig2,fig3,temp)


###Make dataframes for analyses

#Make diet dataframe
mat.d<-m[m$Morph_diet!="",c(1:3,6)] #get morph ID info
mat.d<-mat.d[mat.d$Lake%in%t_m.d$Lake,] #only keep lakes where everything is IDed
mat.d<-as.data.frame(table(mat.d$Morph_diet,mat.d$Lake),stringsAsFactors=F) #get frequency info
mat.d<-spread(mat.d,Var2,Freq) #turn into contingency table

mat.d<-data.frame(t(mat.d)) #transpose
colnames(mat.d)<-mat.d[1,]; mat.d<-mat.d[-1,]
mat.d$Lake<-rownames(mat.d); rownames(mat.d)<-NULL #make "Lake" column
mat.d<-distinct(dplyr::left_join(mat.d,m,by="Lake")[,1:7]) #join mat.d and m to get species names
mat.d<-data.frame(t(mat.d)) #transpose again
colnames(mat.d)<-mat.d[6,]; mat.d<-mat.d[-6,] #make lake names colnames and remove lake names row
mat.d[6,]<-substr(mat.d[6,],1,3); mat.d[6,1]<-"Genus" #get genus names and rename row
rownames(mat.d)<-mat.d$Var1; mat.d$Var1<-NULL #make ecotype names rownames and remove ecotype column

mat.d.cor<-mat.d[1:5,mat.d[6,]%in%c(NA,"Pro","Cor")] #Coregonus only and remove Genus row
mat.d.sal<-mat.d[1:5,mat.d[6,]%in%c(NA,"Sal")] #Salvelinus only and remove Genus row

mat.d.cor[,1:ncol(mat.d.cor)]<-lapply(mat.d.cor[,1:ncol(mat.d.cor)],as.numeric) #make dataframe numeric
mat.d.cor<-mat.d.cor[rowSums(mat.d.cor)>0,] #remove omnivore and piscivore, because none in Coregonus
mat.d.cor.2<-mat.d.cor[,colSums(mat.d.cor)==2] #only keep 2 morph lakes
mat.d.cor.2<-mat.d.cor.2[rowSums(mat.d.cor.2)>0,]

mat.d.cor.3<-mat.d.cor[,colSums(mat.d.cor)==3] #only keep 3 morph lakes
mat.d.cor.3<-mat.d.cor.3[rowSums(mat.d.cor.3)>0,]

#only 3 examples for 4 morphs

mat.d.sal[,1:ncol(mat.d.sal)]<-lapply(mat.d.sal[,1:ncol(mat.d.sal)],as.numeric) #make dataframe numeric
mat.d.sal.2<-mat.d.sal[,colSums(mat.d.sal)==2] #only keep 2 morph lakes
mat.d.sal.3<-mat.d.sal[,colSums(mat.d.sal)==3]



#Make habitat dataframe
mat.z<-m[m$Morph_zone!="",c(1:3,5)] #morph ID info
mat.z<-mat.z[mat.z$Lake%in%t_m.z$Lake,] #only keep lakes where everything is IDed
mat.z<-as.data.frame(table(mat.z$Morph_zone,mat.z$Lake),stringsAsFactors=F)
mat.z<-spread(mat.z,Var2,Freq) #make contingency table

mat.z<-data.frame(t(mat.z))
colnames(mat.z)<-mat.z[1,]; mat.z<-mat.z[-1,]
mat.z$Lake<-rownames(mat.z) #set lake names as rownames
rownames(mat.z)<-NULL

mat.z<-distinct(left_join(mat.z,m,by="Lake")[,1:7]) #join with m to get species names
mat.z<-data.frame(t(mat.z))
colnames(mat.z)<-mat.z[6,] #make lake names colnames
mat.z<-mat.z[-6,] #remove lake names row
mat.z[6,]<-substr(mat.z[6,],1,3) #get genus names
mat.z[6,1]<-"Genus"
rownames(mat.z)<-mat.z$Var1; mat.z$Var1<-NULL
#seperate by Genus
mat.z.cor<-mat.z[1:5,mat.z[6,]%in%c(NA,"Pro","Cor")]
mat.z.sal<-mat.z[1:5,mat.z[6,]%in%c(NA,"Sal")]
#convert to numeric
mat.z.cor[,1:ncol(mat.z.cor)]<-lapply(mat.z.cor[,1:ncol(mat.z.cor)],as.numeric)
mat.z.cor.2<-mat.z.cor[,colSums(mat.z.cor)==2] #make two species matrix
mat.z.cor.2<-mat.z.cor.2[rowSums(mat.z.cor.2)>0,] #remove ecotypes with no observations
mat.z.cor.3<-mat.z.cor[,colSums(mat.z.cor)==3]
mat.z.cor.3<-mat.z.cor.3[rowSums(mat.z.cor.3)>0,]
#convert to numeric
mat.z.sal[,1:ncol(mat.z.sal)]<-lapply(mat.z.sal[,1:ncol(mat.z.sal)],as.numeric)
mat.z.sal.2<-mat.z.sal[,colSums(mat.z.sal)==2] #make two species matrix
mat.z.sal.3<-mat.z.sal[,colSums(mat.z.sal)==3]
mat.z.sal.3<-mat.z.sal.3[rowSums(mat.z.sal.3)>0,]
#8 dataframes, Salvelinus + Coregonus, diet + habitat, 2 morph + all morph


#Looking for less deviation from expected than occurs by chance
#for that, need row and column sums to differ less from expected than would occur by chance


simulate.chisq(mat.d.cor.2)
simulate.chisq(mat.d.cor.3)
simulate.chisq(mat.d.cor)

simulate.chisq(mat.d.sal.2)
simulate.chisq(mat.d.sal.3)
simulate.chisq(mat.d.sal)

simulate.chisq(mat.z.cor.2)
simulate.chisq(mat.z.cor.3)
simulate.chisq(mat.z.cor)
simulate.chisq(mat.z.sal.2)
simulate.chisq(mat.z.sal.3)
simulate.chisq(mat.z.sal)

#list of lakes included in niche repeatability analysis
lakes.list.eco<-unique(c(colnames(mat.z.cor),colnames(mat.z.sal),colnames(mat.d.cor),colnames(mat.d.sal)))

#random tables have many more "2" columns than expected

#Test for association between different morphs using presence-absence matrices



pr.ab.tests(mat.d.cor)
pr.ab.tests(mat.z.cor)
pr.ab.tests(mat.d.sal)
pr.ab.tests(mat.z.sal)


#repeatability of phenotypes

morphs$GRc_e<-as.numeric(morphs$GRc_e); morphs$GRc_n<-as.numeric(morphs$GRc_n)
morphs$GRc_mean<-as.numeric(morphs$GRc_mean)
morphs.GR<-morphs[!is.na(morphs$GRc_n)&!is.na(morphs$GRc_e)&!is.na(morphs$GRc_mean),]


#phenotype repeatability - gill rakers
phenotype_repeat(morphs.GR,"GRc_mean","GRc_n","Cor",T)

#convert fork etc lengths to total length with fish base conversion
morphs$length_e<-as.numeric(morphs$length_e); morphs$length_n<-as.numeric(morphs$length_n)
morphs$length_mean<-as.numeric(morphs$length_mean)
morphs$length_mean[morphs$length_units=="cm"]<-morphs$length_mean[morphs$length_units=="cm"]*10
morphs$length_e[morphs$length_units=="cm"]<-morphs$length_e[morphs$length_units=="cm"]*10
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

###Fix typo
morphs.fl[morphs.fl$Lake=="Muddus",c("Total_length","Paper_ID")]
morphs.fl[morphs.fl$Lake=="Muddus"&morphs.fl$Paper_ID%in%c("40B","55B"),
          c("Total_length")]<-
  morphs.fl[morphs.fl$Lake=="Muddus"&morphs.fl$Paper_ID%in%c("40B","55B"),
            c("Total_length")]*10


phenotype_repeat(morphs.fl,"Total_length","length_n","Sal",T)
phenotype_repeat(morphs.fl,"Total_length","length_n","Cor",T)

morphs.fl.pr.S<-phenotype_repeat(morphs.fl,"Total_length","length_n","Sal",F)
morphs.fl.pr.S$Lake<-factor(morphs.fl.pr.S$Lake,
                            levels=rev(unique(arrange(morphs.fl.pr.S,Lacustrine_morphs,trait)$Lake)))
morphs.fl.pr.C<-phenotype_repeat(morphs.fl,"Total_length","length_n","Cor",F)
morphs.fl.pr.C$Lake<-factor(morphs.fl.pr.C$Lake,
                            levels=rev(unique(arrange(morphs.fl.pr.C,Lacustrine_morphs,trait)$Lake)))


A<-ggplot(data=morphs.fl.pr.S,
          aes(y=Lake,x=trait,shape=Species))+
  geom_point(size=4,colour=col.sal,stroke=2)+
  scale_shape_manual(values=c(0,1,2,5))+
  theme_classic()+xlab("Total Length")+ylab("Lake")+
  theme(panel.grid.major.y=element_line(size=.1,color="grey80"))
A
B<-ggplot(data=morphs.fl.pr.C,
          aes(y=Lake,x=trait,shape=Species))+
  geom_point(size=4,color=col.cor,stroke=2)+
  scale_shape_manual(values=c(0,1,2,5))+
  theme_classic()+xlab("Total Length")+ylab("Lake")+
  
  theme(panel.grid.major.y=element_line(size=.1,color="grey80"))
B




ggarrange(A,B,ncol=2,labels=c("A","B"))

morphs.GR.S.pr<-phenotype_repeat(morphs.GR,"GRc_mean","GRc_n","Sal",F)
morphs.GR.S.pr$Lake<-factor(morphs.GR.S.pr$Lake,
                            levels=rev(unique(arrange(morphs.GR.S.pr,Lacustrine_morphs,trait)$Lake)))
morphs.GR.C.pr<-phenotype_repeat(morphs.GR,"GRc_mean","GRc_n","Cor",F)
morphs.GR.C.pr$Lake<-factor(morphs.GR.C.pr$Lake,
                            levels=rev(unique(arrange(morphs.GR.C.pr,Lacustrine_morphs,trait)$Lake)))



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
setdiff(lakes.list.trait,lakes.list.eco)

#start habitat and diet histograms

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
A

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
GR.pb<-subset(GR.pb,GR.pb$Lake%in%z.bp$Lake)
GR.pb<-GR.pb[GR.pb$Lake!="Nipigon",]
GR.pb

temp1<-as.data.frame(pivot_wider(GR.pb,names_from=Lake,values_from=all_of("GRc_mean"))) #means in columns with each lake as a column name
ad.test(temp1[,6:ncol(temp1)])

A<-ggplot(data=GR.pb,
          aes(y=fct_reorder(Lake,Lacustrine_morphs,mean,.desc=T),
              x=GRc_mean,colour=Species,shape=Morph_diet))+
  geom_point(size=4)+
  theme_classic()+xlab("Gill Raker Count")+ylab("Lake")+
  scale_color_manual(values=cols.cor)+
  theme(panel.grid.major.y=element_line(size=.1,color="grey80"))#+
#theme(legend.position=c(.9,.95))+
#ggtitle("Coregonus")


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






z.z.lp<-subset(z.z,z.z$Morphs=="lit\npel")
morphs.GR$GRc_mean
GR.lp<-phenotype_repeat(morphs.GR,"GRc_mean","Cor",F)
GR.lp<-subset(GR.lp,GR.lp$Lake%in%z.z.lp$Lake)[,1:7]
temp1<-as.data.frame(pivot_wider(GR.lp,names_from=Lake,values_from=all_of("trait"))) #means in columns with each lake as a column name
ad.test(temp1[,6:ncol(temp1)])

# GR.lp<-phenotype_repeat(morphs.GR,"GRc_mean","Sal",F)
# GR.lp<-subset(GR.lp,GR.lp$Lake%in%z.z.lp$Lake)
# temp1<-as.data.frame(pivot_wider(GR.lp,names_from=Lake,values_from=all_of("GRc_mean"))) #means in columns with each lake as a column name
# ad.test(temp1[,6:ncol(temp1)])


FL.lp<-phenotype_repeat(morphs.FL,"length_mean","Sal",F)
FL.lp<-subset(FL.lp,FL.lp$Lake%in%z.z.lp$Lake)[,1:7]
temp1<-as.data.frame(pivot_wider(FL.lp,names_from=Lake,values_from=all_of("trait"))) #means in columns with each lake as a column name
ad.test(temp1[,6:ncol(temp1)])


TL.lp<-phenotype_repeat(morphs.TL,"length_mean","Cor",F)
TL.lp<-subset(TL.lp,TL.lp$Lake%in%z.z.lp$Lake)
temp1<-as.data.frame(pivot_wider(TL.lp,names_from=Lake,values_from=all_of("length_mean"))) #means in columns with each lake as a column name
ad.test(temp1[,6:ncol(temp1)])



B<-ggplot(data=GR.lp,
          aes(y=fct_reorder(Lake,Lacustrine_morphs,mean,.desc=T),
              x=GRc_mean,colour=Species,shape=Morph_zone))+
  geom_point(size=4)+
  theme_classic()+xlab("Gill Raker Count")+ylab("Lake")+
  scale_color_manual(values=cols.cor)+
  theme(panel.grid.major.y=element_line(size=.1,color="grey80"))#+
ggarrange(A,B,labels=c("A","B"))

FL.lp<-phenotype_repeat(morphs.FL,"length_mean","Sal",F)
FL.lp<-subset(FL.lp,FL.lp$Lake%in%z.z.lp$Lake)
FL.lp<-FL.lp[FL.lp$Lake!="Nipigon",]
#GR.pb<-GR.pb[GR.pb$Species=="Coregonus_lavaretus",]
temp1<-as.data.frame(pivot_wider(GR.pb,names_from=Lake,values_from=all_of("GRc_mean"))) #means in columns with each lake as a column name
ad.test(temp1[,6:ncol(temp1)])

ggplot(data=FL.lp,
       aes(y=fct_reorder(Lake,Lacustrine_morphs,mean,.desc=T),
           x=length_mean,colour=Species,shape=Morph_zone))+
  geom_point(size=4)+
  theme_classic()+xlab("Fork Length")+ylab("Lake")+
  scale_color_manual(values=cols.cor)+
  theme(panel.grid.major.y=element_line(size=.1,color="grey80"))#+




FL.pb<-phenotype_repeat(morphs.FL,"length_mean","Sal",F)
FL.pb<-subset(FL.pb,FL.pb$Lake%in%z.bp$Lake)
#GR.pb<-GR.pb[GR.pb$Species=="Coregonus_lavaretus",]
trait_col<-"length_mean"
temp1<-as.data.frame(pivot_wider(FL.pb,names_from=Lake,values_from=all_of("length_mean"))) #means in columns with each lake as a column name
ad.test(temp1[,6:ncol(temp1)])

temp<-data.frame(table(FL.pb$Lake))
FL.pb<-FL.pb[FL.pb$Lake%in%temp[temp$Freq>1,"Var1"],]

A<-ggplot(data=FL.pb,
          aes(y=fct_reorder(Lake,Lacustrine_morphs,mean,.desc=T),
              x=length_mean,colour=Species,shape=Morph_diet))+
  geom_point(size=4)+
  theme_classic()+xlab("Fork Length")+ylab("Lake")+
  scale_color_manual(values=cols.sal)+
  scale_shape_manual(values=c(17,15))+
  theme(panel.grid.major.y=element_line(size=.1,color="grey80"))+
  theme(legend.position="none")
nlevels(as.factor(morphs$Paper_ID))

FL.pb<-phenotype_repeat(morphs.fl,"length_mean","length_n","Sal",F)
FL.pb<-subset(FL.pb,FL.pb$Morph_diet%in%c("benthivore","planktivore"))
#GR.pb<-GR.pb[GR.pb$Species=="Coregonus_lavaretus",]
trait_col<-"length_mean"
temp1<-as.data.frame(pivot_wider(FL.pb,names_from=Lake,values_from=all_of("trait"))) #means in columns with each lake as a column name
ad.test(temp1[,7:ncol(temp1)])

temp<-data.frame(table(FL.pb$Lake))
FL.pb<-FL.pb[FL.pb$Lake%in%temp[temp$Freq>1,"Var1"],]

B<-ggplot(data=FL.pb,
          aes(y=fct_reorder(Lake,Lacustrine_morphs,mean,.desc=T),
              x=length_mean,colour=Species,shape=Morph_diet))+
  geom_point(size=4)+
  scale_shape_manual(values=c(17,16))+
  theme_classic()+xlab("Fork Length")+ylab("Lake")+
  scale_color_manual(values=cols.sal)+
  theme(panel.grid.major.y=element_line(size=.1,color="grey80"))+
  theme(legend.position="none")


FL.pb<-phenotype_repeat(morphs.FL,"length_mean","Sal",F)
FL.pb<-subset(FL.pb,FL.pb$Morph_diet%in%c("piscivore","planktivore"))
#GR.pb<-GR.pb[GR.pb$Species=="Coregonus_lavaretus",]
trait_col<-"length_mean"
temp1<-as.data.frame(pivot_wider(FL.pb,names_from=Lake,values_from=all_of("length_mean"))) #means in columns with each lake as a column name
ad.test(temp1[,6:ncol(temp1)])

temp<-data.frame(table(FL.pb$Lake))
FL.pb<-FL.pb[FL.pb$Lake%in%temp[temp$Freq>1,"Var1"],]

C<-ggplot(data=FL.pb,
          aes(y=fct_reorder(Lake,Lacustrine_morphs,mean,.desc=T),
              x=length_mean,colour=Species,shape=Morph_diet))+
  geom_point(size=4)+
  theme_classic()+xlab("Fork Length")+ylab("Lake")+
  scale_color_manual(values=cols.sal)+
  scale_shape_manual(values=c(15,16))+
  theme(panel.grid.major.y=element_line(size=.1,color="grey80"))+
  theme(legend.position="none")

grid.arrange(A,B,C,ncol=3)



#GET LEGEND
FL.pb<-phenotype_repeat(morphs.fl,"length_mean","length_n","Sal",F)
FL.pb<-subset(FL.pb,FL.pb$Morph_diet%in%c("piscivore","planktivore","benthivore"))
#GR.pb<-GR.pb[GR.pb$Species=="Coregonus_lavaretus",]
trait_col<-"length_mean"
temp1<-as.data.frame(pivot_wider(FL.pb,names_from=Lake,values_from=all_of("trait"))) #means in columns with each lake as a column name
ad.test(temp1[,7:ncol(temp1)])

temp<-data.frame(table(FL.pb$Lake))
FL.pb<-FL.pb[FL.pb$Lake%in%temp[temp$Freq>1,"Var1"],]

ggplot(data=FL.pb,
       aes(y=fct_reorder(Lake,Lacustrine_morphs,mean,.desc=T),
           x=length_mean,colour=Species,shape=Morph_diet))+
  geom_point(size=4)+
  theme_classic()+xlab("Fork Length")+ylab("Lake")+
  scale_color_manual(values=cols.sal)+
  scale_shape_manual(values=c(17,15,16))+
  theme(panel.grid.major.y=element_line(size=.1,color="grey80"))

grid.arrange(A,B,C,ncol=3)






####estimate repeatability: interaction relative to main effects


#add error type column for total length
morphs.fl$Total_length_etype<-"Variance"

morphs.fl[morphs.fl$Lake=="Vaggatem",c(1:6,106:109,174:176)]


#make dataframe for results
res.tr<-data.frame()
#cycle through each trait
for(tr.i in c("gill raker count","total body length","C13","N15","age")){
  
  #for each trait, assign correct column numbers for mean, sample size, error, error type
  if(tr.i=="gill raker count"){
    tr_cols<-7:10
    tr<-morphs[c(1:6,tr_cols)]}
  if(tr.i=="total body length"){ #use total length, based on adjustments
    tr_cols<-c(174,107,175,176)
    tr<-morphs.fl[c(1:6,tr_cols)]}
  if(tr.i=="C13"){
    tr_cols<-98:101
    tr<-morphs[c(1:6,tr_cols)]}
  if(tr.i=="N15"){
    tr_cols<-102:105
    tr<-morphs[c(1:6,tr_cols)]}
  if(tr.i=="age"){
    tr_cols<-160:163
    tr<-morphs[c(1:6,tr_cols)]}
  
  #cycle through ecological categories
  for(eco.id in c("diet","habitat")){
    
    #rename columns, set as numeric
    colnames(tr)[7:10]<-c("tr.mean","tr.n","tr.e","tr.etype")
    tr$tr.mean<-as.numeric(tr$tr.mean)
    tr$tr.n<-as.numeric(tr$tr.n)
    tr$tr.e<-as.numeric(tr$tr.e)
    
    tr<-tr[!is.na(tr$tr.n),] # remove if no sample size
    tr<-tr[!is.na(tr$tr.e),] # remove if no error
    
    #remove if morph sample size less than 5
    tr<-tr[tr$tr.n>=5,]
    
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



###################summarize # of lakes
#put together metadata for lakes

lakes.list.trait.2<-unique(res.tr$Lake)
lakes.list.all<-unique(c(lakes.list.eco,lakes.list.trait,lakes.list.trait.2))
df.lakes<-df%>%filter(Lake%in%lakes.list.all)
unique(c(lakes.list.trait,lakes.list.trait.2))
lakes.list.eco
#write.csv(df.lakes,"metaLakes.csv")
morphs.lakes<-morphs%>%filter(Lake%in%lakes.list.all)
#write.csv(morphs.lakes,"metaMorphs.csv")
length(unique(morphs.lakes$Paper_ID))

colnames(df.lakes)
df.lakes<-left_join(df.lakes,m%>%select(Lake,Genus)%>%distinct())


lakes_summary<-df.lakes%>%select(Lake,Genus,Lat,Lon,Lacustrine_morphs)%>%
  mutate(Continent=case_when((Lon>(-28)&Lon<45)~"Europe",Lon<(-28)~"N_America",
                             Lon>45~"Asia"))%>%
  mutate(Ecotype=if_else(Lake%in%lakes.list.eco,"x",""),
         Trait=if_else(Lake%in%lakes.list.trait,"x",""),
         EcotypeXTrait=if_else(Lake%in%lakes.list.trait.2,"x",""))%>%
  rename(EcotypeCount=Lacustrine_morphs,Latitude=Lat,Longitude=Lon)
lakes_summary%>%group_by(Continent)%>%
  summarise(count=n())
#write.csv(lakes_summary,"supp_table_lakes.csv",row.names=F)



##############plot lake locations



world <- ne_countries(scale = "medium", returnclass = "sf")
#Restrict to lakes used in study
df<-subset(df,df$Lake%in%lakes.list.all)

inset_coords<-data.frame(x1=c(15,110,5,-8),x2=c(35,121,10,-1),
                         y1=c(66,55,45,56),y2=c(72,58,48,59))
all(inset_coords$x1<inset_coords$x2&inset_coords$y1<inset_coords$y2)

ggplot(data=world) + 
  geom_sf(color="grey50",fill="grey95")+
  geom_point(data=df,aes(x=Lon,y=Lat,shape=morphs_level),
             colour="grey20",stroke=1.5,size=2)+
  scale_shape_manual(values=c(1,0,2))+
  theme(panel.background = element_rect(fill="#bcc4d6"))+
  geom_rect(data=inset_coords,
            mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2),fill=NA,colour="#313D5A",size=1.5)+
  ylim(40,83)+xlim(-163,163)+
  xlab("Longitude")+ylab("Latitude")

#Scandanavia
ggplot(data=world) + 
  geom_sf(color="grey50",fill="grey95")+
  geom_point(data=df,aes(x=Lon,y=Lat,shape=morphs_level),
             colour="grey20",stroke=1.5,size=2)+
  scale_shape_manual(values=c(1,0,2))+
  theme(panel.background = element_rect(fill="#bcc4d6"),legend.position="none")+
  ylim(66,72)+xlim(15,35)+ylab("")+xlab("")


#Iceland
ggplot(data=world) + 
  geom_sf(color="grey50",fill="grey95")+
  geom_point(data=df,aes(x=Lon,y=Lat,shape=morphs_level),
             colour="grey20",stroke=1.5,size=2)+
  scale_shape_manual(values=c(1,0,2))+
  theme(panel.background = element_rect(fill="#bcc4d6"),legend.position="none")+
  ylim(63,67)+xlim(-25,-12)+ylab("")+xlab("")


#Russia
ggplot(data=world) + 
  geom_sf(color="grey50",fill="grey95")+
  geom_point(data=df,aes(x=Lon,y=Lat,shape=morphs_level),
             colour="grey20",stroke=1.5,size=2)+
  scale_shape_manual(values=c(1,0,2))+
  theme(panel.background = element_rect(fill="#bcc4d6"),legend.position="none")+
  ylim(55,58)+xlim(110,121)+ylab("")+xlab("")

#Switzerland
ggplot(data=world) + 
  geom_sf(color="grey50",fill="grey95")+
  geom_point(data=df,aes(x=Lon,y=Lat,shape=morphs_level),
             colour="grey20",stroke=1.5,size=2)+
  scale_shape_manual(values=c(1,0,2))+
  theme(panel.background = element_rect(fill="#bcc4d6"),legend.position="none")+
  ylim(45,48)+xlim(5,10)+ylab("")+xlab("")

#Scotland
ggplot(data=world) + 
  geom_sf(color="grey50",fill="grey95")+
  geom_point(data=df,aes(x=Lon,y=Lat,shape=morphs_level),
             colour="grey20",stroke=1.5,size=2)+
  scale_shape_manual(values=c(1,0,2))+
  theme(panel.background = element_rect(fill="#bcc4d6"),legend.position="none")+
  ylim(56,59)+xlim(-8,-1)+ylab("")+xlab("")




rm(world)














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
plots.res.tr[[16]]


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


#


#figure 4

morphs.fl.pr<-rbind(morphs.fl.pr.S,morphs.fl.pr.C)
morphs.fl.pr$Genus<-t(data.frame(strsplit(morphs.fl.pr$Species,"_")))[,1]
TL.plot<-ggplot(data=morphs.fl.pr,
                aes(y=Lake,x=trait,colour=Genus,shape=Genus))+
  geom_point(size=3,stroke=2)+
  scale_colour_manual(values=c(col.cor,col.sal))+
  scale_shape_manual(values=c(0,2))+
  theme_classic()+xlab("Total Length")+ylab("Lake")+
  theme(panel.grid.major.y=element_line(size=.1,color="grey80"))

TL.plot


ggarrange(TL.plot,
          ggarrange(plots.mean.res.tr[[15]],
                    plots.mean.res.tr[[7]],
                    nrow=2,ncol=1,
                    labels=c("B","C")),
          common.legend=T,legend="right",
          ncol=2,widths=c(1,0.8),
          labels=c("A",""))







#Main text figures

res.gr.cor.h.2$Morph_niche<-factor(res.gr.cor.h.2$Morph_niche,levels=c("pelagic",
                                                                       "littoral/benthic","profundal"))
IJK<-ggplot(data=res.gr.cor.h.2,
            aes(x=Morph_niche,y=tr.mean.2,colour=Genus))+
  geom_point(size=3)+
  geom_line(aes(group=Lake),size=1)+
  ylab("Gill raker count")+xlab("Ecotype - Habitat")+
  theme_classic()+scale_colour_manual(values=c(col.cor))+
  theme(legend.position="none")

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
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
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
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
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
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
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
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
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



fig4<-ggarrange(TL.plot,
          ggarrange(A1,B,C,A2,get_legend(TL.plot),
                    nrow=2,ncol=3,
                    labels=c("B","C","D","E")),
          common.legend=T,legend="none",
          ncol=2,widths=c(1,1.2),
          labels=c("A",""))+ bgcolor("white") 

ggsave("figure4.png",plot=fig4,dpi=300,width=22,height=16,units="cm")


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
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
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
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
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
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
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


fig5<-ggarrange(GR.plot,
          ggarrange(I,J,K,get_legend(GR.plot),
                    nrow=2,ncol=2,
                    labels=c("B","C","D")),
          common.legend=T,legend="none",
          ncol=2,widths=c(1,1),
          labels=c("A",""))+ bgcolor("white") 

ggsave("figure5.png",plot=fig5,dpi=300,width=22,height=16,units="cm")





#


lakes.list.trait<-unique(c(datN.all$Lake,lakes.list.trait))
lakes.list<-unique(c(lakes.list.eco,lakes.list.trait))





df.2<-df%>%filter(df$Lake%in%lakes.list)%>%
  left_join(m,by="Lake")%>%
  select(Lake,Genus,Species,Lacustrine_morphs,Country,Lat,Lon,Depth_max,Area,morphs_level)%>%
  distinct()
table(df.2$Genus,df.2$morphs_level)
table(df.2$Genus,df.2$Lacustrine_morphs)

df.2%>%group_by(Genus)%>%
  summarise(Depth=range(Depth_max,na.rm=T),
            Area=range(Area,na.rm=T),
            Lat=range(Lat,na.rm=T),
            Lon=range(Lon,na.rm=T))

range(df$Area,na.rm=T)

distinct(df[df$Area>10000,])

#write.csv(df.2,file="lake_summary_table.csv")


morphs.papers<-subset(morphs,morphs$Lake%in%lakes.list)
length(unique(morphs$Paper_ID))
getwd()
m2<-read.csv("../../data/fewer_than_five.csv")
m2<-subset(m2,m2$Fewer_than_five==T)
m2[m2$ï..Lake%in%lakes.list,]
