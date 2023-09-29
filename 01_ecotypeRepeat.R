###Run analyses evaluating repeatability of ecotype assemblages by habitat/diet categories
##Includes both analyses and figures 2 and 3

###Read in packages
pkgs<-c("tidyverse","ggplot2","coxed","gridExtra","matrixcalc","forcats",
        "vegan","rgeos","rnaturalearth", "rnaturalearthdata","sf","bipartite",
        "kSamples","ggpubr","metaDigitise","metafor","ape","rfishbase")

lapply(pkgs,library,character.only=T); rm(pkgs)
source("00_functions_morphs_GEB.R")

###Colours
col.cor<-"#8cb5a7"
col.sal<-"#b58c93"
cols.sal<-c("#241c1d","#5b464a","#917076","#c4a3a9")
cols.cor<-c("#2a3632","#546d64","#709186","#a3c4b9")


###Read in data

df<-read.csv("lakes.csv",stringsAsFactors=F)
m<-read.csv("ecotypes_sum.csv",stringsAsFactors=F)
morphs<-read.csv("ecotypes.csv",stringsAsFactors=F,fileEncoding="latin1")


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
m.d<-m%>%filter(Morph_diet!="")%>%
  select(Lake,Species,Morph_diet)

t_m.d<-data.frame(addmargins(table(m.d$Lake,m.d$Morph_diet)))
colnames(t_m.d)<-c("Lake","Ecotype","Freq")
temp<-t_m.d[t_m.d$Ecotype=="Sum",]
temp<-merge(temp,df,by="Lake")%>%select(Lake,Freq,Lacustrine_morphs)
#get vector of lakes where number of morphs with ID is the number in the lake
temp<-temp[temp$Freq==temp$Lacustrine_morphs,]
t_m.d<-t_m.d[t_m.d$Lake%in%temp$Lake,]
t_m.d<-t_m.d[t_m.d$Ecotype!="Sum",]
temp<-temp[order(temp$Lacustrine_morphs),]
#reorder factor levels for Lake by number of morphs
t_m.d$Lake<-factor(t_m.d$Lake,levels=temp$Lake)
t_m.d$Lake<-fct_rev(t_m.d$Lake)
t_m.d<-merge(t_m.d,temp%>%select(-Freq),by="Lake")

t_m.d.p<-t_m.d
t_m.d.p$Ecotype<-factor(t_m.d.p$Ecotype,levels=c("benthivore","planktivore","piscivore",
                                                 "generalist.invertivore",
                                                 "generalist invertivore","omnivore"))


temp<-data.frame(table(m$Lake,m$Species))
temp<-temp[temp$Freq!=0,]
colnames(temp)<-c("Lake","Species","Freq2")
t_m.d.p<-merge(t_m.d.p,temp,by="Lake")
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
t_m.d.cor<-subset(t_m.d.p,t_m.d.p$Genus=="Cor")
t_m.d.sal<-subset(t_m.d.p,t_m.d.p$Genus=="Sal")

t_m.d.cor.2<-t_m.d.cor[t_m.d.cor$Lacustrine_morphs==2,]
t_m.d.sal.2<-t_m.d.sal[t_m.d.sal$Lacustrine_morphs==2,]



#Morph ID habitat
m.z<-m%>%filter(Morph_zone!="")%>%
  select(Lake,Species,Morph_zone)

t_m.z<-data.frame(addmargins(table(m.z$Lake,m.z$Morph_zone)))
colnames(t_m.z)<-c("Lake","Ecotype","Freq")
temp<-t_m.z[t_m.z$Ecotype=="Sum",]
temp<-merge(temp,df,by="Lake")%>%select(Lake,Freq,Lacustrine_morphs)

#get vector of lakes where number of morphs with ID is the number in the lake
temp<-temp[temp$Freq==temp$Lacustrine_morphs,]
t_m.z<-t_m.z[t_m.z$Lake%in%temp$Lake,]
t_m.z<-t_m.z[t_m.z$Ecotype!="Sum",]
temp<-temp[order(temp$Lacustrine_morphs),]

#reorder factor levels for Lake by number of morphs
t_m.z$Lake<-factor(t_m.z$Lake,levels=temp$Lake)
t_m.z$Lake<-fct_rev(t_m.z$Lake)
t_m.z<-merge(t_m.z,temp%>%select(-Freq),by="Lake")

t_m.z.p<-t_m.z
t_m.z.p$Ecotype<-factor(t_m.z.p$Ecotype,levels=c("littoral.benthic","littoral/benthic",
                                                 "pelagic","profundal","generalist","shallow"))
t_m.z.p$Ecotype[t_m.z.p$Ecotype=="littoral/benthic"]<-"littoral.benthic"


temp<-data.frame(table(m$Lake,m$Species))
temp<-temp[temp$Freq!=0,]
colnames(temp)<-c("Lake","Species","Freq2")
t_m.z.p<-merge(t_m.z.p,temp,by="Lake")


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
#ggarrange(e2,e1,ncol=2,nrow=1,labels=c("A","B"),align="v")
#ggarrange(f2,f1,ncol=2,nrow=1,labels=c("A","B"),align="v")


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
#ggarrange(c1,d1,c2,d2,ncol=2,nrow=2,labels=c("A","B","C","D"),align="v")

#clean up environment
rm(t_m.d.cor,t_m.d.cor.2,t_m.d.cor.pa,
   t_m.z.cor,t_m.z.cor.2,t_m.z.cor.pa,
   t_m.d.sal,t_m.d.sal.2,t_m.d.sal.pa,
   t_m.z.sal,t_m.z.sal.2,t_m.z.sal.pa,
   t_m.d.p,t_m.z.p,m.d,m.z,
   a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,fig2,fig3,temp)


###Make dataframes for analyses

#Make diet dataframe
mat.d<-m%>%filter(Morph_diet!="")%>%select(Lake,Species,Morph_diet)
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

mat.d.cor<-mat.d[1:5,mat.d[6,]%in%c(NA,"Cor")] #Coregonus only and remove Genus row
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
mat.z<-m%>%filter(Morph_zone!="")%>%select(Lake,Species,Morph_zone) #morph ID info
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
mat.z.cor<-mat.z[1:5,mat.z[6,]%in%c(NA,"Cor")]
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

