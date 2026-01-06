#################################################################
#         		Supplement S6: Analysis script										
#################################################################
# Dahl et al 2022. Long-term warming-induced trophic downgrading in the soil microbial food web 

# ********************************************************************************************************************

# ********************************************************************************************************************
#Figure 1: NMDS Quantification & taxomoic responses

# load packages
ls(all.names=TRUE)
rm(list=ls(all=TRUE))

library(vegan)
library(ggplot2)
library(qdapTools)
library(reshape2)

### Load the data:
data<-read.table("forhot_200K_trim_GO.txt", sep="\t", header=T, quote=NULL, comment='')
meta <- read.csv("meta.csv",header = TRUE, sep = ";")
dataV <-read.table("viruses_forhot.txt", sep="\t", header=T, comment='')

# ********************************************************************************************************************

### Figure 1: NMDS on rRNA families:

### aggregate on family level
data_fam=data[,c(17,1:8)] 
data_fam <- aggregate(.~fam, FUN="sum", data=data_fam) 
data_fam <- data_fam[data_fam$fam!="",] ### remove unclassified
sum(data_fam[,2:8])/sum(data[,1:8])*100 # 45.1% of all data 

### transform
var <- meta[,c(1, 7:22)] 
data_t<- as.data.frame(t(data_fam[,2:9]))
data_t <- data_t[ order(match(row.names(data_t), var$ID)), ]

### NMDS plot
NMDS=metaMDS(data_t,k=2,trymax=1000, distance = "bray") 
NMDS$stress #stress: 0.04533972 

### Get sample scores
samples2d <-as.data.frame(NMDS$points)

### Get meta data
place <- meta[,c(1,3)]
place$ID <- as.character(place$ID)
Treatment <- meta[,c(1,4)]
Treatment$ID <- as.character(Treatment$ID)
RNA <- meta[,c(1,6)] 
RNA$ID <- as.character(RNA$ID)
samples2d$place  <- lookup(rownames(samples2d), place[1:2])
samples2d$placecol  <- lookup(rownames(samples2d), place[1:2])
samples2d$placecol <-as.character(samples2d$placecol)
samples2d$placecol[which(samples2d$placecol == "LTW")]  <- "hollow"
samples2d$treat  <- lookup(rownames(samples2d), Treatment[1:2])
samples2d$RNA  <- lookup(rownames(samples2d), RNA[1:2])
rownames(samples2d) <- meta$samples


### Plot
ggplot(samples2d, aes(x=MDS1, y=MDS2, color=factor(treat),shape=factor(place))) + 
  geom_point(size=sqrt(samples2d$RNA)/sqrt(max(samples2d$RNA))*15,
             aes(fill=factor(treat), alpha = as.factor(place)))+
  geom_point(size=sqrt(samples2d$RNA)/sqrt(max(samples2d$RNA))*15) +    
  geom_text(aes(label=row.names(samples2d)),size=3, position = position_nudge(y = -0.02)) + 
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  scale_shape_manual(values=c(21,24)) +
  scale_alpha_manual(values=c(1,0)) +
  scale_colour_manual(values = c("indianred3","lightskyblue4"))+
  scale_fill_manual(values =alpha(c("indianred3","lightskyblue4")))+ 
  theme_bw()

NMDS.plot

# ********************************************************************************************************************

### Figure 1: NMDS on viral families:

### aggregate on family level
data_fam=dataV[,c(7,11:26)] 
data_fam <- aggregate(.~family, FUN="sum", data=data_fam) 
data_fam <- data_fam[data_fam$family!="f__unclassified_LCA",] # rm unclassified 
data_fam <- data_fam[data_fam$family!="f__unclassified",] # rm unclassified 
sum(data_fam[,2:17])/sum(dataV[,11:26])*100 # 68.42852 of all data 

### transform
var <- meta[,c(1, 7:22)] 
var$ID<-c(colnames(data_fam[,c(2:9)]))
data_t<- as.data.frame(t(data_fam[,2:9]))
data_t <- data_t[ order(match(row.names(data_t), var$ID)), ]

### NMDS plot
NMDS=metaMDS(data_t,k=2,trymax=1000, distance = "bray") 
stressplot(NMDS)
NMDS$stress #0.05950881

#### get sample scores
samples2d <-as.data.frame(NMDS$points)
rownames(samples2d) <- meta$samples

### Get meta data
place <- meta[,c(1,3)]
place$ID <- as.character(meta$samples)
Treatment <- meta[,c(1,4)]
Treatment$ID <- as.character(meta$samples)
RNA <- meta[,c(1,6)] 
RNA$ID <- as.character(meta$samples)
samples2d$place  <- lookup(rownames(samples2d), place[1:2])
samples2d$placecol  <- lookup(rownames(samples2d), place[1:2])
samples2d$placecol <-as.character(samples2d$placecol)
samples2d$placecol[which(samples2d$placecol == "LTW")]  <- "hollow"
samples2d$treat  <- lookup(rownames(samples2d), Treatment[1:2])
samples2d$RNA  <- lookup(rownames(samples2d), RNA[1:2])

### Plot
ggplot(samples2d, aes(x=MDS1, y=MDS2, color=factor(treat),shape=factor(place))) + 
  geom_point(size=sqrt(samples2d$RNA)/sqrt(max(samples2d$RNA))*15,
             aes(fill=factor(treat), alpha = as.factor(place)))+
  geom_point(size=sqrt(samples2d$RNA)/sqrt(max(samples2d$RNA))*15) +    
  geom_text(aes(label=row.names(samples2d)),size=3, position = position_nudge(y = -0.02)) + # label=row.names(samples2d) or samples2d$date
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  scale_shape_manual(values=c(21,24)) +
  scale_alpha_manual(values=c(1,0)) +
  scale_colour_manual(values = c("indianred3","lightskyblue4"))+
  scale_fill_manual(values =alpha(c("indianred3","lightskyblue4")))+ 
  theme_bw()

# ********************************************************************************************************************

### Figure 1: Quantification viruses

### load non-normalized mRNA tax only table
data_tax_ONLY<-read.table("mRNA_only_tax_allMerged_LCA100_final.tsv", sep="\t", header=T, quote=NULL, comment='') ###LOAD THE data
newColNames <- c("general_type","domain_kingdom","phylum","class","order","family","genus","species") ###split the taxString column into single columns
newCols <- colsplit(data_tax_ONLY$taxString, "::", newColNames)
data_sep <- cbind(newCols, data_tax_ONLY)
data_tax_ONLY <- data_sep ###this data_tax_ONLY table contains all taxa (species level)
rm(data_sep)
rm(newColNames)
rm(newCols)
data_tax_ONLY <- data_tax_ONLY[which(data_tax_ONLY$general_type !="unclassified_LCA"),] 

### quantify (LTW data only; aggregate viruses)
df_raw_V <- data_tax_ONLY[which(data_tax_ONLY$general_type =="_Viruses"),]
V_851nt_DW <- as.numeric(as.vector(quanT["V_851nt_DW",])) #make a vector
df_raw_V[,c(10:25)] <- sweep(df_raw_V[,c(10:25)],2,V_851nt_DW, FUN="*")
df_V <- df_raw_V[,c(1,10:17)]
df_V <- aggregate(.~general_type, FUN="sum", data=df_V) 
rownames(df_V) <- df_V[,1]
df_V <- df_V[,-1]

dfV_t <- as.data.frame(t(df_V))
dfV_t$site <- meta$category
dfV_t$x <- meta$samples
names(dfV_t) <- c('vira', 'site', 'x')
### plot individual bars 
df.mean = dfV_t %>% 
  group_by(site) %>% 
  mutate(ymean = mean(vira))

ggplot(dfV_t, aes(x, vira, fill=site)) +
  geom_col() +
  geom_errorbar(data=df.mean, aes(x, ymax = ymean, ymin = ymean),
                size=0.5, linetype = "longdash", inherit.aes = F, width = 1)+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### plot
df_V$LTWA <- rowMeans(df_V[,1:4])
df_V$LTWE <- rowMeans(df_V[,5:8])
df_V_ <- df_V[,9:10]
df_V_$dom <- "viruses"
df_m <- melt(df_V_, id=c("dom"))
colnames(df_m)

ggplot() + geom_bar(aes(y = value, x = variable, fill = dom), data = df_m,stat="identity") +
  theme_bw() +theme(axis.text.x= element_text(angle=45, hjust = 1.3, vjust = 1.2), legend.title = element_blank())+
  xlab("") +ylab("Proportion of sequences")

# ********************************************************************************************************************

### Significant differences?
stats_df <- as.data.frame(t(df_V[,1:8]))

t.test(`_Viruses` ~ cat , data = stats_df, var.equal = TRUE)

cordf <- cbind(dfV_t[1], meta$Cmic..µg.g.DW.)
names(cordf) <- c('vira', 'micC')

# correlation between micC and vira
cor.test(cordfr$vira, cordfr$micC, method = c("pearson"))

# removing outlier
cordfr <- cordf[-2,]
# ********************************************************************************************************************

### Figure 1: total microbial C per g soil

ls(all.names=TRUE)
rm(list=ls(all=TRUE))

meta <- read.csv("meta.csv",header = TRUE, sep = ";")
rownames(meta) <- meta$ID
meta_ <- meta[,-1:-5]

meta_t <- as.data.frame(t(meta_))
miC <- meta_t[12,1:8]

dfV_t <- as.data.frame(t(miC))
dfV_t$id <- meta$ID
dfV_t$site <- meta$category
dfV_t$x <- meta$samples
names(dfV_t) <- c('miC', 'id', 'site', 'x')

### plot individual bars 
df.mean = dfV_t %>% 
  group_by(site) %>% 
  mutate(ymean = mean(miC))

ggplot(dfV_t, aes(x, miC, fill=site)) +
  geom_col() +
  geom_errorbar(data=df.mean, aes(x, ymax = ymean, ymin = ymean),
                size=0.5, linetype = "longdash", inherit.aes = F, width = 1)+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# ********************************************************************************************************************

### Significant differences?
stats_df <- meta[1:8,]
t.test(Cmic..µg.g.DW. ~ cat , data = stats_df, var.equal = TRUE)

# ********************************************************************************************************************
# Figure 1: Taxonomic overview - customized 

ls(all.names=TRUE)
rm(list=ls(all=TRUE))

library(ggplot2)
library(qdapTools)
library(reshape2)
library(RColorBrewer)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

### Load the data:
data<-read.table("forhot_200K_trim_GO.txt", sep="\t", header=T, quote=NULL, comment='')
meta <- read.csv("STAMP_meta.csv",header = TRUE, sep = ";")
# ********************************************************************************************************************

### Bactria phyla 
data_bac <- data[which(data$dom =="Bacteria"),] 
data_bac <- aggregate(.~phy, FUN="sum", data=data_bac[,c(1:8,14)]) 
data_bac <- data_bac[data_bac$phy!="",] 

row.names(data_bac) <- as.character(data_bac$phy)
data_bac <- data_bac[,-1]
data_bac$tax <- rownames(data_bac)
data_bac$sum <- rowSums(data_bac[1:8])

### select abundant bacterial phyla
for (i in 1:nrow(data_bac)){
  data_bac$tax[i] <- if (data_bac$sum[i]<50000) {"xBac_Other"} else {data_bac$tax[i]} #or 6e+11 100000
}
data_bac <- aggregate(.~tax, FUN="sum", data=data_bac)
data_bac <-data_bac[,-10]
names(data_bac) <-gsub("X", "", names(data_bac))

### Proteobacteria (pull them out and merge with the rest)
data_pro <- data[which(data$phy =="Proteobacteria"),] 
data_pro<-data_pro[,c(-9:-14,-16:-21)]
data_pro <- aggregate(.~class, FUN="sum", data=data_pro) 
data_pro$tax <-data_pro$class

data_pro<-data_pro[,c(10,2:9)]
names(data_pro) <-gsub("X", "", names(data_pro))
data_bac <- rbind(data_bac,data_pro)

data_bac <- data_bac[which(data_bac$tax!="Proteobacteria"),] 
data_bac$tax <- gsub("Gammaproteobacteria", "xBac_Other",data_bac$tax)
data_bac$tax <- gsub("Ca. Tenderia class incertae sedis", "xBac_Other",data_bac$tax)
data_bac$tax <- gsub("Ca. Thiobios class incertae sedis", "xBac_Other",data_bac$tax)
data_bac$tax <- gsub("Zetaproteobacteria", "xBac_Other",data_bac$tax)
data_bac[7,1]<-"xBac_Other"

data_bac$tax <-paste("Bac_",data_bac$tax)
data_bac$tax <- gsub(" ", "",data_bac$tax)
data_bac <- aggregate(.~tax, FUN="sum", data=data_bac)

### Fungi
data_fun <- data[which(data$king =="Fungi"),] 
data_fun <-  data_fun[,c(-9:-14,-16:-21)]
data_fun <- aggregate(.~class, FUN="sum", data=data_fun) 
data_fun <- data_fun[data_fun$class!="",] ### Fungal reads not assigned to class level are removed

row.names(data_fun) <- as.character(data_fun$class)
data_fun <- data_fun[,-1]
data_fun$tax <- rownames(data_fun)
data_fun$sum <- rowSums(data_fun[1:8])

### select abundant fungal classes
for (i in 1:nrow(data_fun)){
  data_fun$tax[i] <- if (data_fun$sum[i]<15000) {"xFun_Other"} else {data_fun$tax[i]} 
}
data_fun <- aggregate(.~tax, FUN="sum", data=data_fun)
data_fun <-data_fun[,-10]
data_fun$tax <-paste("EuF_",data_fun$tax)
names(data_fun) <-gsub("X", "", names(data_fun))
data_fun$tax <- gsub(" ", "",data_fun$tax)

### protist

dataeu <- data[which(data$dom =="Eukaryota"),] 
data_prot <- dataeu[which(dataeu$king!="Fungi"),]
data_prot <- data_prot[which(data_prot$king!="Plantae"),]
data_prot <- data_prot[which(data_prot$king!="Metazoa (Animalia)"),]

data_prot$phy <- as.character(data_prot$phy)
data_prot$supK <- as.character(data_prot$supK)
for(i in 1:length(data_prot$phy)){
  if(grepl("Amoebozoa", data_prot$supK[i], perl=T) ==T) {data_prot$phy[i] <- data_prot$supK[i]} else {data_prot$phy[i] <- data_prot$phy[i]}
  if(grepl("Excavata", data_prot$supK[i], perl=T) ==T) {data_prot$phy[i] <- data_prot$supK[i]} else {data_prot$phy[i] <- data_prot$phy[i]}
  if(grepl("Hacrobia", data_prot$supK[i], perl=T) ==T) {data_prot$phy[i] <- data_prot$supK[i]} else {data_prot$phy[i] <- data_prot$phy[i]}
} 
prot <- data_prot[,c(1:8,14)]
data_prot <- aggregate(.~phy, FUN="sum", data=prot) 

data_prot <- data_prot[data_prot$phy!="",] ### reads not assigned to phy level are removed

row.names(data_prot) <- as.character(data_prot$phy)
data_prot <- data_prot[,-1]
data_prot$tax <- rownames(data_prot)
data_prot$sum <- rowSums(data_prot[1:8])

### select abundant protist phyla
for (i in 1:nrow(data_prot)){
  data_prot$tax[i] <- if (data_prot$sum[i]<5000) {"xPro_Other"} else {data_prot$tax[i]} 
} 
data_prot <- aggregate(.~tax, FUN="sum", data=data_prot)
data_prot <-data_prot[,-10]
data_prot$tax <-paste("EuP_",data_prot$tax)
names(data_prot) <-gsub("X", "", names(data_prot))
data_prot$tax <-gsub("EuP_ Retaria", "EuP_ Foraminifera", data_prot$tax)
data_prot$tax <- gsub(" ", "",data_prot$tax)

### metazoa
data_met <- data[which(data$king =="Metazoa (Animalia)"),] 
data_met <-  data_met[,c(-9:-13,-15:-21)]### minus other tax levels 
data_met <- aggregate(.~phy, FUN="sum", data=data_met) 

data_met <- data_met[data_met$phy!="",] ### reads not assigned to phy level are removed

row.names(data_met) <- as.character(data_met$phy)
data_met <- data_met[,-1]
data_met$tax <- rownames(data_met)
data_met$sum <- rowSums(data_met[1:8])

### select abundant metazoa phyla
for (i in 1:nrow(data_met)){
  data_met$tax[i] <- if (data_met$sum[i]<5000) {"xmet_Other"} else {data_met$tax[i]} 
}
data_met <- aggregate(.~tax, FUN="sum", data=data_met)
data_met <-data_met[,-10]
data_met$tax <-paste("EuM_",data_met$tax)
names(data_met) <-gsub("X", "", names(data_met))
data_met$tax <- gsub(" ", "",data_met$tax)

### Archeae
data_ar <- data[which(data$dom =="Archaea"),] 
data_ar <-  data_ar[,c(-9:-10,-12:-21)] 
data_ar <- aggregate(.~dom, FUN="sum", data=data_ar) 
data_ar$tax <-data_ar$dom
data_ar <- data_ar[,c(10,2:9)]
names(data_ar) <-gsub("X", "", names(data_ar))
data_ar$tax <- gsub(" ", "",data_ar$tax)


### Combine the data frames
data_cust <- rbind(data_bac, data_fun)
data_cust <- rbind(data_cust, data_prot)
data_cust <- rbind(data_cust, data_met)
data_cust <- rbind(data_cust, data_ar)
data_cust$tax <- gsub("EuP_xPro_Other", "other_Eukaryota",data_cust$tax)
data_cust$tax <- gsub("EuM_xmet_Other", "other_Eukaryota",data_cust$tax)
data_cust <- aggregate(.~tax, FUN="sum", data=data_cust)

df_m <- melt(data_cust, id=c("tax"))
df_m$value <- as.numeric(df_m$value)

### look-up the sample category
matchlist <- meta[,c(1,5)] 
matchlist$category <- as.character(matchlist$category)
matchlist$ID <- as.character(matchlist$ID)
df_m$cat <- NA
df_m$cat <- lookup(df_m$variable, matchlist[1:2])

#### reorder the dataframe
df_m$variable <- factor(df_m$variable, levels = unique(df_m$variable[order(df_m$cat)]))

set.seed(2)
mycols <- length(unique(df_m$tax))
p <- ggplot() + geom_bar(aes(y = value, x = cat, fill = tax), data = df_m, stat="identity", position = "fill") +
  theme_bw() +theme(axis.text.x= element_text(angle=45, hjust = 1.3, vjust = 1.2), legend.title = element_blank())+
  xlab("") +ylab("Proportion of sequences") +scale_fill_manual(values = getPalette(mycols))
p

# ********************************************************************************************************************

### Significant differences?
df <- data_cust
rownames(df) <- df$tax
df <- df[,-1]
df <- rbind(df, df[15,]+df[16,])
df_stat <- as.data.frame(t(df))
df_stat$cat <- c("LTWA","LTWA","LTWA","LTWA","LTWE","LTWE","LTWE","LTWE")

t.test(Archaea ~ cat , data = df_stat, var.equal = TRUE)
t.test(Bac_Alphaproteobacteria ~ cat , data = df_stat, var.equal = TRUE)
t.test(Bac_Betaproteobacteria ~ cat , data = df_stat, var.equal = TRUE)
t.test(Bac_Deltaproteobacteria ~ cat , data = df_stat, var.equal = TRUE)
t.test(Bac_Actinobacteria ~ cat , data = df_stat, var.equal = TRUE)
t.test(Bac_Acidobacteria ~ cat , data = df_stat, var.equal = TRUE)
t.test(Bac_Planctomycetes ~ cat , data = df_stat, var.equal = TRUE)
t.test(Bac_Verrucomicrobia ~ cat , data = df_stat, var.equal = TRUE)
t.test(Bac_Chloroflexi ~ cat , data = df_stat, var.equal = TRUE)
t.test(Bac_xBac_Other ~ cat , data = df_stat, var.equal = TRUE)
t.test(EuP_Amoebozoa  ~ cat , data = df_stat, var.equal = TRUE)
t.test(EuP_Foraminifera ~ cat , data = df_stat, var.equal = TRUE)
t.test(EuP_Cercozoa ~ cat , data = df_stat, var.equal = TRUE)
t.test(EuF_Glomeromycetes ~ cat , data = df_stat, var.equal = TRUE)
t.test(EuF_Archaeorhizomycetes ~ cat , data = df_stat, var.equal = TRUE)
t.test(EuF_Agaricomycetes ~ cat , data = df_stat, var.equal = TRUE)
t.test(EuF_xFun_Other ~ cat , data = df_stat, var.equal = TRUE)
t.test(EuM_Arthropoda ~ cat , data = df_stat, var.equal = TRUE)
t.test(EuM_Nematoda ~ cat , data = df_stat, var.equal = TRUE)
t.test(other_Eukaryota ~ cat , data = df_stat, var.equal = TRUE)
t.test(EuM_Arthropoda1 ~ cat , data = df_stat, var.equal = TRUE) #Metazoa combined, p-value = 0.5822

# ********************************************************************************************************************
# Figure 2: Microbial food-web response plot

ls(all.names=TRUE)
rm(list=ls(all=TRUE))

library(ggplot2)
library(reshape2)

### Load the data:
data<-read.table("forhot_200K_trim_GO.txt", sep="\t", header=T, quote=NULL, comment='')

# ********************************************************************************************************************

### now select and aggregate on different levels:

### Archaea ****************************************
archaea_t <- data[which(data$dom == 'Archaea'),]
archaea_t <- archaea_t[,c(1:8,11)]
archaea_t <- aggregate(.~dom, FUN="sum", data=archaea_t)
archaea_t$tax <- archaea_t$dom
archaea_t$color <- archaea_t$dom
archaea_t <- archaea_t[,-1]
### *******************************************

rest <- data[which(data$dom != 'Archaea'),]

### Fungi ****************************************
fungi_t <- data[which(rest$king == 'Fungi'),]
fungi_t <- fungi_t[,c(1:8,14)]
fungi_t <- aggregate(.~phy, FUN="sum", data=fungi_t)
fungi_t$tax <- "Fungi_(excl_AM)"
fungi_t[10,10] <- "AM_Fungi"
fungi_t <- fungi_t[,-1]
fungi_t <- aggregate(.~tax, FUN="sum", data=fungi_t)
fungi_t$color <- "Fungi"
fungi_t <- fungi_t[,c(2:9,1,10)]
### *******************************************

THE_TABLE <- rbind(archaea_t,fungi_t)
rest <- rest[which(rest$king != 'Fungi'),]

### Bdellovibrionales ****************************************
bdello_t <- rest[which(rest$ord == 'Bdellovibrionales'),]
bdello_t <- bdello_t[,c(1:8,16)]
bdello_t <- aggregate(.~ord, FUN="sum", data=bdello_t)
bdello_t$tax <- bdello_t$ord
bdello_t$color <- "Bacterivorous_Bacteria"
bdello_t <- bdello_t[,-1]
### *******************************************

THE_TABLE <- rbind(THE_TABLE,bdello_t)
rest <- rest[which(rest$ord != 'Bdellovibrionales'),]

### Myxococcales ****************************************
myxo_t <- rest[which(rest$ord == 'Myxococcales'),]
myxo_t <- myxo_t[,c(1:8,17)]
myxo_t <- aggregate(.~fam, FUN="sum", data=myxo_t)
myxo_t$tax <- myxo_t$fam
myxo_t[1,10] <- "uncl. Myxococcales"
myxo_t[4,10] <- "other_Myxococcales"
myxo_t[5,10] <- "other_Myxococcales"
myxo_t[10,10] <- "other_Myxococcales"
myxo_t <- myxo_t[,-1]
myxo_t <- aggregate(.~tax, FUN="sum", data=myxo_t)
myxo_t$color <- "Bacterivorous_Bacteria"
myxo_t <- myxo_t[,c(2:9,1,10)]
### *******************************************

THE_TABLE <- rbind(THE_TABLE,myxo_t)
rest <- rest[which(rest$ord != 'Myxococcales'),]

### (food) Bacteria ****************************************
Bacteria_t <- rest[which(rest$dom == 'Bacteria'),]
Bacteria_t <- Bacteria_t[,c(1:8,11)]
Bacteria_t <- aggregate(.~dom, FUN="sum", data=Bacteria_t)
Bacteria_t$tax <- Bacteria_t$dom
Bacteria_t$color <- "Bacteria"
Bacteria_t <- Bacteria_t[,-1]
## *******************************************

THE_TABLE <- rbind(THE_TABLE,Bacteria_t)
rest <- rest[which(rest$dom != 'Bacteria'),]

table(rest$dom)
rest <- rest[which(rest$dom != 'Eukaryota (Chloroplast)'),]
rest <- rest[which(rest$dom != 'Eukaryota (Mitochondrion)'),]
rest <- rest[which(rest$dom != ''),]

# Nematoda ****************************************
nematoda_t <- rest[which(rest$phy == 'Nematoda'),]
### Nematode analysis according to the following two webpages, and based on Yeates 1993
### http://nemaplex.ucdavis.edu/Ecology/bioindicators.htm
### http://nemaplex.ucdavis.edu/Ecology/feeding_habits.htm
### final assignment by Tim Urich

library(readxl)
feedingmode<-read_xlsx("C:/Users/dahlm/Dropbox/PostDoc/forhot/manus_forhot/Rscripts_Andrea/_rRNA_soellinger_scripts_and_results/nematoda_feeding_mode_summary_per_row.xlsx", 
                       sheet = "feeding_mode")
nematoda_t <- nematoda_t[,c(1:17)]
nematoda_t$tax <- feedingmode$feeding_mode
THE_nematoda_table <- nematoda_t #maybe needed later
nematoda_t <- nematoda_t[,c(1:8,18)]
nematoda_t <- aggregate(.~tax, FUN="sum", data=nematoda_t)
nematoda_t[2,1]<- "Bacterial-feeding"
nematoda_t <- aggregate(.~tax, FUN="sum", data=nematoda_t)
nematoda_t$color <- "Nematoda"
nematoda_t <- nematoda_t[,c(2:9,1,10)]
### *******************************************

THE_TABLE <- rbind(THE_TABLE,nematoda_t)
rest <- rest[which(rest$phy != 'Nematoda'),]

### Entognatha ****************************************
Entognatha_t <- rest[which(rest$class == 'Entognatha'),] 
Entognatha_t <- Entognatha_t[,c(1:8,15)]
Entognatha_t <- aggregate(.~class, FUN="sum", data=Entognatha_t)
Entognatha_t$tax <- "Collembola_and_Protura"
Entognatha_t$color <- "Arthropoda"
Entognatha_t <- Entognatha_t[,-1]
### *******************************************

THE_TABLE <- rbind(THE_TABLE,Entognatha_t)
rest <- rest[which(rest$class != 'Entognatha'),]

### Haplotaxida ****************************************
Haplotaxida_t <- rest[which(rest$ord == 'Haplotaxida'),] 
Haplotaxida_t <- Haplotaxida_t[,c(1:8,16)]
Haplotaxida_t <- aggregate(.~ord, FUN="sum", data=Haplotaxida_t)
Haplotaxida_t$tax <- "Haplotaxida (earth worms)"
Haplotaxida_t$color <- "Metazoa"
Haplotaxida_t <- Haplotaxida_t[,-1]
### *******************************************

THE_TABLE <- rbind(THE_TABLE,Haplotaxida_t)
rest <- rest[which(rest$ord != 'Haplotaxida'),]

### Rotifera ****************************************
Rotifera_t <- rest[which(rest$phy == 'Rotifera'),] 
Rotifera_t <- Rotifera_t[,c(1:8,14)]
Rotifera_t <- aggregate(.~phy, FUN="sum", data=Rotifera_t)
Rotifera_t$tax <- Rotifera_t$phy
Rotifera_t$color <- "Metazoa"
Rotifera_t <- Rotifera_t[,-1]
### *******************************************

THE_TABLE <- rbind(THE_TABLE,Rotifera_t)
rest <- rest[which(rest$phy != 'Rotifera'),]

### Tardigrada ****************************************
Tardigrada_t <- rest[which(rest$phy == 'Tardigrada'),] 
Tardigrada_t <- Tardigrada_t[,c(1:8,14)]
Tardigrada_t <- aggregate(.~phy, FUN="sum", data=Tardigrada_t)
Tardigrada_t$tax <- Tardigrada_t$phy
Tardigrada_t$color <- "Metazoa"
Tardigrada_t <- Tardigrada_t[,-1]
### *******************************************

THE_TABLE <- rbind(THE_TABLE,Tardigrada_t)
rest <- rest[which(rest$phy != 'Tardigrada'),]

### Porifera ****************************************
Porifera_t <- rest[which(rest$phy == 'Porifera'),] 
Porifera_t <- Porifera_t[,c(1:8,14)]
Porifera_t <- aggregate(.~phy, FUN="sum", data=Porifera_t)
Porifera_t$tax <- Porifera_t$phy
Porifera_t$color <- "Metazoa"
Porifera_t <- Porifera_t[,-1]
### *******************************************

THE_TABLE <- rbind(THE_TABLE,Porifera_t)
rest <- rest[which(rest$phy != 'Porifera'),]

### Platyhelminthes ****************************************
Platyhelminthes_t <- rest[which(rest$phy == 'Platyhelminthes'),] 
Platyhelminthes_t <- Platyhelminthes_t[,c(1:8,14)]
Platyhelminthes_t <- aggregate(.~phy, FUN="sum", data=Platyhelminthes_t)
Platyhelminthes_t$tax <- Platyhelminthes_t$phy
Platyhelminthes_t$color <- "Metazoa"
Platyhelminthes_t <- Platyhelminthes_t[,-1]
### *******************************************************

THE_TABLE <- rbind(THE_TABLE,Platyhelminthes_t)
rest <- rest[which(rest$phy != 'Platyhelminthes'),]

### plants ****************************************
plant_t <- rest[which(rest$supK == 'Archaeplastida'),] 
plant_t <- plant_t[,c(1:8,13)]
plant_t <- aggregate(.~king, FUN="sum", data=plant_t)
### all Archaeplastida are simply removed from the dataset as it is mainly Plantea
### *******************************************************

rest <- rest[which(rest$supK != 'Archaeplastida'),]

### Insecta ****************************************
Insecta_t <- rest[which(rest$class == 'Insecta'),] 
Insecta_t <- Insecta_t[,c(1:8,15)]
Insecta_t <- aggregate(.~class, FUN="sum", data=Insecta_t)
Insecta_t$tax <- Insecta_t$class
Insecta_t$color <- "Arthropoda"
Insecta_t <- Insecta_t[,-1] 
### *******************************************************

THE_TABLE <- rbind(THE_TABLE,Insecta_t)
rest <- rest[which(rest$class != 'Insecta'),]

### Arthropoda ****************************************
Arthropoda_t <- rest[which(rest$phy == 'Arthropoda'),] 
Arthropoda_t <- Arthropoda_t[,c(1:8,14)]
Arthropoda_t <- aggregate(.~phy, FUN="sum", data=Arthropoda_t)
Arthropoda_t$tax <- "other Arthropoda"
Arthropoda_t$color <- "Arthropoda"
Arthropoda_t <- Arthropoda_t[,-1]
### *******************************************************

THE_TABLE <- rbind(THE_TABLE,Arthropoda_t)
rest <- rest[which(rest$phy != 'Arthropoda'),]

### Metazoa (Animalia) ****************************************
Metazoa_t <- rest[which(rest$king == 'Metazoa (Animalia)'),] 
Metazoa_t <- Metazoa_t[,c(1:8,14)]
Metazoa_t <- aggregate(.~phy, FUN="sum", data=Metazoa_t) #just to see which are left
Metazoa_t <- rest[which(rest$king == 'Metazoa (Animalia)'),] 
Metazoa_t <- Metazoa_t[,c(1:8,13)]
Metazoa_t <- aggregate(.~king, FUN="sum", data=Metazoa_t) 
Metazoa_t$tax <- "other Metazoa"
Metazoa_t$color <- "Metazoa"
Metazoa_t <- Metazoa_t[,-1]
### *******************************************************

THE_TABLE <- rbind(THE_TABLE,Metazoa_t)
rest <- rest[which(rest$king != 'Metazoa (Animalia)'),]

colnames(rest)
table(rest$dom)
table(rest$supK)
table(rest$king)

### Protists:incl Tim Urich's and Stefan Geisen's input **************************************
### Amoebozoa  ****************************************
Amoebozoa_t <- rest[which(rest$supK == 'Amoebozoa'),] 
Amoebozoa_t <- Amoebozoa_t[,c(1:8,13)]
Amoebozoa_t <- aggregate(.~king, FUN="sum", data=Amoebozoa_t) #just to see which are left
table(Amoebozoa_t$king)

Amoebozoa_tx <- rest[which(rest$supK == 'Amoebozoa'),] 
Amoebozoa_tx <- Amoebozoa_tx[which(Amoebozoa_tx$king == 'Amoebozoa (kingdom)'),] 
Amoebozoa_tx <- Amoebozoa_tx[,c(1:8,15)]
Amoebozoa_tx <- aggregate(.~class, FUN="sum", data=Amoebozoa_tx)

Myxogastria_t <- Amoebozoa_t[6,]
Myxogastria_t$tax <- "Myxogastria"
Myxogastria_t$color <- "Protist"
Myxogastria_t <- Myxogastria_t[,-1]

Amoebozoa_t2 <- rest[which(rest$supK == 'Amoebozoa'),]
Amoebozoa_t2 <- Amoebozoa_t2[which(Amoebozoa_t2$king != 'Myxogastria'),]
Amoebozoa_t2 <- Amoebozoa_t2[which(Amoebozoa_t2$king == 'Amoebozoa (kingdom)'),]
table(Amoebozoa_t2$king)
table(Amoebozoa_t2$phy)
Amoebozoa_t2 <- Amoebozoa_t2[,c(1:8,14)]
Amoebozoa_t2 <- aggregate(.~phy, FUN="sum", data=Amoebozoa_t2)

Amoebozoa_tx <- Amoebozoa_t2
Amoebozoa_tx$meanLTW <- rowMeans(Amoebozoa_tx[,2:9])

Amoebozoa_t2$tax <- Amoebozoa_t2$phy
Amoebozoa_t2[1,10] <- "other_and_uncl_Amoebozoa_(kingdom)"
Amoebozoa_t2[3,10] <- "other_and_uncl_Amoebozoa_(kingdom)"
Amoebozoa_t2[7,10] <- "other_and_uncl_Amoebozoa_(kingdom)"
Amoebozoa_t2[9,10] <- "other_and_uncl_Amoebozoa_(kingdom)"
Amoebozoa_t2[10,10] <- "other_and_uncl_Amoebozoa_(kingdom)"

Amoebozoa_t2 <- Amoebozoa_t2[,-1]
Amoebozoa_t2 <- aggregate(.~tax, FUN="sum", data=Amoebozoa_t2)
Amoebozoa_t2$color <- "Protist-Micropredator"
Amoebozoa_t2 <- Amoebozoa_t2[,c(2:9,1,10)]
Amoebozoa_t2[6,10] <- "Protist"

Amoebozoa_t <- rest[which(rest$supK == 'Amoebozoa'),] 
Amoebozoa_t <- Amoebozoa_t[,c(1:8,13)]
Amoebozoa_t <- aggregate(.~king, FUN="sum", data=Amoebozoa_t)
Amoebozoa_t$tax <- c("other_Amoebozoa", "Amoebozoa (kingdom)","other_Amoebozoa",
                     "other_Amoebozoa","other_Amoebozoa","Myxogastria")
Amoebozoa_t <- Amoebozoa_t[,-1]
Amoebozoa_t <- aggregate(.~tax, FUN="sum", data=Amoebozoa_t)
Amoebozoa_t$color <- "Protist"
Amoebozoa_t <- Amoebozoa_t[,c(2:9,1,10)]
Amoebozoa_t <- Amoebozoa_t[-1,] #otherwise they are 2 times present

Amoebozoa_t <- rbind(Amoebozoa_t,Amoebozoa_t2)
### *******************************************************

THE_TABLE <- rbind(THE_TABLE,Amoebozoa_t)
rest <- rest[which(rest$supK != 'Amoebozoa'),]

### Excavata ****************************************
Excavata_t <- rest[which(rest$supK == 'Excavata'),] 
Excavata_t <- Excavata_t[,c(1:8,13)]
Excavata_t <- aggregate(.~king, FUN="sum", data=Excavata_t) #just to see which are left
table(Excavata_t$king)

Excavata_t2 <- rest[which(rest$supK == 'Excavata'),] 
Excavata_t2 <- Excavata_t2[which(Excavata_t2$king == 'Discoba'),] 
Excavata_t2 <- Excavata_t2[,c(1:8,14)]
Excavata_t2 <- aggregate(.~phy, FUN="sum", data=Excavata_t2) #just to see which are left

Excavata_t <- rest[which(rest$supK == 'Excavata'),] 
Excavata_t <- Excavata_t[,c(1:8,12)]
Excavata_t <- aggregate(.~supK, FUN="sum", data=Excavata_t)
Excavata_t$tax <- "Excavata"
Excavata_t$color <- "Protist"
Excavata_t <- Excavata_t[,-1]
# *******************************************************

THE_TABLE <- rbind(THE_TABLE,Excavata_t)
rest <- rest[which(rest$supK != 'Excavata'),]

### Alveolata ****************************************
Alveolata_t <- rest[which(rest$king == 'Alveolata'),] 
Alveolata_t <- Alveolata_t[,c(1:8,14)]
Alveolata_t <- aggregate(.~phy, FUN="sum", data=Alveolata_t)

Alve_t1 <- Alveolata_t[3,]
Alve_t1$tax <- Alve_t1$phy
Alve_t1$color <- "Protist"
Alve_t1 <- Alve_t1[,-1]

Alve_t2 <- Alveolata_t[c(1,2,4,6:9),]
Alve_t2$tax <- "other_Alveolata"
Alve_t2 <- Alve_t2[,-1]
Alve_t2 <- aggregate(.~tax, FUN="sum", data=Alve_t2)
Alve_t2$color <- "Protist"
Alve_t2 <- Alve_t2[,c(2:9,1,10)]

Alve_t3 <- rest[which(rest$phy == 'Ciliophora'),]
Alve_t3 <- Alve_t3[,c(1:8,15)]
Alve_t3 <- aggregate(.~class, FUN="sum", data=Alve_t3)
Alve_t3$tax <- Alve_t3$class
Alve_t3[1,10] <- "other_Ciliophora"
Alve_t3[2,10] <- "other_Ciliophora"
Alve_t3[4,10] <- "other_Ciliophora"
Alve_t3[6,10] <- "other_Ciliophora"
Alve_t3[8,10] <- "other_Ciliophora"
Alve_t3[9,10] <- "other_Ciliophora"
Alve_t3[10,10] <- "other_Ciliophora"
Alve_t3 <- Alve_t3[,-1]
Alve_t3 <- aggregate(.~tax, FUN="sum", data=Alve_t3)
Alve_t3$color <- c("Protist-Micropredator","Protist-Micropredator","Protist-Micropredator"
                   ,"Protist","Protist-Micropredator")
Alve_t3 <- Alve_t3[,c(2:9,1,10)]

Alveolata_t <- rbind(Alve_t1,Alve_t2,Alve_t3)
### *******************************************************

THE_TABLE <- rbind(THE_TABLE,Alveolata_t)
rest <- rest[which(rest$king != 'Alveolata'),]

### Rhizaria ****************************************
Rhizaria_t <- rest[which(rest$king == 'Rhizaria'),] 
Rhizaria_t <- Rhizaria_t[,c(1:8,14)]
Rhizaria_t <- aggregate(.~phy, FUN="sum", data=Rhizaria_t)

Rhizaria_tx <- rest[which(rest$phy == 'Cercozoa'),] 
Rhizaria_tx <- Rhizaria_tx[,c(1:8,15)]
Rhizaria_tx <- aggregate(.~class, FUN="sum", data=Rhizaria_tx)
Rhizaria_tx$tax <- "Cercozoa"
Rhizaria_tx[13,10] <- "Phytomyxea"
Rhizaria_tx <- Rhizaria_tx[,-1]
Rhizaria_tx <- aggregate(.~tax, FUN="sum", data=Rhizaria_tx)
Rhizaria_tx$color <- c("Protist-Micropredator","Protist")
Rhizaria_tx <- Rhizaria_tx[,c(2:9,1,10)]

Rhizaria_t <- rest[which(rest$king == 'Rhizaria'),] 
Rhizaria_t <- Rhizaria_t[,c(1:8,14)]
Rhizaria_t <- aggregate(.~phy, FUN="sum", data=Rhizaria_t)

Rhizaria_tx2 <- Rhizaria_t[c(1,3),]
Rhizaria_tx2$tax <- c("other_Rhizaria","Foraminifera")
Rhizaria_tx2$color <- c("Protist","Protist-Micropredator")
Rhizaria_tx2 <- Rhizaria_tx2[,-1]

Rhizaria_t <- rbind(Rhizaria_tx,Rhizaria_tx2)

### *******************************************************

THE_TABLE <- rbind(THE_TABLE,Rhizaria_t)
rest <- rest[which(rest$king != 'Rhizaria'),]

### Hacrobia ****************************************
Hacrobia_t <- rest[which(rest$supK == 'Hacrobia'),] 
Hacrobia_t <- Hacrobia_t[,c(1:8,13)]
Hacrobia_t <- aggregate(.~king, FUN="sum", data=Hacrobia_t)

Hacrobia_t <- rest[which(rest$supK == 'Hacrobia'),] 
Hacrobia_t <- Hacrobia_t[,c(1:8,12)]
Hacrobia_t <- aggregate(.~supK, FUN="sum", data=Hacrobia_t)
Hacrobia_t$tax <- Hacrobia_t$supK
Hacrobia_t$color <- "Protist"
Hacrobia_t <- Hacrobia_t[,-1]

### *******************************************************

THE_TABLE <- rbind(THE_TABLE,Hacrobia_t)
rest <- rest[which(rest$supK != 'Hacrobia'),]

### Stramenopiles ****************************************
Stramenopiles_t <- rest[which(rest$king == 'Stramenopiles'),] 
Stramenopiles_t <- Stramenopiles_t[,c(1:8,14)]
Stramenopiles_t <- aggregate(.~phy, FUN="sum", data=Stramenopiles_t)

Stramenopiles_t <- rest[which(rest$king == 'Stramenopiles'),] 
Stramenopiles_t <- Stramenopiles_t[,c(1:8,13)]
Stramenopiles_t <- aggregate(.~king, FUN="sum", data=Stramenopiles_t)
Stramenopiles_t$tax <- Stramenopiles_t$king
Stramenopiles_t$color <- "Protist"
Stramenopiles_t <- Stramenopiles_t[,-1]
# *******************************************************

THE_TABLE <- rbind(THE_TABLE,Stramenopiles_t)
rest <- rest[which(rest$king != 'Stramenopiles'),]

### maybe_the_rest ****************************************
maybe_the_rest_t <- rest
maybe_the_rest_t <- maybe_the_rest_t[,c(1:8,12)]
maybe_the_rest_t <- aggregate(.~supK, FUN="sum", data=maybe_the_rest_t)

maybe_the_rest_t2 <- rest
maybe_the_rest_t2 <- maybe_the_rest_t2[which(maybe_the_rest_t2$supK == 'Opisthokonta'),]
maybe_the_rest_t2 <- maybe_the_rest_t2[,c(1:8,13)]
maybe_the_rest_t2 <- aggregate(.~king, FUN="sum", data=maybe_the_rest_t2)

maybe_the_rest_t <- rest
maybe_the_rest_t <- maybe_the_rest_t[,c(1:8,12)]
maybe_the_rest_t <- aggregate(.~supK, FUN="sum", data=maybe_the_rest_t)

maybe_the_rest_t$tax <- c("unclassified_Eukaryota", "other_Protists","other_Protists",
                          "other_Protists","other_Protists","other_Protists",
                          "other_and_uncl_Opisthokonta","other_Protists","other_Protists","other_Protists")

maybe_the_rest_t <- maybe_the_rest_t[,-1]
maybe_the_rest_t <- aggregate(.~tax, FUN="sum", data=maybe_the_rest_t)
maybe_the_rest_t$color <- c("other", "Protist","other")
maybe_the_rest_t <- maybe_the_rest_t[,c(2:9,1,10)]
# *******************************************************

THE_TABLE <- rbind(THE_TABLE,maybe_the_rest_t)

### a quick test if the data make sense
#plant_t_ <- plant_t
#plant_t_$tax <- plant_t_$king
#plant_t_$color <- "other"
#plant_t_ <- plant_t_[,-1]
#test_t <- rbind(THE_TABLE,plant_t_)
#lala<- colSums(test_t[,1:8]) ### LTWA1: 195 302 reads
#blabla <- colSums(data[,1:8]) ### LTWA1: 198 694 reads
#198694-195302 ### = 3392; 3259 are not assigned on root2 level --> 
#3392-3259 ###133
#78+1+54 ###133 these are the data assigned to mitochondria, chloroplast and main genome
### --> no reads are missing! 
### *******************************************************

### normalize the rRNA data ####

###cleanup
rm(list= ls()[!(ls() %in% c("THE_TABLE"))])
df <- THE_TABLE
###normalization
yyy <- colSums(df[,1:8]) ###get lib sizes
df_n <- df
df_n[,c(-9,-10)] = sweep(df_n[,c(-9,-10)],2,yyy, `/`)
df_n[,c(-9,-10)] = apply(df_n[,c(-9,-10)],2,function(x){x*1000}) #transform to permille
### **************************************************************************

### add the viral data #### (they are normalized to total mRNA reads)
viralT <-read.table("viruses_forhot.txt", sep="\t", header=T, comment='')
colnames(viralT)
viralT$tax_ <- paste(viralT$host_simp,"::",viralT$nuclAc,"::",viralT$tax)
viralT <- viralT[,c(11:26,33)]
viralT <- aggregate(.~tax_, FUN="sum", data=viralT)

###select the top 3 eukaryotic and top 3 prokaryotic
viralT$sum <- rowSums(viralT[2:9])
for (i in 1:length(viralT$sum)){
  viralT$level[i] <- if (viralT$sum[i]<0.1) {"xOther"} else {viralT$tax_[i]} 
}
viralT <- viralT[,c(-1,-18)]
viralT <- aggregate(.~level, FUN="sum", data=viralT)

viralT$tax <- viralT$level
viralT$sum <- rowSums(viralT[2:9])

viralT[6,18]  <- "other_and_unclassified_viruses"
viralT[8,18]  <- "other_and_unclassified_viruses"
viralT[9,18]  <- "other_and_unclassified_viruses"
viralT[10,18]  <- "other_and_unclassified_viruses"
viralT[11,18]  <- "other_and_unclassified_viruses"

viralT <- viralT[,c(-1,-19)]
viralT <- aggregate(.~tax, FUN="sum", data=viralT)
viralT$color <- "virus"
viralT <- viralT[,c(2:9,1,18)]

###finaly change colnames of THE_TABLE (i.e. df_n) and merge it with viralT
colnames(THE_TABLE) <- colnames(viralT)
colnames(df_n) <- colnames(viralT)

THE_TABLE_n <- rbind(df_n,viralT)

### the calculations ####

### calculate mean and log2 etc.
THE_TABLE_n$mean_LTWA <- rowMeans(THE_TABLE_n[,1:4])
THE_TABLE_n$mean_LTWE <- rowMeans(THE_TABLE_n[,5:8])
THE_TABLE_n$mean_LTW <- rowMeans(THE_TABLE_n[,1:8])

THE_TABLE_n$xfold_inc_E <- THE_TABLE_n$mean_LTWE / THE_TABLE_n$mean_LTWA
THE_TABLE_n$log2fold <- log2(THE_TABLE_n$xfold_inc_E)

THE_TABLE_n$tax_ <- paste(THE_TABLE_n$color,"::",THE_TABLE_n$tax)
THE_TABLE_n$size_sqrtPi <- sqrt(THE_TABLE_n$mean_LTW/pi)

### manually edit the table --> export --> edit --> import ####
# write.table(THE_TABLE_n, "THE_TABLE_responder.txt", sep="\t")
write.csv(THE_TABLE_n, "THE_TABLE_responder_MBD.csv")
# import the edited table:
THE_TABLE_n_modif <- read_xlsx("THE_TABLE_responder_modif_MBD.xlsx", sheet = "R_input")
colnames(THE_TABLE_n_modif)
THE_TABLE_n_modif_ <- THE_TABLE_n_modif[which(THE_TABLE_n_modif$include == 'keep'),]
THE_TABLE_n_modif_$tax_ <- paste(THE_TABLE_n_modif_$new_color,"::",
                                 THE_TABLE_n_modif_$new_tax)

# the plot ####
ggplot(THE_TABLE_n_modif_, aes(x=reorder(tax_,new_order), y=log2fold, size=size_sqrtPi,color=new_color))+
  geom_point()+
  geom_hline(yintercept=0,linetype="dashed")+
  guides(size=guide_legend(title="Abundance\n(sqrt(meanAbundance/pi))"),
         color=guide_legend(title="taxonomic group"))+
  #theme_classic()+
  theme(axis.text.x= element_text(angle=90,size=10,colour="black", hjust = 1, vjust = 0.25))

#graph2pdf(file="foodweb_responder.pdf", width=10, height=7)

# add dummy for correct size legend
str(THE_TABLE_n_modif_)
THE_TABLE_n_modif_[37,16] <- "ten"
THE_TABLE_n_modif_[37,17] <- 10
THE_TABLE_n_modif_[37,15] <- 1
THE_TABLE_n_modif_[37,17] <- sqrt(THE_TABLE_n_modif_[37,17]/pi)
THE_TABLE_n_modif_[37,18] <- 0.9

THE_TABLE_n_modif_[38,16] <- "one"
THE_TABLE_n_modif_[38,17] <- 1
THE_TABLE_n_modif_[38,15] <- 1
THE_TABLE_n_modif_[38,17] <- sqrt(THE_TABLE_n_modif_[38,17]/pi)
THE_TABLE_n_modif_[38,18] <- 0.5

THE_TABLE_n_modif_[39,16] <- "zero.one"
THE_TABLE_n_modif_[39,17] <- 0.1
THE_TABLE_n_modif_[39,15] <- 1
THE_TABLE_n_modif_[39,17] <- sqrt(THE_TABLE_n_modif_[39,17]/pi)
THE_TABLE_n_modif_[39,18] <- 0.1

ggplot(THE_TABLE_n_modif_, aes(x=reorder(tax_,new_order), y=log2fold, size=size_sqrtPi,color=new_color))+
  geom_point()+
  geom_hline(yintercept=0,linetype="dashed")+
  guides(size=guide_legend(title="Abundance\n(sqrt(meanAbundance/pi))"),
         color=guide_legend(title="taxonomic group"))+
  theme(axis.text.x= element_text(angle=90,size=10,colour="black", hjust = 1, vjust = 0.25))

# ********************************************************************************************************************

### significant differences? ####
sig_t <- as.data.frame(THE_TABLE_n_modif_[1:36,c(1:8,16)])
sig_t <- as.data.frame(THE_TABLE_n[1:66,c(1:8,16)])
rownames(sig_t) <- sig_t$tax_
sig_t <- sig_t[,-9]

sig_t <- as.data.frame(t(sig_t))
sig_t$group <- c("LTWA","LTWA","LTWA","LTWA","LTWE","LTWE","LTWE","LTWE")

# t_test loop (warning in the end is only because group~group can't be calculated)

taxList <- colnames(sig_t)
for (i in 1:length(taxList)){
  t_test <- t.test(sig_t[,i] ~ group , data = sig_t, var.equal = TRUE)
  print(taxList[i])
  print(t_test)
}

# ***********************************************************************************************************************************************************************
### food web bar charts ####

### input = THE_TABLE_n_responderTable_modif_AS_v1.xlsx sheet: R_input
rm(list=ls(all=TRUE))
THE_TABLE_n_modif <- read_xlsx("THE_TABLE_responder_modif.xlsx", sheet = "R_input")

### Protist :: Amoebozoa (kingdom) NEEDS to be removes as it would be otherwise twice in the table!!
THE_TABLE_n_modif <- THE_TABLE_n_modif[which(THE_TABLE_n_modif$tax_ != 'Protist :: Amoebozoa (kingdom)'),] ###see excel file

# check out the foodweb assignments
table(THE_TABLE_n_modif$new_color)

# Metazoa need to be devided into Metazoa and Arthropoda and AM Fungi need to be tagged
df <- THE_TABLE_n_modif
df$web <- "web"
for(i in 1:length(df$tax)){
  if (grepl("Arthropoda", df$color[i], perl=T) == T){df$web[i] <- "Arthropoda"} else {df$web[i] <- df$new_color[i]}
  if (grepl("AM_Fungi", df$tax[i], perl=T) == T){df$web[i] <- "AM_Fungi"} else {df$web[i] <- df$web[i]}
  if (grepl("Collembola_and_Protura", df$tax[i], perl=T) == T){df$web[i] <- "Collembola_and_Protura"} else {df$web[i] <- df$web[i]}
  if (grepl("Rotifera", df$tax[i], perl=T) == T){df$web[i] <- "Rotifera"} else {df$web[i] <- df$web[i]}
}

# sum up the foodweb data
df <- df[,c(1:8,23)]
df <- aggregate(.~web, FUN="sum", data=df)
#warnings() #these warnings (e.g. Unknown or uninitialised column: 'color'.) can be ignored! if they come
colSums(df[,2:9])

#get a vector with bacteria
bact <- df[4,2:9]
bact_t <- t(bact)
bact <- bact_t[,1]

# divide everything by bacteria
df_ <- df
df_[1:14,c(-1)] = sweep(df_[1:14,c(-1)],2,bact, `/`)

# get mean and SD for LTWA and LTWE samples
df_$meanLTWA <- rowMeans(df_[,2:5])
df_$SD_LTWA <- apply(df_[,2:5], 1, sd)

df_$meanLTWE <- rowMeans(df_[,6:9])
df_$SD_LTWE <- apply(df_[,6:9], 1, sd)


t1 <- df_[,c(1,10,11)]
t1$temp <- "LTWA"
colnames(t1) <- c("web","mean","SD","temp")

t2 <- df_[,c(1,12,13)]
t2$temp <- "LTWE"
colnames(t2) <- c("web","mean","SD","temp")

df <- rbind(t1,t2)

# foodweb bars V2: 
ggplot(data=df, aes(x=web, y=mean, fill=temp)) +
  geom_bar(stat="identity", position=position_dodge())

#graph2pdf(file="foodweb_bars_v2.pdf", width=10, height=5) 
#in illustrator: Bact were made 40x40 and the rest was also set to 40x40
#viruses are not correct here!

# excl bacteria
df_exclB <- df[which(df$web != 'Bacteria'),]
ggplot(data=df_exclB, aes(x=web, y=mean, fill=temp)) +
  geom_bar(stat="identity", position=position_dodge())

# ********************************************************************************************************************

### food web boxes to match bar charts ####

rm(list=ls(all=TRUE))
THE_TABLE_n_modif <- read_xlsx("THE_TABLE_responder_modif.xlsx", sheet = "R_input")

# check out the foodweb assignments
table(THE_TABLE_n_modif$new_color)

# Split Metazoa into Metazoa and Arthropoda and identify AM Fungi 
df <- THE_TABLE_n_modif
df$web <- "web"
for(i in 1:length(df$tax)){
  if (grepl("Arthropoda", df$color[i], perl=T) == T){df$web[i] <- "Arthropoda"} else {df$web[i] <- df$new_color[i]}
  if (grepl("AM_Fungi", df$tax[i], perl=T) == T){df$web[i] <- "AM_Fungi"} else {df$web[i] <- df$web[i]}
  if (grepl("Collembola_and_Protura", df$tax[i], perl=T) == T){df$web[i] <- "Collembola_and_Protura"} else {df$web[i] <- df$web[i]}
  if (grepl("Rotifera", df$tax[i], perl=T) == T){df$web[i] <- "Rotifera"} else {df$web[i] <- df$web[i]}
}

# sum up the foodweb data
df <- df[,c(1:8,23)]
df <- aggregate(.~web, FUN="sum", data=df)
colSums(df[,2:9])

#get a vector with bacteria
bact <- df[4,2:9]
bact_t <- t(bact)
bact <- bact_t[,1]

# divide everything by bacteria
df_ <- df
df_[1:14,c(-1)] = sweep(df_[1:14,c(-1)],2,bact, `/`)

# get mean and SD for LTWA and LTWE samples
df_$meanLTWA <- rowMeans(df_[,2:5])
df_$SD_LTWA <- apply(df_[,2:5], 1, sd)

df_$meanLTWE <- rowMeans(df_[,6:9])
df_$SD_LTWE <- apply(df_[,6:9], 1, sd)

df_$LTWA_mean_SD <- df_$meanLTWA + df_$SD_LTWA
df_$LTWE_mean_SD <- df_$meanLTWE + df_$SD_LTWE

#BUT I want a figure with the overall mean:
df_$meanLTWAE <- rowMeans(df_[,2:9])
df_$SD_LTWAE <- apply(df_[,2:9], 1, sd)
df_$LTWAE_mean_SD <- df_$meanLTWAE + df_$SD_LTWAE

#and I need dummies for the scale
df_[15,] <- 0.1
df_[16,] <- 0.01
df_[17,] <- 0.001

t1 <- df_[,c(1,10,11,14)]
t1$temp <- "LTWA"
colnames(t1) <- c("web","mean","SD","mean_plus_SD","temp")

t2 <- df_[,c(1,12,13,15)]
t2$temp <- "LTWE"
colnames(t2) <- c("web","mean","SD","mean_plus_SD","temp")

t3 <- df_[,c(1,16:18)]
t3$temp <- "LTWA_and_E"
colnames(t3) <- c("web","mean","SD","mean_plus_SD","temp")

df <- rbind(t1,t2,t3)

#unscaled figure
p1<-ggplot(df, aes(x=temp, y=web, size=mean_plus_SD,color=web))+
  geom_point(shape=15)
p1

p2<-ggplot(df, aes(x=temp, y=web, size=mean_plus_SD))+
  geom_point(shape=15)
p2

p3<- p2+ geom_point(data=df, aes(x=temp, y=web, size=mean,color=web),shape=15)
p3


#scale by size area
p1<-ggplot(df, aes(x=temp, y=web, size=mean_plus_SD,color=web))+
  geom_point(shape=15)+
  scale_size_area()
p1

p2<-ggplot(df, aes(x=temp, y=web, size=mean_plus_SD))+
  geom_point(shape=15)+
  scale_size_area()
p2

p3<- p2+ geom_point(data=df, aes(x=temp, y=web, size=mean,color=web),shape=15)
p3

# ********************************************************************************************************************

### the SD deviation is so hard to see, thus only the mean is used in the Figure

### significant diff. in foodweb boxes? ####

df_2 <- df_[1:14,1:9]

#transform table
rownames(df_2) <- df_2$web
df_2 <- df_2[,-1]

df_t <- as.data.frame(t(df_2))

df_t$group <- c("LTWA","LTWA","LTWA","LTWA","LTWE","LTWE","LTWE","LTWE")

# t-test 
df_t <- df_t[,-4]
taxList <- colnames(df_t)
for (i in 1:length(taxList)){
  t_test <- t.test(df_t[,i] ~ group ,var.equal = TRUE, data = df_t)
  print(taxList[i])
  print(t_test)
}

# ********************************************************************************************************************
# Supplement figure S4: Schematic overview of community profile 

# load packages
library(ggalluvial)
library(ggrepel)
library(reshape2)

data <- read.csv("forhot_200K_trim.csv", header = T, row.names = 1)

#### Domain overview
data_dom <- aggregate(.~dom, FUN="sum", data=data[, c(1:8,11)]) 
data_dom$sum <- rowSums(data_dom[2:9])
row.names(data_dom) <- as.character(data_dom$dom)
data_dom <-data_dom[c(-1,-5,-6),c(-1,-10)] #  Eukaryota (Chloroplast), Eukaryota (Mitochondrion) and unclassifed

### showing the fungal fraction in domain plot
dataeu <- data[which(data$dom =="Eukaryota"),] 
data_fun <- dataeu[which(dataeu$king =="Fungi"),]
data_fun <- aggregate(.~king, FUN="sum", data=data_fun[, c(1:8,13)])
row.names(data_fun) <- data_fun$king
data_fun <- data_fun[,c(-1)]
data_fun <- t(data_fun)

data_dom <- t(data_dom)
data_dom <- cbind(data_dom, data_fun)
data_dom <- as.data.frame(data_dom)
data_dom$OtherEukaryota <- data_dom$Eukaryota - data_dom$Fungi
data_dom <- data_dom[,c(-3)]

data_dom <- as.data.frame(t(data_dom))
data_dom$dom <- row.names(data_dom)

sum_dom <- sum(data_dom[,1:8])
sum_dom/sum_all*100 #98.27% of data

### Bacerial overview on phyla-level ###
#--------------------------------------#
databac <- data[which(data$dom =="Bacteria"),] 
data_bac <- aggregate(.~phy, FUN="sum", data=databac[, c(1:8,14)]) 
data_bac$sum <- rowSums(data_bac[2:9])
uncl_bac <- sum(data_bac$sum[which(data_bac$phy =="")])/sum_dom*100 #  7.770458
data_bac <- data_bac[data_bac$phy!="",]
row.names(data_bac) <- as.character(data_bac$phy)
data_bac <- data.frame(sum=rowSums(data_bac[,2:9]), data_bac[,2:9])

data_bac$tax <-row.names(data_bac)
for (i in 1:length(data_bac$sum)){
  data_bac$tax[i] <- if (data_bac$sum[i]<50000) {"xBac_Other"} else {data_bac$tax[i]} 
}
data_bac <- aggregate(.~tax, FUN="sum", data=data_bac)
data_bac <-data_bac[,-2]

### resolve proteobacteria
databac <- data[which(data$phy =="Proteobacteria"),] 
data_pro <- aggregate(.~class, FUN="sum", data=databac[, c(1:8,15)]) 
data_pro$tax <-data_pro$class
data_pro$tax[1]<-"xBac_Other"
data_bac <- rbind(data_bac,data_pro[,c(10,2:9)])
data_bac <- data_bac[which(data_bac$tax!="Proteobacteria"),] #d
data_bac$tax <- gsub("Betaproteobacteria", "xBac_Other",data_bac$tax)
data_bac$tax <- gsub("Gammaproteobacteria", "xBac_Other",data_bac$tax)
data_bac$tax <- gsub("Zetaproteobacteria", "xBac_Other",data_bac$tax)
data_bac$tax <- gsub("Ca. Tenderia class incertae sedis", "xBac_Other",data_bac$tax)
data_bac$tax <- gsub("Ca. Thiobios class incertae sedis", "xBac_Other",data_bac$tax)

data_bac$tax <- gsub(" ", "",data_bac$tax)
data_bac <- aggregate(.~tax, FUN="sum", data=data_bac)
data_bac$cat <- 'Bacteria'
data_bac$count <- 1:nrow(data_bac) 

### Fungi overview on class-level ###
#-----------------------------------#
data_fun <- data[which(data$king =="Fungi"),] 
data_fun <- aggregate(.~class, FUN="sum", data=data_fun[, c(1:8,15)]) 
data_fun$sum <- rowSums(data_fun[2:9])
data_fun <- data_fun[data_fun$class!="",]
row.names(data_fun) <- as.character(data_fun$class)
data_fun <- data.frame(sum=rowSums(data_fun[,2:9]), data_fun[,2:9])
data_fun$tax <-row.names(data_fun)

for (i in 1:length(data_fun$sum)){
  data_fun$tax[i] <- if (data_fun$sum[i]<3500) {"xFun_Other"} else {data_fun$tax[i]} 
}
data_fun <- aggregate(.~tax, FUN="sum", data=data_fun)
data_fun <-data_fun[,-2]
data_fun$cat <- 'Fungi'
data_fun$count <- 1:nrow(data_fun) 

### Protist overview on kingdom-level ###
#---------------------------------------#
dataeu <- data[which(data$dom =="Eukaryota"),] 
data_prot <- dataeu[which(dataeu$king!="Fungi"),]
data_prot <- data_prot[which(data_prot$king!="Plantae"),]
data_prot <- data_prot[which(data_prot$king!="Metazoa (Animalia)"),]
data_prot$phy <- as.character(data_prot$phy)
data_prot$supK <- as.character(data_prot$supK)
prot <- data_prot[,c(1:8,13)]
data_prot <- aggregate(.~king, FUN="sum", data=prot) 
data_prot$sum <- rowSums(data_prot[2:9])
data_prot <- data_prot[data_prot$king!="",]
row.names(data_prot) <- as.character(data_prot$king)
data_prot <- data.frame(sum=rowSums(data_prot[,2:9]), data_prot[,2:9])
data_prot$tax <-row.names(data_prot)
for (i in 1:length(data_prot$sum)){
  data_prot$tax[i] <- if (data_prot$sum[i]<300) {"xPro_Other"} else {data_prot$tax[i]} 
}
data_prot <- aggregate(.~tax, FUN="sum", data=data_prot)
data_prot <-data_prot[,-2]
data_prot$cat <- 'protists'
data_prot$count <- 1:nrow(data_prot) 

### Metazoa overview on phyla-level ###
#-------------------------------------#
data_met <- data[which(data$king =="Metazoa (Animalia)"),] 
data_met <- aggregate(.~phy, FUN="sum", data=data_met[, c(1:8,14)]) 
data_met$sum <- rowSums(data_met[2:9])
data_met <- data_met[data_met$phy!="",]
row.names(data_met) <- as.character(data_met$phy)
data_met <- data.frame(sum=rowSums(data_met[,2:9]), data_met[,2:9])
data_met$tax <-row.names(data_met)

for (i in 1:length(data_met$sum)){
  data_met$tax[i] <- if (data_met$sum[i]<150) {"xmet_Other"} else {data_met$tax[i]} 
}
data_met <- aggregate(.~tax, FUN="sum", data=data_met)
data_met <-data_met[,-2]
data_met$cat <- 'Metazoa'
data_met$count <- 1:nrow(data_met) 

### Archeae overview on kingdom-level ###
#--------------------------------------#
data_ar <- data[which(data$dom =="Archaea"),] 
data_ar <- aggregate(.~king, FUN="sum", data=data_ar[, c(1:8,13)]) 
data_ar$sum <- rowSums(data_ar[2:9])
data_ar <- data_ar[data_ar$king!="",]
row.names(data_ar) <- as.character(data_ar$king)
data_ar <- data.frame(sum=rowSums(data_ar[,2:9]), data_ar[,2:9])
data_ar$tax <-row.names(data_ar)
data_ar <-data_ar[,-1]
data_ar$cat <- 'Archaea'
data_ar$count <- 1:nrow(data_ar) 

### Combined overview - Run AFTER the generation of the individual plots, skip down to 'plotting' to complete ###
#-------------------------------------------------------------------------------------------------------------#
data_ar <- data[which(data$dom =="Archaea"),]
data_ar <- aggregate(.~dom, FUN="sum", data=data_ar[, c(1:8,11)]) 
names(data_ar)[1] <- 'tax'

data_met <- data[which(data$king =="Metazoa (Animalia)"),] 
data_met <- aggregate(.~king, FUN="sum", data=data_met[, c(1:8,13)]) 
names(data_met)[1] <- 'tax'

dataeu <- data[which(data$dom =="Eukaryota"),] 
data_prot <- dataeu[which(dataeu$king!="Fungi"),]
data_prot <- data_prot[which(data_prot$king!="Plantae"),]
data_prot <- data_prot[which(data_prot$king!="Metazoa (Animalia)"),]
data_prot <- aggregate(.~dom, FUN="sum", data=data_prot[, c(1:8,11)]) 
data_prot$dom <- 'Protist'
names(data_prot)[1] <- 'tax'

data_fun <- data[which(data$king =="Fungi"),] 
data_fun <- aggregate(.~king, FUN="sum", data=data_fun[, c(1:8,13)])
names(data_fun)[1] <- 'tax'

databac <- data[which(data$dom =="Bacteria"),] 
data_bac <- aggregate(.~dom, FUN="sum", data=databac[, c(1:8,11)]) 
names(data_bac)[1] <- 'tax'

data_plants <- data[which(data$king=="Plantae"),]
data_plants <- aggregate(.~king, FUN="sum", data=data_plants[, c(1:8,13)])
names(data_plants)[1] <- 'tax'

data_other <- data[which(data$dom!="Eukaryota"),]
data_other <- data_other[which(data_other$dom!="Bacteria"),]
data_other <- data_other[which(data_other$dom!="Archaea"),]
data_other <- aggregate(.~root, FUN="sum", data=data_other[, c(1:9)])
names(data_other)[1] <- 'tax'

#### Combine the data frames ###
all <- rbind(data_ar,data_met)
all <- rbind(all, data_prot)
all <- rbind(all, data_fun)
all <- rbind(all, data_bac)
all <- rbind(all, data_plants)
all <- rbind(all, data_other)
all$cat <- 'all'

### Plotting   

### Go through one by one 
df_m <- melt(data_bac, id=c("tax",'cat', 'count'))
df_m <- melt(data_fun, id=c("tax",'cat', 'count'))
df_m <- melt(data_met, id=c("tax",'cat', 'count'))
df_m <- melt(data_prot, id=c("tax",'cat', 'count'))
df_m <- melt(data_ar, id=c("tax",'cat', 'count'))
df_m <- melt(all, id=c("tax",'cat'))

# sum across samples
df_sum <- aggregate(value~tax, FUN="sum", data=df_m)
df_sum$cat <- lookup(df_sum$tax, unique(df_m[1:2]),df_sum$cat)

### Transform to relative abundance
sum(df_sum$value)/sum_dom 
df_sum$pct <- df_sum$value / sum(df_sum$value)*100

### Sorting procentage and making labels for plot
df_sum <- df_sum[order(df_sum$pct),]
df_sum$sort <- paste(df_sum$tax, "(", round(df_sum$pct,1), "%)", sep = " ")
df_sum$sort <- factor(df_sum$sort, levels = df_sum$sort)

### Introduce dummy-variable to get nice curves and spacers for sankey plot
df_sum2 <- df_sum
df_sum2$pct <- 10
df_sum2$cat <- 'dummy'
df_sum <- rbind(df_sum, df_sum2)

### Sankey diagram ###
#--------------------#
p <- ggplot(df_sum,
            aes(y = pct, axis1 = cat, axis2 = sort)) +
  geom_flow(aes(fill = cat), width = 1/10, alpha = 1) +
  geom_stratum(width = 1/10, fill = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("cat", "tax"), expand = c(.05, .05)) +
  scale_fill_manual(values =  c("blue", "red"), )+
  ggtitle("") +theme_minimal()+  theme(legend.position = 'none', panel.grid = element_blank())
p 

# ********************************************************************************************************************
### Supplement Figure S4: Quantification rRNA data

ls(all.names=TRUE)
rm(list=ls(all=TRUE))

### Data
data<-read.table("forhot_200K_trim_GO.txt", sep="\t", header=T, quote=NULL, comment='')
meta <- read.csv("meta.csv",header = TRUE, sep = ";")

### get quant. factors 
quanT<-read.table("quantific_table.txt", sep="\t", header=T, quote=NULL, comment='',row.names = 1) # this tab contains various quantification values (WW, DW, vaious protein length, rRNA)
quanT <- quanT[,c(1:8)]
q18 <- as.numeric(as.vector(quanT["Eukary_1800_rRNA",]))
q16 <- as.numeric(as.vector(quanT["Prokary_1500_rRNA",]))
s <- c("60759","60773","60774","60775","60776","60777","60778","60779")

### make dataframe for 16s and 18S
quanty <- as.data.frame(s)
Q18S <- quanty
Q18S$q <- q18
Q16S <- quanty
Q16S$q <- q16

### quantified data
qdata <- data
qdata$dom <- as.character(qdata$dom)
sample_order <- names(qdata[1:8]) ## get the order of the samples 

# Sort the 'calculation factor' data frame to match the sample data frame
Q16S<- Q16S[order(match(Q16S$s, sample_order)), ] 
Q18S<- Q18S[order(match(Q18S$s, sample_order)), ]

### Quantify
for (i in 1:8){
  qdata[[i]] <- ifelse(qdata$dom =="Bacteria", qdata[[i]] * Q16S$q[i], qdata[[i]])
  qdata[[i]]<- ifelse(qdata$dom =="Archaea", qdata[[i]] * Q16S$q[i], qdata[[i]])
  qdata[[i]] <- ifelse(qdata$dom =="Eukaryota", qdata[[i]] * Q18S$q[i], qdata[[i]])
}
qdata[] <- lapply(qdata, function(x) type.convert(as.character(x)))

####  Quantified domain overview 
data_dom <- aggregate(.~dom, FUN="sum", data=qdata[c(1:8,11)]) 
data_dom <-  data_dom[c(2:4),] 

### showing the fungal fraction in domain plot
dataeu <- qdata[which(data$dom =="Eukaryota"),] 
data_fun <- dataeu[which(dataeu$king =="Fungi"),]
data_fun <- aggregate(.~king, FUN="sum", data=data_fun[c(1:8,13)])
row.names(data_fun) <- data_fun$king
data_fun <- data_fun[,c(-1)]
row.names(data_dom) <- data_dom$dom
data_dom <- data_dom[,c(-1)]
data_dom <- rbind(data_dom, data_fun)
data_dom <- t(data_dom)
data_dom <- as.data.frame(data_dom)
data_dom$OtherEukaryota <- data_dom$Eukaryota - data_dom$Fungi
data_dom <- data_dom[,c(-3)]
# --> --> stats below

### plot individual bars
dfV_t <- data_dom
dfV_t$id <- meta$ID
dfV_t$site <- meta$category
dfV_t$x <- meta$samples

df_m <- melt(dfV_t, id=c("x", 'id', 'site'))
ggplot(df_m, aes(x,value, fill=variable)) +
  geom_col() +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# plot summary for LTW-A and LTW-E
data_dom <- as.data.frame(t(data_dom))
data_dom$dom <- row.names(data_dom)
data_dom_ <- data_dom
data_dom_$LTWA <- rowMeans(data_dom_[,1:4])
data_dom_$LTWE <- rowMeans(data_dom_[,5:8])
data_dom_ <- data_dom_[,9:11]
df_m <- melt(data_dom_, id=c("dom"))
colnames(df_m)

p1 <- ggplot() + geom_bar(aes(y = value, x = variable, fill = dom), data = df_m,stat="identity") +
  theme_bw() +theme(axis.text.x= element_text(angle=45, hjust = 1.3, vjust = 1.2), legend.title = element_blank())+
  xlab("") +ylab("Proportion of sequences") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p1

# ********************************************************************************************************************
### Significant differences?
stats_df <- as.data.frame(data_dom)

t.test(Archaea ~ cat , data = stats_df, var.equal = TRUE)
t.test(Bacteria ~ cat , data = stats_df, var.equal = TRUE)
t.test(Fungi ~ cat , data = stats_df, var.equal = TRUE)
t.test(OtherEukaryota ~ cat , data = stats_df, var.equal = TRUE)

#**************************************************************************************************************************************
### Supplement S5 Correlation analysis ###

#-------------- NOT in function ------------#
'%!in%' <- function(x,y)!('%in%'(x,y))
#-------------------------------------------#
options(scipen=999) # avoid scientific numbers

## variables
var_basic <- meta[,c(19,20,25, 32:35 )] #"C.N.mic"  "DOC..µg.g.DW."  "soil.pH" "Ntot..mg.g.DW." "C.N.soil" "BD_0.5"    "BD_5.10" 
## NMDS axes
var_basic <- cbind(site=as.factor(meta$temperature), var_basic)
var_basic <- cbind(NMDS1=samples2d$MDS1, var_basic)
var_basic <- cbind(NMDS2=samples2d$MDS2, var_basic)

#### predatory myxobacteria
Hali <- data_fam[which(data_fam$fam == 'Haliangiaceae'),]
Sanda <- data_fam[which(data_fam$fam == 'Sandaracinaceae'),]
Bdello <- data_fam[which(data_fam$fam == 'Bacteriovoracaceae'),]

myxo <-rbind(Hali, Sanda, Bdello)
myxo <- as.data.frame(t(myxo))
myxo$myxo <- rowSums(myxo)
names(myxo) <- lapply(myxo[1, ], as.character)
myxo <- myxo[-1,] 

var_basic <- cbind(var_basic, myxo$myxo)

### correlation NMDS1 and NMDS2 (exchange)
combs <- do.call(cbind.data.frame,
                 lapply("NMDS1", rbind, combn(names(var_basic)[names(var_basic) %!in% c("NMDS1", "site")], 1)))
combs

res.l <- lapply(combs, function(x) {
  `attr<-`(cor(var_basic[,x[1]], var_basic[,x[2]]),
           "what", {
             paste0(x[1], ", ", paste(x[2]))})
})
res <- setNames(unlist(res.l), sapply(res.l, attr, "what"))
res

# get p values run cor.text
p.l <- lapply(combs, function(x) {
  `attr<-`(cor.test(var_basic[,x[1]], var_basic[,x[2]])$p.value,
           "what", {
             paste0(x[1], ", ", paste(x[2]))})
})
pres <- setNames(unlist(p.l), sapply(res.l, attr, "what"))
pres

## myxo test removing outlier 
cor.test(var_basic$`myxo$myxo`[-6], var_basic$BD_5.10[-6], method = "pearson")

## plotting correlations for supplement S5

#  Function to add correlation coefficients 
#*************************************************************************
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  Cor <- abs(cor(x, y)) # Remove abs function if desired
  txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
  if(missing(cex.cor)) {
    cex.cor <- 0.4 / strwidth(txt)
  }
  text(0.5, 0.5, txt,
       cex = 1 + cex.cor * Cor) # Resize the text by level of correlation
}
#*************************************************************************
# functuion for liniaer regression line 
# made to draw regression line instead of lowess smooth line
reg <- function(x, y, ...) {
  points(x,y, ...)
  abline(lm(y~x)) 
}
#*************************************************************************

# Plotting the correlation matrix
pairs(var_basic[,c(1:2, 4:11)], 
      upper.panel = panel.cor,    # Correlation panel
      lower.panel = reg, #  linier regression lines
      labels = colnames(var_basic[,c(1:2, 4:11)]),  # Variable names
      pch = 21,                 # Pch symbol
      bg = rainbow(2)[var_basic$site],  # Background color of the symbol (pch 21 to 25)
      col = rainbow(2)[var_basic$site], # Border color of the symbol
      main = "",    # Title of the plot
      row1attop = TRUE,         # If FALSE, changes the direction of the diagonal
      gap = 1,                  # Distance between subplots
      cex.labels = NULL,        # Size of the diagonal text
      font.labels = 1)

### Significant response to temperature
var_basic <- meta[,c(3,4)] 

adonis2(data_t ~temperature, var_basic, permutations = 999, method = "bray", by = NULL)

