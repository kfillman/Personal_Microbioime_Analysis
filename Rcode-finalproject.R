setwd("~/MICE5035/Final_Project")

library("Maaslin2")
library("vegan")
library("car")
library("ape")
library("vegan")
library("dplyr")
library("hagis")
library("ggplot2")
library ("ggpubr")
library("multcompView")
library("labdsv")


data = read.csv("taxa_table_L7.txt", header = TRUE, sep = "\t", row.names = 1)
map= read.csv("meta.txt", header = TRUE, sep = "\t", row.names = 1)
bray = read.delim('bray_curtis_taxa_table_L7_final.txt', row.names = 1)
alpha= read.delim('alpha-diversity.txt', row.names = 1)

jacca <- read.delim('binary_jaccard_taxa_table_L7_final.txt', row.names = 1)

data1 <- data[, -ncol(data)] #Remove last column 
data1[] <- lapply(data1, as.numeric) # Convert as numeric

data2=t(data1) #Transpose the taxa

taxa_relab<- decostand(data2, "total") * 100 # Convert to relative abundance 


common.rownames <- intersect(rownames(map),rownames(taxa_relab)) # Check the column names 

taxa_relab <- taxa_relab[common.rownames,]


# Fit Maasslin package  

fit_data <- Maaslin2(
  taxa_relab, 
  map,
  min_prevalence = 0,
  output = 'FrozenMasslin',
  fixed_effects = c('num_frozen_foods'),
  standardize = FALSE,
)


fit_data <- Maaslin2(
  taxa_relab, 
  map,
  min_prevalence = 0,
  output = 'AlMasslin',
  fixed_effects = c('serv_alcohol'),
  standardize = FALSE,
)



###beta diversity Frozen ##################


pc <- cmdscale(bray)


plot(pc[,1],pc[,2],pch=3)



map$num_frozen_foods <- factor(map$num_frozen_foods,levels=c('1','2','3','4','5','6', '7', '8', '9', '10'))

colors <- c("red", "blue", "green", "yellow", "orange", 
            "purple", "brown", "pink", "gray", "black")

plot(pc[,1],pc[,2],pch=c(1, 2, 3, 4, 5, 6)[map$num_frozen_foods])


plot(pc[,1], pc[,2], pch=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)[map$num_frozen_foods], col=colors[map$num_frozen_foods])
legend('topleft', levels(map$num_frozen_foods), pch=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), col=colors, cex=.75)


dataEllipse(x=pc[,1], y=pc[,2], groups=map$num_frozen_foods,plot.points=FALSE,levels=0.40,robust=TRUE,col=colors,segments=100)



################Alpha Frozen ####################


boxplot(alpha$shannon ~ map$num_frozen_foods,las=2, col=colors, xlab='')



####ANCOM-BC analysis#####


library("ANCOMBC")


smd = S4Vectors::DataFrame(map)
assays = S4Vectors::SimpleList(counts = t(taxa_relab))

tse = TreeSummarizedExperiment::TreeSummarizedExperiment(assays=assays, colData = smd)

output.generation.bmi <- ancombc2(data = tse, assay_name = "counts", 
                                  fix_formula = "num_frozen_foods", lib_cut = 0,
                                  verbose = TRUE)


####ANCOM-BC analysis Alcohol#####


library('ANCOMBC')
library('TreeSummarizedExperiment')

smd = S4Vectors::DataFrame(map)
assays = S4Vectors::SimpleList(counts = t(taxa_relab))

tse = TreeSummarizedExperiment::TreeSummarizedExperiment(assays=assays, colData = smd)

output.generation.bmi <- ancombc2(data = tse, assay_name = "counts", 
                                  fix_formula = "serv_alcohol_Cat", lib_cut = 0,
                                  verbose = TRUE)
