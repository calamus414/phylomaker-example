library("V.PhyloMaker2")
library("nlme")
#create tree by V.PhyloMaker2
specieslist<-read.csv(url('https://raw.githubusercontent.com/calamus414/phylomaker-example/main/climbingplant_species.csv'))
climbing_plant_phylo <- phylo.maker(sp.list = specieslist, tree = GBOTB.extended.TPL, nodes = nodes.info.1.TPL, scenarios = "S3")

is.ultrametric(climbing_plant_phylo$scenario.3)#check each branch length have >0, which mean the tree is  dichotomous 
is.rooted(climbing_plant_phylo$scenario.3)
plot(climbing_plant_phylo$scenario.3)#visualize the tree
#the basal node is not polytomy although it looks like that


#load the trait data
trait<-read.csv(url('https://raw.githubusercontent.com/calamus414/phylomaker-example/main/trait.csv'),row.names = 1)
sunnpq<-trait$sunnpq
index90<-trait$index90
#Give names to the vectors
names(sunnpq)<-rownames(trait)
names(index90)<-rownames(trait)
sort(climbing_plant_phylo$scenario.3$tip.label)==sort(rownames(trait))#make sure that the species name for trait is the same as in tree

#PIC
sunnpqcontrast<-pic(sunnpq,climbing_plant_phylo$scenario.3)
index90contrast<-pic(index90,climbing_plant_phylo$scenario.3)
summary(lm(index90contrast~sunnpqcontrast-1))#p value is 0.03632, slope is -3.984 


#PGLS
bmcorrclimb<-corBrownian(phy = climbing_plant_phylo$scenario.3)
sunnpq90<-gls(index90~sunnpq,data = trait,correlation = bmcorrclimb)
summary(sunnpq90)#p value is 0.0372 slope is -3.964459

#i forget to sort the trait data into the right order, so PIC and PGLS get different result
trait<-trait[climbing_plant_phylo$scenario.3$tip.label,]

#PIC
sunnpqcontrast<-pic(sunnpq,climbing_plant_phylo$scenario.3)
index90contrast<-pic(index90,climbing_plant_phylo$scenario.3)
summary(lm(index90contrast~sunnpqcontrast-1))#p value is 0.03632, slope is -3.984

#PGLS
bmcorrclimb<-corBrownian(phy = climbing_plant_phylo$scenario.3)
sunnpq90<-gls(index90~sunnpq,data = trait,correlation = bmcorrclimb)
summary(sunnpq90)#result change and PGLS have the same result with PIC after sorting the trait data


#example that PIC an PGLS have same result from https://github.com/simjoly/CourseComparativeMethods
##seedplant data is from https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.1456, and reupload to my github

seedplantstree <- read.nexus("https://raw.githubusercontent.com/calamus414/phylomaker-example/main/phylo.tre")
seedplantsdata <- read.csv2("https://raw.githubusercontent.com/calamus414/phylomaker-example/main/trait_github.csv",sep=',')
# Remove species for which we don't have complete data
seedplantsdata <- na.omit(seedplantsdata)
# Remove species in the tree that are not in the data matrix
species.to.exclude <- seedplantstree$tip.label[!(seedplantstree$tip.label %in% 
                                                   seedplantsdata$Code)]
seedplantstree <- drop.tip(seedplantstree,species.to.exclude)
rm(species.to.exclude)

plot(seedplantstree)

rownames(seedplantsdata) <- seedplantsdata$Code
seedplantsdata <- seedplantsdata[,-1]#remove species code

seedplantsdata <- seedplantsdata[seedplantstree$tip.label,]#sort the dataframe
seedplantsdata$maxH<-as.numeric(seedplantsdata$maxH)
seedplantsdata$Wd<-as.numeric(seedplantsdata$Wd)
seedplantsdata$Sm<-as.numeric(seedplantsdata$Sm)
seedplantsdata$Shade<-as.numeric(seedplantsdata$Shade)
seedplantsdata$N<-as.numeric(seedplantsdata$N)

Wd <- seedplantsdata$Wd
Shade <- seedplantsdata$Shade
Sm <- seedplantsdata$Sm
N <- seedplantsdata$N
#Give names to the vectors
names(Wd) <- names(Shade) <- names(Sm) <- names(N) <- row.names(seedplantsdata)

#PIC
#default in argument scale is true, so it is not necessary
Wd.contrast <- pic(Wd,seedplantstree,scaled=TRUE)
Shade.contrast <- pic(Shade,seedplantstree,scaled=TRUE)

RegressShade.pic <- lm(Shade.contrast~Wd.contrast -1)
summary(RegressShade.pic)#p-value is 0.01273 slope is 4.361

#PGLS
bm.corr <- corBrownian(phy=seedplantstree)
shade.bm.pgls <- gls(Shade ~ Wd, data = seedplantsdata, correlation = bm.corr)
summary(shade.bm.pgls)#p-value is 0.0127 slope is 4.361



