library(tidyverse)

###1. Get our data.
####dataset one. Need to do some data cleaning/filtering before we can use it.
nsf<-read.csv("nsf0607_spec.csv")
nsfnet<-filter(nsf, method=="net")
table(nsfnet$site)

nsfnet2<-unite(nsfnet, genus.species, genus, species)
nsfnet.no2007<-droplevels(filter(nsfnet2, year!=2007 & site_type!="reference" & site_type!="matrix"))
length(unique(nsfnet.no2007$site))

###dataset two. some data cleaning/filtering before we can use it.
cig<-read.csv("cig1115_spec.csv")
cig<-droplevels(filter(cig, round!="NULL" & year=="2013"))
cig2<-unite(cig, genus.species, genus, species)
cig3<-unite(cig2, plant_genus_sp, plant_genus, plant_species)
table(cig$site)
unique(cig3$site)

write.csv(cig3, "cig3 for dataclub.csv")

##2. Format data for analysis

##filter to get records for just one site in one dataset
site<-filter(cig3, site=="DR")

##each row in "site" corresponds to one insect specimen. The column "genus.species" indicates
#the species of the insect. The column "plant_genus_sp" indicates the flowers species that the specimen
#was visiting when the specimen was caught. We want to know how many visits each flower species
#at the site recieved from each insect species. 

interactions<-site%>%
  group_by(genus.species, plant_genus_sp)%>%
  summarize(visits=length(uniqueID))

##the below code turns "interactions" into an insect species by plant species matrix
matrix<-spread(interactions, genus.species, as.numeric(visits))
View(matrix)

##Next need to do some formatting to get the matrix into shape for the analysis.

matrix[is.na(matrix)]<-0
matrix2<-data.matrix(matrix)
matrix3<-matrix2[,-1]
row.names(matrix3)<-c(as.matrix(matrix[,1]))
View(matrix3)


#### 3. Used published code (adapted from publicly available function "nestednodf") to get species-level nestedness values for each plant and pollinator species
###at our focal site.


#This first bit sorts our matrix by decreasing row and column fill (number of bees a plant interacts with, number of plants a bee interacts with)
comm<-matrix3

bin.comm <- ifelse(comm > 0, 1, 0)
rfill <- rowSums(bin.comm)
cfill <- colSums(bin.comm)

rgrad <- rowSums(comm)
cgrad <- colSums(comm)

rorder <- order(rfill, rgrad, decreasing = TRUE)
corder <- order(cfill, cgrad, decreasing = TRUE)
 
comm <- comm[rorder, corder]
rfill <- rfill[rorder]
cfill <- cfill[corder]

#This bit gets some info from our interaction matrix that will be used in the calculation of the WNODF values
#(number of rows and columns, and the "fill" of the matrix). 

nr <- NROW(comm)
nc <- NCOL(comm)
fill <- sum(rfill)/prod(dim(comm))

#we then use those values to create two vectors of 0s, each the length of the number of 
#unique pairs of rows and columns, respectively. The loops below will replace these zeros with the pairwise
##nestedness values of each pair of rows and each pair of columns, respectively. So "N.paired.rows" will
##ultimately become an object that stores pairwise nestedness of rows(aka pairwise nestedness of plant species), 
#and the same for "N.paired.cols" for columns (aka pairwise nestedness of insect species) 
N.paired.rows <- numeric(nr * (nr - 1)/2)
N.paired.cols <- numeric(nc * (nc - 1)/2)


###CALCULATE NESTEDNESS OF ALL PAIRS OF ROWS###
#The line of code below is just an empty data frame in which to store nestedness values of each pair of rows i and j, 
#plus the total abundance of focal row species in the network. 

thegoods.rows<-data.frame(value=NA, row.species.i=NA, row.species.j=NA, row.species.i.abundance=NA)

##this loops through each row (indexed by "i"), 
#and calculates the proportion of interactions of each subsequent row "j"
#that are nested within the interactions of the focal row ("i"). 

counter <- 0
for (i in 1:(nr - 1)) {
  #get the ith row. starts with the first row in the matrix.
  first <- comm[i, ]
  
  for (j in (i + 1):nr) {
    counter <- counter + 1
    
    #if the jth row has higher or equal fill to the ith row OR if the fill of either row is 0, move on to the next value of j
    ##(the nestedness value is 0)
    if (rfill[i] <= rfill[j] || any(rfill[c(i, j)] == 0)) 
      next
    
    #otherwise, get the jth row and do the calculation. j starts with the first row AFTER row i.
      second <- comm[j, ]
      #this is the actual calculation.
      N.paired.rows[counter] <-
      sum(first - second > 0 & second > 0)/sum(second > 0)
      
      #store that info in my empty data frame "thegoods.rows"
      thegoods.rows<-rbind(
        thegoods.rows, 
        c(N.paired.rows[counter], names(comm[,1][i]), names((comm[,1][j])), rgrad[rorder][i])
        )
  
  }
  
}

#get rid of the empty first row of our dataframe 
thegoods.rows<-thegoods.rows[-1,]

####CALCULATE NESTEDNESS OF ALL PAIRS OF COLUMNS
##the same as above, but for each pair of columns.

counter <- 0
thegoods.cols<-data.frame(value=NA, column.species.i=NA, column.species.j=NA, col.species.i.abundance=NA)

for (i in 1:(nc - 1)) {
  first <- comm[, i]
  for (j in (i + 1):nc) {
    counter <- counter + 1
    if (cfill[i] <= cfill[j] || any(cfill[c(i, j)] == 0)) 
      next
      second <- comm[, j]
        N.paired.cols[counter] <- sum(first - second > 0 & second > 0)/sum(second > 0)
        
        thegoods.cols<-rbind(
          thegoods.cols, 
          c(N.paired.cols[counter], names(comm[1,][i]), names((comm[1,][j])), cgrad[corder][i])
          )
     
  }
}

thegoods.cols<-thegoods.cols[-1,]


#### 4.
####so we now have all of the pairwise nestedness values for each pollinator (columns) and plant (rows).
### to come up with a nestedness value for each species, let's take the mean of their pairwise values.

thegoods.rows$value<-as.numeric(thegoods.rows$value)
thegoods.cols$value<-as.numeric(thegoods.cols$value)


rmeans<-thegoods.rows%>%
  group_by(row.species.i, row.species.i.abundance)%>%
  summarize(species.wnodf=mean(value))

rmeans$row.species.i.abundance<-as.numeric(rmeans$row.species.i.abundance)


cmeans<-thegoods.cols%>%
  group_by(column.species.i, col.species.i.abundance)%>%
  summarize(species.wnodf=mean(value))

cmeans$col.species.i.abundance<-as.numeric(cmeans$col.species.i.abundance)


###5. visualizations
#I can now look at the distribution of species-level nestedness values, and the plot species' nestedness against
##their abundance.

hist(cmeans$species.wnodf)
hist(rmeans$species.wnodf)

ggplot(rmeans, aes(x=row.species.i.abundance,y=species.wnodf ))+
  geom_point()

ggplot(cmeans, aes(x=col.species.i.abundance,y=species.wnodf ))+
  geom_point()


