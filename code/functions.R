# This is a set of outdated but working functions to parse data
# Custom functions were created to build metadata off of 
# sample ID names. 

# There are a few  functions 
# 1. Not in
# 2. make_MA_metadata - build metadata from sample ID Names 
# 3. matround - for better rounding 
# 3. scale_reads - scaling the data
# 4. make_MA2_metadata - build metadata from combined replicates.

#----------------------------------------------------------------
# 1. Function for "not in" a vector
'%!in%' <- function(x,y)!('%in%'(x,y))


# 2. make_MA_metadata, which wil build metadata from sample names. 

## This function adds the categorical metadata to a dataframe based on the sample name
# IMPORTANT!!!  The dataframe MUST have column named "names"
# EXAMPLE: 20221014_FSP_128_WS 19932
make_FSP_metadata <- function(dataframe){ 
  
  # Create new columnes based on information in the sample name
  dataframe$year <- substr(dataframe$names, 1,4) # 1-4th letters = year sampled
  dataframe$month <- substr(dataframe$names, 5,6) # 5-6th letter = month sampled
  dataframe$day <- substr(dataframe$names, 7,8) # 7-8th letter = day sampled
  
  dataframe$project <- substr(dataframe$names, 10,12) # 10-11th letter = project
  
  dataframe$pond <- substr(dataframe$names,14,16) # 13-14 Pond Number info
  dataframe$fraction <- substr(dataframe$names, 18, 18) # 19 where in water column sample was taken
  # Samples were taken by dunking the bottle into the surface of the pond
  dataframe$depth <- substr(dataframe$names, 19, 19) # 19 where in water column sample was taken
  
  # Fraction
  dataframe$fraction <- ifelse(dataframe$fraction == "F", "Free", 
                               ifelse(dataframe$fraction == "P", "Particle", 
                                      ifelse(dataframe$fraction == "W","Whole", "Control")))
  dataframe$fraction <- as.factor(dataframe$fraction)
  dataframe$fraction <- factor(dataframe$fraction,levels = c("Particle", "Free", "Whole", "Control"))
  
  # Month
  dataframe$month <- ifelse(dataframe$month == "08", "August", 
                            ifelse(dataframe$month == "09", "September", 
                                   ifelse(dataframe$month == "10", "October", "Control")))
  dataframe$month <- as.factor(dataframe$month)
  dataframe$month <- factor(dataframe$month,levels = c("August", "September", "October", "Control"))
  
  # depth
  dataframe$depth <- ifelse(dataframe$depth == "S", "Surface", "Control")
  dataframe$depth <- as.factor(dataframe$depth)
  dataframe$depth <- factor(dataframe$depth,levels = c("Surface", "Control"))
  
  # year
  dataframe$year <- ifelse(dataframe$year == "2022", "2022", "2022")
  dataframe$year <- as.factor(dataframe$year)
  
  # Return the data
  return(dataframe)
  
}





# 
### scale_reads function
#################################################################################### 2
# Function to scale reads 
# Modified from code written by Michelle Berry, available at http://deneflab.github.io/MicrobeMiseq/ 
# Scales reads by 
# 1) taking proportions
# 2) multiplying by a given library size of n
# 3) rounding 
# Default for n is the minimum sample size in your library
# Default for round is floor

# 3. matround for better rounding 
matround <- function(x){trunc(x+0.5)}


# 4. scale reads function 
scale_reads <- function(physeq, n = min(sample_sums(physeq)), round = "round") {
  
  # transform counts to n
  physeq.scale <- transform_sample_counts(physeq, function(x) {(n * x/sum(x))})
  
  # Pick the rounding functions
  if (round == "floor"){
    otu_table(physeq.scale) <- floor(otu_table(physeq.scale))
  } else if (round == "round"){
    otu_table(physeq.scale) <- round(otu_table(physeq.scale))
  } else if (round == "matround"){
    otu_table(physeq.scale) <- matround(otu_table(physeq.scale))
  }
  
  # Prune taxa and return new phyloseq object
  physeq.scale <- prune_taxa(taxa_sums(physeq.scale) > 0, physeq.scale)
  return(physeq.scale)
}



######
#####
######

# 5. making metadata from replicates. 

make_MA2_metadata <- function(dataframe){ 
  
  # Create new columnes based on information in the sample name
  dataframe$year <- substr(dataframe$names, 2,5) # 1-4th letters = year sampled
  dataframe$month <- substr(dataframe$names, 6,7) # 5-6th letter = month sampled
  dataframe$day <- substr(dataframe$names, 8,9) # 7-8th letter = day sampled
  
  dataframe$project <- substr(dataframe$names, 11,12) # 10-11th letter = project
  
  dataframe$station <- substr(dataframe$names,14,15) # 13-14 station info
  dataframe$depth <- substr(dataframe$names, 16, 17) # 15th = where in water column sample was taken
  dataframe$replicate <- substr(dataframe$names, 17, 17) # 16th = Which sample replicate 
  dataframe$fraction <- substr(dataframe$names, 18,18) # 17th = filter fraction (e.g. particle, whole, free)
  
  # Fraction
  dataframe$fraction <- ifelse(dataframe$fraction == "F", "Free", 
                               ifelse(dataframe$fraction == "P", "Particle", 
                                      ifelse(dataframe$fraction == "W","Whole", "Control")))
  dataframe$fraction <- as.factor(dataframe$fraction)
  dataframe$fraction <- factor(dataframe$fraction,levels = c("Particle", "Free", "Whole", "Control"))
  
  # Month
  dataframe$month <- ifelse(dataframe$month == "06", "June", 
                            ifelse(dataframe$month == "10", "October", "Control"))
  dataframe$month <- as.factor(dataframe$month)
  dataframe$month <- factor(dataframe$month,levels = c("June", "October", "Control"))
  
  # Station
  dataframe$station <- ifelse(dataframe$station == "AB", "Aransas Bay", 
                              ifelse(dataframe$station == "CE", "Copano East", 
                                     ifelse(dataframe$station == "CW","Copano West", 
                                            ifelse(dataframe$station == "MB", "Mesquite Bay",
                                                   ifelse(dataframe$station == "SC", "Shipping Channel",
                                                          "Control")))))
  dataframe$station <- as.factor(dataframe$station)
  dataframe$station <- factor(dataframe$station,levels = c("Copano West", "Copano East", "Mesquite Bay", 
                                                           "Aransas Bay", "Shipping Channel", "Control"))
  
  # Return the data
  return(dataframe)
  
}


