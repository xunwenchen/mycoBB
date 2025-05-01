
# Simulation of randomly removing species/ASVs from the network and measure species/ASVs remained. ----
# These lines are a bit repetitive (not DRY) as they are not presented in functions, yet.
# These lines are time-consuming


#### ~~~~~~~~~~~~~~~~~~~ ####
#### Rhizosphere - Cd ----

# set conditions and thresholds
rho_cutoff <- 0.8

c_cutoff <- 5 

pl_ur_rzs_cd <- otu_table(pl_un_rarefied %>% subset_samples(Compartment =='Rhizosphere' & Group == 'Cd'))

dim(pl_ur_rzs_cd)
# column sum of otutab
colSums(pl_ur_rzs_cd) 

# Check if all column sums are equal
if (length(unique(colSums(pl_ur_rzs_cd))) == 1) {
  print("normalized OTU table")
} else {
  print("unnormalized OTU table")
}

# rarefy the OTU table using phyloseq package
pl_ur_rzs_cd_rare <- rarefy_even_depth(pl_ur_rzs_cd, rngseed = 123,
                                       sample.size = min(colSums(pl_ur_rzs_cd)))

if (length(unique(colSums(pl_ur_rzs_cd_rare))) == 1) {
  print("normalized OTU table")
} else {
  print("unnormalized OTU table")
}


pl_ur_rzs_cd_rare[is.na(pl_ur_rzs_cd_rare)] <- 0


counts <- rowSums(pl_ur_rzs_cd_rare > 0) # sum the non-zero values in each row (each otu)
hist(counts, breaks = 50) # plot histogram of counts

dim(pl_ur_rzs_cd_rare) # check df size

pl_ur_rzs_cd_rare <- pl_ur_rzs_cd_rare[counts >= c_cutoff, ] # keep the otus with at least c_cutoff counts
dim(pl_ur_rzs_cd_rare) # check df size again



comm <- t(pl_ur_rzs_cd_rare) # transpose df for later analysis
comm <- as(otu_table(comm), "matrix") 
rowSums(pl_ur_rzs_cd_rare)
depth <- min(colSums(pl_ur_rzs_cd))

if(length(depth) == 1){
  sp.ra <- colMeans(comm)/unique(min(colSums(pl_ur_rzs_cd)))  
} else {
  print("Error: The value of 'length' should be 1. Please check length(depth). The otu table is probably not normalized")
}

if(sum(sp.ra) == 1){
  print('OK to proceed')
} else { print('Sum of relative abundances should be 1. It should equal to 1, but still use depth of sub-otutab to calculate relative abundance of each species.')
}  


cormt = matrix(0, ncol(comm), ncol(comm)) # create a matrix with filled with zeros in required dimension
dim(cormt) # check size

for (i in 1:ncol(comm)){
  for (j in i:ncol(comm)){
    speciesi <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,i] > 0, comm[k,i], ifelse(comm[k,j] > 0, 0.01, NA))
    })
    speciesj <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,j] > 0, comm[k,j], ifelse(comm[k,i] > 0, 0.01, NA))
    })
    corij <- cor(log(speciesi)[!is.na(speciesi)],log(speciesj)[!is.na(speciesj)])
    
    cormt[i,j] <- cormt[j,i] <- corij
    
  }}

# plot histogram of counts of cormt, ignore 0 values
hist(cormt[cormt != 0], breaks = 50) 

cormt2 <- cormt*(abs(cormt) >= rho_cutoff)  # only keep links above the cutoff point
# dimension of matrix will not change, but the values will be changed to 0 if the absolute value is less than rho_cutoff

# plot histogram of counts of cormt2, and ignore 0 values
hist(cormt2[cormt2 != 0], breaks = 50)

sum(is.na(cormt2)) # check how many NAs
cormt2[is.na(cormt2)] <- 0 # assign 0 to NA values
sum(is.na(cormt2)) # check if there are any NA values left

diag(cormt2) <- 0 # remove self links
hist(cormt2[cormt2 != 0], breaks = 50) # plot histogram of counts of cormt2, and ignore 0 values again

sum(abs(cormt2) > 0) / 2  # this should be the number of links. 
sum(colSums(abs(cormt2)) > 0)  # node number: number of species with at least one linkage with others.

network.raw <- cormt2[colSums(abs(cormt2)) > 0, colSums(abs(cormt2)) > 0]
sp.ra2 <- sp.ra[colSums(abs(cormt2)) > 0]

dim(network.raw)
length(sp.ra)
length(sp.ra2)
# sum(row.names(network.raw) == names(sp.ra2))  # check if matched, but it is strange that network.raw has no row names. 

## node removal simulation 
#input network matrix, percentage of randomly removed species, and ra of all species
#return the proportion of species remained

rm.p.list = seq(0.05, 0.2, by = 0.05)

Weighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = T, nperm = 100)

Unweighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = F, nperm = 100)

rmsimu_rs_rhs_cd <- data.frame(Proportion.removed = rep(seq(0.05, 1, by = 0.05), 2), 
                 rbind(Weighted.simu, Unweighted.simu), 
                 weighted = rep(c("weighted", "unweighted"), each = 20))



#### Rhizosphere - M+Cd ----

# set conditions and threlsholds
rho_cutoff <- 0.8

c_cutoff <- 5 

pl_ur_rzs_M_cd <- otu_table(pl_un_rarefied %>% subset_samples(Compartment =='Rhizosphere' & Group == 'M+Cd'))

dim(pl_ur_rzs_M_cd)
# column sum of otutab
colSums(pl_ur_rzs_M_cd) # it is a normalized OTU table, so the sum of each column is the same, i.e., 3000 in this example. 
# return 'normalized OTU table' with the sum of each column is 3000, otherwise return 'unnormalized OTU table'

# Check if all column sums are equal
if (length(unique(colSums(pl_ur_rzs_M_cd))) == 1) {
  print("normalized OTU table")
} else {
  print("unnormalized OTU table")
}

# rarefy the OTU table using phyloseq package
pl_ur_rzs_M_cd_rare <- rarefy_even_depth(pl_ur_rzs_M_cd, rngseed = 123,
                                         sample.size = min(colSums(pl_ur_rzs_M_cd)))

if (length(unique(colSums(pl_ur_rzs_M_cd_rare))) == 1) {
  print("normalized OTU table")
} else {
  print("unnormalized OTU table")
}


pl_ur_rzs_M_cd_rare[is.na(pl_ur_rzs_M_cd_rare)] <- 0


counts <- rowSums(pl_ur_rzs_M_cd_rare > 0) # sum the non-zero values in each row (each otu)
hist(counts, breaks = 50) # plot histogram of counts

dim(pl_ur_rzs_M_cd_rare) # check df size

pl_ur_rzs_M_cd_rare <- pl_ur_rzs_M_cd_rare[counts >= c_cutoff, ] # keep the otus with at least c_cutoff counts
dim(pl_ur_rzs_M_cd_rare) # check df size again



comm <- t(pl_ur_rzs_M_cd_rare) # transpose df for later analysis
comm <- as(otu_table(comm), "matrix") 
rowSums(pl_ur_rzs_M_cd_rare)
depth <- min(colSums(pl_ur_rzs_M_cd))

if(length(depth) == 1){
  sp.ra <- colMeans(comm)/unique(min(colSums(pl_ur_rzs_M_cd)))  
} else {
  print("Error: The value of 'length' should be 1. Please check length(depth). The otutab_0 is probably not normalized")
}

if(sum(sp.ra) == 1){
  print('OK to proceed')
} else { print('Sum of relative abundances should be 1. It should equal to 1, but still use depth of sub-otutab to calculate relative abundance of each species.')
}  


cormt = matrix(0, ncol(comm), ncol(comm)) # create a matrix with filled with zeros in required dimension
dim(cormt) # check size

for (i in 1:ncol(comm)){
  for (j in i:ncol(comm)){
    speciesi <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,i] > 0, comm[k,i], ifelse(comm[k,j] > 0, 0.01, NA))
    })
    speciesj <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,j] > 0, comm[k,j], ifelse(comm[k,i] > 0, 0.01, NA))
    })
    corij <- cor(log(speciesi)[!is.na(speciesi)],log(speciesj)[!is.na(speciesj)])
    
    cormt[i,j] <- cormt[j,i] <- corij
    
  }}

# plot histogram of counts of cormt, ignore 0 values
hist(cormt[cormt != 0], breaks = 50) 

cormt2 <- cormt*(abs(cormt) >= rho_cutoff)  # only keep links above the cutoff point
# dimension of matrix will not change, but the values will be changed to 0 if the absolute value is less than rho_cutoff

# plot histogram of counts of cormt2, and ignore 0 values
hist(cormt2[cormt2 != 0], breaks = 50)

sum(is.na(cormt2)) # check how many NAs
cormt2[is.na(cormt2)] <- 0 # assign 0 to NA values
sum(is.na(cormt2)) # check if there are any NA values left

diag(cormt2) <- 0 # remove self links
hist(cormt2[cormt2 != 0], breaks = 50) # plot histogram of counts of cormt2, and ignore 0 values again

sum(abs(cormt2) > 0) / 2  # this should be the number of links. 
sum(colSums(abs(cormt2)) > 0)  # node number: number of species with at least one linkage with others.

network.raw <- cormt2[colSums(abs(cormt2)) > 0, colSums(abs(cormt2)) > 0]
sp.ra2 <- sp.ra[colSums(abs(cormt2)) > 0]

dim(network.raw)
length(sp.ra)
length(sp.ra2)
# sum(row.names(network.raw) == names(sp.ra2))  # check if matched, but it is strange that network.raw has no row names. 

## node removal simulation 
#input network matrix, percentage of randomly removed species, and ra of all species
#return the proportion of species remained

Weighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = T, nperm = 100)

Unweighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = F, nperm = 100)

rmsimu_rs_rhs_M_cd <- data.frame(Proportion.removed = rep(seq(0.05, 1, by = 0.05), 2), 
                                 rbind(Weighted.simu, Unweighted.simu), 
                                 weighted = rep(c("weighted", "unweighted"), each = 20))




#### Rhizosphere - Ctrl ----
# set conditions and thresholds
rho_cutoff <- 0.8

c_cutoff <- 5 

sub <- otu_table(pl_un_rarefied %>% subset_samples(Compartment =='Rhizosphere' & Group == 'Ctrl'))
sample_names(sub)
dim(sub)

# column sum of otutab
colSums(sub) 


# rarefy the OTU table using phyloseq package
sub_rare <- rarefy_even_depth(sub, rngseed = 123,
                              sample.size = min(colSums(sub)))

if (length(unique(colSums(sub_rare))) == 1) {
  print("normalized OTU table")
} else {
  print("unnormalized OTU table")
}


sub_rare[is.na(sub_rare)] <- 0


counts <- rowSums(sub_rare > 0) # sum the non-zero values in each row (each otu)
hist(counts, breaks = 50) # plot histogram of counts

dim(sub_rare) # check df size

sub_rare <- sub_rare[counts >= c_cutoff, ] # keep the otus with at least c_cutoff counts
dim(sub_rare) # check df size again



comm <- t(sub_rare) # transpose df for later analysis
comm <- as(otu_table(comm), "matrix") 
rowSums(sub_rare)
depth <- min(colSums(sub_rare))

if(length(depth) == 1){
  sp.ra <- colMeans(comm)/unique(min(colSums(sub)))  
} else {
  print("Error: The value of 'length' should be 1. Please check length(depth). The otutab_0 is probably not normalized")
}

if(sum(sp.ra) == 1){
  print('OK to proceed')
} else { print('Sum of relative abundances should be 1. It should equal to 1, but still use depth of sub-otutab to calculate relative abundance of each species.')
}  


cormt = matrix(0, ncol(comm), ncol(comm)) # create a matrix with filled with zeros in required dimension
dim(cormt) # check size

for (i in 1:ncol(comm)){
  for (j in i:ncol(comm)){
    speciesi <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,i] > 0, comm[k,i], ifelse(comm[k,j] > 0, 0.01, NA))
    })
    speciesj <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,j] > 0, comm[k,j], ifelse(comm[k,i] > 0, 0.01, NA))
    })
    corij <- cor(log(speciesi)[!is.na(speciesi)],log(speciesj)[!is.na(speciesj)])
    
    cormt[i,j] <- cormt[j,i] <- corij
    
  }}

# plot histogram of counts of cormt, ignore 0 values
hist(cormt[cormt != 0], breaks = 50) 

cormt2 <- cormt*(abs(cormt) >= rho_cutoff)  # only keep links above the cutoff point
# dimension of matrix will not change, but the values will be changed to 0 if the absolute value is less than rho_cutoff

# plot histogram of counts of cormt2, and ignore 0 values
hist(cormt2[cormt2 != 0], breaks = 50)

sum(is.na(cormt2)) # check how many NAs
cormt2[is.na(cormt2)] <- 0 # assign 0 to NA values
sum(is.na(cormt2)) # check if there are any NA values left

diag(cormt2) <- 0 # remove self links
hist(cormt2[cormt2 != 0], breaks = 50) # plot histogram of counts of cormt2, and ignore 0 values again

sum(abs(cormt2) > 0) / 2  # this should be the number of links. 
sum(colSums(abs(cormt2)) > 0)  # node number: number of species with at least one linkage with others.

network.raw <- cormt2[colSums(abs(cormt2)) > 0, colSums(abs(cormt2)) > 0]
sp.ra2 <- sp.ra[colSums(abs(cormt2)) > 0]

dim(network.raw)
length(sp.ra)
length(sp.ra2)
# sum(row.names(network.raw) == names(sp.ra2))  # check if matched, but it is strange that network.raw has no row names. 

## node removal simulation 
#input network matrix, percentage of randomly removed species, and ra of all species
#return the proportion of species remained

Weighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = T, nperm = 100)

Unweighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = F, nperm = 100)

rmsimu_rs_rhs_ctrl <- data.frame(Proportion.removed = rep(seq(0.05, 1, by = 0.05), 2), 
                                rbind(Weighted.simu, Unweighted.simu), 
                                weighted = rep(c("weighted", "unweighted"), each = 20))


#### Rhizosphere - M ----
# set conditions and thresholds
rho_cutoff <- 0.8

c_cutoff <- 5 

sub <- otu_table(pl_un_rarefied %>% subset_samples(Compartment =='Rhizosphere' & Group == 'M'))
sample_names(sub)
dim(sub)

# column sum of otutab
colSums(sub) 


# rarefy the OTU table using phyloseq package
sub_rare <- rarefy_even_depth(sub, rngseed = 123,
                              sample.size = min(colSums(sub)))

if (length(unique(colSums(sub_rare))) == 1) {
  print("normalized OTU table")
} else {
  print("unnormalized OTU table")
}


sub_rare[is.na(sub_rare)] <- 0


counts <- rowSums(sub_rare > 0) # sum the non-zero values in each row (each otu)
hist(counts, breaks = 50) # plot histogram of counts

dim(sub_rare) # check df size

sub_rare <- sub_rare[counts >= c_cutoff, ] # keep the otus with at least c_cutoff counts
dim(sub_rare) # check df size again



comm <- t(sub_rare) # transpose df for later analysis
comm <- as(otu_table(comm), "matrix") 
rowSums(sub_rare)
depth <- min(colSums(sub_rare))

if(length(depth) == 1){
  sp.ra <- colMeans(comm)/unique(min(colSums(sub)))  
} else {
  print("Error: The value of 'length' should be 1. Please check length(depth). The otutab_0 is probably not normalized")
}

if(sum(sp.ra) == 1){
  print('OK to proceed')
} else { print('Sum of relative abundances should be 1. It should equal to 1, but still use depth of sub-otutab to calculate relative abundance of each species.')
}  


cormt = matrix(0, ncol(comm), ncol(comm)) # create a matrix with filled with zeros in required dimension
dim(cormt) # check size

for (i in 1:ncol(comm)){
  for (j in i:ncol(comm)){
    speciesi <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,i] > 0, comm[k,i], ifelse(comm[k,j] > 0, 0.01, NA))
    })
    speciesj <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,j] > 0, comm[k,j], ifelse(comm[k,i] > 0, 0.01, NA))
    })
    corij <- cor(log(speciesi)[!is.na(speciesi)],log(speciesj)[!is.na(speciesj)])
    
    cormt[i,j] <- cormt[j,i] <- corij
    
  }}

# plot histogram of counts of cormt, ignore 0 values
hist(cormt[cormt != 0], breaks = 50) 

cormt2 <- cormt*(abs(cormt) >= rho_cutoff)  # only keep links above the cutoff point
# dimension of matrix will not change, but the values will be changed to 0 if the absolute value is less than rho_cutoff

# plot histogram of counts of cormt2, and ignore 0 values
hist(cormt2[cormt2 != 0], breaks = 50)

sum(is.na(cormt2)) # check how many NAs
cormt2[is.na(cormt2)] <- 0 # assign 0 to NA values
sum(is.na(cormt2)) # check if there are any NA values left

diag(cormt2) <- 0 # remove self links
hist(cormt2[cormt2 != 0], breaks = 50) # plot histogram of counts of cormt2, and ignore 0 values again

sum(abs(cormt2) > 0) / 2  # this should be the number of links. 
sum(colSums(abs(cormt2)) > 0)  # node number: number of species with at least one linkage with others.

network.raw <- cormt2[colSums(abs(cormt2)) > 0, colSums(abs(cormt2)) > 0]
sp.ra2 <- sp.ra[colSums(abs(cormt2)) > 0]

dim(network.raw)
length(sp.ra)
length(sp.ra2)
# sum(row.names(network.raw) == names(sp.ra2))  # check if matched, but it is strange that network.raw has no row names. 

## node removal simulation 
#input network matrix, percentage of randomly removed species, and ra of all species
#return the proportion of species remained

Weighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = T, nperm = 100)

Unweighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = F, nperm = 100)

rmsimu_rs_rhs_M <- data.frame(Proportion.removed = rep(seq(0.05, 1, by = 0.05), 2), 
                                 rbind(Weighted.simu, Unweighted.simu), 
                                 weighted = rep(c("weighted", "unweighted"), each = 20))

#### ~~~~~~~~~~~~~~~~~~~ ####
#### Hyphosphere - Cd ----
# set conditions and threlsholds
rho_cutoff <- 0.8

c_cutoff <- 5 

pl_ur_hps_cd <- otu_table(pl_un_rarefied %>% subset_samples(Compartment =='Hyphosphere' & Group == 'Cd'))

dim(pl_ur_hps_cd)
# column sum of otutab
colSums(pl_ur_hps_cd) # it is a normalized OTU table, so the sum of each column is the same, i.e., 3000 in this example. 
# return 'normalized OTU table' with the sum of each column is 3000, otherwise return 'unnormalized OTU table'

# Check if all column sums are equal
if (length(unique(colSums(pl_ur_hps_cd))) == 1) {
  print("normalized OTU table")
} else {
  print("unnormalized OTU table")
}

# rarefy the OTU table using phyloseq package
pl_ur_hps_cd_rare <- rarefy_even_depth(pl_ur_hps_cd, rngseed = 123,
                                       sample.size = min(colSums(pl_ur_hps_cd)))

if (length(unique(colSums(pl_ur_hps_cd_rare))) == 1) {
  print("normalized OTU table")
} else {
  print("unnormalized OTU table")
}


pl_ur_hps_cd_rare[is.na(pl_ur_hps_cd_rare)] <- 0


counts <- rowSums(pl_ur_hps_cd_rare > 0) # sum the non-zero values in each row (each otu)
hist(counts, breaks = 50) # plot histogram of counts

dim(pl_ur_hps_cd_rare) # check df size

pl_ur_hps_cd_rare <- pl_ur_hps_cd_rare[counts >= c_cutoff, ] # keep the otus with at least c_cutoff counts
dim(pl_ur_hps_cd_rare) # check df size again



comm <- t(pl_ur_hps_cd_rare) # transpose df for later analysis
comm <- as(otu_table(comm), "matrix") 
rowSums(pl_ur_hps_cd_rare)
depth <- min(colSums(pl_ur_hps_cd_rare))

if(length(depth) == 1){
  sp.ra <- colMeans(comm)/unique(min(colSums(pl_ur_hps_cd)))  
} else {
  print("Error: The value of 'length' should be 1. Please check length(depth). The otutab_0 is probably not normalized")
}

if(sum(sp.ra) == 1){
  print('OK to proceed')
} else { print('Sum of relative abundances should be 1. It should equal to 1, but still use depth of sub-otutab to calculate relative abundance of each species.')
}  

cormt = matrix(0, ncol(comm), ncol(comm)) # create a matrix with filled with zeros in required dimension
dim(cormt) # check size

for (i in 1:ncol(comm)){
  for (j in i:ncol(comm)){
    speciesi <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,i] > 0, comm[k,i], ifelse(comm[k,j] > 0, 0.01, NA))
    })
    speciesj <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,j] > 0, comm[k,j], ifelse(comm[k,i] > 0, 0.01, NA))
    })
    corij <- cor(log(speciesi)[!is.na(speciesi)],log(speciesj)[!is.na(speciesj)])
    
    cormt[i,j] <- cormt[j,i] <- corij
    
  }}

# plot histogram of counts of cormt, ignore 0 values
hist(cormt[cormt != 0], breaks = 50) 

cormt2 <- cormt*(abs(cormt) >= rho_cutoff)  # only keep links above the cutoff point
# dimension of matrix will not change, but the values will be changed to 0 if the absolute value is less than rho_cutoff

# plot histogram of counts of cormt2, and ignore 0 values
hist(cormt2[cormt2 != 0], breaks = 50)

sum(is.na(cormt2)) # check how many NAs
cormt2[is.na(cormt2)] <- 0 # assign 0 to NA values
sum(is.na(cormt2)) # check if there are any NA values left

diag(cormt2) <- 0 # remove self links
hist(cormt2[cormt2 != 0], breaks = 50) # plot histogram of counts of cormt2, and ignore 0 values again

sum(abs(cormt2) > 0) / 2  # this should be the number of links. 
sum(colSums(abs(cormt2)) > 0)  # node number: number of species with at least one linkage with others.

network.raw <- cormt2[colSums(abs(cormt2)) > 0, colSums(abs(cormt2)) > 0]
sp.ra2 <- sp.ra[colSums(abs(cormt2)) > 0]

dim(network.raw)
length(sp.ra)
length(sp.ra2)
# sum(row.names(network.raw) == names(sp.ra2))  # check if matched, but it is strange that network.raw has no row names. 

## node removal simulation 
#input network matrix, percentage of randomly removed species, and ra of all species
#return the proportion of species remained

Weighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = T, nperm = 100)

Unweighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = F, nperm = 100)

rmsimu_rs_hps_cd <- data.frame(Proportion.removed = rep(seq(0.05, 1, by = 0.05), 2), 
                               rbind(Weighted.simu, Unweighted.simu), 
                               weighted = rep(c("weighted", "unweighted"), each = 20))


#### Hyphosphere - M+cd ----
# set conditions and threlsholds
rho_cutoff <- 0.8

c_cutoff <- 5 

sub <- otu_table(pl_un_rarefied %>% subset_samples(Compartment =='Hyphosphere' & Group == 'M+Cd'))
sample_names(sub)
dim(sub)

# column sum of otutab
colSums(sub) 


# rarefy the OTU table using phyloseq package
sub_rare <- rarefy_even_depth(sub, rngseed = 123,
                              sample.size = min(colSums(sub)))

if (length(unique(colSums(sub_rare))) == 1) {
  print("normalized OTU table")
} else {
  print("unnormalized OTU table")
}


sub_rare[is.na(sub_rare)] <- 0


counts <- rowSums(sub_rare > 0) # sum the non-zero values in each row (each otu)
hist(counts, breaks = 50) # plot histogram of counts

dim(sub_rare) # check df size

sub_rare <- sub_rare[counts >= c_cutoff, ] # keep the otus with at least c_cutoff counts
dim(sub_rare) # check df size again



comm <- t(sub_rare) # transpose df for later analysis
comm <- as(otu_table(comm), "matrix") 
rowSums(sub_rare)
depth <- min(colSums(sub_rare))

if(length(depth) == 1){
  sp.ra <- colMeans(comm)/unique(min(colSums(sub)))  
} else {
  print("Error: The value of 'length' should be 1. Please check length(depth). The otutab_0 is probably not normalized")
}

if(sum(sp.ra) == 1){
  print('OK to proceed')
} else { print('Sum of relative abundances should be 1. It should equal to 1, but still use depth of sub-otutab to calculate relative abundance of each species.')
}  


cormt = matrix(0, ncol(comm), ncol(comm)) # create a matrix with filled with zeros in required dimension
dim(cormt) # check size

for (i in 1:ncol(comm)){
  for (j in i:ncol(comm)){
    speciesi <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,i] > 0, comm[k,i], ifelse(comm[k,j] > 0, 0.01, NA))
    })
    speciesj <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,j] > 0, comm[k,j], ifelse(comm[k,i] > 0, 0.01, NA))
    })
    corij <- cor(log(speciesi)[!is.na(speciesi)],log(speciesj)[!is.na(speciesj)])
    
    cormt[i,j] <- cormt[j,i] <- corij
    
  }}

# plot histogram of counts of cormt, ignore 0 values
hist(cormt[cormt != 0], breaks = 50) 

cormt2 <- cormt*(abs(cormt) >= rho_cutoff)  # only keep links above the cutoff point
# dimension of matrix will not change, but the values will be changed to 0 if the absolute value is less than rho_cutoff

# plot histogram of counts of cormt2, and ignore 0 values
hist(cormt2[cormt2 != 0], breaks = 50)

sum(is.na(cormt2)) # check how many NAs
cormt2[is.na(cormt2)] <- 0 # assign 0 to NA values
sum(is.na(cormt2)) # check if there are any NA values left

diag(cormt2) <- 0 # remove self links
hist(cormt2[cormt2 != 0], breaks = 50) # plot histogram of counts of cormt2, and ignore 0 values again

sum(abs(cormt2) > 0) / 2  # this should be the number of links. 
sum(colSums(abs(cormt2)) > 0)  # node number: number of species with at least one linkage with others.

network.raw <- cormt2[colSums(abs(cormt2)) > 0, colSums(abs(cormt2)) > 0]
sp.ra2 <- sp.ra[colSums(abs(cormt2)) > 0]

dim(network.raw)
length(sp.ra)
length(sp.ra2)
# sum(row.names(network.raw) == names(sp.ra2))  # check if matched, but it is strange that network.raw has no row names. 

## node removal simulation 
#input network matrix, percentage of randomly removed species, and ra of all species
#return the proportion of species remained

Weighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = T, nperm = 100)

Unweighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = F, nperm = 100)

rmsimu_rs_hps_M_cd <- data.frame(Proportion.removed = rep(seq(0.05, 1, by = 0.05), 2), 
                                 rbind(Weighted.simu, Unweighted.simu), 
                                 weighted = rep(c("weighted", "unweighted"), each = 20))



#### Hyphosphere - Ctrl ----
# set conditions and thresholds
rho_cutoff <- 0.8

c_cutoff <- 5 

sub <- otu_table(pl_un_rarefied %>% subset_samples(Compartment =='Hyphosphere' & Group == 'Ctrl'))
sample_names(sub)
dim(sub)

# column sum of otutab
colSums(sub) 


# rarefy the OTU table using phyloseq package
sub_rare <- rarefy_even_depth(sub, rngseed = 123,
                              sample.size = min(colSums(sub)))

if (length(unique(colSums(sub_rare))) == 1) {
  print("normalized OTU table")
} else {
  print("unnormalized OTU table")
}


sub_rare[is.na(sub_rare)] <- 0


counts <- rowSums(sub_rare > 0) # sum the non-zero values in each row (each otu)
hist(counts, breaks = 50) # plot histogram of counts

dim(sub_rare) # check df size

sub_rare <- sub_rare[counts >= c_cutoff, ] # keep the otus with at least c_cutoff counts
dim(sub_rare) # check df size again



comm <- t(sub_rare) # transpose df for later analysis
comm <- as(otu_table(comm), "matrix") 
rowSums(sub_rare)
depth <- min(colSums(sub_rare))

if(length(depth) == 1){
  sp.ra <- colMeans(comm)/unique(min(colSums(sub)))  
} else {
  print("Error: The value of 'length' should be 1. Please check length(depth). The otutab_0 is probably not normalized")
}

if(sum(sp.ra) == 1){
  print('OK to proceed')
} else { print('Sum of relative abundances should be 1. It should equal to 1, but still use depth of sub-otutab to calculate relative abundance of each species.')
}  


cormt = matrix(0, ncol(comm), ncol(comm)) # create a matrix with filled with zeros in required dimension
dim(cormt) # check size

for (i in 1:ncol(comm)){
  for (j in i:ncol(comm)){
    speciesi <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,i] > 0, comm[k,i], ifelse(comm[k,j] > 0, 0.01, NA))
    })
    speciesj <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,j] > 0, comm[k,j], ifelse(comm[k,i] > 0, 0.01, NA))
    })
    corij <- cor(log(speciesi)[!is.na(speciesi)],log(speciesj)[!is.na(speciesj)])
    
    cormt[i,j] <- cormt[j,i] <- corij
    
  }}

# plot histogram of counts of cormt, ignore 0 values
hist(cormt[cormt != 0], breaks = 50) 

cormt2 <- cormt*(abs(cormt) >= rho_cutoff)  # only keep links above the cutoff point
# dimension of matrix will not change, but the values will be changed to 0 if the absolute value is less than rho_cutoff

# plot histogram of counts of cormt2, and ignore 0 values
hist(cormt2[cormt2 != 0], breaks = 50)

sum(is.na(cormt2)) # check how many NAs
cormt2[is.na(cormt2)] <- 0 # assign 0 to NA values
sum(is.na(cormt2)) # check if there are any NA values left

diag(cormt2) <- 0 # remove self links
hist(cormt2[cormt2 != 0], breaks = 50) # plot histogram of counts of cormt2, and ignore 0 values again

sum(abs(cormt2) > 0) / 2  # this should be the number of links. 
sum(colSums(abs(cormt2)) > 0)  # node number: number of species with at least one linkage with others.

network.raw <- cormt2[colSums(abs(cormt2)) > 0, colSums(abs(cormt2)) > 0]
sp.ra2 <- sp.ra[colSums(abs(cormt2)) > 0]

dim(network.raw)
length(sp.ra)
length(sp.ra2)
# sum(row.names(network.raw) == names(sp.ra2))  # check if matched, but it is strange that network.raw has no row names. 

## node removal simulation 
#input network matrix, percentage of randomly removed species, and ra of all species
#return the proportion of species remained

Weighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = T, nperm = 100)

Unweighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = F, nperm = 100)

rmsimu_rs_hps_ctrl <- data.frame(Proportion.removed = rep(seq(0.05, 1, by = 0.05), 2), 
                              rbind(Weighted.simu, Unweighted.simu), 
                              weighted = rep(c("weighted", "unweighted"), each = 20))
#### Hyphosphere - M ----
# set conditions and thresholds
rho_cutoff <- 0.8

c_cutoff <- 5 

sub <- otu_table(pl_un_rarefied %>% subset_samples(Compartment =='Hyphosphere' & Group == 'M'))
sample_names(sub)
dim(sub)

# column sum of otutab
colSums(sub) 


# rarefy the OTU table using phyloseq package
sub_rare <- rarefy_even_depth(sub, rngseed = 123,
                              sample.size = min(colSums(sub)))

if (length(unique(colSums(sub_rare))) == 1) {
  print("normalized OTU table")
} else {
  print("unnormalized OTU table")
}


sub_rare[is.na(sub_rare)] <- 0


counts <- rowSums(sub_rare > 0) # sum the non-zero values in each row (each otu)
hist(counts, breaks = 50) # plot histogram of counts

dim(sub_rare) # check df size

sub_rare <- sub_rare[counts >= c_cutoff, ] # keep the otus with at least c_cutoff counts
dim(sub_rare) # check df size again



comm <- t(sub_rare) # transpose df for later analysis
comm <- as(otu_table(comm), "matrix") 
rowSums(sub_rare)
depth <- min(colSums(sub_rare))

if(length(depth) == 1){
  sp.ra <- colMeans(comm)/unique(min(colSums(sub)))  
} else {
  print("Error: The value of 'length' should be 1. Please check length(depth). The otutab_0 is probably not normalized")
}

if(sum(sp.ra) == 1){
  print('OK to proceed')
} else { print('Sum of relative abundances should be 1. It should equal to 1, but still use depth of sub-otutab to calculate relative abundance of each species.')
}  


cormt = matrix(0, ncol(comm), ncol(comm)) # create a matrix with filled with zeros in required dimension
dim(cormt) # check size

for (i in 1:ncol(comm)){
  for (j in i:ncol(comm)){
    speciesi <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,i] > 0, comm[k,i], ifelse(comm[k,j] > 0, 0.01, NA))
    })
    speciesj <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,j] > 0, comm[k,j], ifelse(comm[k,i] > 0, 0.01, NA))
    })
    corij <- cor(log(speciesi)[!is.na(speciesi)],log(speciesj)[!is.na(speciesj)])
    
    cormt[i,j] <- cormt[j,i] <- corij
    
  }}

# plot histogram of counts of cormt, ignore 0 values
hist(cormt[cormt != 0], breaks = 50) 

cormt2 <- cormt*(abs(cormt) >= rho_cutoff)  # only keep links above the cutoff point
# dimension of matrix will not change, but the values will be changed to 0 if the absolute value is less than rho_cutoff

# plot histogram of counts of cormt2, and ignore 0 values
hist(cormt2[cormt2 != 0], breaks = 50)

sum(is.na(cormt2)) # check how many NAs
cormt2[is.na(cormt2)] <- 0 # assign 0 to NA values
sum(is.na(cormt2)) # check if there are any NA values left

diag(cormt2) <- 0 # remove self links
hist(cormt2[cormt2 != 0], breaks = 50) # plot histogram of counts of cormt2, and ignore 0 values again

sum(abs(cormt2) > 0) / 2  # this should be the number of links. 
sum(colSums(abs(cormt2)) > 0)  # node number: number of species with at least one linkage with others.

network.raw <- cormt2[colSums(abs(cormt2)) > 0, colSums(abs(cormt2)) > 0]
sp.ra2 <- sp.ra[colSums(abs(cormt2)) > 0]

dim(network.raw)
length(sp.ra)
length(sp.ra2)
# sum(row.names(network.raw) == names(sp.ra2))  # check if matched, but it is strange that network.raw has no row names. 

## node removal simulation 
#input network matrix, percentage of randomly removed species, and ra of all species
#return the proportion of species remained

Weighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = T, nperm = 100)

Unweighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = F, nperm = 100)

rmsimu_rs_hps_M <- data.frame(Proportion.removed = rep(seq(0.05, 1, by = 0.05), 2), 
                                 rbind(Weighted.simu, Unweighted.simu), 
                                 weighted = rep(c("weighted", "unweighted"), each = 20))

#### ~~~~~~~~~~~~~~~~~~~ ####
#### Hyphae - M ----
# set conditions and thresholds
rho_cutoff <- 0.8

c_cutoff <- 5 

sub <- otu_table(pl_un_rarefied %>% subset_samples(Compartment =='Hyphae' & Group == 'M'))
sample_names(sub)
dim(sub)

# column sum of otutab
colSums(sub) 


# rarefy the OTU table using phyloseq package
sub_rare <- rarefy_even_depth(sub, rngseed = 123,
                              sample.size = min(colSums(sub)))

if (length(unique(colSums(sub_rare))) == 1) {
  print("normalized OTU table")
} else {
  print("unnormalized OTU table")
}


sub_rare[is.na(sub_rare)] <- 0


counts <- rowSums(sub_rare > 0) # sum the non-zero values in each row (each otu)
hist(counts, breaks = 50) # plot histogram of counts

dim(sub_rare) # check df size

sub_rare <- sub_rare[counts >= c_cutoff, ] # keep the otus with at least c_cutoff counts
dim(sub_rare) # check df size again



comm <- t(sub_rare) # transpose df for later analysis
comm <- as(otu_table(comm), "matrix") 
rowSums(sub_rare)
depth <- min(colSums(sub_rare))

if(length(depth) == 1){
  sp.ra <- colMeans(comm)/unique(min(colSums(sub)))  
} else {
  print("Error: The value of 'length' should be 1. Please check length(depth). The otutab_0 is probably not normalized")
}

if(sum(sp.ra) == 1){
  print('OK to proceed')
} else { print('Sum of relative abundances should be 1. It should equal to 1, but still use depth of sub-otutab to calculate relative abundance of each species.')
}  


cormt = matrix(0, ncol(comm), ncol(comm)) # create a matrix with filled with zeros in required dimension
dim(cormt) # check size

for (i in 1:ncol(comm)){
  for (j in i:ncol(comm)){
    speciesi <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,i] > 0, comm[k,i], ifelse(comm[k,j] > 0, 0.01, NA))
    })
    speciesj <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,j] > 0, comm[k,j], ifelse(comm[k,i] > 0, 0.01, NA))
    })
    corij <- cor(log(speciesi)[!is.na(speciesi)],log(speciesj)[!is.na(speciesj)])
    
    cormt[i,j] <- cormt[j,i] <- corij
    
  }}

# plot histogram of counts of cormt, ignore 0 values
hist(cormt[cormt != 0], breaks = 50) 

cormt2 <- cormt*(abs(cormt) >= rho_cutoff)  # only keep links above the cutoff point
# dimension of matrix will not change, but the values will be changed to 0 if the absolute value is less than rho_cutoff

# plot histogram of counts of cormt2, and ignore 0 values
hist(cormt2[cormt2 != 0], breaks = 50)

sum(is.na(cormt2)) # check how many NAs
cormt2[is.na(cormt2)] <- 0 # assign 0 to NA values
sum(is.na(cormt2)) # check if there are any NA values left

diag(cormt2) <- 0 # remove self links
hist(cormt2[cormt2 != 0], breaks = 50) # plot histogram of counts of cormt2, and ignore 0 values again

sum(abs(cormt2) > 0) / 2  # this should be the number of links. 
sum(colSums(abs(cormt2)) > 0)  # node number: number of species with at least one linkage with others.

network.raw <- cormt2[colSums(abs(cormt2)) > 0, colSums(abs(cormt2)) > 0]
sp.ra2 <- sp.ra[colSums(abs(cormt2)) > 0]

dim(network.raw)
length(sp.ra)
length(sp.ra2)
# sum(row.names(network.raw) == names(sp.ra2))  # check if matched, but it is strange that network.raw has no row names. 

## node removal simulation 
#input network matrix, percentage of randomly removed species, and ra of all species
#return the proportion of species remained

Weighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = T, nperm = 100)

Unweighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = F, nperm = 100)

rmsimu_rs_hp_M <- data.frame(Proportion.removed = rep(seq(0.05, 1, by = 0.05), 2), 
                                 rbind(Weighted.simu, Unweighted.simu), 
                                 weighted = rep(c("weighted", "unweighted"), each = 20))

#### Hyphae - M+Cd ----
# set conditions and thresholds
rho_cutoff <- 0.8

c_cutoff <- 5 

sub <- otu_table(pl_un_rarefied %>% subset_samples(Compartment =='Hyphae' & Group == 'M+Cd'))
sample_names(sub)
dim(sub)

# column sum of otutab
colSums(sub) 


# rarefy the OTU table using phyloseq package
sub_rare <- rarefy_even_depth(sub, rngseed = 123,
                              sample.size = min(colSums(sub)))

if (length(unique(colSums(sub_rare))) == 1) {
  print("normalized OTU table")
} else {
  print("unnormalized OTU table")
}


sub_rare[is.na(sub_rare)] <- 0


counts <- rowSums(sub_rare > 0) # sum the non-zero values in each row (each otu)
hist(counts, breaks = 50) # plot histogram of counts

dim(sub_rare) # check df size

sub_rare <- sub_rare[counts >= c_cutoff, ] # keep the otus with at least c_cutoff counts
dim(sub_rare) # check df size again



comm <- t(sub_rare) # transpose df for later analysis
comm <- as(otu_table(comm), "matrix") 
rowSums(sub_rare)
depth <- min(colSums(sub_rare))

if(length(depth) == 1){
  sp.ra <- colMeans(comm)/unique(min(colSums(sub)))  
} else {
  print("Error: The value of 'length' should be 1. Please check length(depth). The otutab_0 is probably not normalized")
}

if(sum(sp.ra) == 1){
  print('OK to proceed')
} else { print('Sum of relative abundances should be 1. It should equal to 1, but still use depth of sub-otutab to calculate relative abundance of each species.')
}  


cormt = matrix(0, ncol(comm), ncol(comm)) # create a matrix with filled with zeros in required dimension
dim(cormt) # check size

for (i in 1:ncol(comm)){
  for (j in i:ncol(comm)){
    speciesi <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,i] > 0, comm[k,i], ifelse(comm[k,j] > 0, 0.01, NA))
    })
    speciesj <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,j] > 0, comm[k,j], ifelse(comm[k,i] > 0, 0.01, NA))
    })
    corij <- cor(log(speciesi)[!is.na(speciesi)],log(speciesj)[!is.na(speciesj)])
    
    cormt[i,j] <- cormt[j,i] <- corij
    
  }}

# plot histogram of counts of cormt, ignore 0 values
hist(cormt[cormt != 0], breaks = 50) 

cormt2 <- cormt*(abs(cormt) >= rho_cutoff)  # only keep links above the cutoff point
# dimension of matrix will not change, but the values will be changed to 0 if the absolute value is less than rho_cutoff

# plot histogram of counts of cormt2, and ignore 0 values
hist(cormt2[cormt2 != 0], breaks = 50)

sum(is.na(cormt2)) # check how many NAs
cormt2[is.na(cormt2)] <- 0 # assign 0 to NA values
sum(is.na(cormt2)) # check if there are any NA values left

diag(cormt2) <- 0 # remove self links
hist(cormt2[cormt2 != 0], breaks = 50) # plot histogram of counts of cormt2, and ignore 0 values again

sum(abs(cormt2) > 0) / 2  # this should be the number of links. 
sum(colSums(abs(cormt2)) > 0)  # node number: number of species with at least one linkage with others.

network.raw <- cormt2[colSums(abs(cormt2)) > 0, colSums(abs(cormt2)) > 0]
sp.ra2 <- sp.ra[colSums(abs(cormt2)) > 0]

dim(network.raw)
length(sp.ra)
length(sp.ra2)
# sum(row.names(network.raw) == names(sp.ra2))  # check if matched, but it is strange that network.raw has no row names. 

## node removal simulation 
#input network matrix, percentage of randomly removed species, and ra of all species
#return the proportion of species remained

Weighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = T, nperm = 100)

Unweighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = F, nperm = 100)

rmsimu_rs_hp_M_cd <- data.frame(Proportion.removed = rep(seq(0.05, 1, by = 0.05), 2), 
                             rbind(Weighted.simu, Unweighted.simu), 
                             weighted = rep(c("weighted", "unweighted"), each = 20))








#### ~~~~~~~~~~~~~~~~~~~ ####
#### Bulk soil - Cd ----
# set conditions and thresholds
rho_cutoff <- 0.8

c_cutoff <- 5 

sub <- otu_table(pl_un_rarefied %>% subset_samples(Compartment =='Bulk soil' & Group == 'Cd'))
sample_names(sub)
dim(sub)

# column sum of otutab
colSums(sub) 


# rarefy the OTU table using phyloseq package
sub_rare <- rarefy_even_depth(sub, rngseed = 123,
                              sample.size = min(colSums(sub)))

if (length(unique(colSums(sub_rare))) == 1) {
  print("normalized OTU table")
} else {
  print("unnormalized OTU table")
}


sub_rare[is.na(sub_rare)] <- 0


counts <- rowSums(sub_rare > 0) # sum the non-zero values in each row (each otu)
hist(counts, breaks = 50) # plot histogram of counts

dim(sub_rare) # check df size

sub_rare <- sub_rare[counts >= c_cutoff, ] # keep the otus with at least c_cutoff counts
dim(sub_rare) # check df size again



comm <- t(sub_rare) # transpose df for later analysis
comm <- as(otu_table(comm), "matrix") 
rowSums(sub_rare)
depth <- min(colSums(sub_rare))

if(length(depth) == 1){
  sp.ra <- colMeans(comm)/unique(min(colSums(sub)))  
} else {
  print("Error: The value of 'length' should be 1. Please check length(depth). The otutab_0 is probably not normalized")
}

if(sum(sp.ra) == 1){
  print('OK to proceed')
} else { print('Sum of relative abundances should be 1. It should equal to 1, but still use depth of sub-otutab to calculate relative abundance of each species.')
}  


cormt = matrix(0, ncol(comm), ncol(comm)) # create a matrix with filled with zeros in required dimension
dim(cormt) # check size

for (i in 1:ncol(comm)){
  for (j in i:ncol(comm)){
    speciesi <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,i] > 0, comm[k,i], ifelse(comm[k,j] > 0, 0.01, NA))
    })
    speciesj <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,j] > 0, comm[k,j], ifelse(comm[k,i] > 0, 0.01, NA))
    })
    corij <- cor(log(speciesi)[!is.na(speciesi)],log(speciesj)[!is.na(speciesj)])
    
    cormt[i,j] <- cormt[j,i] <- corij
    
  }}

# plot histogram of counts of cormt, ignore 0 values
hist(cormt[cormt != 0], breaks = 50) 

cormt2 <- cormt*(abs(cormt) >= rho_cutoff)  # only keep links above the cutoff point
# dimension of matrix will not change, but the values will be changed to 0 if the absolute value is less than rho_cutoff

# plot histogram of counts of cormt2, and ignore 0 values
hist(cormt2[cormt2 != 0], breaks = 50)

sum(is.na(cormt2)) # check how many NAs
cormt2[is.na(cormt2)] <- 0 # assign 0 to NA values
sum(is.na(cormt2)) # check if there are any NA values left

diag(cormt2) <- 0 # remove self links
hist(cormt2[cormt2 != 0], breaks = 50) # plot histogram of counts of cormt2, and ignore 0 values again

sum(abs(cormt2) > 0) / 2  # this should be the number of links. 
sum(colSums(abs(cormt2)) > 0)  # node number: number of species with at least one linkage with others.

network.raw <- cormt2[colSums(abs(cormt2)) > 0, colSums(abs(cormt2)) > 0]
sp.ra2 <- sp.ra[colSums(abs(cormt2)) > 0]

dim(network.raw)
length(sp.ra)
length(sp.ra2)
# sum(row.names(network.raw) == names(sp.ra2))  # check if matched, but it is strange that network.raw has no row names. 

## node removal simulation 
#input network matrix, percentage of randomly removed species, and ra of all species
#return the proportion of species remained

Weighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = T, nperm = 100)

Unweighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = F, nperm = 100)

rmsimu_rs_bs_cd <- data.frame(Proportion.removed = rep(seq(0.05, 1, by = 0.05), 2), 
                                rbind(Weighted.simu, Unweighted.simu), 
                                weighted = rep(c("weighted", "unweighted"), each = 20))

#### Bulk soil - M+Cd ----
# set conditions and thresholds
rho_cutoff <- 0.8

c_cutoff <- 5 

sub <- otu_table(pl_un_rarefied %>% subset_samples(Compartment =='Bulk soil' & Group == 'M+Cd'))
sample_names(sub)
dim(sub)

# column sum of otutab
colSums(sub) 


# rarefy the OTU table using phyloseq package
sub_rare <- rarefy_even_depth(sub, rngseed = 123,
                              sample.size = min(colSums(sub)))

if (length(unique(colSums(sub_rare))) == 1) {
  print("normalized OTU table")
} else {
  print("unnormalized OTU table")
}


sub_rare[is.na(sub_rare)] <- 0


counts <- rowSums(sub_rare > 0) # sum the non-zero values in each row (each otu)
hist(counts, breaks = 50) # plot histogram of counts

dim(sub_rare) # check df size

sub_rare <- sub_rare[counts >= c_cutoff, ] # keep the otus with at least c_cutoff counts
dim(sub_rare) # check df size again



comm <- t(sub_rare) # transpose df for later analysis
comm <- as(otu_table(comm), "matrix") 
rowSums(sub_rare)
depth <- min(colSums(sub_rare))

if(length(depth) == 1){
  sp.ra <- colMeans(comm)/unique(min(colSums(sub)))  
} else {
  print("Error: The value of 'length' should be 1. Please check length(depth). The otutab_0 is probably not normalized")
}

if(sum(sp.ra) == 1){
  print('OK to proceed')
} else { print('Sum of relative abundances should be 1. It should equal to 1, but still use depth of sub-otutab to calculate relative abundance of each species.')
}  


cormt = matrix(0, ncol(comm), ncol(comm)) # create a matrix with filled with zeros in required dimension
dim(cormt) # check size

for (i in 1:ncol(comm)){
  for (j in i:ncol(comm)){
    speciesi <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,i] > 0, comm[k,i], ifelse(comm[k,j] > 0, 0.01, NA))
    })
    speciesj <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,j] > 0, comm[k,j], ifelse(comm[k,i] > 0, 0.01, NA))
    })
    corij <- cor(log(speciesi)[!is.na(speciesi)],log(speciesj)[!is.na(speciesj)])
    
    cormt[i,j] <- cormt[j,i] <- corij
    
  }}

# plot histogram of counts of cormt, ignore 0 values
hist(cormt[cormt != 0], breaks = 50) 

cormt2 <- cormt*(abs(cormt) >= rho_cutoff)  # only keep links above the cutoff point
# dimension of matrix will not change, but the values will be changed to 0 if the absolute value is less than rho_cutoff

# plot histogram of counts of cormt2, and ignore 0 values
hist(cormt2[cormt2 != 0], breaks = 50)

sum(is.na(cormt2)) # check how many NAs
cormt2[is.na(cormt2)] <- 0 # assign 0 to NA values
sum(is.na(cormt2)) # check if there are any NA values left

diag(cormt2) <- 0 # remove self links
hist(cormt2[cormt2 != 0], breaks = 50) # plot histogram of counts of cormt2, and ignore 0 values again

sum(abs(cormt2) > 0) / 2  # this should be the number of links. 
sum(colSums(abs(cormt2)) > 0)  # node number: number of species with at least one linkage with others.

network.raw <- cormt2[colSums(abs(cormt2)) > 0, colSums(abs(cormt2)) > 0]
sp.ra2 <- sp.ra[colSums(abs(cormt2)) > 0]

dim(network.raw)
length(sp.ra)
length(sp.ra2)
# sum(row.names(network.raw) == names(sp.ra2))  # check if matched, but it is strange that network.raw has no row names. 

## node removal simulation 
#input network matrix, percentage of randomly removed species, and ra of all species
#return the proportion of species remained

Weighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = T, nperm = 100)

Unweighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = F, nperm = 100)

rmsimu_rs_bs_M_cd <- data.frame(Proportion.removed = rep(seq(0.05, 1, by = 0.05), 2), 
                              rbind(Weighted.simu, Unweighted.simu), 
                              weighted = rep(c("weighted", "unweighted"), each = 20))

#### Bulk soil - Ctrl ----
# set conditions and thresholds
rho_cutoff <- 0.8

c_cutoff <- 5 

sub <- otu_table(pl_un_rarefied %>% subset_samples(Compartment =='Bulk soil' & Group == 'Ctrl'))
sample_names(sub)
dim(sub)

# column sum of otutab
colSums(sub) 


# rarefy the OTU table using phyloseq package
sub_rare <- rarefy_even_depth(sub, rngseed = 123,
                              sample.size = min(colSums(sub)))

if (length(unique(colSums(sub_rare))) == 1) {
  print("normalized OTU table")
} else {
  print("unnormalized OTU table")
}


sub_rare[is.na(sub_rare)] <- 0


counts <- rowSums(sub_rare > 0) # sum the non-zero values in each row (each otu)
hist(counts, breaks = 50) # plot histogram of counts

dim(sub_rare) # check df size

sub_rare <- sub_rare[counts >= c_cutoff, ] # keep the otus with at least c_cutoff counts
dim(sub_rare) # check df size again



comm <- t(sub_rare) # transpose df for later analysis
comm <- as(otu_table(comm), "matrix") 
rowSums(sub_rare)
depth <- min(colSums(sub_rare))

if(length(depth) == 1){
  sp.ra <- colMeans(comm)/unique(min(colSums(sub)))  
} else {
  print("Error: The value of 'length' should be 1. Please check length(depth). The otutab_0 is probably not normalized")
}

if(sum(sp.ra) == 1){
  print('OK to proceed')
} else { print('Sum of relative abundances should be 1. It should equal to 1, but still use depth of sub-otutab to calculate relative abundance of each species.')
}  


cormt = matrix(0, ncol(comm), ncol(comm)) # create a matrix with filled with zeros in required dimension
dim(cormt) # check size

for (i in 1:ncol(comm)){
  for (j in i:ncol(comm)){
    speciesi <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,i] > 0, comm[k,i], ifelse(comm[k,j] > 0, 0.01, NA))
    })
    speciesj <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,j] > 0, comm[k,j], ifelse(comm[k,i] > 0, 0.01, NA))
    })
    corij <- cor(log(speciesi)[!is.na(speciesi)],log(speciesj)[!is.na(speciesj)])
    
    cormt[i,j] <- cormt[j,i] <- corij
    
  }}

# plot histogram of counts of cormt, ignore 0 values
hist(cormt[cormt != 0], breaks = 50) 

cormt2 <- cormt*(abs(cormt) >= rho_cutoff)  # only keep links above the cutoff point
# dimension of matrix will not change, but the values will be changed to 0 if the absolute value is less than rho_cutoff

# plot histogram of counts of cormt2, and ignore 0 values
hist(cormt2[cormt2 != 0], breaks = 50)

sum(is.na(cormt2)) # check how many NAs
cormt2[is.na(cormt2)] <- 0 # assign 0 to NA values
sum(is.na(cormt2)) # check if there are any NA values left

diag(cormt2) <- 0 # remove self links
hist(cormt2[cormt2 != 0], breaks = 50) # plot histogram of counts of cormt2, and ignore 0 values again

sum(abs(cormt2) > 0) / 2  # this should be the number of links. 
sum(colSums(abs(cormt2)) > 0)  # node number: number of species with at least one linkage with others.

network.raw <- cormt2[colSums(abs(cormt2)) > 0, colSums(abs(cormt2)) > 0]
sp.ra2 <- sp.ra[colSums(abs(cormt2)) > 0]

dim(network.raw)
length(sp.ra)
length(sp.ra2)
# sum(row.names(network.raw) == names(sp.ra2))  # check if matched, but it is strange that network.raw has no row names. 

## node removal simulation 
#input network matrix, percentage of randomly removed species, and ra of all species
#return the proportion of species remained

Weighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = T, nperm = 100)

Unweighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = F, nperm = 100)

rmsimu_rs_bs_ctrl <- data.frame(Proportion.removed = rep(seq(0.05, 1, by = 0.05), 2), 
                              rbind(Weighted.simu, Unweighted.simu), 
                              weighted = rep(c("weighted", "unweighted"), each = 20))

#### Bulk soil - M ----
# set conditions and thresholds
rho_cutoff <- 0.8

c_cutoff <- 5 

sub <- otu_table(pl_un_rarefied %>% subset_samples(Compartment =='Bulk soil' & Group == 'M'))
sample_names(sub)
dim(sub)

# column sum of otutab
colSums(sub) 


# rarefy the OTU table using phyloseq package
sub_rare <- rarefy_even_depth(sub, rngseed = 123,
                              sample.size = min(colSums(sub)))

if (length(unique(colSums(sub_rare))) == 1) {
  print("normalized OTU table")
} else {
  print("unnormalized OTU table")
}


sub_rare[is.na(sub_rare)] <- 0


counts <- rowSums(sub_rare > 0) # sum the non-zero values in each row (each otu)
hist(counts, breaks = 50) # plot histogram of counts

dim(sub_rare) # check df size

sub_rare <- sub_rare[counts >= c_cutoff, ] # keep the otus with at least c_cutoff counts
dim(sub_rare) # check df size again



comm <- t(sub_rare) # transpose df for later analysis
comm <- as(otu_table(comm), "matrix") 
rowSums(sub_rare)
depth <- min(colSums(sub_rare))

if(length(depth) == 1){
  sp.ra <- colMeans(comm)/unique(min(colSums(sub)))  
} else {
  print("Error: The value of 'length' should be 1. Please check length(depth). The otutab_0 is probably not normalized")
}

if(sum(sp.ra) == 1){
  print('OK to proceed')
} else { print('Sum of relative abundances should be 1. It should equal to 1, but still use depth of sub-otutab to calculate relative abundance of each species.')
}  


cormt = matrix(0, ncol(comm), ncol(comm)) # create a matrix with filled with zeros in required dimension
dim(cormt) # check size

for (i in 1:ncol(comm)){
  for (j in i:ncol(comm)){
    speciesi <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,i] > 0, comm[k,i], ifelse(comm[k,j] > 0, 0.01, NA))
    })
    speciesj <- sapply(1:nrow(comm),function(k){
      ifelse(comm[k,j] > 0, comm[k,j], ifelse(comm[k,i] > 0, 0.01, NA))
    })
    corij <- cor(log(speciesi)[!is.na(speciesi)],log(speciesj)[!is.na(speciesj)])
    
    cormt[i,j] <- cormt[j,i] <- corij
    
  }}

# plot histogram of counts of cormt, ignore 0 values
hist(cormt[cormt != 0], breaks = 50) 

cormt2 <- cormt*(abs(cormt) >= rho_cutoff)  # only keep links above the cutoff point
# dimension of matrix will not change, but the values will be changed to 0 if the absolute value is less than rho_cutoff

# plot histogram of counts of cormt2, and ignore 0 values
hist(cormt2[cormt2 != 0], breaks = 50)

sum(is.na(cormt2)) # check how many NAs
cormt2[is.na(cormt2)] <- 0 # assign 0 to NA values
sum(is.na(cormt2)) # check if there are any NA values left

diag(cormt2) <- 0 # remove self links
hist(cormt2[cormt2 != 0], breaks = 50) # plot histogram of counts of cormt2, and ignore 0 values again

sum(abs(cormt2) > 0) / 2  # this should be the number of links. 
sum(colSums(abs(cormt2)) > 0)  # node number: number of species with at least one linkage with others.

network.raw <- cormt2[colSums(abs(cormt2)) > 0, colSums(abs(cormt2)) > 0]
sp.ra2 <- sp.ra[colSums(abs(cormt2)) > 0]

dim(network.raw)
length(sp.ra)
length(sp.ra2)
# sum(row.names(network.raw) == names(sp.ra2))  # check if matched, but it is strange that network.raw has no row names. 

## node removal simulation 
#input network matrix, percentage of randomly removed species, and ra of all species
#return the proportion of species remained

Weighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = T, nperm = 100)

Unweighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = F, nperm = 100)

rmsimu_rs_bs_M <- data.frame(Proportion.removed = rep(seq(0.05, 1, by = 0.05), 2), 
                              rbind(Weighted.simu, Unweighted.simu), 
                              weighted = rep(c("weighted", "unweighted"), each = 20))

#### ~~~~~~~~~~~~~~~~~~~ ####
# ~~ Combine results ----

# add group name to data then combine
# Rhizosphere
rmsimu_rs_rhs_cd$compartment <- 'Rhizosphere'
rmsimu_rs_rhs_cd$group <- 'Cd'
rmsimu_rs_rhs_M_cd$compartment <- 'Rhizosphere'
rmsimu_rs_rhs_M_cd$group <- 'M+Cd'
rmsimu_rs_rhs_ctrl$compartment <- 'Rhizosphere'
rmsimu_rs_rhs_ctrl$group <- 'Ctrl'
rmsimu_rs_rhs_M$compartment <- 'Rhizosphere'
rmsimu_rs_rhs_M$group <- 'M'

# Hyphosphere
rmsimu_rs_hps_cd$compartment <- 'Hyphosphere'
rmsimu_rs_hps_cd$group <- 'Cd'
rmsimu_rs_hps_M_cd$compartment <- 'Hyphosphere'
rmsimu_rs_hps_M_cd$group <- 'M+Cd'
rmsimu_rs_hps_ctrl$compartment <- 'Hyphosphere'
rmsimu_rs_hps_ctrl$group <- 'Ctrl'
rmsimu_rs_hps_M$compartment <- 'Hyphosphere'
rmsimu_rs_hps_M$group <- 'M'

# Hyphae
rmsimu_rs_hp_M$compartment <- 'Hyphae'
rmsimu_rs_hp_M$group <- 'M'
rmsimu_rs_hp_M_cd$compartment <- 'Hyphae'
rmsimu_rs_hp_M_cd$group <- 'M+Cd'

# Bulk soil
rmsimu_rs_bs_cd$compartment <- 'Bulk soil'
rmsimu_rs_bs_cd$group <- 'Cd'
rmsimu_rs_bs_M_cd$compartment <- 'Bulk soil'
rmsimu_rs_bs_M_cd$group <- 'M+Cd'
rmsimu_rs_bs_ctrl$compartment <- 'Bulk soil'
rmsimu_rs_bs_ctrl$group <- 'Ctrl'
rmsimu_rs_bs_M$compartment <- 'Bulk soil'
rmsimu_rs_bs_M$group <- 'M'

# combine results of rhizosphere, hyphosphere, hyphae, and bulk soil
rmsimu_rs <- rbind(rmsimu_rs_rhs_ctrl, rmsimu_rs_rhs_cd,
                   rmsimu_rs_rhs_M, rmsimu_rs_rhs_M_cd, 
                   
                   rmsimu_rs_hps_ctrl, rmsimu_rs_hps_cd,
                   rmsimu_rs_hps_M, rmsimu_rs_hps_M_cd,
                   
                   rmsimu_rs_hp_M, rmsimu_rs_hp_M_cd,
                   
                   rmsimu_rs_bs_ctrl, rmsimu_rs_bs_cd,
                   rmsimu_rs_bs_M, rmsimu_rs_bs_M_cd)

# set factor levels
rmsimu_rs$group <- factor(rmsimu_rs$group, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))
levels(rmsimu_rs$group)

# ~~ Plotting ----
# plot line and pointrange using weighted data, compartment Rhizosphere, y = remain.mean, x = Proportion.removed, color = group
# ~~~~ weighted ----
p_simu_rhs <- ggplot(rmsimu_rs[rmsimu_rs$weighted=="weighted" & rmsimu_rs$compartment == 'Rhizosphere',], aes(x = Proportion.removed, y = remain.mean, group = group, color = group)) + 
  geom_line()+
  geom_pointrange(aes(ymin = remain.mean - remain.sd, 
                      ymax = remain.mean + remain.sd),
                  size=0.05)+
  ggtitle("Rhizosphere (weighted)")+
  scale_color_manual(values = c("#dfc27d","#a6611a", "#80cdc1", "#018571"))+
  xlab("Ratio of node removed")+
  ylab("Ratio of node remained")+
  theme_bw()+
  # remove legend title and move legend within plot
  theme_bw()+
  theme(legend.position = "none")

# plot compartment Hydrosphere
p_simu_hps <- ggplot(rmsimu_rs[rmsimu_rs$weighted=="weighted" & rmsimu_rs$compartment == 'Hyphosphere',], aes(x = Proportion.removed, y = remain.mean, group = group, color = group)) + 
  geom_line()+
  geom_pointrange(aes(ymin = remain.mean - remain.sd, 
                      ymax = remain.mean + remain.sd),
                  size=0.05)+
  ggtitle("Hyphosphere (weighted)")+
  scale_color_manual(values = c("#dfc27d","#a6611a", "#80cdc1", "#018571"))+
  xlab("Ratio of node removed")+
  ylab("Ratio of node remained")+
  theme_bw()+
  theme(legend.position = "none")

p_simu_bs <- ggplot(rmsimu_rs[rmsimu_rs$weighted=="weighted" & rmsimu_rs$compartment == 'Bulk soil',], aes(x = Proportion.removed, y = remain.mean, group = group, color = group)) + 
  geom_line()+
  geom_pointrange(aes(ymin = remain.mean - remain.sd, 
                      ymax = remain.mean + remain.sd),
                  size=0.05)+
  ggtitle("Bulk soil (weighted)")+
  scale_color_manual(values = c("#dfc27d","#a6611a", "#80cdc1", "#018571"))+
  xlab("Ratio of node removed")+
  ylab("Ratio of node remained")+
  theme(legend.title = element_blank(),
        legend.position = c(0.75, 0.72),
        legend.key.height = unit(0.8, "lines"),
        legend.key.width = unit(1, "lines")) 

  

p_simu_hp <- ggplot(rmsimu_rs[rmsimu_rs$weighted=="weighted" & rmsimu_rs$compartment == 'Hyphae',], aes(x = Proportion.removed, y = remain.mean, group = group, color = group)) + 
  geom_line()+
  geom_pointrange(aes(ymin = remain.mean - remain.sd, 
                      ymax = remain.mean + remain.sd),
                  size=0.05)+
  ggtitle("Hyphae (weighted)")+
  scale_color_manual(values = c("#80cdc1", "#018571"))+
  xlab("Ratio of node removed")+
  ylab("Ratio of node remained")+
  theme_bw()+
  theme(legend.position = "none")


# combine plot using ggarrange and save the plot as jpg and pdf using ggsave(); add labels a and b.
p_simu_w <- ggarrange(p_simu_bs, p_simu_rhs, 
                    p_simu_hps, p_simu_hp, 
                    ncol = 2, nrow = 2, 
                    labels = c('a', 'b', 'c', 'd'))
p_simu_w
# ~~~~ Extended Data Fig. 2. Ratio of node remained after random removal of nodes unweighted ----
if(save){
  ggsave("out/ExtDataFig. 2. node_rm_simu_w.jpg", p_simu_w, width = 88*1.5, height = 66*2, units = "mm")
  ggsave("out/ExtDataFig. 2. node_rm_simu_w.pdf", p_simu_w, width = 88*1.5, height = 66*2, units = "mm")
}

# ~~~~ unweighted ----
# plot the unweighted data for Rhizosphere 
p_simu_rhs_uw <- ggplot(rmsimu_rs[rmsimu_rs$weighted=="unweighted" & rmsimu_rs$compartment == 'Rhizosphere',], aes(x = Proportion.removed, y = remain.mean, group = group, color = group)) + 
  geom_line()+
  geom_pointrange(aes(ymin = remain.mean - remain.sd, 
                      ymax = remain.mean + remain.sd),
                  size=0.05)+
  ggtitle("Rhizosphere (u.w.)")+
  scale_color_manual(values = c("#dfc27d","#a6611a", "#80cdc1", "#018571"))+
  xlab("Ratio of node removed")+
  ylab("Ratio of node remained")+
  theme_bw()+
  # remove legend title and move legend within plot
  theme_bw()+
  theme(legend.position = "none")

# plot compartment Hydrosphere
p_simu_hps_uw <- ggplot(rmsimu_rs[rmsimu_rs$weighted=="unweighted" & rmsimu_rs$compartment == 'Hyphosphere',], aes(x = Proportion.removed, y = remain.mean, group = group, color = group)) + 
  geom_line()+
  geom_pointrange(aes(ymin = remain.mean - remain.sd, 
                      ymax = remain.mean + remain.sd),
                  size=0.05)+
  ggtitle("Hyphosphere (u.w.)")+
  scale_color_manual(values = c("#dfc27d","#a6611a", "#80cdc1", "#018571"))+
  xlab("Ratio of node removed")+
  ylab("Ratio of node remained")+
  theme_bw()+
  theme(legend.position = "none")

p_simu_bs_uw <- ggplot(rmsimu_rs[rmsimu_rs$weighted=="unweighted" & rmsimu_rs$compartment == 'Bulk soil',], aes(x = Proportion.removed, y = remain.mean, group = group, color = group)) + 
  geom_line()+
  geom_pointrange(aes(ymin = remain.mean - remain.sd, 
                      ymax = remain.mean + remain.sd),
                  size=0.05)+
  ggtitle("Bulk soil (u.w.)")+
  scale_color_manual(values = c("#dfc27d","#a6611a", "#80cdc1", "#018571"))+
  xlab("Ratio of node removed")+
  ylab("Ratio of node remained")+
  theme(legend.title = element_blank(),
        legend.position = c(0.75, 0.72),
        legend.key.height = unit(0.8, "lines"),
        legend.key.width = unit(1, "lines"))



p_simu_hp_uw <- ggplot(rmsimu_rs[rmsimu_rs$weighted=="unweighted" & rmsimu_rs$compartment == 'Hyphae',], aes(x = Proportion.removed, y = remain.mean, group = group, color = group)) + 
  geom_line()+
  geom_pointrange(aes(ymin = remain.mean - remain.sd, 
                      ymax = remain.mean + remain.sd),
                  size=0.05)+
  ggtitle("Hyphae (u.w.)")+
  scale_color_manual(values = c("#80cdc1", "#018571"))+
  xlab("Ratio of node removed")+
  ylab("Ratio of node remained")+
  theme_bw()+
  theme(legend.position = "none")


# combine plot using ggarrange and save the plot as jpg and pdf using ggsave(); add labels a and b.
# ~~~~ Extended Data Fig. 2. Ratio of node remained after random removal of nodes unweighted ----
p_simu_uw <- ggarrange(p_simu_bs_uw, p_simu_rhs_uw, 
                    p_simu_hps_uw, p_simu_hp_uw, 
                    ncol = 2, nrow = 2, 
                    labels = c('e', 'f', 'g', 'h'))
p_simu_uw
if(save){
  ggsave("out/ExtDataFig. 2. node_rm_simu_uw.jpg", p_simu_uw, width = 88*1.5, height = 66*2, units = "mm")
  ggsave("out//ExtDataFig. 2. node_rm_simu_uw.pdf", p_simu_uw, width = 88*1.5, height = 66*2, units = "mm")
}

# ~~ End ~~~~~~~~~~~~~~ ----
