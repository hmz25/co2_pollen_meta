# Meta-analysis on the effects of increased CO2 on pollen and reproduction
# Study authors: Hannah Zonnevylle, Allison Kozak, Dan Katz
# Summer 2023 - Spring 2026


# set up work environment -------------------------------------------------

library(googlesheets4)
library(metafor)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(janitor)
library(forcats)
library(patchwork)
library(orchaRd)
library(dmetar) #install dmetar: https://dmetar.protectlab.org/
library(rotl)
library(ape)
library(ggcorrplot)
library(maps)


# devtools::install_github("daniel1noble/orchaRd", ref = "main", force = TRUE)
pak::pak("daniel1noble/orchaRd")
pacman::p_load(devtools, tidyverse, metafor, patchwork, R.rsp, orchaRd, emmeans,
               ape, phytools, flextable)


# read in + clean data ----------------------------------------------------
data_url <- "https://docs.google.com/spreadsheets/d/1Xlvh1YfJ3H5yebCseq1KC_i0rE4HW86keOMqpLNKZrs/edit?usp=sharing"

p_raw <- googlesheets4::read_sheet(data_url, sheet ="data", .name_repair = "universal")

#filter data without SD + mean, calculate effect sizes 
p <- p_raw %>% 
  filter(!is.na(eCO2.SD)) %>% 
  filter(!is.na(eCO2.mean)) %>% 
  mutate(lnR = round(log(eCO2.mean) - log(aCO2.mean), 3),
         lnR = case_when(measurement.type == "start of reproduction" ~ lnR * -1,  #ES should be flipped for this one since an earlier start date 
                         measurement.type != "start of reproduction" ~ lnR), #means more pollen exposure
         vlnR = (eCO2.SD^2 / (eCO2.n * eCO2.mean^2)) +
           (aCO2.SD^2 / (aCO2.n * aCO2.mean^2)),
         sdlnR = sqrt(vlnR),
         ES_ratio = exp(lnR),
         ES_sd = exp(sdlnR)
  ) 

unique(p_raw$measurement.type)

tabyl(p_raw, measurement.type)

tabyl(p_raw, Experiment.Type)

# summary stats -----------------------------------------------------------

#number of papers
length(unique(p$paper.index))

#number of unique observations
nrow(p)

#number of different species 
length(unique(p$species))

#how many observations for each species
p %>% 
  group_by(species) %>%
  summarize(n = n())

#mean response ratio
mean(p$lnR[!is.infinite(p$lnR)], na.rm = T)

#how many studies were grown with neighbors 
unique(p$Experiment.Type)
p |> 
  filter(Experiment.Type != "FACE") |> 
  count(neighbors)

#how many studies were included in each response type
p |> 
  distinct(paper.index, measurement.type) |> 
  count(measurement.type)

#how many studies were from FACE experiments
length(p$Experiment.Type[p$Experiment.Type == "FACE" & !is.na(p$aCO2.n.indv.plants)])

#how many studies were not from FACE experiments
length(p$Experiment.Type[p$Experiment.Type != "FACE"])

# p_indiv <- p %>% filter(!is.na(aCO2.n.indv.plants))
# p_area <- p %>% filter(is.na(aCO2.n.indv.plants))
# p_area <- p %>%  filter(Experiment.Type == "FACE")
# p_area %>% group_by(study.name) %>% 
#   summarize(n = n())

#calculate n, mean response ratio, and confidence interval for each response type
p |>
  filter(!is.infinite(lnR)) |>
  group_by(measurement.type) |>
  summarize(
    n_obs = n(),
    n_papers = n_distinct(paper.index),
    mean_lnR = mean(lnR, na.rm = TRUE),
    ci_lower = round(t.test(lnR, conf.level = 0.95)$conf.int[1], 3),
    ci_upper = round(t.test(lnR, conf.level = 0.95)$conf.int[2], 3)
  )

# p |>
#   filter(!is.infinite(lnR)) |>
#   group_by(measurement.type) |>
#   summarize(
#     n_obs = n(),
#     n_papers = n_distinct(paper.index),
#     mean_lnR = mean(lnR, na.rm = TRUE),
#     ci_lower = mean(lnR) - mean(sdlnR),
#     ci_upper = mean(lnR) + mean(sdlnR)
#   )

# data viz ----------------------------------------------------------------

#most basic results
ggplot(p, aes(x = lnR)) + geom_histogram() + theme_bw() + 
  geom_vline(xintercept = 0, lty = 2, color = "red", lwd = 1.3) +
  geom_vline(xintercept = mean(p$lnR[!is.infinite(p$lnR)], na.rm = T), lty = 2, color = "blue", lwd = 1.3)

#back transformed mean of lnR
exp(mean(p$lnR[!is.infinite(p$lnR)], na.rm = T))

p %>% 
  ggplot(aes(x = lnR)) + geom_histogram() + theme_bw() + facet_wrap(~Experiment.Type)

p %>% #group_by(Experiment.Type) %>% 
  filter(measurement.type == "allergenicity") %>% 
  filter(!is.infinite(lnR)) %>% 
  summarize(mean_lnR = mean(lnR, na.rm = TRUE),
            mean_r = exp(mean_lnR))

p %>% #group_by(Experiment.Type) %>% 
  filter(measurement.type == "amount of reproductive tissue") %>% 
  ggplot(aes(x = lnR)) + geom_histogram() + theme_bw() + 
  geom_vline(xintercept = 0, lty = 2, color = "red", lwd = 1.3) +
  geom_vline(xintercept = mean(p$lnR[!is.infinite(p$lnR) & p$measurement.type == "amount of reproductive tissue"], na.rm = T), lty = 2, color = "blue", lwd = 1.3)

p %>% 
  #filter(Experiment.Type == "FACE") %>%  
  mutate(study_obs = as.numeric(as.factor(paper.index))) %>% 
  ggplot(aes(y = study_obs, xmin = lnR - sdlnR, x = lnR, xmax = lnR + sdlnR,
             col = Growth.Form)) +  #Experiment.Type wind.pollinated Growth.Form photosynthesis.type Country
  geom_pointrange(alpha = 0.5)  + facet_wrap(~measurement.type, scales = "free_x") + 
  geom_vline(xintercept = 0, lty = 2) + theme_bw()

length(unique(p$study.name))

names(p)


# examining measurement types ------------------------------------------------

unique(p$measurement.type)

#allergenicity
p %>% 
  filter(measurement.type == "allergenicity") %>% 
  mutate(study_obs = as.numeric(as.factor(paper.index)),
         row_n = row_number(),
         study_obs_n = paste(study.name, row_n)) %>% 
  ggplot(aes(y = study_obs_n, xmin = lnR - sdlnR, x = lnR, xmax = lnR + sdlnR,
             col = species)) +  #Experiment.Type wind.pollinated Growth.Form photosynthesis.type Country
  geom_pointrange(alpha = 0.5)  + facet_wrap(~measurement.type, scales = "free") + 
  geom_vline(xintercept = 0, lty = 2) + ggthemes::theme_few() 

#pollen size
p %>% 
  filter(measurement.type == "pollen size") %>% 
  mutate(study_obs = as.numeric(as.factor(paper.index)),
         row_n = row_number(),
         study_obs_n = paste(study.name, row_n)) %>% 
  ggplot(aes(y = study_obs_n, xmin = lnR - sdlnR, x = lnR, xmax = lnR + sdlnR, col = species)) +  #Experiment.Type wind.pollinated Growth.Form photosynthesis.type Country
  geom_pointrange(alpha = 0.5)  + facet_wrap(~measurement.type, scales = "free") + 
  geom_vline(xintercept = 0, lty = 2) + ggthemes::theme_few() 

#reproductive tissue production
p %>% 
  filter(measurement.type == "amount of reproductive tissue") %>% 
  mutate(study_obs = as.numeric(as.factor(paper.index)),
         row_n = row_number(),
         study_obs_n = paste(study.name, row_n)) %>% 
  ggplot(aes(y = study_obs_n, xmin = lnR - sdlnR, x = lnR, xmax = lnR + sdlnR, color = Experiment.Type)) +  #Experiment.Type wind.pollinated Growth.Form photosynthesis.type Country
  geom_pointrange(alpha = 0.5)  + facet_wrap(~measurement.type, scales = "free") + 
  geom_vline(xintercept = 0, lty = 2) + ggthemes::theme_few() 

p %>% 
  filter(measurement.type == "amount of reproductive tissue") %>% 
  mutate(study_obs = as.numeric(as.factor(paper.index)),
         row_n = row_number(),
         study_obs_n = paste(study.name, row_n)) %>% 
  ggplot(aes(y = study_obs_n, xmin = lnR - sdlnR, x = lnR, xmax = lnR + sdlnR)) +  #Experiment.Type wind.pollinated Growth.Form photosynthesis.type Country
  geom_pointrange(alpha = 0.5)  + facet_wrap(~measurement.type, scales = "free") + 
  geom_vline(xintercept = 0, lty = 2) + ggthemes::theme_few() 

p |> 
  filter(measurement.type == "amount of reproductive tissue") |> 
  group_by(study.name, species) |> 
  summarize(num_study = n()) |> 
  View()

p %>% 
  filter(measurement.type == "amount of reproductive tissue") %>% 
  mutate(study_obs = as.numeric(as.factor(paper.index)),
         row_n = row_number(),
         study_obs_n = paste(study.name, row_n)) %>% 
  ggplot(aes(x = lnR)) +  #Experiment.Type wind.pollinated Growth.Form photosynthesis.type Country
  
  geom_histogram() + theme_bw() + 
  geom_vline(xintercept = 0, lty = 2, color = "red", lwd = 1.3) +
  geom_vline(xintercept = mean(p$lnR[!is.infinite(p$lnR)], na.rm = T), lty = 2, color = "blue", lwd = 1.3)

#phenology: reproductive duration and start of reproduction
p %>% 
  filter(measurement.type == "reproductive duration" | measurement.type == "start of reproduction") %>% 
  mutate(study_obs = as.numeric(as.factor(paper.index)),
         row_n = row_number(),
         study_obs_n = paste(study.name, row_n)) %>% 
  ggplot(aes(y = study_obs_n, xmin = lnR - sdlnR, x = lnR, xmax = lnR + sdlnR)) +  
  #Experiment.Type wind.pollinated Growth.Form photosynthesis.type Country
  geom_pointrange(alpha = 0.5)  + facet_wrap(~measurement.type, scales = "free") + 
  geom_vline(xintercept = 0, lty = 2) + ggthemes::theme_few() 

#reproductive allocation
p %>% 
  filter(measurement.type == "reproductive allocation") %>% 
  mutate(study_obs = as.numeric(as.factor(paper.index)),
         row_n = row_number(),
         study_obs_n = paste(study.name, row_n)) %>% 
  ggplot(aes(y = study_obs_n, xmin = lnR - sdlnR, x = lnR, xmax = lnR + sdlnR, col = species)) +  
  #Experiment.Type wind.pollinated Growth.Form photosynthesis.type Country
  geom_pointrange(alpha = 0.5)  + facet_wrap(~measurement.type, scales = "free") + 
  geom_vline(xintercept = 0, lty = 2) + ggthemes::theme_few() 

#percent flowering
p %>% 
  filter(measurement.type == "percent flowering") %>% 
  mutate(study_obs = as.numeric(as.factor(paper.index)),
         row_n = row_number(),
         study_obs_n = paste(study.name, row_n)) %>% 
  ggplot(aes(y = study.name, xmin = lnR - sdlnR, x = lnR, xmax = lnR + sdlnR)) +  
  #Experiment.Type wind.pollinated Growth.Form photosynthesis.type Country
  geom_pointrange(alpha = 0.5)  + facet_wrap(~measurement.type, scales = "free") + 
  geom_vline(xintercept = 0, lty = 2) + ggthemes::theme_few() 


# meta-analysis of different measurement types ------


## calculate effect sizes ---------------------------------

#clean data
p_subset <- p_raw |>  
  filter(!is.na(eCO2.SD)) |> 
  filter(!is.na(eCO2.mean)) |>  
  filter(eCO2.mean != 0) |> 
  filter(aCO2.mean != 0) |>  
  filter(eCO2.SD != 0) |> 
  filter(aCO2.SD != 0) |>  
  filter(!is.na(eCO2.n)) |>  
  # filter(Experiment.Type != "FACE") |> #include this line for sensitivity analysis to exclude FACE studies
  mutate(lnR = round(log(eCO2.mean) - log(aCO2.mean), 3), #calculate effect size
         lnR = case_when(measurement.type == "start of reproduction" ~ lnR * -1,  #ES should be flipped for this one since an earlier start date 
                         measurement.type != "start of reproduction" ~ lnR), #means more pollen exposure
         vlnR = eCO2.SD^2 / (eCO2.n * eCO2.mean^2) +
           aCO2.SD^2 / (aCO2.n * aCO2.mean^2), #calculate variance
         sdlnR = sqrt(vlnR), #calculate standard deviation
         ES_ratio = exp(lnR), #back transform effect size
         ES_sd = exp(sdlnR), #back transform standard deviation
         obs = 1:n(), #add observation number
         dif_co2 = eCO2 - aCO2,
         study_n = as.numeric(as.factor(study.name)),
         authors = sub("\\d.*", "", study.name),
         effect_id = seq.int(n()), #add column for effect size ID
         Experiment.Type = ifelse(Experiment.Type == "open chamber", "FACE", Experiment.Type) #group OTC and FACE
         ) |> 
  select(study.name, paper.index, Country, species, Growth.Form:eCO2.n, measurement.type:effect_id) #select relevant columns

# unique(p_subset$study.name)

#calculate effect sizes and corresponding sampling variances
dat <- escalc(measure="ROM", m1i = eCO2.mean, m2i = aCO2.mean, sd1i = eCO2.SD, sd2i = aCO2.SD, 
              n1i = eCO2.n, n2i = aCO2.n,
              data = p_subset) #ROM = ratio of means log(m1i/m2i), symmetric around 0
#yi = effect size 
#vi = sampling variances

#names(dat)
# str(dat)


## generate phylogeny ------------------------------------------------------

#check species in data
length(unique(dat$species)) #95 unique species names, if no misspelling

#find Open Tree Taxonomy (OTT) IDs for each species
taxa <- tnrs_match_names(names = unique(dat$species), context = "Land plants")

#Austrostipa spp., Leontodon-Hypochaeris, Vulpia spp are not matched

#all from Hovenden et al 2007, sent original data 
#Leontodon taraxacoides = 25, Hypochaeris radicata = 11, Hypochaeris glabra = 2
#Vulpia bromoides only one listed in data
dat$species <- ifelse(dat$species =="Austrostipa spp.","Austrostipa mollis", dat$species) #mollis was only species in original data
dat$species <- ifelse(dat$species == "Leontodon-Hypochaeris", "Leontodon taraxacoides", dat$species) #Leontodon taraxacoides is more abundant
dat$species <-  ifelse(dat$species == "Vulpia spp.", "Vulpia bromoides", dat$species)

#try again
taxa <- tnrs_match_names(names = unique(dat$species), context = "Land plants") #no warning :-)
length(taxa$unique_name) #95 unique species names - all matched
tabyl(taxa,approximate_match) #10 approximate match

taxa %>% filter(approximate_match == "TRUE") 
#leptorhynchos sqamatus = misspelled --> Leptorhynchos squamatus
#austrodanthonia spp. should be Rytidosperma caespitosum based on data
#hypoachaeris radicata = misspelled --> Hypochaeris radicata
#leontodon taraxaxoides = misspelled --> Leontodon taraxacoides
#geranium retrotsum = misspelled --> Geranium retrorsum
#stipagrostis hirtiglumis = typo --> Stipagrostis hirtigluma
#styloanthes capitata = misspelled --> Stylosanthes capitata
#anthrozanthum odoratum = misspelled --> Anthoxanthum odoratum
#linum cathaticum = misspelled --> Linum catharticum
#agrotis capillaris = misspelled --> Agrostis capillaris

#fix the typos
dat$species <- ifelse(dat$species=="Leptorhynchos sqamatus","Leptorhynchos squamatus", dat$species)
dat$species <- ifelse(dat$species=="Austrodanthonia spp.","Rytidosperma caespitosum", dat$species)
dat$species <- ifelse(dat$species=="Hypoachaeris radicata","Leontodon taraxacoides", dat$species)
dat$species <- ifelse(dat$species=="Leontodon taraxaxoides","Hypochaeris radicata", dat$species)
dat$species <- ifelse(dat$species=="Geranium retrotsum","Geranium retrorsum", dat$species)
dat$species <- ifelse(dat$species=="Stipagrostis hirtiglumis","Stipagrostis hirtigluma", dat$species)
dat$species <- ifelse(dat$species=="Styloanthes capitata","Stylosanthes capitata", dat$species)
dat$species <- ifelse(dat$species=="Anthrozanthum odoratum","Anthoxanthum odoratum", dat$species)
dat$species <- ifelse(dat$species=="Linum cathaticum","Linum catharticum", dat$species)
dat$species <- ifelse(dat$species=="Agrotis capillaris","Agrostis capillaris", dat$species)

taxa <- tnrs_match_names(names = unique(dat$species), context = "Land plants")
taxa %>% filter(approximate_match == "TRUE") #all good

#any flag
tabyl(taxa,flags) 
#1 extinct_inherited, incertae_sedis
#1 hybrid
#1 incertae_sedis
#1 incertae_sedis_inherited
#7 sibling_higher 

# filter(taxa, flags == "incertae_sedis")$ott_id
# inspect(taxa, ott_id = "816902")

#any synonym
tabyl(taxa,is_synonym) #11
taxa %>% filter(is_synonym == "TRUE") # this is problematic, causing problems, which were resolved below

#check whether all otts occur in the synthetic tree
ott_in_tree <- ott_id(taxa)[is_in_tree(ott_id(taxa))]  # all good
length(ott_in_tree) #89 

#trim tree corresponding to our taxa
my.tr <- tol_induced_subtree(ott_ids=ott_in_tree)
my.tr 
#87 tips

#check how many unique species now in data
length(unique(dat$species)) #91 

#some species in the data set (search_string) were mapped as the same species in the synthetic tree (unique_name or tip.label)?

#find the duplicated species
taxa[duplicated(taxa$unique_name) | duplicated(taxa$unique_name, fromLast=T),] # 4 duplicates need to resolve; https://stackoverflow.com/questions/16905425/find-duplicate-values-in-r

# # remove duplicates to see whether the number matches
# taxa %>%
#   distinct(unique_name, .keep_all = F)  |> 
#   nrow() #89 - the number seems correct

#fix the duplicates
dat$species <- ifelse(dat$species=="Austrodanthonia caespitosa","Rytidosperma caespitosum", dat$species)
dat$species <- ifelse(dat$species=="rytidosperma caespitosum","Rytidosperma caespitosum", dat$species)
dat$species <- ifelse(dat$species=="Setaria lutescens","Setaria pumila", dat$species)
dat$species <- ifelse(dat$species=="Setaria pumila","Setaria pumila", dat$species)

taxa <- tnrs_match_names(names = unique(dat$species), context = "Land plants")

#check if duplicates resolved
taxa[duplicated(taxa$unique_name) | duplicated(taxa$unique_name,fromLast=T),]

#check whether the number of species match well
length(unique(dat$species)) == length(taxa$unique_name)

#the tip labels contain OTTs, which means they will not perfectly match the species names in our dataset or the taxon map that we created earlier,
my.tr$tip.label[1:2] # remove the extra information from the tip labels later; with the IDs removed, we can use our taxon map to replace the tip labels in the tree with the species names from dataset 
my.tr$tip.label <- strip_ott_ids(my.tr$tip.label, remove_underscores = TRUE)

#test whether a tree is binary
is.binary(my.tr)  
# my.tr <- multi2di(my.tr, random = T)  # set.seed(2023) set a seed to resolve politomies at random, and obtain similar results; resolve polytomies at random

# write.tree(my.tr, file = "my.tr.tre")  #save the tree
# plot
plot.phylo(my.tr, cex = 0.6, label.offset = 0.1, no.margin = T)


# *NOTE:* underscores within species names on tree tip labals are added
# automatically tree <- read.tree(file='plot_cooked_fish_MA.tre') #if you need
# to read in the tree tree$tip.label <- gsub('_',' ', tree$tip.label) #get rid
# of the underscores tree$node.label <- NULL #you can delete internal node labels


# remove ott ID from species name to match it to the data set
# my.tr$tip.label <- strip_ott_ids(my.tr$tip.label, remove_underscores = TRUE)


#decapitalize species names to match with the search string names in taxa
dat <- dat %>% mutate(search_string = tolower(species))  

#align data
dat <- left_join(dat, dplyr::select(taxa, search_string, unique_name, ott_id), by = "search_string")  

#create the variables of phy and sp
dat <- dat %>% mutate(spp = search_string, phy = unique_name)

dat <- dat[dat$unique_name %in% my.tr$tip.label, ] 

#check whether tips in the tree match well with species in dataset 
intersect(as.character(my.tr$tip.label), unique(dat$phy)) |>  
  length() #87

# setdiff(taxa$unique_name,my.tr$tip.label)
#Pelargonium x hortorum and Trifolium montanum are in taxa$unique_name but not in my.tr$tip.label

#calculate phylogenetic correlation
# my.tr$tip.label <- gsub('_',' ', my.tr$tip.label) # _ is automatic after re-uploading

# roughly approximate branch lengths using default method (Grafen's method with power = 1)
my.tr <- compute.brlen(my.tr, method = "Grafen", power = 1)
#my.tr0.5 <- compute.brlen(my.tr, method = "Grafen", power = 1)
#my.tr2 <- compute.brlen(my.tr, method = "Grafen", power = 2)
# test whether a tree is ultrametric, based on the distances from each tip to the root
is.ultrametric(my.tr)

# compute phylogenetic correlation matrix of the melatonin assuming it evolves under a Brownian model (Felsenstein 1985, Martins and Hansen 1997); https://rdrr.io/cran/ape/man/corBrownian.html
# https://www.journals.uchicago.edu/doi/full/10.1086/303327
tr_matrix <- vcv.phylo(my.tr, model = "Brownian", corr = T)   
#tr_matrix0.5 <- vcv.phylo(my.tr0.5, model = "Brownian", corr = T)   
#tr_matrix2 <- vcv.phylo(my.tr2, model = "Brownian", corr = T)  

# visualize correlation
ggcorrplot::ggcorrplot(tr_matrix, sig.level=0.05, lab_size = 4.5, p.mat = NULL, insig = c("pch", "blank"), pch = 1, pch.col = "black", pch.cex =1, tl.cex = 14) +
  theme(axis.text.x = element_text(size = 10, margin=margin(-2,0,0,0)),
        axis.text.y = element_text(size = 10, margin=margin(0,-2,0,0)),
        panel.grid.minor = element_line(size=10)) + 
  geom_tile(fill="white") +
  geom_tile(height=0.8, width=0.8) +
  scale_fill_gradient2(low = "#E69F00", mid = "white", high = "#56B4E9", midpoint = 0.5, breaks=c(0, 1), limit=c(0, 1)) + 
  labs(fill = "Correlation")

# plot tree
plot.phylo(vcv2phylo(tr_matrix), cex = 0.6, label.offset = 0.1, no.margin = T)


## fit model ---------------------------------------------------------------

# res <- rma.mv(yi, vi, mods = ~ measurement.type - 1,
#               random = ~ 1| factor(study_n),
#               data=dat)

# res <- rma.mv(yi, vi, mods = ~ measurement.type - 1,
#               random = list(~1 | effect_id, ~1 | factor(study_n), ~1 | species),
#               data=dat)

res <- rma.mv(yi, vi, mods = ~ measurement.type - 1,
              random = list(~1 | effect_id, ~1 | factor(study_n)),
              data=dat)
#mods = measurement.type as fixed effect, -1 removes the intercept
#random effect of study_n 
#removing intercept = tests whether mean outcome is 0 for moderators, keeping intercept tests whether mean outcome is the same for moderators 

summary(res)
confint(res)
profile(res)
forest(res)

res_phylo <- rma.mv(yi, vi, mods = ~ measurement.type - 1,
                    random = list(~1 | spp, ~1 | phy, ~1 | study_n, ~1 | effect_id), 
                    method = "ML", test = "t", data = dat, sparse = T, 
                    R = list(phy=tr_matrix))

summary(res_phylo)

res_test <- rma.mv(yi, vi, mods = ~ measurement.type - 1 + spp,
                    random = list(~1 | study_n, ~1 | effect_id),
                    method = "ML", test = "t", data = dat, sparse = T)

summary(res_test)

# variance_decomposition <- function(m){
#   n <- m$k
#   vector.inv.var <- 1/(diag(m$V))
#   sum.inv.var <- sum(vector.inv.var)
#   sum.sq.inv.var <- (sum.inv.var)^2
#   vector.inv.var.sq <- 1/(diag(m$V)^2)
#   sum.inv.var.sq <- sum(vector.inv.var.sq)
#   num <- (n-1)*sum.inv.var
#   den <- sum.sq.inv.var - sum.inv.var.sq
#   est.samp.var <- num/den
#   if(length(m$sigma2) > 2) stop("Cannot handle more than three levels.")
#   total_var <- (sum(m$sigma2)+est.samp.var)/100
#   Variance <- c(est.samp.var, m$sigma2)/total_var
#   names(Variance) <- c("Level1", m$s.names)
#   Variance
# }
# 
# variance_decomposition(res)
# #0% variance for level 1 (sampling variance variance of observed effect size... sampling error)
# #36.4% for level 2 (within study variance...correlation btwn effect sizes in the same study)
# #63.1% for level 3 (between study variance)

## model selection -------------------------------------------------

#fit null model
null_mod <- rma.mv(yi, vi, mods = ~ measurement.type - 1,
                   method = "ML", test = "t", data = dat, sparse = TRUE )

summary(null_mod)

#add effect size ID as a random effect and compare it to the null model
mod_effectID <- rma.mv(yi, yi, mods = ~ measurement.type - 1,
                       random =  ~1 | effect_id, 
                       method = "ML", test = "t", data = dat, sparse = T)

summary(mod_effectID)

# #compare
# anova(null_mod,mod_effectID)

## summary statistics -----------------------------------------------

#how many different observations and how papers for each response type
dat |>
  filter(!is.infinite(lnR)) |>
  group_by(measurement.type) |>
  summarize(
    n_obs = n(),
    n_papers = n_distinct(paper.index))

#how many different observations for each species 
dat |>
  filter(!is.infinite(lnR)) |>
  group_by(species) |>
  summarize(
    n_obs = n()) |> 
  View()

#how many different authors
unique(dat$authors)

#create data frame of results to calculate percent change for each response type
res_df <- data.frame(
  estimate = round(exp(res$b),3),
  perc_change = round(((exp(res$b)-1)*100),3), #perc change = (exp(lnRR) - 1)*100
  ci.lb = round(exp(res$ci.lb),3),
  ci.ub = round(exp(res$ci.ub),3),
  pval  = round((res$pval),5)
)

res_df

## visualizations for measurement type analysis ------------------------------

# a figure for how response varies according to response type
data.frame(row.names(res$b), res$b[1:7], res$ci.lb, res$ci.ub) %>% 
  mutate( response_type = gsub(x = row.names.res.b., pattern = "measurement.type", replacement = "")) %>% 
  # filter(response_type != "percent flowering") %>% 
  # filter(response_type != "start of reproduction") %>% 
  # filter(response_type != "pollen size") %>% 
  # filter(response_type != "reproductive allocation") %>% 
  ggplot(aes(y = fct_rev(response_type), x = res.b.1.7., xmin = res.ci.lb, xmax = res.ci.ub)) + geom_errorbarh(height = 0.2) + geom_point() +
  ggthemes::theme_few(base_size = 16) +
  xlab("effect size (logRR)") + 
  geom_vline(xintercept = 0, lty = 2) + 
  ylab("response to CO2 enrichment")


#with back transformation to response ratio
data.frame(row.names(res$b), res$b[1:7], res$ci.lb, res$ci.ub) %>% 
  mutate( response_type = gsub(x = row.names.res.b., pattern = "measurement.type", replacement = "")) %>% 
  mutate(
    res.b.1.7. = case_when(response_type == "start of reproduction" ~ res.b.1.7. * -1,  #ES should be flipped for this one since an earlier start date
                           response_type != "start of reproduction" ~ res.b.1.7.),
    res.ci.lb = case_when(response_type == "start of reproduction" ~ res.ci.lb * -1,  #ES should be flipped for this one since an earlier start date
                          response_type != "start of reproduction" ~ res.ci.lb),
    res.ci.ub = case_when(response_type == "start of reproduction" ~ res.ci.ub * -1,  #ES should be flipped for this one since an earlier start date
                          response_type != "start of reproduction" ~ res.ci.ub)) %>%
  ggplot(aes(y = response_type, x = exp(res.b.1.7.), xmin = exp(res.ci.lb), xmax = exp(res.ci.ub))) + geom_errorbarh(height = 0.2) + geom_point() +
  ggthemes::theme_few(base_size = 16) +
  xlab("effect size (response ratio)") + geom_vline(xintercept = 1, lty = 2) + ylab(expression(response~to~CO[2]~enrichment)) 

#orchard plot visualization for how effect size varies by response type
p_measure_type <- orchaRd::orchard_plot(res, mod = "measurement.type", xlab = "lnRR", 
                      group = "study_n", cb = TRUE, k = TRUE, twig.size = 0.5, 
                      trunk.size = 0.4, branch.size = 2, angle = 0, N = res$aCO2.n, alpha = 0.5)

# ggsave("measurement_type_plot.png", p_measure_type, dpi = 300)

#forest plots 
forest(dat$yi, dat$vi,
       xlim=c(-2.5,3.5),        ### adjust horizontal plot region limits
       order="obs",             ### order by size of yi
       slab=NA, annotate=FALSE, ### remove study labels and annotations
       efac=0,                  ### remove vertical bars at end of CIs
       pch=19,                  ### changing point symbol to filled circle
       col="gray40",            ### change color of points/CIs
       psize=2,                 ### increase point size
       cex.lab=1, cex.axis=1,   ### increase size of x-axis title/labels
       lty=c("solid","blank"))  ### remove horizontal line at top of plot
#addpoly(res, mlab="", cex=1)

forest(res, atransf=exp, at=log(c(0.05, 0.25, 1, 4)), xlim=c(-5,6),
       #ilab=cbind(tpos, tneg, cpos, cneg), ilab.xpos=c(-9.5,-8,-6,-4.5),
       cex=0.75, header="Author(s) and Year", mlab="", shade=TRUE)

## publication bias ---------------------------------------------------------

#funnel plots to detect publication bias
funnel(res, main="Standard Error")

funnel(res, yaxis = "sei")


#meta regression with sampling variance and standard error added to fixed effects
# res_pubbias <- rma.mv(yi, vi, mods = ~ measurement.type - 1 + vi + standard_error,
#               random = ~ 1| factor(study_n),
#               data=dat)

#create standard error variable
dat$sei <- sqrt(dat$vi)

#precision-effect test
res_PET <- rma.mv(yi, vi, 
                  mods = ~ measurement.type + sei - 1,
                  random = list(~1 | effect_id, ~1 | factor(study_n)),
                  data = dat)

summary(res_PET)

#sampling variance test
res_PEESE <- rma.mv(yi, vi, 
                    mods = ~ measurement.type + vi - 1,
                    random = list(~1 | effect_id, ~1 | factor(study_n)),
                    data = dat)

summary(res_PEESE)

summary(allergen_pb)

# meta-analysis of studies measuring reproductive tissue ------------------

## fit models ----------------------------------------------------

# res <- rma.mv(yi, vi, subset = (measurement.type == "amount of reproductive tissue"), mods = ~  factor(Growth.Form) + factor(N2.Fixing) + factor(wind.pollinated) +
#                 factor(neighbors) + 
#                 factor(Experiment.Type),
#               random = ~ 1| factor(study_n),
#               data = dat)

# res <- rma.mv(yi, vi, subset = (measurement.type == "amount of reproductive tissue"), mods = ~  factor(Growth.Form) + factor(N2.Fixing) + factor(wind.pollinated) +
#                 factor(neighbors) +
#                 factor(Experiment.Type),
#               random = list(~1 | effect_id, ~1 | factor(study_n), ~1 | species),
#               data = dat)

res <- rma.mv(yi, vi, subset = (measurement.type == "amount of reproductive tissue"),
              mods = ~  factor(Growth.Form) + factor(N2.Fixing) + factor(wind.pollinated) +
                factor(neighbors) + 
                factor(Experiment.Type),
              random = list(~1 | spp, ~1 | phy, ~1 | study_n, ~1 | effect_id), 
              method = "ML", test = "t", data = dat, sparse = T, 
              R = list(phy=tr_matrix))

summary(res)

# res_test <- rma.mv(yi, vi, subset = (measurement.type == "amount of reproductive tissue"),
#               mods = ~  factor(Growth.Form) + factor(N2.Fixing) + factor(wind.pollinated) +
#                 factor(neighbors) + 
#                 factor(Experiment.Type),
#               random = list(~1 | study_n, ~1 | effect_id), 
#               method = "ML", test = "t", data = dat, sparse = T)
# 
# summary(res_test)

## orchard plots -------------------------------------------------------

p_growthform <- orchaRd::orchard_plot(res, mod = "Growth.Form", xlab = "lnRR", group = "study_n", 
                      transfm = "none", twig.size = 0.5, trunk.size = .8, branch.size = 2, angle = 0) +
  ggtitle("growth form")

p_nfix <- orchaRd::orchard_plot(res, mod = "N2.Fixing", xlab = "lnRR", group = "study_n", 
                      transfm = "none", twig.size = 0.5, trunk.size = .8, branch.size = 2, angle = 0) + 
  ggtitle("nitrogen-fixing capacity")

p_wind <- orchaRd::orchard_plot(res, mod = "wind.pollinated", xlab = "lnRR", group = "study_n", 
                      transfm = "none", twig.size = 0.5, trunk.size = .8, branch.size = 2, angle = 0) + 
  ggtitle("anemophilous")

p_experiment <- orchaRd::orchard_plot(res, mod = "Experiment.Type", xlab = "lnRR", group = "study_n", 
                                      transfm = "none", twig.size = 0.5, trunk.size = .8, branch.size = 2, angle = 0) + 
  ggtitle("experiment type")

orchaRd::orchard_plot(res, mod = "Experiment.Type", xlab = "lnRR", group = "study_n", 
                      transfm = "none", twig.size = 0.5, trunk.size = .8, branch.size = 2, angle = 0)

orchaRd::orchard_plot(res, mod = "neighbors", xlab = "lnRR", group = "study_n", 
                      transfm = "none", twig.size = 0.5, trunk.size = .8, branch.size = 2, angle = 0)

p_growthform/p_nfix/p_wind/p_experiment

# back trasnformed box plot for how response varies according to plant type
data.frame(row.names(res$b), res$b, res$ci.lb, res$ci.ub) %>% 
  mutate( response_type = gsub(x = row.names.res.b., pattern = "measurement.type", replacement = "")) %>% 
  mutate( response_type = gsub(x = response_type, pattern = "factor(", replacement = "", fixed = TRUE)) %>% 
  mutate( response_type = gsub(x = response_type, pattern = "Growth.Form)", replacement = "", fixed = TRUE)) %>% 
  mutate( response_type = gsub(x = response_type, pattern = ")yes", replacement = "", fixed = TRUE)) %>% 
  ggplot(aes(y = response_type, x = exp(res.b), xmin = exp(res.ci.lb), xmax = exp(res.ci.ub))) + geom_errorbarh(height = 0.2) + geom_point() +
  ggthemes::theme_few(base_size = 16) +
  xlab("effect size (RR)") + geom_vline(xintercept = 1, lty = 2) + ylab(expression(plant~type)) 

## forest plot ---------------------------------------------------------------
forest(res)

forest(dat$yi, dat$vi,
       xlim=c(-2.5,3.5),        ### adjust horizontal plot region limits
       order="obs",             ### order by size of yi
       slab=NA, annotate=FALSE, ### remove study labels and annotations
       efac=0,                  ### remove vertical bars at end of CIs
       pch=19,                  ### changing point symbol to filled circle
       col="gray40",            ### change color of points/CIs
       psize=2,                 ### increase point size
       cex.lab=1, cex.axis=1,   ### increase size of x-axis title/labels
       lty=c("solid","blank"))  ### remove horizontal line at top of plot
#addpoly(res, mlab="", cex=1)

forest(res, atransf=exp, at=log(c(0.05, 0.25, 1, 4)), xlim=c(-5,6),
       #ilab=cbind(tpos, tneg, cpos, cneg), ilab.xpos=c(-9.5,-8,-6,-4.5),
       cex=0.75, header="Author(s) and Year", mlab="", shade=TRUE)

## publication bias ---------------------------------------------------------

#funnel plot
funnel(res, main="Standard Error")

funnel(res, yaxis = "seinv")

#run multilevel meta regression with sampling variance added to fixed effects
res <- rma.mv(yi, vi, mods = ~  factor(Growth.Form) + factor(N2.Fixing) + factor(wind.pollinated) +
                factor(neighbors) + 
                factor(Experiment.Type) + vi,
              random = list(~1 | effect_id, ~1 | factor(study_n), ~1 | species),
              data = dat)

summary(res)


# phylogenetic sensitivity analysis ---------------------------------------

# res <- rma.mv(yi, vi, mods = ~ measurement.type - 1,
#               random = list(~1 | effect_id, ~1 | factor(study_n), ~1 | species, ~1 | phylo),
#               R = list(phylo = cor),
#               control = list(optimizer = "optim"),
#               data=dat)

# # construct phylogenetic tree matrix for use as random factor
# species <- unique(p_subset$species) # list of unique species in meta-analysis
# species <- as.character(species) # change to character object
# taxa <- tnrs_match_names(species) |> 
#   drop_na(unique_name) #dropped unmatched species, maybe revisit this later 
# 
# #filter to only IDs that are in the synthetic tree
# taxa <- taxa |>
#   mutate(in_tree = map_lgl(ott_id, function(id) {
#     tryCatch({
#       info <- rotl::taxonomy_taxon_info(id, include_lineage = FALSE)
#       flags <- info[[1]]$flags
#       # Exclude taxa flagged as problematic for the synthetic tree
#       !any(flags %in% c("pruned_ott_id", "incertae_sedis", 
#                         "major_rank_conflict", "unplaced",
#                         "extinct", "hybrid", "viral"))
#     }, error = function(e) FALSE)
#   }))
# 
# #see which species were dropped
# dropped <- taxa |> filter(!in_tree)
# cat("Dropped species:\n")
# print(dropped$unique_name)
# 
# #build tree with only valid IDs
# taxa_clean <- taxa |> filter(in_tree)
# tree <- rotl::tol_induced_subtree(taxa_clean$ott_id)
# 
# plot.phylo(tree)
# 
# tree$tip.label <-
#   strip_ott_ids(tree$tip.label, remove_underscores = TRUE) # change ids to the names from dataset
# 
# #calculate correlations between all species (cor)
# tree2 <- compute.brlen(tree)
# cor <- vcv(tree2, cor = T)
# 
# # sort(rownames(cor))
# # sort(unique(dat$phylo))
# 
# setdiff(unique(dat$phylo), rownames(cor))  # these are the problem species
# setdiff(rownames(cor), unique(dat$phylo))  # in tree but not in data (less critical)
# 
# 
# #filter to only species present in the phylogenetic cor matrix
# p_subset <- p_subset |>
#   filter(species %in% rownames(cor))
# 
# #add phylo random factor
# p_subset$phylo <- p_subset$species


# create map for studies --------------------------------------------------

#count observations per location
study_location_summary <- p_raw %>%
  group_by(lat, lon, paper.index, Experiment.Type) %>%
  summarise(n_obs = n(), .groups = "drop")

#get world map
world_map <- map_data("world")

# world_map <- map_data("world") %>% 
#   filter(! long > 180) #remove countries with longitude >180 to make equal projection-like map without artifacts

#create map
ggplot() +
  geom_polygon(data = world_map,
               aes(x = long, y = lat, group = group),
               fill = "grey", color = "white", linewidth = 0.1) +
  geom_point(data = study_location_summary,
             aes(x = lon, y = lat, size = n_obs,
             col = Experiment.Type), alpha = 0.5) +
  scale_size_continuous(name = "number of observations", breaks = seq(0, 100, by = 25)) +
  scale_color_discrete(name = "experiment type") +
  theme_void() +
  coord_fixed(1.3)

ggsave("study_map.png", dpi = 300)



