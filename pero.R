#install phyloseq, ANCOMBC
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("phyloseq")
# BiocManager::install("ANCOMBC")

library(vegan)
library(phyloseq)
library(here)
library(janitor)
library(lubridate)
library(ggtext) #to add stress value
# library(RColorBrewer)
# library(colorBlindness)
library(tidyverse)

#clear environment
rm(list = ls())

#######################################################################
####### Load and Prep the Pero data - from 2019 field captures ########
#######################################################################

# #load data
# pero <- read.csv(here("MN_pero_raw_labeled_v10.8.19.csv"), na.strings=c(""," ","NA"))
# #pull ecto parasite columns if desired
# pero_ecto <- pero %>% dplyr::select(c(ear_tag_number, feces_number, fleas, botflies, botscars, wounds,
#                                   ticks_head, ticks_neck, tick_taken_number))

#### THE FOLLOWING CODE has been commented out and instead just pull the .rds file - MARCH 9 2023 ####

# # load and clean data
# pero <- read.csv(here("MN_pero_raw_labeled_v10.8.19.csv"), na.strings=c(""," ","NA")) %>%
#   clean_names() %>%
#   filter(species=="pspp") %>% #filter out non-Pero captures
#   select(!species) %>%
#   filter(capture!="WR") %>% #remove WR entries (duplicate fecal samples eg "026-2" were combined prior to DNA extraction)
#   drop_na(feces_number) %>% #drop entries with no feces sample
#   dplyr::select(!c(id, date_data_entry, data_enterer,
#                    field_site_detail, trap_number, animals_per_trap,
#                    weather, processing_time, handler, data_recorder,
#                    dna_number, blood_number, eye_bled, dead, dead_in_trap, rigor_mortis,
#                    vd_eye_number, notes)) %>% #drop unneeded columns
#   dplyr::select(!c(fleas, botflies, botscars, wounds,
#                    ticks_head, ticks_neck, tick_taken_number)) %>% #drop ecto parasite data
#   mutate(processing_date = as_date(processing_date, format="%m/%d/%Y"),
#          month = lubridate::month(processing_date, label=TRUE, abbr = FALSE),
#          location = as.factor(location),
#          landscape_type = as.factor(case_when(location=="CCESR" ~ "Agricultural",
#                                               location=="IBSL" ~ "Undeveloped")),
#          field_site = as.character(field_site),
#          site_type = as.factor(ifelse(grepl("perid", field_site), "Synanthropic",
#                             ifelse(grepl("maintenance", field_site), "Synanthropic", "Forest"))),
#          capture = as.factor(capture),
#          FDNA = as.character(paste0("F", feces_number)),
#          feces_number = as.numeric(feces_number),
#          pelage = as.factor(pelage),
#          sex = as.factor(sex),
#          testes = as.factor(testes),
#          perforate = as.factor(perforate),
#          nipples = as.factor(nipples),
#          lactating = as.factor(lactating),
#          pregnant = as.factor(pregnant)) %>% #reformat columns
#   mutate(reproductive = case_when(testes == "1" ~ "1",
#                                   perforate == "1" | nipples == "1" | lactating == "1" | pregnant == "1" ~ "1",
#                                   testes == "0" ~ "0",
#                                   perforate == "0" & nipples == "0" & lactating == "0" & pregnant == "0" ~ "0")) %>%
#   rows_patch(tibble(FDNA=c("F070", "F124", "F136"), body_mass=c(24,21,17))) %>% #patch in missing body_mass values (based on field estimates from NOTES column)
#   dplyr::select(!c(processing_date, field_site)) %>% #remove unneeded columns
#   unite(loc_site, landscape_type, site_type, remove=FALSE) %>%
#   mutate(loc_site = as.factor(loc_site)) %>%
#   mutate(ear_tag_number = ifelse(is.na(ear_tag_number), feces_number, ear_tag_number)) %>% #give all the IBSL mice an 'ear tag'
#   select(location, landscape_type, loc_site, site_type, month,
#          ear_tag_number, capture, FDNA, feces_number, body_mass, body_length,
#          sex, pelage, reproductive) %>% #reorder columns (and drop unneeded)
#   arrange(FDNA) #arrange numerically by FDNA#
# 
# 
# #check to make sure nothing is weird
# # unique(pero$processing_date)
# # unique(pero$location)
# # unique(pero$field_site)
# # unique(pero$pelage)
# # unique(pero$sex)
# # unique(pero$testes)
# # unique(pero$perforate)
# # unique(pero$nipples)
# # unique(pero$lactating)
# # unique(pero$pregnant)
# 
# #read in Nanopore barcode IDs
# nanopore <- read.csv(here("sampleID_to_barcode.csv")) %>%
#   separate(sample, c(NA, "feces_number")) %>% #separate fecal sample ID (F-###) into two columns, just keep ###
#   mutate(feces_number = as.numeric(feces_number)) %>%
#   arrange(feces_number)
# 
# # #confirm that every row in 'nanopore' matches a row in 'pero'
# # check <- anti_join(nanopore, pero, by="feces_number") #if this is 0, you're good
# 
# #join pero and nanopore
# pero <- pero %>% left_join(nanopore, by="feces_number") %>%
#   select(!feces_number) #drop unneeded column
# 
# #separate out samples not sequenced
# pero_notsequenced <- pero %>% filter(is.na(date_barcode))
# 
# #filter pero down to only samples sequenced
# pero <- pero %>% drop_na(date_barcode) %>%
#   arrange(FDNA)
# 
# # #Restore CC and IBSL meteo data from the rdata file
# # meteo_summary <- readRDS(file = "meteo_summary.rds") %>%
# #   mutate(location = as.factor(location))
# # #join to pero data
# # pero <- pero %>% left_join(meteo_summary, by=c("location", "month")) %>%
# #   arrange(FDNA) #order one more time, just to be safe
# 
# 
# ######## WRITE PERO160 data to rds ############
# # Save to a rdata file
# saveRDS(pero, file = here("pero_data_04.19.23.rds"))
# #03.22.23 version has only location words, peridomestic, meteo data
# #04.19.23 version has landscape words, synanthropic, and no meteo data
# 
# #Restore Nanopore-sampled pero160 data from the rdata file
# pero_data <- readRDS(file = "pero_data_04.19.23.rds")
# 
# #Create a short version of 160metadata
# # #if meteo data is included
# # pero_datashort <- pero_data %>%
# #   select(location, landscape_type, site_type, loc_site, month, 
# #          FDNA, pelage, sex, avg_hitemp, avg_lotemp, tot_precip, prev_precip)
# pero_datashort <- pero_data %>%
#   select(location, landscape_type, site_type, loc_site, month,
#          FDNA, ear_tag_number, sex, body_mass, pelage, reproductive)
# 
# #### END 160 samples - BEGIN Pero140 data (no recaptures)
# 
# #subset pero_data to see only animals with multiple entries
# recapped <- pero_data %>%
#   filter(location=="CCESR") %>%
#   group_by(ear_tag_number) %>%
#   arrange(ear_tag_number) %>%
#   mutate(ncaps = length(ear_tag_number)) %>%
#   relocate(ncaps, .after=ear_tag_number) %>%
#   filter(ncaps>1) %>%
#   ungroup()
# #every recapped animal was caught in July, so we'll keep those entries and remove others
# remove <- recapped %>% filter(month != "July")
# removeFDNAs.v <- remove$FDNA
# 
# #REMOVE DUPLICATE CAPTURES OF ANIMALS
# pero140_data <- pero_data %>%
#   filter(!FDNA %in% removeFDNAs.v)
# 
# ### WRITE Pero140 to RDS
# saveRDS(pero140_data, file="pero140_data_01.05.24.rds") #landscape words, synanthropic, and no meteo data

#Restore Nanopore-sampled pero140 data from the rdata file
pero140_data <- readRDS(file = "pero140_data_01.05.24.rds")

#Create a short version of metadata
# #with meteo data
# pero140_datashort <- pero140_data %>%
#   select(location, landscape_type, site_type, loc_site, month, 
#          FDNA, pelage, sex, avg_hitemp, avg_lotemp, tot_precip, prev_precip)
pero140_datashort <- pero140_data %>%
  select(location, landscape_type, site_type, loc_site, month,
         FDNA, ear_tag_number, sex, body_mass, pelage, reproductive)

#Just FDNA and site
samp_site <- pero140_data %>% select(FDNA, loc_site)
tag_loc <- pero140_datashort %>% select(location, landscape_type, site_type, month, FDNA)

###########################################################################
####### Load and Prep the Emu abundance data - from UMN MSI output ########
##### Create (count) abundance table per sample to the SPECIES LEVEL ######
###########################################################################

### NON-THRESHOLDED ABUNDANCE DATA - this is the 'raw' output from EMU - relabundance can be as low as it wants
### AS OF FEB 27 2023 - THESE ARE THE EMU DATA I'M USING

############# NON thresholded ABUNDANCE DATA #######################
################ TAXONOMIC INFO REMOVED ############################
# #read in the emu output from each site, remove the taxonomic info columns
# CC_peri_june_abund <- read.table(file = here("C_p_june_unfilt_emu-combined-tax_id-counts.tsv"), sep = '\t', header = TRUE)  %>%
#   select(!c(species, genus, family, order, class, phylum, superkingdom))
# CC_peri_july_abund <- read.table(file = here("C_p_july_unfilt_emu-combined-tax_id-counts.tsv"), sep = '\t', header = TRUE) %>%
#   select(!c(species, genus, family, order, class, phylum, superkingdom))
# CC_peri_aug_abund <- read.table(file = here("C_p_aug_unfilt_emu-combined-tax_id-counts.tsv"), sep = '\t', header = TRUE) %>%
#   select(!c(species, genus, family, order, class, phylum, superkingdom))
# 
# CC_forest_june_abund <- read.table(file = here("C_f_june_unfilt_emu-combined-tax_id-counts.tsv"), sep = '\t', header = TRUE) %>%
#   select(!c(species, genus, family, order, class, phylum, superkingdom))
# CC_forest_july_abund <- read.table(file = here("C_f_july_unfilt_emu-combined-tax_id-counts.tsv"), sep = '\t', header = TRUE) %>%
#   select(!c(species, genus, family, order, class, phylum, superkingdom))
# CC_forest_aug_abund <- read.table(file = here("C_f_aug_unfilt_emu-combined-tax_id-counts.tsv"), sep = '\t', header = TRUE) %>%
#   select(!c(species, genus, family, order, class, phylum, superkingdom))
# 
# IBSL_peri_abund <- read.table(file = here("I_p_unfilt_emu-combined-tax_id-counts.tsv"), sep = '\t', header = TRUE) %>%
#   select(!c(species, genus, family, order, class, phylum, superkingdom))
# IBSL_forest_abund <- read.table(file = here("I_f_unfilt_emu-combined-tax_id-counts.tsv"), sep = '\t', header = TRUE) %>%
#   select(!c(species, genus, family, order, class, phylum, superkingdom))
# 
# #stick all these suckers together
# abundance <- full_join(CC_peri_june_abund, CC_peri_july_abund, by="tax_id")
# abundance <- full_join(abundance, CC_peri_aug_abund, by="tax_id")
# abundance <- full_join(abundance, CC_forest_june_abund, by="tax_id")
# abundance <- full_join(abundance, CC_forest_july_abund, by="tax_id")
# abundance <- full_join(abundance, CC_forest_aug_abund, by="tax_id")
# abundance <- full_join(abundance, IBSL_forest_abund, by="tax_id")
# abundance <- full_join(abundance, IBSL_peri_abund, by="tax_id")
# 
# #input 0 for count=NA
# abundance <- abundance %>% replace(is.na(.), 0) %>%
#   dplyr::rename(taxID = tax_id) %>%
#   filter(!taxID=="unassigned") #remove 'unassigned' ID
# 
# ########### this little bit of code will cut down the sample ID columns ###############
# ###### https://community.rstudio.com/t/shorten-column-names-in-data-frame/143129 ######
# 
# long_names <- colnames(abundance)
# 
# short_names <- lapply(
#   X = strsplit(x = long_names, split = "_"),
#   FUN = function(x) paste0(head(x, n = 1)))
# 
# colnames(abundance) <- short_names
# ########################################################################################
# 
# # Regarding "count" values in abundance table: PLEASE NOTE!! (from Emu Github documentation):
# # "Note: Estimated read counts are based on likelihood probabilities and therefore may not be integer values."
# 
# #change counts to integers, remove duplicate FDNAs, arrange FDNA numerically
# abundance <- abundance %>%
#   pivot_longer(-taxID, names_to= "FDNA", values_to = "count") %>% #pivot longer, ignore the taxID column and pivot FDNA and counts
#   mutate(count = round(count, digits=0)) %>% #round "counts" to nearest whole number
#   filter(!FDNA %in% removeFDNAs.v) %>%
#   arrange(FDNA) %>% #order by FDNA#
#   pivot_wider(names_from = FDNA, values_from = count)
# 
# #### ABUNDANCE is a df where: SAMPLES are COLUMNS, TAXA are ROWS, TAXID is a COLUMN
# #Save abundance table to a rdata file
# saveRDS(abundance, file = here("abundance_140_03.22.23.rds"))

#Restore abundance from the rdata file
#rows are taxa, columns are samples
abundance <- readRDS(file = "abundance_140_03.22.23.rds") 

# ### THIS IS THE 160 - includes duplicates - abundance data
# abundance <- readRDS(file = "abundance_03.01.23.rds")

################ CREATE WORKING VERSIONS OF ABUNDANCE ######################

#phyloseq needs a not-trimmed (no thresholding) version of abundance table
#rows are taxa, columns are samples - in numerical order, rowname is taxID
abundance_phylo <- abundance %>%
  column_to_rownames(var="taxID")

#to go into vegan, abundance data needs to be a df where 
#taxa are the columns, samples are the rows, sampleIDs are row names
abundance_vegan <- abundance %>% 
  column_to_rownames(var="taxID") %>%
  t() %>% #transpose the matrix
  as.data.frame()

#make a long version of abundance for ggplot and others
abundance_long <- abundance %>%
  pivot_longer(-taxID, names_to = "FDNA", values_to = "count")


#########################################################################
# #OTUs per sample
# abundance %>% pivot_longer(-taxID, names_to = "FDNA", values_to="count") %>%
#   filter(count != 0) %>%
#   group_by(FDNA) %>%
#   summarise(n_OTU = length(count)) %>%
#   summarise(sd = sd(n_OTU)) #summarise(mean = mean(n_OTU)) 
#   #summarise(max = max(n_OTU)) #summarise(min = min(n_OTU))
#########################################################################


######### RAREFY ABUNDANCE (count) DATA ############

#how many should we rarefy to? 74517 is minimum number of sequences in a sample
#THANKS RIFFOMONAS PROJECT!
min_n_seq <- abundance_long %>% group_by(FDNA) %>%
  summarise(n_seq = sum(count)) %>%
  summarise(min = min(n_seq)) %>%
  pull(min)
 
# #mean and std dev of sequence length (also in the nanoq summary plots)
# seq_length_summary <- abundance_long %>% group_by(FDNA) %>%
#   summarise(n_seq = sum(count)) %>%
#   inner_join(pero_datashort, by="FDNA") %>%
#   group_by(loc_site) %>%
#   summarise(mean = mean(n_seq),
#             sd = sd(n_seq))

# #to see number of sequences per sample
# abundance %>%
#   rownames_to_column(var="taxID") %>%
#   pivot_longer(-taxID, names_to="FDNA", values_to = "count") %>%
#   group_by(FDNA) %>%
#   summarise(N = sum(count)) %>%
#   arrange(N) %>% print(n=20)

#### VISUALIZE RAREFACTION CURVES ###########

# #visualize the rarefaction curve in base R
# rarecurve(abundance_vegan, step=1000)
# abline(v=74517, col="blue")

# #pull the data from rarecurve() **THANKS RIFFOMONAS PROJECT!
# rarecurve_data <- rarecurve(abundance_vegan, step=1000)
# #clean it up, plot rarecurve() data in ggplot **THANKS RIFFOMONAS PROJECT!
# png(here("rarecurves140.png"), width=1000, height=600)
# map_dfr(rarecurve_data, bind_rows) %>%
#   bind_cols(FDNA = rownames(abundance_vegan),.) %>%
#   pivot_longer(-FDNA, names_to = "n_seq", values_to = "species_obs") %>%
#   left_join(samp_site, by="FDNA") %>%
#   drop_na() %>%
#   mutate(n_seq = as.numeric(str_replace(n_seq, "N", ""))) %>% #change the N### to just # of sequences
#   ggplot(aes(x=n_seq, y=species_obs, group=FDNA)) +
#   geom_line() +
#   geom_vline(xintercept = min_n_seq, color="red") +
#   labs(x="Number of Sequences", y="Number of Bacterial Species Observed") +
#   theme(axis.text = element_text(size=15),
#         axis.title = element_text(size=18))
# dev.off()

####### CREATE RAREFIED ABUNDANCE DATA (counts of all obs species, rarefied to min_n_seq ) #############
################### CALCULATE ALPHA DIVERSITY METRICS ON RAREFIED DATA #################################

####### NOTE: the EMU method DOES NOT give singleton or doubleton counts... 
###### even with the non-thresholded EMU output, lowest read count for a sample is 10 reads
## there will be warnings about using thresholded data because of this - can't do anything about it
  ## that's as few counts as the EMU algorithm is going to give ¯\_(ツ)_/¯
## as a result, there are some alpha diversity indices (related to richness) that cannot be estimated
  ## ACE groups reads into rare (less than 10 reads) and not rare (more than 10 reads)
  ## Chao looks at singletons
  ### so these richness estimators will not work with EMU data

#following Riffomonas project: "How to rarefy community data in R with vegan and tidyverse"
#https://youtu.be/_OEdFjc1D9I?si=e7dRDXoP1MAHf9ry starts 11:22-12:34
#using vegan::rrarefy() (instead of phyloseq::rarefy_even_depth) 
  #to resample the count data at min_n_seq, calculate the alpha diversity measures, 
  #save, and repeat many times
  #then average the alpha measures per sample from the many iterations (this is now 'rarefaction')

###### Estimate alpha diversity metrics using rarefaction #######
##### using vegan:rrarefy() to resample the count numbers #######

#commented out on Dec 18 2023
# #do a lil set.seed so I can reproduce these results
# set.seed(2111994)
# 
# #initialize a list to store results
# rare_alpha.ls <- list()
# 
# #for 100 iterations, rarefy count data, calculate alpha diversity, save, and repeat
# for(i in 1:100){
#   
#   #rarefy data (one iteration)
#   #output in vegan format: taxa are the columns, samples are the rows, sampleIDs are row names
#   rareit_vegan <- rrarefy(abundance_vegan, sample=min_n_seq) #min_n_seq = 74,517
#   
#   #convert rareit_vegan to phyloseq abundance table
#   #rows are taxa, columns are samples - in numerical order, rowname is taxID
#   rareit_phylo <- rareit_vegan %>% t()
#   #convert phyloseq to otu table
#   rareit_phylo_otu <- otu_table(rareit_phylo, taxa_are_rows=TRUE)
#   
#   #use phyloseq::estimate_richness() to get Richness, Shannon, Simpson diversity
#   #calculate Shannon evenness manually
#   #add column for iteration number, just to keep track
#   out_it <- estimate_richness(rareit_phylo_otu, measures = c("Observed", "Shannon", "Simpson")) %>%
#     mutate(Evenness = Shannon/log(Observed)) %>% #calculate evenness (Shannon/ln(Richness))
#     # mutate(iteration = i) %>%
#     rownames_to_column(var="FDNA")
#   
#   #save each iteration as an object in a list
#   rare_alpha.ls[[i]] <- out_it
#   
# }
# 
# #collapse the rareit_output list (one list item per iteration) down to a df
# #counts of all taxa for all samples, across 100 iterations
# #really big, n_iterations * 140 samples number of rows
# rare_alpha.df <- do.call(rbind, rare_alpha.ls) 
# #save it
# saveRDS(rare_alpha.df, file="rare_alpha.df_v12.18.2023.RDS")
#load it
rare_alpha.df <- readRDS(file = "rare_alpha.df_v12.18.2023.RDS") 

#group by sample, average diversity indices for a given FDNA sample across all 100 iterations
#join summary data to pero_datashort (metadata)
rare_alpha.avg <- rare_alpha.df %>% group_by(FDNA) %>% 
  summarise(rare_obs = mean(Observed),
            rare_shan = mean(Shannon),
            rare_simp = mean(Simpson),
            rare_even = mean(Evenness)) %>%
  inner_join(pero140_datashort, by="FDNA") %>%
  mutate(month = as.character(month)) %>%
  mutate(month.n = case_when(month == "June" ~ 1,
                             month == "July" ~ 2,
                             month == "August" ~ 3)) %>%
  mutate(month.n = factor(month.n, levels=c("2", "1", "3")))

#summary table of alpha diversity metrics
rare_alpha_month <- rare_alpha.avg %>% 
  group_by(landscape_type, site_type, month) %>%
  summarise(n = length(FDNA),
            mean_obs = mean(rare_obs),
            sd_obs = sd(rare_obs),
            mean_shan = mean(rare_shan),
            sd_shan = sd(rare_shan),
            mean_simp = mean(rare_simp),
            sd_simp = sd(rare_simp),
            mean_even = mean(rare_even),
            sd_even = sd(rare_even)) %>%
  ungroup() %>%
  mutate_if(is.numeric, round, digits=2) 

#alpha metrics for Cedar Creek across months
rare_alpha_CC <- rare_alpha.avg %>% 
  filter(landscape_type=="Agricultural") %>%
  group_by(landscape_type, site_type) %>%
  summarise(n = length(FDNA),
            mean_obs = mean(rare_obs),
            sd_obs = sd(rare_obs),
            mean_shan = mean(rare_shan),
            sd_shan = sd(rare_shan),
            mean_simp = mean(rare_simp),
            sd_simp = sd(rare_simp),
            mean_even = mean(rare_even),
            sd_even = sd(rare_even)) %>%
  ungroup() %>%
  mutate_if(is.numeric, round, digits=2) %>%
  mutate(month = "Summer") %>%
  relocate(month, .after=site_type)

rare_alpha_summary <- rbind(rare_alpha_month, rare_alpha_CC) %>%
  arrange(landscape_type, site_type, month)

#alpha diversity metric summary table - MANUSCRIPT FIGURE S3
alphatable <- rare_alpha_summary %>%
  mutate(obs = paste(mean_obs, "\u00B1", sd_obs, sep=" "),
         shan = paste(mean_shan, "\u00B1", sd_shan, sep=" "),
         simp = paste(mean_simp, "\u00B1", sd_simp, sep=" "),
         even = paste(mean_even, "\u00B1", sd_even, sep=" ")) %>%
  select(!c(mean_obs, sd_obs, mean_shan, sd_shan, mean_simp, sd_simp, mean_even, sd_even)) %>%
  rename(Landscape = landscape_type) %>% 
  rename(Habitat = site_type) %>%
  rename(Month = month) %>% 
  rename(N = n) %>%
  rename('Observed Richness \u2020' = obs) %>%
  rename("Shannon Diversity \u2021" = shan) %>%
  rename("Simpson Diversity \u00a7" = simp) %>% 
  rename("Evenness \u00b6" = even)

#save as .csv
write.csv(alphatable, here("Table_S3_v2.csv"))

# library(kableExtra)
# #pretty summary table - .png VERSION
# rare_alpha_summary %>%
#   mutate(obs = paste(mean_obs, "\u00B1", sd_obs, sep=" "),
#          shan = paste(mean_shan, "\u00B1", sd_shan, sep=" "),
#          simp = paste(mean_simp, "\u00B1", sd_simp, sep=" "),
#          even = paste(mean_even, "\u00B1", sd_even, sep=" ")) %>%
#   select(!c(mean_obs, sd_obs, mean_shan, sd_shan, mean_simp, sd_simp, mean_even, sd_even)) %>%
#   kbl(col.names = c("Landscape", "Habitat", "Month", "N",
#                     paste("Observed Richness", footnote_marker_alphabet(1)), 
#                     paste("Shannon Diversity", footnote_marker_alphabet(2)),
#                     paste("Simpson Diversity", footnote_marker_alphabet(3)), 
#                     paste("Evenness", footnote_marker_alphabet(4))),
#                     escape = FALSE, 
#       align="llllcccc") %>%
#   kable_styling(full_width=FALSE) %>%
#   row_spec(0, bold=T, color="black", background="#DAD7D7") %>%
#   row_spec(1:3, italic=TRUE, background="#F8F6F6") %>%
#   row_spec(4, bold=TRUE) %>%
#   row_spec(5:7, italic=TRUE, background="#F8F6F6") %>%
#   row_spec(8:10, bold=TRUE) %>%
#   footnote(alphabet=c("Estimated number of observed species as identified by the expectation-maximization algorithm (Emu)", 
#                       "Shannon diversity index gives equal weight to common and rare species", 
#                       "Simpson diversity index gives higher weight to common species when calculating diversity", 
#                       "Evenness measures the relative abundance of different taxa (Evenness = Shannon / Richness)"))

# rare_alpha_summary %>%
#   mutate(obs = paste(mean_obs, "\u00B1", sd_obs, sep=" "),
#          shan = paste(mean_shan, "\u00B1", sd_shan, sep=" "),
#          simp = paste(mean_simp, "\u00B1", sd_simp, sep=" "),
#          even = paste(mean_even, "\u00B1", sd_even, sep=" ")) %>%
#   select(!c(mean_obs, sd_obs, mean_shan, sd_shan, mean_simp, sd_simp, mean_even, sd_even)) %>%
#   kbl(col.names = c("Landscape", "Habitat", "Month", "N",
#                     paste("Observed Richness \u2020"), 
#                     paste("Shannon Diversity \u2021"),
#                     paste("Simpson Diversity \u00a7"), 
#                     paste("Evenness \u00b6")),
#       escape = FALSE, 
#       align="llllcccc") %>%
#   kable_styling(full_width=FALSE) %>%
#   row_spec(0, bold=T, color="black", background="#DAD7D7") %>%
#   row_spec(1:3, italic=TRUE, background="#F8F6F6") %>%
#   row_spec(4, bold=TRUE) %>%
#   row_spec(5:7, italic=TRUE, background="#F8F6F6") %>%
#   row_spec(8:10, bold=TRUE)


########## BUILD A MODEL #################

#cedar creek by month

# cc <- rare_alpha.avg %>% filter(location=="CCESR") %>%
#   mutate(month=as.character(month)) %>%
#   mutate(month = factor(month, levels=c("June", "July", "August")))
#   
# mod <- lm(Observed ~ site_type + month, data=cc)
# summary(mod)
# 
# mod <- lm(Shannon ~ site_type + month, data=cc)
# summary(mod)
# 
# mod <- lm(Simpson ~ site_type + month, data=cc)
# summary(mod)
# 
# mod <- lm(Evenness ~ site_type + month, data=cc)
# summary(mod)


########################################

#### before modeling: ######

# #visualize response distribution
# hist(rare_alpha.avg$Observed) #normalish
# hist(rare_alpha.avg$Shannon) #slight left skew
# hist(rare_alpha.avg$Simpson) #left skew
# hist(rare_alpha.avg$Evenness) #slight left skew

### MODELING ###

###-----all of these models have been run with Dec 18 2023 version of alpha_rare

#### RICHNESS MODEL - LINEAR REGRESSION ####

# 1. plot the full model
rich.mod <- lm(rare_obs ~ landscape_type + site_type + landscape_type:site_type + 
                 sex + reproductive + body_mass + month.n, 
               data=rare_alpha.avg)
# #model diagnostics
# # make a null model and compare full to null by AIC (difference of 2-4 is meaningful, lower AIC wins)
# null.mod <- lm(rare_obs ~ 1, data=rare_alpha.avg)
# AIC(null.mod, rich.mod)
# # 2. check out the summary
rich.summ <- summary(rich.mod)
# anova(rich.mod)
#to get 95% CI for model params
richCI <- as.data.frame(confint(rich.mod, level=0.95)) %>% #gives CI for all params, can also specify which you want
  rownames_to_column(var = "param")
rich <- as.data.frame(rich.summ$coefficients) %>% rownames_to_column(var="param") %>% 
  left_join(richCI, by="param") %>%
  mutate(mod = "Richness") %>% mutate_if(is.numeric, round, digits=3)

write.csv(rich, here("mod.rich.csv")) #save this to make a table for the supplement

# # 3. look at the residuals vs fitted plots (want nice clouds, look weird patterns!)
# #### USE YOUR EYEBALLS FIRST! If the residual plots are fucked up, then go for the formal diagnostics
# #base R plots (hit return to see them all)
# plot(rich.mod) #diagnostic plots
# #or use 'performance' package
# library(performance)
# check_model(rich.mod) #diagnostic plots
# 
# # 4. IF THINGS ARE UGLY - then run some formal diagnostics:
# #LINEARITY
# 
# #INDEPENDENCE (of data and of residuals) 
# 
# #MULTICOLLINEARITY
# 
# library(car)
# vif(rich.mod) #variance inflation factor (if parameters are correlated)
# #under 2 is great, under 4 is good enough (Zuur et al 2009)
# 
# #NORMALITY (of residuals)
# 
# #plot a histogram of your response variable and numeric independent vars (should be normal or normalish)
# #Q-Q plot also addresses this
# #use Shapiro-Wilkes as well if things seem wonky
# 
# #EQUIVARIANCE (homoscedasticity)
# 
# library(lmtest)
# bptest(rich.mod) #for regressions for homoscedasticity (use if your data is goofy, ignore elsewise)
# #Levine test (use in place of bp) for glm or mixed-models
# 
# ####### JASMINE has good explanations with more detail in her code for Lucie and Mathilde


#### SHANNON MODEL - LINEAR REGRESSION ####

shan.mod <- lm(rare_shan ~ landscape_type + site_type + landscape_type:site_type + 
                 sex + reproductive + body_mass + month.n, 
               data=rare_alpha.avg)
summary(shan.mod)
shan.summ <- summary(shan.mod)
#to get 95% CI for model params
shanCI <- as.data.frame(confint(shan.mod, level=0.95)) %>% #gives CI for all params, can also specify which you want
  rownames_to_column(var = "param")
shan <- as.data.frame(shan.summ$coefficients) %>% rownames_to_column(var="param") %>%
  left_join(shanCI, by="param") %>%
  mutate(mod = "Shannon") %>% mutate_if(is.numeric, round, digits=3)
write.csv(shan, here("mod.shan.csv"))

# #model diagnostics
# plot(shan.mod)
# # vif(shan.mod) #variance inflation factor (if parameters are correlated) <4 is okay, <2 is great
# bptest(shan.mod) #for homoskedasticity, this is 0.051 - not a amazing but maybe okay, plots looked okay
# check_model(shan.mod)
# 
# null.mod <- lm(Shannon ~ 1, data=rare_alpha.avg)
# AIC(null.mod, shan.mod)


#### Beta regression for proportional continuous data (between 0-1) like Simpson and Evennness
## this works well with slightly skewed (ie non-normal) data, lack of homoscedasticity

library(betareg) 

#### SIMPSON MODEL - BETA REGRESSION ####

simp.mod <- betareg(rare_simp ~ landscape_type + site_type + landscape_type:site_type + 
                      sex + reproductive + body_mass + month.n, 
                    link="logit",
                    data=rare_alpha.avg)
simp.summ <- summary(simp.mod)
#to get 95% CI for model params
simpCI <- as.data.frame(confint(simp.mod, level=0.95)) %>% #gives CI for all params, can also specify which you want
  rownames_to_column(var = "param")
simp <- as.data.frame(simp.summ$coefficients$mean) %>% rownames_to_column(var="param") %>% 
  left_join(simpCI, by="param") %>%
  mutate(mod = "Simpson") %>% mutate_if(is.numeric, round, digits=3)
write.csv(simp, here("mod.simp.csv"))

# #model diagnostics
# plot(simp.mod)
# # vif(simp.mod) #variance inflation factor (if parameters are correlated) <4 is okay, <2 is great
# check_model(simp.mod)
# 
# null.mod <- lm(Simpson ~ 1, data=rare_alpha.avg)
# AIC(null.mod, simp.mod)

#### EVENNESS MODEL - BETA REGRESSION ####

even.mod <- betareg(rare_even ~ landscape_type + site_type + landscape_type:site_type + 
                      sex + reproductive + body_mass + month.n, 
                    link="logit",
                    data=rare_alpha.avg)
summary(even.mod)
even.summ <- summary(even.mod)
#to get 95% CI for model params
evenCI <- as.data.frame(confint(even.mod, level=0.95)) %>% #gives CI for all params, can also specify which you want
  rownames_to_column(var = "param")
even <- as.data.frame(even.summ$coefficients$mean) %>% rownames_to_column(var="param") %>%
  left_join(evenCI, by="param") %>%
  mutate(mod = "Evenness") %>% mutate_if(is.numeric, round, digits=3)
write.csv(even, here("mod.even.csv"))

# #model diagnostics
# plot(even.mod)
# # vif(even.mod) #variance inflation factor (if parameters are correlated) <4 is okay, <2 is great
# bptest(even.mod) #for homoskedasticity
# check_model(even.mod)
# 
# null.mod <- lm(Evenness ~ 1, data=rare_alpha.avg)
# AIC(null.mod, even.mod)


######### visualize ###########

########### WORKING HERE - to show pvalues on plots
#### but can't do it with facets... so figure that out

# library(ggh4x) #this might need to be loaded without phyloseq running
# ggplot(aes(y=rare_obs, x=interaction(site_type, landscape_type), 
#                fill=loc_site, alpha=0.7)) +
#   scale_x_discrete(guide = guide_axis_nested(), name="location_type") 
#############################################################

#plot observed

#https://www.r-bloggers.com/2017/06/add-p-values-and-significance-levels-to-ggplots/
library(ggpubr)
# Perform pairwise comparisons
compare_means(rare_obs ~ loc_site,  data = rare_alpha.avg)
# Visualize: Specify the comparisons you want
my_comparisons <- list( c("Agricultural_Forest", "Undeveloped_Forest"), 
                        c("Agricultural_Synanthropic", "Undeveloped_Forest"), 
                        c("Undeveloped_Forest", "Undeveloped_Synanthropic") )

# obs.plot <- rare_alpha.avg %>%
#   ggplot(aes(y=rare_obs, x=site_type, fill=loc_site, alpha=0.7)) +
#   facet_wrap(~landscape_type) +
#   geom_violin() +
#   theme_light() +
#   scale_fill_manual(values=c("#A67326", "#F0C907",
#                              "#2c7c94", "#a6d0c8")) +
#   stat_summary(fun = "mean",
#                geom = "crossbar", 
#                width = 0.5,
#                colour = "black") +
#   geom_jitter(width=0.1, alpha=0.5) +
#   theme(legend.position = "none",
#         strip.text = element_text(size=18, color="black"),
#         strip.background=element_rect(fill="#d3d3d3"),
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.title.y = element_text(size=22),
#         axis.text.y = element_text(size=14),
#         plot.margin = margin(t = 10, r = 25, b = 0, l = 10, unit = "pt")) +
#   labs(y = "Observed Species Richness")

obs.plot <- rare_alpha.avg %>%
  ggplot(aes(y=rare_obs, x=loc_site, fill=loc_site, alpha=0.7)) +
  geom_violin() +
  geom_jitter(width=0.1, alpha=0.5) + #add the jittered points
  theme_half_open() +
  scale_fill_manual(values=c("#A67326", "#F0C907",
                             "#2c7c94", "#a6d0c8")) +
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = "black") +
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value (option: label="p.signif" )
  stat_compare_means(label.y = 440, size=5) +  # Add global p-value
  theme(legend.position = "none",       
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=22),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=8),
        plot.margin = margin(t = 10, r = 25, b = 0, l = 10, unit = "pt")) +
  labs(y = "Observed Species Richness")
#to edit the size of the pvalue text
#https://stackoverflow.com/questions/48550525/ggpubr-change-font-size-of-stat-compare-means-kruskal-wallis-p-values
library(gginnards)
obs.plot$layers[[which_layers(obs.plot, "GeomSignif")]]$aes_params$textsize <- 5


#plot Shannon

# Perform pairwise comparisons
compare_means(rare_shan ~ loc_site,  data = rare_alpha.avg)
# Visualize: Specify the comparisons you want
my_comparisons <- list( c("Agricultural_Forest", "Undeveloped_Forest"), 
                        c("Agricultural_Synanthropic", "Undeveloped_Forest"), 
                        c("Undeveloped_Forest", "Undeveloped_Synanthropic") )

# shan.plot <- rare_alpha.avg %>%
#   ggplot(aes(y=rare_shan, x=site_type, fill=loc_site, alpha=0.7)) +
#   facet_wrap(~landscape_type) +
#   geom_violin() +
#   theme_light() +
#     scale_fill_manual(values=c("#A67326", "#F0C907",
#                                "#2c7c94", "#a6d0c8")) +
#   stat_summary(fun = "mean",
#                geom = "crossbar", 
#                width = 0.5,
#                colour = "black") +
#   geom_jitter(width=0.1, alpha=0.5) +
#     theme(legend.position = "none",
#           strip.text = element_text(size=18, color="black"),
#           strip.background=element_rect(fill="#d3d3d3"),
#           axis.title.x = element_blank(),
#           axis.text.x = element_blank(),
#           axis.ticks.x = element_blank(),
#           axis.title.y = element_text(size=22),
#           axis.text.y = element_text(size=14),
#           plot.margin = margin(t = 10, r = 10, b = 0, l = 25, unit = "pt")) +
#   labs(y = "Shannon Diversity Index")

shan.plot <- rare_alpha.avg %>%
  ggplot(aes(y=rare_shan, x=loc_site, fill=loc_site, alpha=0.7)) +
  geom_violin() +
  geom_jitter(width=0.1, alpha=0.5) + #add the jittered points
  theme_half_open() +
  scale_fill_manual(values=c("#A67326", "#F0C907",
                             "#2c7c94", "#a6d0c8")) +
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = "black") +
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value (option: label="p.signif" )
  stat_compare_means(label.y = 5.25, size=5) +  # Add global p-value
  theme(legend.position = "none",       
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=22),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=8),
        plot.margin = margin(t = 10, r = 10, b = 0, l = 25, unit = "pt")) +
  labs(y = "Shannon Diversity Index")
#to edit the size of the pvalue text
shan.plot$layers[[which_layers(shan.plot, "GeomSignif")]]$aes_params$textsize <- 5


# c("#2c7c94", "#a6d0c8", 
#   "#c4aa23", "#fbe45b")
# 
# c("#0060A2", "#61CFFA",
#   "#A67326", "#F0C907")

#plot Simpson

# Perform pairwise comparisons
compare_means(rare_simp ~ loc_site,  data = rare_alpha.avg)
# Visualize: Specify the comparisons you want
my_comparisons <- list( c("Agricultural_Forest", "Undeveloped_Forest"), 
                        c("Agricultural_Synanthropic", "Undeveloped_Forest"), 
                        c("Undeveloped_Forest", "Undeveloped_Synanthropic") )

# simp.plot <- rare_alpha.avg %>%
#   ggplot(aes(y=rare_simp, x=site_type, fill=loc_site, alpha=0.7)) +
#   facet_wrap(~landscape_type) +
#   geom_violin() +
#   theme_light() +
#   scale_fill_manual(values=c("#A67326", "#F0C907",
#                              "#2c7c94", "#a6d0c8")) +
#   stat_summary(fun = "mean",
#                geom = "crossbar", 
#                width = 0.5,
#                colour = "black") +
#   geom_jitter(width=0.1, alpha=0.5) +
#   theme(legend.position = "none",
#         strip.text = element_text(size=18, color="black"),
#         strip.background=element_rect(fill="#d3d3d3"),
#         axis.title.x = element_blank(),
#         axis.title.y = element_text(size=22),
#         axis.text.y = element_text(size=14),
#         axis.text.x = element_text(size=18),
#         plot.margin = margin(t = 25, r = 25, b = 10, l = 10, unit = "pt")) +
#   labs(y = "Simpson Diversity Index")

simp.plot <- rare_alpha.avg %>%
  ggplot(aes(y=rare_simp, x=loc_site, fill=loc_site, alpha=0.7)) +
  geom_violin() +
  geom_jitter(width=0.1, alpha=0.5) + #add the jittered points
  theme_half_open() +
  scale_fill_manual(values=c("#A67326", "#F0C907",
                             "#2c7c94", "#a6d0c8")) +
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = "black") +
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value (option: label="p.signif" )
  stat_compare_means(label.y = 1.25, size=5) +  # Add global p-value
  theme(legend.position = "none",       
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=22),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=8),
        plot.margin = margin(t = 25, r = 25, b = 10, l = 10, unit = "pt")) +
  labs(y = "Simpson Diversity Index")
#to edit the size of the pvalue text
simp.plot$layers[[which_layers(simp.plot, "GeomSignif")]]$aes_params$textsize <- 5


#plot Evennness

# Perform pairwise comparisons
compare_means(rare_even ~ loc_site,  data = rare_alpha.avg)
# Visualize: Specify the comparisons you want
my_comparisons <- list( c("Agricultural_Forest", "Undeveloped_Forest"), 
                        c("Agricultural_Synanthropic", "Undeveloped_Forest"), 
                        c("Undeveloped_Forest", "Undeveloped_Synanthropic") )

# even.plot <- rare_alpha.avg %>%
#   ggplot(aes(y=rare_even, x=site_type, fill=loc_site, alpha=0.7)) +
#   facet_wrap(~landscape_type) +
#   geom_violin() +
#   theme_light() +
#   scale_fill_manual(values=c("#A67326", "#F0C907",
#                              "#2c7c94", "#a6d0c8")) +
#   stat_summary(fun = "mean",
#                geom = "crossbar", 
#                width = 0.5,
#                colour = "black") +
#   geom_jitter(width=0.1, alpha=0.5) +
#   theme(legend.position = "none",
#         strip.text = element_text(size=18, color="black"),
#         strip.background=element_rect(fill="#d3d3d3"),        
#         axis.title.x = element_blank(),
#         axis.title.y = element_text(size=22),
#         axis.text.y = element_text(size=14),
#         axis.text.x = element_text(size=18),
#         plot.margin = margin(t = 25, r = 10, b = 10, l = 25, unit = "pt")) +
#   labs(y = "Shannon Species Evenness")

even.plot <- rare_alpha.avg %>%
  ggplot(aes(y=rare_even, x=loc_site, fill=loc_site, alpha=0.7)) +
  geom_violin() +
  geom_jitter(width=0.1, alpha=0.5) + #add the jittered points
  theme_half_open() +
  scale_fill_manual(values=c("#A67326", "#F0C907",
                             "#2c7c94", "#a6d0c8")) +
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = "black") +
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value (option: label="p.signif" )
  stat_compare_means(label.y = 0.94, size=5) +  # Add global p-value
  theme(legend.position = "none",       
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=22),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=8),
        plot.margin = margin(t = 25, r = 10, b = 10, l = 25, unit = "pt")) +
  labs(y = "Shannon Evenness Index")
#to edit the size of the pvalue text
even.plot$layers[[which_layers(even.plot, "GeomSignif")]]$aes_params$textsize <- 5

######### COMPOSITE PLOT (publication-ready using cowplot) ##########

library(cowplot)

png(here("alphadiv_plotv2.png"), width=1100, height=800)
plot_grid(obs.plot, shan.plot, simp.plot, even.plot, 
          labels = "AUTO",label_size = 28)
dev.off()


# cvdPlot(plot)

####### (simpler) PLOT IT #############

# library(ggbeeswarm)
# 
# rare_alpha.avg %>% 
#   select(rare_obs, rare_shan, rare_simp, rare_Even, loc_site) %>%
#   pivot_longer(-loc_site, values_to = "est", names_to="stat") %>%
#   mutate(stat = factor(stat, levels=c("rare_obs", "rare_shan", "rare_simp", "rare_even"))) %>%
#   ggplot(aes(x=loc_site, y=est, fill=loc_site)) +
#   geom_boxplot() +
#   geom_beeswarm(alpha=0.3) +
#   facet_wrap(~stat, scales = "free_y", ncol=4)


################################### END ALPHA DIVERSITY #####################################








###########################################################################################
######################## BETA DIVERSITY - DISTANCE MATRIX #################################
###########################################################################################

### THIS WILL COMPUTE a distance matrix between all the samples based on the full microbiome community
  ### samples are rarefied to least numerous sample (74 517) before computing the distance matrix

# #avgdist needs as input: dataframe with columns as OTU and rows as samples (sampleID in row name)
# #confirm that FDNA numbers are already in numerical order
# abundance_dist <- abundance_vegan %>% #rows are mice, columns are OTUs
#   avgdist(sample=74517, iterations=100, dmethod="bray") #calculate the BrayCurtis dist - rarefy and repeat 100 iterations
# 
# # Save to a rdata file
# saveRDS(abundance_dist, file = here("abundance_dist_140_03.22.23.rds"))

# Restore from the rdata file
abundance_dist <- readRDS(file = "abundance_dist_140_03.22.23.rds")

#stress plot (scree plot) (check the k values! get STRESS down to < 0.1)
# library(goeveg) #for dimcheckMNDS() function
# dimcheckMDS(abundance_vegan, distance="bray")

################# NMDS ######################

#run NMDS on distance matrix (FDNA IN NUMERICAL ORDER!)
set.seed(19940211) #best solution not repeated after 20 tries >> IS THIS A PROBLEM??
nmds_results <- metaMDS(abundance_dist, k=4) 

#create tibble of NMDS values, add FDNA numbers, metadata
pero_nmds <- nmds_results %>%
  scores() %>% 
  as_tibble(rownames="FDNA") %>%
  full_join(pero140_datashort, by="FDNA") %>%
  ungroup()

#generate centroids for each loc_site
centroid <- pero_nmds %>% group_by(loc_site) %>%
  summarise(NMDS1=mean(NMDS1), NMDS2=mean(NMDS2))

# c("#a65852", "#c4aa23",
#   "#a6d0c8", "#fbe45b",
#   "#2c7c94", "#A0999A")
# 
# c("#0060A2", "#61CFFA",
#   "#A67326", "#F0C907")

#plot NMDS, centroids, and ellipses
png(here("community_NMDS_k4.png"), width=700, height = 600)
pero_nmds %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=loc_site)) +
  geom_point(size = 2) +
  scale_color_manual(name="Landscape-Habitat Type", 
                       labels=c("Agricultural-Forest", "Agricultural-Synanthropic", 
                                "Undeveloped-Forest", "Undeveloped-Synanthropic"),
                     values=c("#A67326", "#F0C907",
                              "#2c7c94", "#a6d0c8")) +
  stat_ellipse(show.legend = FALSE, alpha=0.75, linewidth=1) +
  geom_point(data=centroid, size=5, shape=21, color="black", alpha=0.75,
             aes(fill=loc_site),
             show.legend = FALSE) +
  scale_fill_manual(values=c("#A67326", "#F0C907",
                             "#2c7c94", "#a6d0c8")) +
  theme_light() +
  theme(legend.position = "bottom",
        legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        axis.text = element_text(size=12),
        axis.title = element_text(size=15),
        panel.background = element_rect(fill='white')) +
  guides(color=guide_legend(nrow=2)) +
  annotate(geom = "label", x=0.52, y=0.52, size = 7,
           label = paste("Stress:", round(nmds_results$stress, digits = 3)))
dev.off()

# #Colorblind friendliness
# cvdPlot(plot)
# 
# #other colors
# values=c("#A82433", "#EA99A3",
#          "#0060A2", "#61CFFA")


######### Reviewer asked to see NMDS in 2 vs 3 vs 4 dimensions ###########

###### four dimensions to fit in a 3-panel figure
nmds_k4 <- pero_nmds %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=loc_site)) +
  geom_point(size = 2) +
  scale_color_manual(name="Landscape-Habitat Type", 
                     labels=c("Agricultural-Forest", "Agricultural-Synanthropic", 
                              "Undeveloped-Forest", "Undeveloped-Synanthropic"),
                     values=c("#A67326", "#F0C907",
                              "#2c7c94", "#a6d0c8")) +
  stat_ellipse(show.legend = FALSE, alpha=0.75, linewidth=1) +
  geom_point(data=centroid, size=5, shape=21, color="black", alpha=0.75,
             aes(fill=loc_site),
             show.legend = FALSE) +
  xlim(-0.52, 0.61) + ylim(-0.58, 0.55) +
  scale_fill_manual(values=c("#A67326", "#F0C907",
                             "#2c7c94", "#a6d0c8")) +
  theme_light() +
  labs(title="NMDS k=4") +
  theme(legend.position = "bottom",
        plot.title = element_text(size=24),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        panel.background = element_rect(fill='white'),
        plot.margin = margin(t = 0, r = 25, b = 0, l = 25, unit = "pt")) +
  guides(color=guide_legend(nrow=2)) +
  annotate(geom = "label", x=0.42, y=0.52, size = 7,
           label = paste("Stress:", round(nmds_results$stress, digits = 3)))

#####two dimensions
#run NMDS on distance matrix (FDNA IN NUMERICAL ORDER!)
set.seed(19940211) #best solution not repeated after 20 tries >> IS THIS A PROBLEM??
nmds_results2 <- metaMDS(abundance_dist, k=2) 

#create tibble of NMDS values, add FDNA numbers, metadata
pero_nmds2 <- nmds_results2 %>%
  scores() %>% 
  as_tibble(rownames="FDNA") %>%
  full_join(pero140_datashort, by="FDNA") %>%
  ungroup()

#generate centroids for each loc_site
centroid2 <- pero_nmds2 %>% group_by(loc_site) %>%
  summarise(NMDS1=mean(NMDS1), NMDS2=mean(NMDS2))

#plot NMDS, centroids, and ellipses
# png(here("community_NMDS_k2.png"), width=700, height = 600)
nmds_k2 <- pero_nmds2 %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=loc_site)) +
  geom_point(size = 2) +
  scale_color_manual(name="Landscape-Habitat Type", 
                     labels=c("Agricultural-Forest", "Agricultural-Synanthropic", 
                              "Undeveloped-Forest", "Undeveloped-Synanthropic"),
                     values=c("#A67326", "#F0C907",
                              "#2c7c94", "#a6d0c8")) +
  stat_ellipse(show.legend = FALSE, alpha=0.75, linewidth=1) +
  geom_point(data=centroid, size=5, shape=21, color="black", alpha=0.75,
             aes(fill=loc_site),
             show.legend = FALSE) +
  xlim(-0.52, 0.61) + ylim(-0.58, 0.55) +
  scale_fill_manual(values=c("#A67326", "#F0C907",
                             "#2c7c94", "#a6d0c8")) +
  theme_light() +
  labs(title="NMDS k=2") +
  theme(legend.position = "bottom",
        plot.title = element_text(size=24),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        panel.background = element_rect(fill='white'),
        plot.margin = margin(t = 0, r = 25, b = 0, l = 25, unit = "pt")) +
  guides(color=guide_legend(nrow=2)) +
  annotate(geom = "label", x=0.42, y=0.52, size = 7,
           label = paste("Stress:", round(nmds_results2$stress, digits = 3)))
# dev.off()


#####three dimensions
#run NMDS on distance matrix (FDNA IN NUMERICAL ORDER!)
set.seed(19940211) #best solution not repeated after 20 tries >> IS THIS A PROBLEM??
nmds_results3 <- metaMDS(abundance_dist, k=3) 

#create tibble of NMDS values, add FDNA numbers, metadata
pero_nmds3 <- nmds_results3 %>%
  scores() %>% 
  as_tibble(rownames="FDNA") %>%
  full_join(pero140_datashort, by="FDNA") %>%
  ungroup()

#generate centroids for each loc_site
centroid3 <- pero_nmds3 %>% group_by(loc_site) %>%
  summarise(NMDS1=mean(NMDS1), NMDS2=mean(NMDS2))

#plot NMDS, centroids, and ellipses
# png(here("community_NMDS_k3.png"), width=700, height = 600)
nmds_k3 <- pero_nmds3 %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=loc_site)) +
  geom_point(size = 2) +
  scale_color_manual(name="Landscape-Habitat Type", 
                     labels=c("Agricultural-Forest", "Agricultural-Synanthropic", 
                              "Undeveloped-Forest", "Undeveloped-Synanthropic"),
                     values=c("#A67326", "#F0C907",
                              "#2c7c94", "#a6d0c8")) +
  stat_ellipse(show.legend = FALSE, alpha=0.75, linewidth=1) +
  geom_point(data=centroid, size=5, shape=21, color="black", alpha=0.75,
             aes(fill=loc_site),
             show.legend = FALSE) +
  xlim(-0.52, 0.61) + ylim(-0.58, 0.55) +
  scale_fill_manual(values=c("#A67326", "#F0C907",
                             "#2c7c94", "#a6d0c8")) +
  theme_light() +
  labs(title="NMDS k=3") +
  theme(legend.position = "bottom",
        plot.title = element_text(size=24),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        panel.background = element_rect(fill='white'),
        plot.margin = margin(t = 0, r = 25, b = 0, l = 25, unit = "pt")) +
  guides(color=guide_legend(nrow=2)) +
  annotate(geom = "label", x=0.42, y=0.52, size = 7,
           label = paste("Stress:", round(nmds_results3$stress, digits = 3)))
# dev.off()


#cowplot to see k=2 vs k=3 vs k=4
png(here("compare_NMDS.png"), width=1800, height = 600)
plot_grid(nmds_k2, nmds_k3, nmds_k4, ncol = 3, 
          labels = "AUTO", label_size = 36)
dev.off()

##################################################################
# # code from: https://www.rpubs.com/RGrieger/545184
# pero.spp.fit <- envfit(nmds_results, abundance_vegan, permutations=999)
# ordiplot(nmds_results, type="n")
# plot(pero.spp.fit, p.max=0.001, col="black")
# 
# 
# pero_datashort_matrix <- pero_datashort %>% column_to_rownames(var="FDNA")
# pero.envfit <- envfit(nmds_results, pero_datashort_matrix, permutations=999)
# head(pero.envfit)
# ordiplot(nmds_results, type="n")
# plot(pero.envfit, p.max=0.001, cex=0.7, col="black")
####################################################################

####### INTERESTING: JUNE communities at CCESR are different than July/Aug
# #plot CCESR data by month
# pero_nmds %>%
#   filter(location=="CCESR") %>%
#   ggplot(aes(x=NMDS1, y=NMDS2, color=month)) +
#   geom_point() +
#   stat_ellipse(show.legend = FALSE)



################# STATISTICS ################

# #using adonis for a single variable (compare the centroids)
# #compare CCESR and IBSL
# adonis2(as.dist(abundance_dist)~pero_nmds$location) #first thing needs to be a dist matrix
# 
# #compare forest vs peridomestic
# adonis2(as.dist(abundance_dist)~pero_nmds$site_type) #first thing needs to be a dist matrix

#compare location and site_type
adonisres <- adonis2(as.dist(abundance_dist) ~ pero_nmds$landscape_type + pero_nmds$site_type + 
                       pero_nmds$sex + pero_nmds$reproductive +
                       pero_nmds$body_mass + pero_nmds$month, 
        by="margin")
##NO DUPLICATES IN THE DATA

rownames(adonisres) <- c("Landscape", "Habitat", "Sex", "Reproductive Status",
                         "Body Mass", "Sampling Month",
                         "Residual", "Total")
adonisres <- adonisres %>% rownames_to_column(var="Variable")
adonisres <- adonisres %>% mutate_if(is.numeric, round, digits=4)
#pretty kable
# library(kableExtra)
adonisres %>% kbl() %>%
  kable_styling(full_width=FALSE) %>%
row_spec(0, bold=T, color="black", background="#DAD7D7")
  
write.csv(adonisres, here("Table_S4.csv"))


#test for variation in beta dispersion by location
bd <- betadisper(abundance_dist, pero_nmds$location)
# anova(bd) #no difference in dispersion
permutest(bd) #no difference in dispersion

#test for variation in beta dispersion by site_type
bd <- betadisper(abundance_dist, pero_nmds$site_type)
# anova(bd) #no difference in dispersion
permutest(bd) #no difference in dispersion

#test for variation in beta dispersion by loc_type
bd <- betadisper(abundance_dist, pero_nmds$loc_site)
# anova(bd) #no difference in dispersion
permutest(bd) #no difference in dispersion

## anosim - is between group variation > within group variation ?
anosim(abundance_dist, pero_nmds$loc_site)

##################### END STATISTICS ###########################









#################### CREATE (count) ABUNDANCE TABLE at the SPECIES LEVEL #########################
###################### Purpose: to compare to bacterial pathogen data ##########################

############# NON thresholded ABUNDANCE DATA #######################
################### WITH TAXONOMIC INFO ############################
# #read in the emu output from each site, DO NOT remove the taxonomic info columns
# CC_peri_june_abund <- read.table(file = here("C_p_june_unfilt_emu-combined-tax_id-counts.tsv"), sep = '\t', header = TRUE)
# CC_peri_july_abund <- read.table(file = here("C_p_july_unfilt_emu-combined-tax_id-counts.tsv"), sep = '\t', header = TRUE)
# CC_peri_aug_abund <- read.table(file = here("C_p_aug_unfilt_emu-combined-tax_id-counts.tsv"), sep = '\t', header = TRUE)
# 
# CC_forest_june_abund <- read.table(file = here("C_f_june_unfilt_emu-combined-tax_id-counts.tsv"), sep = '\t', header = TRUE)
# CC_forest_july_abund <- read.table(file = here("C_f_july_unfilt_emu-combined-tax_id-counts.tsv"), sep = '\t', header = TRUE)
# CC_forest_aug_abund <- read.table(file = here("C_f_aug_unfilt_emu-combined-tax_id-counts.tsv"), sep = '\t', header = TRUE)
# 
# IBSL_peri_abund <- read.table(file = here("I_p_unfilt_emu-combined-tax_id-counts.tsv"), sep = '\t', header = TRUE)
# IBSL_forest_abund <- read.table(file = here("I_f_unfilt_emu-combined-tax_id-counts.tsv"), sep = '\t', header = TRUE)
# 
# #stick all these suckers together
# abundance_tax <- full_join(CC_peri_june_abund, CC_peri_july_abund,
#                            by=c("tax_id", "species", "genus", "family", "order", "class", "phylum", "superkingdom"))
# abundance_tax <- full_join(abundance_tax, CC_peri_aug_abund,
#                            by=c("tax_id", "species", "genus", "family", "order", "class", "phylum", "superkingdom"))
# abundance_tax <- full_join(abundance_tax, CC_forest_june_abund,
#                            by=c("tax_id", "species", "genus", "family", "order", "class", "phylum", "superkingdom"))
# abundance_tax <- full_join(abundance_tax, CC_forest_july_abund,
#                            by=c("tax_id", "species", "genus", "family", "order", "class", "phylum", "superkingdom"))
# abundance_tax <- full_join(abundance_tax, CC_forest_aug_abund,
#                            by=c("tax_id", "species", "genus", "family", "order", "class", "phylum", "superkingdom"))
# abundance_tax <- full_join(abundance_tax, IBSL_forest_abund,
#                            by=c("tax_id", "species", "genus", "family", "order", "class", "phylum", "superkingdom"))
# abundance_tax <- full_join(abundance_tax, IBSL_peri_abund,
#                            by=c("tax_id", "species", "genus", "family", "order", "class", "phylum", "superkingdom"))
# 
# #input 0 for count=NA
# abundance_tax <- abundance_tax %>% replace(is.na(.), 0) %>%
#   dplyr::rename(taxID = tax_id) %>%
#   filter(!taxID=="unassigned") %>% #remove 'unassigned' ID
#   arrange(taxID)
# 
# ########### this little bit of code will cut down the sample ID columns ###############
# ###### https://community.rstudio.com/t/shorten-column-names-in-data-frame/143129 ######
# 
# long_names <- colnames(abundance_tax)
# 
# short_names <- lapply(
#   X = strsplit(x = long_names, split = "_"),
#   FUN = function(x) paste0(head(x, n = 1)))
# 
# colnames(abundance_tax) <- short_names
# ########################################################################################
# 
# # Save abundance_tax table to a rdata file
# saveRDS(abundance_tax, file = here("abundance_tax_03.01.23.rds"))

# # Restore abundance_tax from the rdata file
# #rows are taxa, columns are samples
# abundance160_tax <- readRDS(file = "abundance_tax_03.01.23.rds")
# ##### ABUNDANCE_TAX HAS ALL 160 SAMPLES, their Emu estimated abundance and the taxonomic info for each taxID
# ##### potential ISSUE: some of the genus field are blank (n=28), they're still in, but maybe remove?

#140 version of abundance_tax (might have entire rows that are 0... idk)
abundance_tax <- readRDS(file = "abundance_tax_03.01.23.rds") %>% #this is abundance_tax with all 160 FDNA samples
  pivot_longer(-c(taxID, species, genus, family, order, class, phylum, superkingdom),
                               names_to = "FDNA", values_to = "count") %>%
  filter(!FDNA %in% removeFDNAs.v) %>% #remove recapture FDNAs (get to the 140)
  pivot_wider(id_cols=c(taxID, species, genus, family, order, class, phylum, superkingdom),
              names_from = "FDNA", values_from = "count")


################## new stuff for V2 - Dec 19 2023 ##########################

###### Rarefy Emu count data for pathogen spp detection #######
##### using vegan:rrarefy() to resample the count numbers #######
##### this is using the 140 data - NO DUPLICATES of recaps ######

# #commented out on 12.19.2023 once I ran it
# #do a lil set.seed so I can reproduce these results
# set.seed(2111994)
# 
# #initialize a list to store results
# rare_counts.ls <- list()
# 
# #for 100 iterations, rarefy count data, save, and repeat
# for(i in 1:100){
# 
#   #rarefy data (one iteration)
#   #running this on abundance_vegan which is the OG abundance140 file 
#     #(same count values as abundance_tax but that has all 160 samples...
#     # and the extra taxonomic info which we don't need yet)
#   #output in vegan format: taxa are the columns, samples are the rows, sampleIDs are row names
#   rareit_vegan <- rrarefy(abundance_vegan, sample=min_n_seq) #min_n_seq = 74,517
# 
#   #save each iteration as an object in a list
#   rare_counts.ls[[i]] <- rareit_vegan
# 
# }
# 
# #collapse the output list (one list item per iteration) down to a df
# #counts of all taxa for all samples, across 100 iterations
# #really big, n_iterations * 140 samples number of rows
# rare_counts.df <- do.call(rbind, rare_counts.ls) %>% 
#   as.data.frame() %>% 
#   rownames_to_column(var="FDNA") %>% #since there are multiple rows named 'FDNA', the iterations get numbered: 'FDNA.1' etc.
#   separate_wider_delim(FDNA, delim=".", names=c("FDNA", "iteration"),
#                        too_few = "align_start") %>% #needed because not all FDNA columns have a ".#" (fills with NA)
#   mutate(iteration = as.numeric(iteration)) %>% #convert iteration to numeric
#   replace_na(list(iteration=0)) #switch NAs to 0s
# #this df has a column for each taxa and columns for FDNA and iteration
# 
# #summarize mean count of each taxa per sample across all iterations
# rare_counts_summary <- rare_counts.df %>% 
#   pivot_longer(-c(FDNA,iteration), names_to="tax", values_to = "count") %>% #long table
#   group_by(FDNA, tax) %>%
#   summarise(rare_count = mean(count),
#             rare_sd = sd(count))
# #output is long format df with FDNA, tax, mean count, and sd of count
# 
# #turn it back to a count matrix
# rare_counts_mat <- rare_counts_summary %>% pivot_wider(id_cols=FDNA, values_from="rare_count", names_from="tax") %>%
#   column_to_rownames(var="FDNA") %>% #taxa are columns, FDNA are row names
#   as.matrix()
# 
# #save it
# saveRDS(rare_counts_mat, file="rare_counts_mat_v12.19.2023.RDS")
#load it
rare_counts_mat <- readRDS(file = "rare_counts_mat_v12.19.2023.RDS")
###### this is now a rarefied count matrix (100 iterations) of the 140 to be used for pathogen detection

#pull tax ID and taxonomic info for all the observed bacteria (from abundance_tax file, which has 160 samples)
taxID_spp <- abundance_tax %>% select(taxID, species)

#join rare_counts_mat with taxID_spp 
#rarefied abundance table - SPECIES LEVEL
#rows are taxa, columns are samples, in numerical order
### ONLY includes the 140 unique mice - NO DUPLICATES for recaps
abundance_spp <- rare_counts_mat %>% t() %>% #transpose so tax are rows
  as.data.frame() %>% rownames_to_column(var="taxID") %>%
  left_join(taxID_spp, by="taxID") %>%
  relocate(species, .after=taxID) %>% select(-taxID) 
#this could be added back in if count needs to be an integer (currently and avg so has decimals)
  # pivot_longer(-species, names_to= "FDNA", values_to = "count") %>%
  # mutate(count = round(count, digits=0)) %>% #round "counts" to nearest whole number
  # pivot_wider(names_from = "FDNA", values_from = "count")

###abundance_spp is the 140 and their avg rarefied counts for all species
    ###df with columns species, and then one column for each FDNA, values are counts
    ###minimum count could be below 10 reads/genus




############ SUBSET ABUNDANCE (count) table to PATHOGENIC SPECIES ONLY #############

# #CODE RUN and COMMENTED OUT on 20 December 2023
# #bacteria on the PHIbase pathogen list (**PLANT-ONLY PATHOGENS REMOVED 3/22/22**)
# phi_pathogens <- read.csv(here("PHIbase_pathogen_species_list.csv"), header=TRUE) %>%
#   separate_wider_delim(species, names=c("genus", "species"), delim=" ",
#                        too_many = "merge") %>%
#   select(!taxID)
# # phi_path_genera <- phi_pathogens %>% group_by(genus) %>% slice(1) %>% select(genus)
# 
# #bacteria on my list (from internet searches and Annie Jahan's Nanopore farm Pero paper)
# my_pathogens <- read.csv(here("bacteria_list.csv"), header=TRUE)
# # my_path_genera <- my_pathogens %>% group_by(genus) %>% slice(1) %>% select(genus)
# 
# # #combine PHIbase and my bacteria list, GENUS LEVEL (remove redundant genera)
# # path_genera <- rbind(phi_path_genera, my_path_genera) %>%
# #   group_by(genus) %>% slice(1)
# # #pull vector of pathogenic genera
# # path_genera.v <- path_genera$genus
# # #save it
# # saveRDS(path_genera.v, file = here("path_genera.v.rds"))
# 
# #combine PHIbase and my list, SPECIES LEVEL
# path_spp <- rbind(phi_pathogens, my_pathogens) %>%
#   arrange(genus, species) %>%
#   distinct(genus, species, .keep_all = TRUE) %>%
#   unite(gspp, genus, species, sep=" ", remove=FALSE)
# write.csv(path_spp, here("Supplementary_Table_S-1_Pathogen_List_v2.csv")) #run on 20 Dec 2023
# path_spp.v <- path_spp$gspp
# 
# #filter abundance path df (essentially abundance matrix) for pathogen genera
# abundance_path <- abundance_spp %>%
#   filter(species %in% path_spp.v) %>%
#   column_to_rownames(var="species")
# 
# # Save abundance_path table to a rdata file
# saveRDS(abundance_path, file = here("abundance_pathsp_12.20.23.rds"))

# Restore abundance_path from the rdata file
#rows are taxa, columns are samples
abundance_path <- readRDS(file = "abundance_pathsp_12.20.23.rds") #has 39 species just like non-rarefied 160 analysis


##########################################################################
#### right here you could re-insert the 'path_genera_NMDS.R' script   ####
##########################################################################

# #summarize number of pathogen species per site
# abundance_path %>%
#   rownames_to_column(var="species") %>%
#   pivot_longer(-species, names_to = "FDNA", values_to = "count") %>%
#   mutate(positive = ifelse(count >= 50, 1, 0)) %>%
#   inner_join(pero_datashort, by="FDNA") %>% #add the metadata
#   group_by(loc_site, species) %>%
#   mutate(n_pos = sum(positive),
#          prev = n_pos / length(FDNA)) %>%
#   mutate(prev_na = ifelse(prev=="0", NA, prev)) %>%
#   slice(1) %>%
#   select(species, loc_site, prev) %>%
#   filter(prev!=0) %>%
#   ungroup() %>% summarise(n_distinct(species)) #18 species
#   group_by(loc_site) %>% summarise(length(species)) #number identified per site



#number of reads per mouse
#number of species per mouse
data <- abundance_path %>%
  rownames_to_column(var="species") %>%
  pivot_longer(-species, names_to = "FDNA", values_to = "n_reads") %>%
  mutate(n_reads = ifelse(n_reads<50, 0, n_reads)) %>%  #change counts 0-50 to 0
  filter(n_reads>0) %>%
  group_by(FDNA) %>% mutate(n_path = n_distinct(species)) %>% ungroup() %>%
  left_join(pero140_datashort, by="FDNA")

n_distinct(data$FDNA) #41 animals had putative pathogens
41/140 #about 29% -- tbh, that's pretty interesting

# Perform pairwise comparisons of n_path by loc_site
compare_means(n_path ~ loc_site,  data = data) #there's nothing significant
# # Visualize: Specify the comparisons you want
# my_comparisons <- list( c("Agricultural_Forest", "Undeveloped_Forest"), 
#                         c("Agricultural_Synanthropic", "Undeveloped_Forest"), 
#                         c("Undeveloped_Forest", "Undeveloped_Synanthropic") )

#plot n pathogens per mouse compare by loc_site -- there's nothing to see
data %>%
  ggplot(aes(y=n_path, x=loc_site, fill=loc_site, alpha=0.7)) +
  geom_violin() +
  geom_jitter(width=0.1, height=0.1, alpha=0.5) + #add the jittered points
  theme_half_open() +
  scale_fill_manual(values=c("#A67326", "#F0C907",
                             "#2c7c94", "#a6d0c8")) +
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = "black") +
  # stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value (option: label="p.signif" )
  stat_compare_means(label.y = 1.8, size=5)  # Add global p-value

####### ermagerd merdeling ########
path_miceIDS <- unique(data$FDNA)

npathmouse <- data %>% group_by(FDNA) %>% slice(1) %>% select(FDNA, n_path) %>%
  mutate(n_path = ifelse(is.na(n_path), 0, n_path))

glmdata <- pero140_datashort %>% 
  mutate(pathfound = ifelse(FDNA %in% path_miceIDS, 1, 0)) %>% #add column for whether mouse had pathbac
  left_join(npathmouse, by="FDNA") %>%
  mutate(landscape_type = factor(landscape_type, levels=c("Undeveloped", "Agricultural")))

#what predicts probability of finding a putative pathogen?
#https://stats.stackexchange.com/questions/561263/poisson-or-binomial-distribution-for-modeling
mod <- glm(pathfound ~ landscape_type + site_type + landscape_type*site_type, 
           family=binomial(link="logit"), data=glmdata)
summary(mod)
plot(mod)

#odds ratio (CI shouldn't cross 1: (0-1) means less likely, (1+) means more likely)
exp(cbind(coef(mod), confint(mod)))


#summary table of mixed model output
#https://cran.r-project.org/web/packages/sjPlot/vignettes/tab_mixed.html
library(sjPlot)
#save it https://stackoverflow.com/questions/67280933/how-to-save-the-output-of-tab-model
tab_model(mod, file="SUPP_pathprobability_glm.doc")


#model diagnostics
#https://rpubs.com/benhorvath/glm_diagnostics

library(statmod)
plot(density(qresid(mod)))
#residuals vs fitted
scatter.smooth(predict(mod, type='response'), rstandard(mod, type='deviance'), col='gray')
#quartile residuals
scatter.smooth(predict(mod, type='response'), qresid(mod), col='gray')

scatter.smooth(predict(mod, type='response'), glmdata$landscape_type, col='hotpink')
scatter.smooth(predict(mod, type='response'), glmdata$site_type, col='hotpink')

qqnorm(qresid(mod)); qqline(qresid(mod))


  



#n unique pathogen species per loc_site
abundance_path %>%
  rownames_to_column(var="species") %>%
  pivot_longer(-species, names_to = "FDNA", values_to = "n_reads") %>%
  mutate(n_reads = ifelse(n_reads<50, 0, n_reads)) %>%  #change counts 0-50 to 0
  filter(n_reads>0) %>%
  left_join(pero140_datashort, by="FDNA") %>%
  group_by(loc_site) %>%
  summarise(rich = n_distinct(species))
#9 pathogens in agforest, 5-6 in others

#total n pathogens identified (doubles okay) per loc_site (basically counting the number of rows)
#the glm above is basically testing this in model form (ag landscape more likely to detect positive)
abundance_path %>%
  rownames_to_column(var="species") %>%
  pivot_longer(-species, names_to = "FDNA", values_to = "n_reads") %>%
  mutate(n_reads = ifelse(n_reads<50, 0, n_reads)) %>%  #change counts 0-50 to 0
  filter(n_reads>0) %>%
  left_join(pero140_datashort, by="FDNA") %>%
  group_by(loc_site) %>%
  summarise(sum = length(FDNA))
#22 and 17 total positives in ag -- 7 and 7 positives in undev

#number, ppn of positive mice per loc_site
glmdata %>% group_by(loc_site) %>%
  summarise(npath = sum(pathfound),
            ntotal = length(FDNA),
            ppn = npath/ntotal)
#number, ppn of positive mice per landscape
glmdata %>% group_by(landscape_type) %>%
  summarise(npath = sum(pathfound),
            ntotal = length(FDNA),
            ppn = npath/ntotal)



#(avg, rarefied) reads per path per mouse, filtered to >=50 reads per pathogen species
checkpath <- abundance_path %>%
  rownames_to_column(var="species") %>%
  pivot_longer(-species, names_to = "FDNA", values_to = "count") %>%
  filter(count>50) %>% arrange(count) %>% #53 total pathogen detections
  # filter(count<200) #37 pathogens under 200 reads
  # filter(count>500) #7 pathogens over 500 reads

#mean path per mouse, per loc_site
abundance_path %>%
  rownames_to_column(var="species") %>%
  pivot_longer(-species, names_to = "FDNA", values_to = "count") %>%
  mutate(present = ifelse(count<50, 0, 1)) %>% #0,1 for pathbac found (if >50 reads)
  left_join(pero140_datashort, by="FDNA") %>%
  group_by(loc_site, FDNA) %>% summarise(n_path = sum(present)) %>% #number pathbac per mouse
  ungroup() %>% group_by(loc_site) %>%
  summarise(mean_path = mean(n_path),
            sd = sd(n_path),
            min = min(n_path),
            max = max(n_path)) #mean path per mouse in each loc_site

#pathogens found and number of mice positive
#adjust to see full list (n=13) or top 5
abundance_path %>%
  rownames_to_column(var="species") %>%
  pivot_longer(-species, names_to = "FDNA", values_to = "count") %>%
  left_join(pero140_datashort, by="FDNA") %>%
  # group_by(loc_site, ear_tag_number, species) %>% #combine across FDNA samples per mouse
  # summarize(count = sum(count)) %>% #stop here for summed read count/bacteria/mouse
  mutate(found = ifelse(count<50, 0, 1)) %>%
  ungroup() %>% group_by(species) %>%
  summarise(n_mice = sum(found)) %>% #number of mice with each pathbac
  arrange(desc(n_mice)) %>%
  # filter(n_mice > 0) #list of 13
  filter(n_mice >3) #top 5


# #number of mice per pathogen, separated by loc_site
# #not terribly interesting but here it is
# library(kableExtra)
# abundance_path %>%
#   rownames_to_column(var="species") %>%
#   pivot_longer(-species, names_to = "FDNA", values_to = "count") %>%
#   inner_join(pero140_datashort, by="FDNA") %>%
#   mutate(found = ifelse(count<50, 0, 1)) %>%
#   ungroup() %>% group_by(loc_site, species) %>%
#   summarise(n_mice = sum(found)) %>% #number of mice with pathbac per loc_site
#   filter(species %in% topfive$species) %>%
#   arrange(species, loc_site) %>%
#   pivot_wider(id_cols=species, names_from = loc_site, values_from = n_mice) %>%
#   mutate(total = rowSums(across(where(is.numeric)))) %>%
#   arrange(desc(total)) %>%
#   kbl(col.names = c("Bacteria Species",
#                     "A-F", "A-S",
#                     "U-F", "U-S",
#                     "Total")) %>%
#   kable_styling(full_width=FALSE) %>%
#   row_spec(0, bold=T, color="black", background="#DAD7D7")


################################################################################
#### right here you could re-insert the 'pathogen_prev_heatmaps.R' script   ####
################################################################################

facet_labs <- as_labeller(c("Agricultural_Forest" = "Agricultural-Forest", 
                            "Agricultural_Synanthropic" = "Agricultural-Synanthropic",
                            "Undeveloped_Forest" = "Undeveloped-Forest", 
                            "Undeveloped_Synanthropic" = "Undeveloped-Synanthropic"))

library(wesanderson)
wes_pal <- wes_palette("Zissou1", type="continuous") #color palette for gradient

abundance_path %>%
  rownames_to_column(var="species") %>%
  pivot_longer(-species, names_to = "FDNA", values_to = "count") %>%
  mutate(count = ifelse(count<50, 0, count)) %>%  #change counts 0-50 to 0
  group_by(species) %>%  mutate(found = sum(count)) %>% filter(found > 0) %>% select(!found) %>% #remove species with no >50 readcounts
  ungroup() %>%
  left_join(pero140_datashort, by="FDNA") %>% #add the metadata
  # select(!FDNA) %>% #drop FDNA column, we'll go off ear_tag so we group together data from recaps
  ggplot(aes(x=FDNA, y=species)) +
  facet_grid(. ~ loc_site, 
             scales="free_x", space="free",
             labeller=facet_labs) +
  geom_tile(aes(fill=log(count))) +
  scale_fill_gradientn(colors = wes_pal, na.value="#f5f5f5", 
                       name="Ln(Read Count)") +
  scale_y_discrete(limits=rev) + #flip the y axis to be alphabetical top to bottom
  theme_light() +
  theme(axis.text.y = element_text(size=12),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size=12),
        legend.text= element_text(size=12),
        legend.position = "right",
        strip.text = element_text(size=12, color="black"),
        strip.background=element_rect(fill="#d3d3d3"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave(here("path_heatmap_v2.eps"),
       plot= last_plot(),
       dpi=300,
       units="px")

###################################################################################















###########################################################################################
### new for revised ms dec 21 2023 :( why am I always in crunchmode right before Christmas?

### CREATE A PHYLOSEQ OBJECT
#following: https://joey711.github.io/phyloseq/import-data.html
#create a tax matrix with taxID as row name
tax_mat <- abundance_tax %>% #140 samples worth of taxIDs
  select("taxID", "species", "genus", "family", "order", "class", "phylum", "superkingdom") %>%
  arrange(taxID) %>%
  column_to_rownames(var="taxID") %>% #row names are taxIDs
  select("superkingdom", "phylum", "class", "order", "family", "genus", "species") %>% #reorder tax to be descending L to R
  as.matrix()

# #convert RAREFIED abundance to matrix, taxID as rowname
# otu_mat2 <- rare_counts_mat %>% #rarefied abundance counts for the 140
#   as.data.frame() %>%
#   rownames_to_column(var="FDNA") %>%
#   pivot_longer(-FDNA, names_to = "taxID", values_to = "count") %>%
#   mutate(count = ifelse(count>50, count, 0)) %>% #remove species with <50 reads
#   pivot_wider(id_cols="taxID", names_from = "FDNA", values_from = "count") %>%
#   column_to_rownames(var="taxID") %>%
#   as.matrix()

#don't rarefy abundance data before ANCOM, prefilter out low abundance
# https://forum.qiime2.org/t/do-ancom-and-lefse-always-require-rarefaction-at-specific-sequencing-depth-feature-counts/18804/3

#convert RAW abundance to matrix, taxID as rowname
#RAW (not rarefied) abundance counts - but counts <50 removed
otu_mat <- abundance_tax %>%
  select(!c("species", "genus", "family", "order", "class", "phylum", "superkingdom")) %>%
  pivot_longer(-taxID, names_to = "FDNA", values_to = "count") %>%
  mutate(count = ifelse(count>50, count, 0)) %>%
  filter(!FDNA %in% removeFDNAs.v) %>% #remove recapture animals (get to n=140)
  pivot_wider(id_cols="taxID", names_from = "FDNA", values_from = "count") %>%
  column_to_rownames(var="taxID") %>%
  as.matrix()

#convert metadata to phyloseq-usable format (FDNA as row name)
pero_datashort_mat <- pero140_datashort %>% 
  mutate(loc_site = fct_relevel(loc_site, c("Undeveloped_Forest", "Undeveloped_Synanthropic",
                                       "Agricultural_Forest", "Agricultural_Synanthropic"))) %>%
  column_to_rownames(var="FDNA")


#transform all to phyloseq objects
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(pero_datashort_mat)


#combine into a single phyloseq object
abundance_phylo <- phyloseq(OTU, TAX, samples)
# abundance_phylo


# #agglomerate by genus
# abund_phylo_genus <- tax_glom(abundance_phylo, taxrank = "genus")
# taxID_genus <- abund_phylo_genus@tax_table %>% as.data.frame() %>% rownames_to_column(var="taxID") %>% select(taxID, genus)

######### ANCOM-BC ##########
#general intro to differential abundance 
  #https://microbiome.github.io/OMA/differential-abundance.html

#the ANCOMBC package for R version 4.1
  #https://rdrr.io/bioc/ANCOMBC/man/ancombc.html

#maybe relevant question about not having a reference group
  #https://github.com/FrederickHuangLin/ANCOMBC/issues/56


#this is documentation for the most recent (Jan 2 2024) version of the package (R versions 4.3)
  #https://bioconductor.org/packages/release/bioc/manuals/ANCOMBC/man/ANCOMBC.pdf
#also the Jan 2024 update:
  #https://rdrr.io/github/FrederickHuangLin/ANCOMBC/man/ancombc.html

#GitHub page
  #https://github.com/FrederickHuangLin/ANCOMBC



# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ANCOMBC")
library(ANCOMBC)

set.seed(2111994)

#ancom-bc2 for multiple pairwise comparisons
output <- ancombc2(data = abundance_phylo, 
                   tax_level = "species",
                   fix_formula = "loc_site + sex + reproductive", 
                   rand_formula = NULL, #no random effects
                   p_adj_method = "holm", #default
                   pseudo_sens = TRUE, #default
                   prv_cut = 0.10, #discard taxon with low prevalence (0.1=default)
                   lib_cut = 0, #do not discard any samples
                   s0_perc = 0.05, #default
                   group = "loc_site", 
                   struc_zero = TRUE, neg_lb = TRUE,
                   alpha = 0.05, 
                   n_cl = 2, #run 2 nodes in parallel (faster computing)
                   verbose = TRUE,
                   global = TRUE, #perform global test
                   pairwise = TRUE, #pairwise directional test
                   dunnet = TRUE, #Dunnett's type of test, even more conservative
                   trend = FALSE, #trend test
                   iter_control = list(tol = 1e-2, max_iter = 20,
                                      verbose = FALSE), #default
                   em_control = list(tol = 1e-5, max_iter = 100), #default
                   lme_control = NULL, #not using a mixed model
                   mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), #default
                   trend_control = NULL) #no trend test

res = output$res[,1:5] %>% rename(species=taxon)
# res_pair = output$res_pair
res_global = output$res_global %>%
  rename(taxID = taxon)


# out <- ancombc(
#   data = abundance_phylo, 
#   formula = "loc_site + sex + reproductive", 
#   p_adj_method = "holm", 
#   prv_cut = 0.1 , #default
#   lib_cut = 0, #default, do not discard any samples
#   group = "loc_site", #character
#   struc_zero = TRUE, 
#   neg_lb = TRUE, #set to TRUE if n per group > 30
#   tol = 1e-5, #default
#   max_iter = 100, #default
#   conserve = TRUE, 
#   alpha = 0.05, #default
#   global = TRUE 
# )
# 
# res <- out$res
# res_global <- out$res_global %>% #pull global comparisons
#   rename(taxID = taxon) %>%
#   left_join(taxID_spp, by="taxID") #join species names by taxID

##########################################################################

#recall the list of 209 pathogenic species
path_spp <- read.csv(here("Supplementary_Table_S-1_Pathogen_List_v2.csv")) #run on 20 Dec 2023
path_spp.v <- path_spp$gspp

#1.3.24!! trying code from the tutorial (Oct 24, 2023)
#https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html

#vector of species with significant differences
sig_taxa <- res_global %>% filter(diff_abn == TRUE) %>% .$taxID

# tab_lfc <- res$lfc %>% rename(taxID = taxon) %>% left_join(taxID_spp, by="taxID")


df_locsite <- res %>% select(-"lfc_(Intercept)") %>%
  filter(species %in% sig_taxa) %>% #only keep the significantly different taxa
  # filter(species %in% path_spp.v) %>%
  rename("A-S" = "lfc_loc_siteAgricultural_Synanthropic",
         "A-F" = "lfc_loc_siteAgricultural_Forest",
         "U-S" = "lfc_loc_siteUndeveloped_Synanthropic")

df_heat = df_locsite %>%
  pivot_longer(cols = -one_of("species"),
               names_to = "loc_site", values_to = "value") %>%
  mutate(value = round(value, 2))
df_heat$species = factor(df_heat$species, levels = sort(sig_taxa))

lo = floor(min(df_heat$value))
up = ceiling(max(df_heat$value))
mid = (lo + up)/2
p_heat = df_heat %>%
  ggplot(aes(x = loc_site, y = species, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(loc_site, species, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, 
       title = "?BETA? for globally significant (putative pathogen species)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
p_heat



#alternatively
#https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/
#boxplots to compare abundance (threshold to at least 50 reads) of pathogen species
  #across loc_sites
abundance_spp %>% pivot_longer(-species, names_to = "FDNA", values_to = "count") %>%
  filter(count>50) %>%
  filter(species %in% path_spp.v) %>% 
  left_join(tag_loc, by="FDNA") %>%
  unite(loc_site, landscape_type, site_type) %>%
ggplot(aes(x = loc_site, y = count)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = species), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ species, scales = "free")
  
  














######################### PHYLA LEVEL - relabundance ################################

#create relative abundance at phylum level

# Restore abundance_tax from the rdata file
#rows are taxa, columns are samples
abundance_tax <- readRDS(file = "abundance_tax_03.01.23.rds")

## ABUNDANCE_TAX has all the species, remove the recapture FDNAs

phy_relabund <- abundance_tax %>%
  select(!c(taxID, species, genus, family, class, order, superkingdom)) %>%
  pivot_longer(-phylum, names_to = "FDNA", values_to = "count") %>%
  filter(!FDNA %in% removeFDNAs.v) %>%
  group_by(FDNA, phylum) %>%
  summarise(phy_count = sum(count)) %>% #sum counts to phylum level per FDNA
  ungroup() %>% group_by(FDNA) %>%
  mutate(total_count = sum(phy_count),
         rel_abund = phy_count/total_count) %>%
  mutate(phylum = case_when(rel_abund <0.01 ~ "Other",
                            rel_abund >= 0.01 ~ phylum)) %>% #group phyla present <1% as 'Other'
  mutate(phylum = as.factor(phylum)) %>%
  mutate(phylum = fct_relevel(phylum, c("Bacteroidetes",
                                        "Deferribacteres",
                                        "Firmicutes",
                                        "Candidatus Melainabacteria",
                                        "Proteobacteria", 
                                        "Other"))) %>%
  left_join(tag_loc, by="FDNA") %>%
  unite(loc_site, landscape_type, site_type, remove=FALSE)

# #check phyla present at 1% or 5% rel abundance
# phy_relabund %>% filter(rel_abund >=  0.01) %>%
#   ggplot(aes(x=FDNA, y=rel_abund, fill=phylum)) +
#   geom_bar(stat="identity", aes(fill=phylum)) +
#   facet_wrap(~loc_site, scales="free_x", nrow=1)

#plot (full) relative abundance
phy_relabund %>%
  ggplot(aes(x=FDNA, y=rel_abund, fill=phylum)) +
  geom_bar(stat="identity", aes(fill=phylum), width=1) +
  scale_fill_manual(labels=c("Bacteroidetes",
                             "Deferribacteres",
                             "Firmicutes",
                             "Melainabacteria",
                             "Proteobacteria", 
                             "Other"),
                    values=c("#a65852", "#c4aa23",
                             "#a6d0c8", "#fbe45b",
                             "#2c7c94", "#A0999A")) + #wesanderson French Dispatch Shadow
  facet_wrap(~loc_site, scales="free_x", nrow=1,
             labeller=facet_labs) +
  labs(y="Relative Abundance",
       fill="Phylum") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=18),
        axis.text.y = element_text(size=15),
        strip.text.x = element_text(size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14),
        panel.background = element_blank())

ggsave("NEWphyrelabundplot.eps", plot=last_plot(), units="px", dpi=600)

#scale fill manual colors: (previous, colorblind-friendly)
# values=c("#D84662", "#FFD3EC",
#          "#68B5E2", "#FDFC86",
#          "#90DB6A", "#A0999A")


#check colorblind friendliness
cvdPlot(plot)

# #boxplot
# phy_relabund %>%
#   filter(rel_abund>0.05) %>%
#   ggplot(aes(x=phylum, y=rel_abund, fill=loc_site)) +
#   geom_boxplot()

#mean sd relabundance by site
# library(kableExtra)
phy_relabund %>%
  filter(rel_abund>0.05) %>%
  group_by(phylum, loc_site) %>%
  summarise(mean = round(mean(rel_abund), digits=4)*100,
            sd = round(sd(rel_abund),digits=4)*100) %>%
  mutate(sd = paste0(sd, "%")) %>%
  mutate(mean_sd = paste(mean, "\u00B1", sd, sep=" ")) %>%
  select(!c(mean, sd)) %>%
  pivot_wider(names_from = phylum, values_from = mean_sd) %>%
  select(loc_site, Firmicutes, Proteobacteria, Bacteroidetes) %>%
  separate(loc_site, c("landscape_type", "site_type"), remove=TRUE) %>%
  kbl(col.names = c("Landscape", "Habitat", 
                    "Firmicutes", "Proteobacteria", "Bacteroidetes")) %>%
  kable_styling(full_width=FALSE)

# #mean sd across all samples
# phy_relabund %>%
#   filter(rel_abund>0.05) %>%
#   group_by(phylum) %>%
#   summarise(mean = round(mean(rel_abund), digits=4)*100,
#             sd = round(sd(rel_abund), digits=4)*100)


################################################################################################################




















##### NanoQ summary stats (for the 140, no recaps!)

##### quick thing first - grab the FDNA ids for recapped animals to remove #####
#Restore Nanopore-sampled pero160 data from the rdata file
pero_data <- readRDS(file = "pero_data_04.19.23.rds")
#subset pero_data to see only animals with multiple entries
recapped <- pero_data %>%
  filter(location=="CCESR") %>%
  group_by(ear_tag_number) %>%
  arrange(ear_tag_number) %>%
  mutate(ncaps = length(ear_tag_number)) %>%
  relocate(ncaps, .after=ear_tag_number) %>%
  filter(ncaps>1) %>%
  ungroup()
#every recapped animal was caught in July, so we'll keep those entries and remove others
remove <- recapped %>% filter(month != "July")
removeFDNAs.v <- remove$FDNA
##################################################

#summaries of Pero140 reads for results text
read_sumstats <- read.table(here("full_FilteredReads_summarystats.txt"), header=TRUE) %>%
  relocate(FDNA, .before=n_reads) %>%
  left_join(tag_loc, by="FDNA") %>%
  filter(!FDNA %in% removeFDNAs.v) #remove this line to keep all 160 sequenced samples (includes recaps)

sum(read_sumstats$n_reads)
mean(read_sumstats$mean_quality)
sd(read_sumstats$mean_quality)
mean(read_sumstats$n_reads)
sd(read_sumstats$n_reads)
min(read_sumstats$n_reads)
max(read_sumstats$n_reads)

#summary table
read_sumstats <- read.table(here("full_FilteredReads_summarystats.txt"), header=TRUE) %>%
  relocate(FDNA, .before=n_reads) %>% 
  left_join(tag_loc, by="FDNA") %>%
  filter(!FDNA %in% removeFDNAs.v) %>% #remove this line to keep all 160 sequenced samples (includes recaps)
  group_by(landscape_type, site_type, month) %>%
  summarise(n_sample = length(FDNA),
            mean_n_reads = round(mean(n_reads)/1000, 2),
            sd_n_reads = round(sd(n_reads)/1000,2),
            mean_bp = round(mean(n_bp)/1000000,2),
            sd_bp = round(sd(n_bp)/1000000,2),
            mean_Q = round(mean(mean_quality),2),
            sd_Q = round(sd(mean_quality),2))

read_sumstats_CC <- read.table(here("full_FilteredReads_summarystats.txt"), header=TRUE) %>%
  relocate(FDNA, .before=n_reads) %>% 
  left_join(tag_loc, by="FDNA") %>%
  filter(location=="CCESR") %>%
  filter(!FDNA %in% removeFDNAs.v) %>% #remove this line to keep all 160 sequenced samples (includes recaps)
  group_by(landscape_type, site_type) %>%
  summarise(n_sample = length(FDNA),
            mean_n_reads = round(mean(n_reads)/1000, 2),
            sd_n_reads = round(sd(n_reads)/1000,2),
            mean_bp = round(mean(n_bp)/1000000,2),
            sd_bp = round(sd(n_bp)/1000000,2),
            mean_Q = round(mean(mean_quality),2),
            sd_Q = round(sd(mean_quality),2)) %>%
  mutate(month = "Summer") %>%
  relocate(month, .after=site_type)

read_sumstats <- rbind(read_sumstats, read_sumstats_CC)
read_sumstats$month <-fct_relevel(read_sumstats$month, c("June", "July", "August", "Summer"))

# #pretty summary table - png version
# library(kableExtra)
# read_sumstats %>%
#   arrange(landscape_type, site_type, month) %>%
#   mutate(n_reads = paste(mean_n_reads, "\u00B1", sd_n_reads, sep=" "),
#          n_bp = paste(mean_bp, "\u00B1", sd_bp, sep=" "),
#          Q = paste(mean_Q, "\u00B1", sd_Q, sep=" ")) %>%
#   ungroup() %>%
#   select(!c(mean_n_reads, sd_n_reads, mean_bp, sd_bp, mean_Q, sd_Q)) %>%
#   kbl(format="html", escape=FALSE,
#       col.names = c("Landscape", "Habitat Type", "Month", "N",
#                     "N reads / sample <br>(thousands of reads)", "N basepairs / sample <br>(Mb)",
#                   "Q Score"), align='lllcccc') %>%
#   kable_styling(full_width=FALSE) %>%
#   row_spec(0, bold=T, color="black", background="#DAD7D7") %>%
#   row_spec(1:3, italic=TRUE, background="#F8F6F6") %>%
#   row_spec(4, bold=TRUE) %>%
#   row_spec(5:7, italic=TRUE, background="#F8F6F6") %>%
#   row_spec(8:10, bold=TRUE)

#pretty summary table - table version for updating in Word
sumstatstable <- read_sumstats %>% 
  arrange(landscape_type, site_type, month) %>%
  mutate(n_reads = paste(mean_n_reads, "\u00B1", sd_n_reads, sep=" "),
         n_bp = paste(mean_bp, "\u00B1", sd_bp, sep=" "),
         Q = paste(mean_Q, "\u00B1", sd_Q, sep=" ")) %>%
  ungroup() %>%
  select(!c(mean_n_reads, sd_n_reads, mean_bp, sd_bp, mean_Q, sd_Q)) %>%
  rename(Landscape = landscape_type) %>%
  rename("Habitat Type" = site_type) %>%
  rename(Month = month) %>%
  rename(N = n_sample) %>%
  rename("N reads / sample (thousands of reads)" = n_reads) %>%
  rename("N basepairs / sample (Mb)" = n_bp) %>%
  rename("Q Score" = Q)
write.csv(sumstatstable, here("Table_1-140.csv"))

##to display scientific notation
# Mean_BP_millions = paste(formatC(mean_bp, format = "e", digits = 2),
#                          "\u00B1", formatC(sd_bp, format = "e", digits = 2), sep=" "),
# Mean_Q_Score = paste(mean_Q, "\u00B1", sd_Q, sep=" ")) %>%






#####---------------- THE FOLLOWING CODE WAS RUN with the 160 data and has not been updated ----------------

########### THE FOLLOWING INCLUDES RECAP EVENTS #################

#summary counts of captures per location, site, month
caps_by_sitemonth <- pero_data %>% 
  group_by(location, site_type, month) %>% 
  summarise(count = length(FDNA))

#summary counts of captures per location, site
caps_by_site <- pero_data %>% group_by(location, site_type) %>% summarise(count = length(feces_number))

############################ END ################################



###### ID RECAPS and their FDNA NUMBERS ###########

#subset pero to see only animals with multiple entries
recapped <- pero %>%
  filter(location=="CCESR") %>% 
  group_by(ear_tag_number) %>%
  arrange(ear_tag_number) %>%
  mutate(ncaps = length(ear_tag_number)) %>%
  relocate(ncaps, .after=ear_tag_number) %>%
  filter(ncaps>1) %>%
  ungroup()

#every recapped animal was caught in July, so we'll keep those entries and remove others
remove <- recapped %>% filter(month != "July")
removeFDNAs.v <- remove$FDNA

#how many animals were recapped?
n_distinct(recapped$ear_tag_number) 
#18 animals were recapped at CCESR

#summary table of recapped animal tag# and associated FDNA#s
recapped_FDNAs <- recapped %>% group_by(ear_tag_number) %>%
  summarise(FDNA_numbers = toString(FDNA)) %>%
  ungroup()
#summary table of recapped animal tag# and associated capture month
recapped_months <- recapped %>% group_by(ear_tag_number) %>%
  summarise(months = toString(month)) %>%
  ungroup()
#summary table of recapped animal tag# and their site_type
recapped_sites <- recapped %>% group_by(ear_tag_number) %>%
  summarise(sites = toString(site_type)) %>%
  ungroup()

#count of individuals per site type
pero %>% group_by(ear_tag_number) %>% 
  slice(1) %>% #keep one entry per eartag (only CC animals)
  ungroup() %>%
  group_by(location, site_type) %>%
  summarise(n=length(ear_tag_number))

#what should I do about these animals? This will make the stats messy to compare across months...
  #option: I could compare their microbiomes from each capture, IF they are not different... uh pick one?
  #option: I could drop all these animals from the analysis because they're annoying - but they have good data
  #option: I could pull these from the analysis at first and then... look at them separately?

########### END ID RECAPS and their FDNA #s ###############

#so what does this mean for our sample size?

n_distinct(pero$FDNA) #160 fecal samples with Nanopore data

pero %>%
  filter(location=="IBSL") %>%
  nrow() #71 unique animals at IBSL (since we didn't process recaps)

pero %>%
  filter(location == "CCESR") %>%
  nrow() #89 CAPTURE EVENTS at CCESR

#but 18 animals were recapped and two of them were recapped 2x
#so 20 of the 89 CAPTURE EVENTS are RECAP EVENTS

#therefore we saw 69 unique animals at CCESR
  #with 20 duplicate captures of animals (16 animals caught 2x, 2 animals caught 3x)

#################### END ########################






############# create individual df for each grouping to pull FDNA #s ###################
################## this was used for identifying groups on MSI #########################
# pero %>% filter(location=="IBSL" & site_type=="peridomestic") %>%
#   select(c(month, location, site_type, FDNA)) %>%
#   write.csv(here("IBSL_perid.csv"), row.names=FALSE)
# 
# pero %>% filter(location=="IBSL" & site_type=="forest") %>%
#   select(c(month, location, site_type, FDNA)) %>%
#   write.csv(here("IBSL_forest.csv"), row.names=FALSE)
# 
# pero %>% filter(location=="CCESR" & site_type=="peridomestic" & month=="June") %>%
#   select(c(month, location, site_type, FDNA)) %>%
#   write.csv(here("CCESR_perid_june.csv"), row.names=FALSE)
# 
# pero %>% filter(location=="CCESR" & site_type=="peridomestic" & month=="July") %>%
#   select(c(month, location, site_type, FDNA)) %>%
#   write.csv(here("CCESR_perid_july.csv"), row.names=FALSE)
# 
# pero %>% filter(location=="CCESR" & site_type=="peridomestic" & month=="August") %>%
#   select(c(month, location, site_type, FDNA)) %>%
#   write.csv(here("CCESR_perid_aug.csv"), row.names=FALSE)
# 
# pero %>% filter(location=="CCESR" & site_type=="forest" & month=="June") %>%
#   select(c(month, location, site_type, FDNA)) %>%
#   write.csv(here("CCESR_forest_june.csv"), row.names=FALSE)
# 
# pero %>% filter(location=="CCESR" & site_type=="forest" & month=="July") %>%
#   select(c(month, location, site_type, FDNA)) %>%
#   write.csv(here("CCESR_forest_july.csv"), row.names=FALSE)
# 
# pero %>% filter(location=="CCESR" & site_type=="forest" & month=="August") %>%
#   select(c(month, location, site_type, FDNA)) %>%
#   write.csv(here("CCESR_forest_aug.csv"), row.names=FALSE)
######################## END write CSV for each group #################################





######### HOW MANY CAPTURE FATALITIES? ###############
# #read in the data
# pero <- read.csv(here("MN_pero_raw_labeled_v10.8.19.csv"), na.strings=c(""," ","NA"))
# #tidy columns
# pero <- clean_names(pero)
# 
# #who died?
# pero <- pero %>% drop_na(species) #477 capture events (not all unique animals)
# pero_dead <- pero %>% filter(dead=="1") #16 dead animals (not just pspp)
# 16/477 #3.35% of all captures ended up dead
######### END CAPTURE FATALITIES ##################