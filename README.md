This repo contains the R code and necessary files to run the analyses for the analysis of the fecal microbiome of rodents trapped in Minnesota, USA in Summer 2019.
These analyses as associated with a manuscript in review: " TITLE " - current preprint version can be found on Authorea: " LINK "

The files in this repo are as follows:

pero.R - the R script for all the analyses

And the following data files called in the R script:
.tsv files containing read count estimates for all sequenced fecal samples and all identified bacteria (based on Emu taxonomic identification algorithm)
C_f_aug_unfilt_emu-combined-tax_id-counts.tsv - read counts for Cedar Creek (agricultural) forest sites in August 2019
C_f_july_unfilt_emu-combined-tax_id-counts.tsv - read counts for Cedar Creek (agricultural) forest sites in July 2019
C_f_june_unfilt_emu-combined-tax_id-counts.tsv - read counts for Cedar Creek (agricultural) forest sites in June 2019
C_p_aug_unfilt_emu-combined-tax_id-counts.tsv - read counts for Cedar Creek (agricultural) synanthropic sites in August 2019
C_p_july_unfilt_emu-combined-tax_id-counts.tsv - read counts for Cedar Creek (agricultural) synanthropic sites in July 2019
C_p_june_unfilt_emu-combined-tax_id-counts.tsv - read counts for Cedar Creek (agricultural) synanthropic sites in June 2019
I_f_unfilt_emu-combined-tax_id-counts.tsv - read counts for Itasca (undeveloped) forest sites in July 2019
I_p_unfilt_emu-combined-tax_id-counts.tsv - read counts for Itasca (undeveloped) synanthropic sites in July 2019

MN_pero_raw_labeled_v10.8.19.csv - the metadata file of all information recorded in the field for all sampled rodents (ie all sequenced fecal samples)

sampleID_to_barcode.csv - file linking ID numbers for fecal samples collected in the field to the Nanopore sequencing run information and specific sample barcode ID

PHIbase_pathogen_species_list.csv - trimmed csv file of bacteria species from the PHIbase 'Pathogens' dataset (http://www.phi-base.org/)
bacteria_list.csv - additional csv file list of zoonotic and foodborne bacterial pathogens

full_FilteredReads_summarystats.txt - file summarizing read information for Nanopore sequencing runs including all 160 sequenced fecal samples (140 samples used in analysis for manuscript)
