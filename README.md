# Minnesota <i>Peromyscus</i> mouse microbiome
This repo contains the R code and necessary files to run the analyses for the analysis of the fecal microbiome of rodents trapped in Minnesota, USA in Summer 2019.

These analyses as associated with the publication: "Microbiome diversity and zoonotic bacterial pathogen prevalence in Peromyscus mice from agricultural landscapes and synanthropic habitat". Janine Mistrick, Evan Kipp, Sarah Weinberg, Collin Adams, Peter Larsen, & Meggan Craft.

This repository has been archived with Zenodo: LINK

---

## The files in this repo are as follows:

- <i>pero.R</i> - the R script for all the analyses


And the following data files called in the R script:

.tsv files containing read count estimates for all sequenced fecal samples and all identified bacteria (based on Emu taxonomic identification algorithm)

- <i>C_f_aug_unfilt_emu-combined-tax_id-counts.tsv</i> - read counts for Cedar Creek (agricultural) forest sites in August 2019
- <i>C_f_july_unfilt_emu-combined-tax_id-counts.tsv</i> - read counts for Cedar Creek (agricultural) forest sites in July 2019
- <i>C_f_june_unfilt_emu-combined-tax_id-counts.tsv</i> - read counts for Cedar Creek (agricultural) forest sites in June 2019
- <i>C_p_aug_unfilt_emu-combined-tax_id-counts.tsv</i> - read counts for Cedar Creek (agricultural) synanthropic sites in August 2019
- <i>C_p_july_unfilt_emu-combined-tax_id-counts.tsv</i> - read counts for Cedar Creek (agricultural) synanthropic sites in July 2019
- <i>C_p_june_unfilt_emu-combined-tax_id-counts.tsv</i> - read counts for Cedar Creek (agricultural) synanthropic sites in June 2019
- <i>I_f_unfilt_emu-combined-tax_id-counts.tsv</i> - read counts for Itasca (undeveloped) forest sites in July 2019
- <i>I_p_unfilt_emu-combined-tax_id-counts.tsv</i> - read counts for Itasca (undeveloped) synanthropic sites in July 2019

Additional data files:

- <i>MN_pero_raw_labeled_v10.8.19.csv</i> - the metadata file of all information recorded in the field for all sampled rodents (ie all sequenced fecal samples)
- <i>sampleID_to_barcode.csv</i> - file linking ID numbers for fecal samples collected in the field to the Nanopore sequencing run information and specific sample barcode ID
- <i>PHIbase_pathogen_species_list.csv</i> - trimmed csv file of bacteria species from the PHIbase 'Pathogens' dataset (http://www.phi-base.org/)
- <i>bacteria_list.csv</i> - additional csv file list of zoonotic and foodborne bacterial pathogens
- <i>full_FilteredReads_summarystats.txt</i> - file summarizing read information for Nanopore sequencing runs including all 160 sequenced fecal samples (140 samples used in analysis for manuscript)
