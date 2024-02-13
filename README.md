# Minnesota <i>Peromyscus</i> mouse microbiome
This repository contains the R code and necessary data files for an analysis of the fecal microbiome of <i>Peromyscus</i> mice trapped in Minnesota, USA in Summer 2019.

These analyses are associated with the publication: "Microbiome diversity and zoonotic bacterial pathogen prevalence in Peromyscus mice from agricultural landscapes and synanthropic habitat". Janine Mistrick, Evan J. Kipp, Sarah I. Weinberg, Collin C. Adams, Peter A. Larsen, & Meggan E. Craft.

This repository has been archived with Zenodo: LINK

---

## The files in this repository are as follows:

<b>R script for all the analyses:</b> <i>pero.R</i> 

In this code, we quantify the alpha and beta diversity of the rodent fecal microbiomes and compare diversity and community composition between landscape and habitat types. We identify the presence of putative pathogenic bacteria and examine landscape- and habitat-level factors influencing the prevalence of potentially pathogenic bacterial species in the sampled rodent fecal microbiomes.


<b>The following data files called in the R script are also included in this repository:</b>

.tsv files containing bacterial read count estimates for all sequenced fecal samples and all identified bacteria (based on Emu taxonomic identification algorithm)

- <i>C_f_aug_unfilt_emu-combined-tax_id-counts.tsv</i> - read counts for Cedar Creek (agricultural) forest sites in August 2019
- <i>C_f_july_unfilt_emu-combined-tax_id-counts.tsv</i> - read counts for Cedar Creek (agricultural) forest sites in July 2019
- <i>C_f_june_unfilt_emu-combined-tax_id-counts.tsv</i> - read counts for Cedar Creek (agricultural) forest sites in June 2019
- <i>C_p_aug_unfilt_emu-combined-tax_id-counts.tsv</i> - read counts for Cedar Creek (agricultural) synanthropic sites in August 2019
- <i>C_p_july_unfilt_emu-combined-tax_id-counts.tsv</i> - read counts for Cedar Creek (agricultural) synanthropic sites in July 2019
- <i>C_p_june_unfilt_emu-combined-tax_id-counts.tsv</i> - read counts for Cedar Creek (agricultural) synanthropic sites in June 2019
- <i>I_f_unfilt_emu-combined-tax_id-counts.tsv</i> - read counts for Itasca (undeveloped) forest sites in July 2019
- <i>I_p_unfilt_emu-combined-tax_id-counts.tsv</i> - read counts for Itasca (undeveloped) synanthropic sites in July 2019

Raw sequence data are accessioned with the NCBI Sequence Read archive: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1068550

Additional data files:

- <i>MN_Pero_field_data_v10.8.19.csv</i> - the metadata file of all information recorded in the field for all sampled rodents (i.e., all sequenced fecal samples)
- <i>sampleID_to_barcode.csv</i> - file linking ID numbers for fecal samples collected in the field to the Nanopore sequencing run information and specific sample barcode ID
- <i>PHIbase_pathogen_species_list.csv</i> - trimmed csv file of bacteria species from the PHIbase 'Pathogens' dataset (http://www.phi-base.org/)
- <i>bacteria_list.csv</i> - additional csv file list of zoonotic and foodborne bacterial pathogens
- <i>full_FilteredReads_summarystats.txt</i> - file summarizing read information for Nanopore sequencing runs including all 160 sequenced fecal samples (140 samples used in analysis for manuscript)
