# Analysis scripts for olfactory GWAMA
Analysis scripts used for further analysis of olfactory meta-GWAS.

## 1. Settings
- **00_SourceFile_smelling_meta.R**: Contains common paths, libraries and helper functions used by other analysis scripts. Variables might need adjustment depending on your system. An overview of used packages and versions is in additional_information/r_package_list.xlsx.

- **00a_create_folders.R**: Create folder structure assumed by the analysis scripts.

- **00b_convert_file_format.R**: Convert GWAMA summary statistic files from GWAS catalog SSF format into the custom format of the analysis scripts. SSF formatted files should be located in folder '01_data_gwcat'.

## 2. Analysis Scripts
### Locus definition:

- **01a_locus_definition.R**: Specifies $\pm$ 250kb regions around genome-wide associated SNPs and merges overlapping loci, keeping the strongest associated variant as index SNP. Is performed on all trait group combinations  to identify loci that simultaneously control multiple odours.

-  **01b_locus_definition_suggestive.R**: Same procedure as 01a_locus_definition.R but with lowered significance threshold to further include suggestive variants. Genome-wide significant loci will be included in this step to prevent the identification of suggestive loci at positions where genome-wide significant loci are located.

- **01c_cleaning_locus_definition_suggestive.R**: Removes loci with genome-wide significance from locus definition with suggestive hits and assigns the suggestive hits a new ID. Used so that hits with and without genome-wide significance can be analysed separately.

- **01d_cleaning_solitary_hits.R**: Removing loci where no SNPs besides the index variant lie above the suggestive threshold. Used to limit the number of false-positive associations. Only applied to the suggestive SNP analysis.

- **01e_attatch_rs_id.R**: Addition of rsID to identified index SNPs. Script has to be run in interactive mode to enable manual selection of alleles in case they are not identical between reference and our SNP.

- **01f_extract_snp_numbers.R**: Extraction of SNP numbers for each locus and update of the locus definition table.

### COJO
- **02a_cojo.R**: Performs COJO for all SNPs within a locus region. Procedure is iterated for all loci in the locus definition file. Statistics are taken from top-associated phenotype.

- **02b_cojo_literature.R**: Perform COJO conditional on known literature variants to check for independent hits in our analysis.

- **02c_sex_conditioning.R**: Conditioning on best associated variant of the opposite sex to check the claim of independent associations between male and female signals.

- **02d_conditioning_chr11.R**: For three regions in relative proximity on chromosome 11, check if the respective index variants keep their significance when conditioning on the index variants of the other two loci.

### Bioinformatic annotation

>[!IMPORTANT] 
>*Note that the in-house annotation pipeline will not be provided.*

- **03a_credible_set_definition.R**: Calculation of ABFs based on Wakefields method for each SNP within a region. If independent SNPs are present in a region, conditional statistics are used.

- **03b_cs_annotation_preparation.R**: Selects SNPs that are part of the 99% credible set of a locus, performs down-lifting of their positions to be compatible with annotation resources and reformats data to the format of our internal annotation pipeline.

- **03c_suggestive_hits_annotation_preparation.R**: Down lifting and reformatting for annotation of suggestive hits. Only the index variant not the credible set will be annotated.

- **03d_find_hre.R**: Checks if annotated genes are in references for AREs and EREs via matching their gene symbols.

- **03e_collect_gene_list.R**: Collect list of all genes near variants of the CS as basis for candidate gene identification.

- **03f_comprehensive_eqtl_coloc.R**: Test for colocalization between odour signals and previously reported eQTLs based on GTEx catalogue.

### Sex interaction

- **04a_sex_interaction.R**: Compares effect sizes between sex stratified analyses and plots the result as a beta-beta-plot. Includes visualisation for the publication in the sections 'Version for publication' and 'Unlabelled version for suggestive hits'.

- **04b_coloc_sex_ia.R**: Performs COLOC between male and female signals for SNPs within loci.

### Mendelian randomization
- **06a_MR1_01_GetInstruments.R**: Instrument selection for MR1 (sex-hormones on olfactory identification) and lookup of the SNP-exposure statistics.

- **06b_MR1_02_SNP_Outcome.R**: Lookup of the SNP-outcome (olfactory perception) statistics for MR1.

- **06c_MR1_03_runMR.R**: Main Script for the MR1.

- **06d_MR1_04_getFDR_plots.R**: FDR correction, plotting of forest and scatter plots and recalculation with alternative MR methods for MR1.

- **06e_MR2_01_GetData.R**: Collect statistics for MR2 forward direction (olfaction on neurodegenerative diseases).

- **06f_MR2_02_MRbySNP.R**: Perform MR2 for single SNPs.

- **06g_MR2_03_MR_IVW.R**: MR2 IVW and scatter plot.

- **06h_MR2r_01_GetInstruments.R**: Collect statistics for reverse analysis of MR2 (neurodegenerative diseases on olfaction).

- **06i_MR2r_02_SNP_Outcome.R**: Lookup of the SNP-outcome (olfactory identification).

- **06j_MR3_MR_IVW.R**: Perform MR2 reverse direction IVW.

- **06k_MR2r_04_getFDR_plots.R**: FDR correction, plotting of forest and scatter plots and recalculation with alternative MR methods for MR2 reverse analysis.

- **06l_MR_summary.R**: Create summary tables for the supplementary material.

## 3. Figures
- **F1_miami_plot.R**: Plot of combined, female and male analysis groups arranged in a Miami plot with annotation of index variants. Annotation data has to be specified manually.

- **Figure 2**: See 04a_sex_interaction.R.

- **F3_coloc_heatmap.R**: Presents COLOC results between sexes in form of a heatmap.

## 4. Supplemental Figures
- **SF1_qq_plots.R**: Q-Q plots with $\lambda$-values of quality filtered SNPs.

- **SF2a_RA_plots.R**: Regional association plots with LD reference based on LIFE-Adult. Requires plot routine from SF2_RegAssocPlot_hg38_EnsemblGenes_LIFE_LD.R.

- **SF2b_RA_plot_LIFE_LD_Region7_all_pairs.R**: Regional association plot for region 7 with manual specification of all pairs of SNPs including the index SNP. Necessary because plink seems not to be calculating the complex LD structure for this region correctly when not specifying all pairs.

- **SF2c_RA_plots_triplets_with_snp_comparison.R**: Regional association plots of combined, female and male analysis group with highlighted top SNPs for each group. Used for visualisation of loci with COLOC support for different variants in males and females.

- **SF3_forest_plots.R**: Forest-plots for index variants are created for all phenotypes where the SNP is genome-wide significant.

- **SF4_power_calculation.R**: Power calculation depending on stick-wise detection rates

## 5. Supplementary Folders
- **additional_data**: Resources used by analysis scripts obtained from manual creation or publicly available resources.

- **data**: Dummy files of formatted meta-data and annotation results to show naming schemes and file structure.

- **helper_scripts**: Functions from in-house analysis pipeline that are called by analysis scripts.

## 6. Disclaimer
THE SOFTWARE AND THE DATA IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

