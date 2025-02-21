# Sequence data from [Crits-Christoph et al (2024)](https://www.cell.com/cell/fulltext/S0092-8674(24)00901-2)

[Table_S1.xlsx](Table_S1.xlsx) is the first supplementary table from [Crits-Christoph et al (2024)](https://www.cell.com/cell/fulltext/S0092-8674(24)00901-2).

[jointWHO_market_annotations.csv](jointWHO_market_annotations.csv) provides additional annotations based on the joint WHO-China report and information about whether sequences are from environmental samples at the Huanan Seafood Market. Specifically, it indicates sequences from Tables 6 and 7 of the the [joint WHO-China report](https://www.who.int/publications/i/item/who-convened-global-study-of-origins-of-sars-cov-2-china-part), which are said by that report to be the only sequences from patients with symptom onset dates before Dec-31-2019. These sequences are annotated by whether the report indicates in Table 7 whether or not they are from a patient linked to the Huanan Seafood Market. Note that I added NMDC60013002-09 and NMDC60013002-06 for consistency with the annotations of the samples in Table S1 of [Crits-Christoph et al, Cell (2024)](https://www.cell.com/cell/fulltext/S0092-8674(24)00901-2). I also indicated which sequences are from environmental samples from the Huanan Market collected by [Liu et et al (2023)](https://www.nature.com/articles/s41586-023-06043-2) on Jan-1-2020.

Run the Jupyter notebook [extract_accessions.ipynb](extract_accessions.ipynb) to extract the list of accessions from Table S1 and add the additional annotations. It creates the following files:
  - [seq_metadata.csv](seq_metadata.csv): all metadata extracted from the table of accessions.
  - [custom_accessions.csv](custom_accessions.csv): metadata for sequences with no accession (custom built)
  - [gisaid_accessions.csv](gisaid_accessions.csv): metadata for sequences with GISAID accessions
  - [ngdc_accessions.csv](ngdc_accessions.csv): metadata for sequences with NGDC accessions
  - [genbank_accessions.csv](genbank_accessions.csv): metadata for sequences with Genbank accessions

The following files were manually downloaded to contain the sequences with the accessions indicated above:
  - [custom_sequences.fa](custom_sequences.fa): downloaded from [here](https://github.com/sars-cov-2-origins/huanan-market-environment/tree/main/sars2_phylogenetics/HSM_sequences) in the GitHub repo of [Crits-Christoph et al, Cell (2024)](https://www.cell.com/cell/fulltext/S0092-8674(24)00901-2).
  - [gisaid_sequences.fa](gisaid_sequences.fa): note my GISAID query turned up 782 sequences not the 783 GISAID accessions given by [Crits-Christoph et al, Cell (2024)](https://www.cell.com/cell/fulltext/S0092-8674(24)00901-2). The accession `EPI_ISL_1143994` for strain *hCoV-19/Germany/BY-RKI-I-013858/2020* is not found on GISAID. **Due to GISAID data sharing restrictions, this file is not tracked in the GitHub repo.**
  - [ngdc_sequences.fa](ngdc_sequences.fa): sequences from NGDC.
  - [genbank_sequences.fa](genbank_sequences.fa): sequences from Genbank.
