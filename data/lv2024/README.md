# Sequence data from [Lv et al (2024)](https://academic.oup.com/ve/article/10/1/veae020/7619252)

[Supplementary Tables 1-3.xlsx](Supplementary Tables 1-3.xlsx) is the supplementary tables from [Lv et al (2024)](https://academic.oup.com/ve/article/10/1/veae020/7619252).

[jointWHO_sequences.csv](jointWHO_sequences.csv) indicates sequences from Tables 6 and 7 of the the [joint WHO-China report](https://www.who.int/publications/i/item/who-convened-global-study-of-origins-of-sars-cov-2-china-part), which are said by that report to be the only sequences from patients with symptom onset dates before Dec-31-2019. The sequences are annotated by whether the report indicates in Table 7 whether or not they are from a patient linked to the Huanan Seafood Market.

The Jupyter notebook [extract_accessions.ipynb](extract_accessions.ipynb) extracts accessions from Supplementary Table 3 in [Supplementary Tables 1-3.xlsx](Supplementary Tables 1-3.xlsx) and adds the annotations in [jointWHO_sequences.csv](jointWHO_sequences.csv). It creates the following files:
  - [all_accessions.csv](all_accessions.csv): all accessions in Supplementary Table 3.
  - [genbank_accessions.csv](genbank_accessions.csv): Genbank accessions in Supplementary Table 3.
  - [gisaid_accessions.csv](gisaid_accessions.csv): GISAID accessions in Supplementary Table 3.

The sequences corresponding to these accessions were then manually downloaded to:
  - [genbank_sequences.fa](genbank_sequences.fa): the Genbank sequences
  - [gisaid_sequences.fa](gisaid_sequences.fa): The GISAID sequences. **Due to GISAID data sharing restrictions, this sequence file is not actually tracked in the repository. You will need to manually generate it yourself from the list of accessions.**
