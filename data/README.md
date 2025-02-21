# Input data

- [reference_sequence.fa](reference_sequence.fa): Wuhan-Hu-1 reference sequence, [Genbank NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/1798174254)

- [auspice_config.json](auspice_config.json): Auspice configuration file.

Comparator sequences that represent relatives of inferred ancestors of SARS-CoV-2:

  - [recCA.fa](recCA.fa): the "recombinant common ancestor" of SARS-CoV-2 inferred by [Pekar et al (2022)](https://www.science.org/doi/10.1126/science.abp8337) as taken from [this Zenodo repository](https://zenodo.org/records/6899613), file `recCA.unmasked.fasta` in the [sarbecovirus_nonrecombinant_region_trees.zip](https://zenodo.org/records/6899613/files/sarbecovirus_nonrecombinant_region_trees.zip) file.

  - [RaTG13.fa](RaTG13.fa): RaTG13, [Genbank MN996532.2](https://www.ncbi.nlm.nih.gov/nuccore/MN996532).

Candidate root sequences are in [./candidate_roots/](candidate_roots):

  - [candidate_roots/lineage-A.fa](candidate_roots/lineage-A.fa): Genbank [LR757995.1](https://www.ncbi.nlm.nih.gov/nuccore/LR757995.1), which is the WH-04 strain of SARS-CoV-2. Matches the lineage A root suggested in the first figure of [Crits-Christoph et al (2024)](https://www.sciencedirect.com/science/article/pii/S0092867424009012).

  - [candidate_roots/lineage-B.fa](candidate_roots/lineage-B.fa):[Genbank NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/1798174254), which is the Wuhan-Hu-1 strain of SARS-CoV-2. Matches the lineage B root from the unconstrained rooting of [Pekar et al (2022)](https://www.science.org/doi/full/10.1126/science.abp8337).

  - [candidate_roots/lineage-A-C18060T.fa](candidate_roots/lineage-A-C18060T.fa): Genbank [MN985325.1](https://www.ncbi.nlm.nih.gov/nuccore/MN985325.1) which, is the "WA1" strain of SARS-CoV-2. Matches the "lineageA plus C18060T" root suggested by [Kumar et al (2021)](https://academic.oup.com/mbe/article/38/8/3046/6257226), who call it proCoV2. This sequence also suggested as a possible root by [Bloom (2021)](https://academic.oup.com/mbe/article/38/12/5211/6353034).

  - [candidate_roots/lineage-A-C29095T.fa](candidate_roots/lineage-A-C29095T.fa): Genbank [MN938384.1](https://www.ncbi.nlm.nih.gov/nucleotide/MN938384.1) which, is the "HKU-SZ-002a" strain of SARS-CoV-2. Matches the "lineageA plus C29095T" possible root suggested by [Lv et al (2024)](https://academic.oup.com/ve/article/10/1/veae020/7619252) and [Bloom (2021)](https://academic.oup.com/mbe/article/38/12/5211/6353034).

- Sequence sets (see READMEs within each subdirectory for details):

  - [./lv2024/](lv2024)

  - [./crits-christoph2024/](crits-christoph2024)
