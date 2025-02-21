import Bio.SeqIO

import pandas as pd


metadata = pd.read_csv(snakemake.input.metadata)

assert not metadata["accession"].str.startswith("NODE_").any()

haplotypes = snakemake.params.haplotypes
mutations_from = snakemake.params.mutations_from

if haplotypes:
    haplotype_records = []
    for seq in Bio.SeqIO.parse(snakemake.input.alignment_inferred, "fasta"):
        if seq.id.startswith("NODE_"):
            continue  # internal node
        record = [seq.id]
        seq = str(seq.seq).upper()
        for haplotype_name, haplotype_d in haplotypes.items():
            sites = haplotype_d["sites"]
            haplotype_values = haplotype_d["values"]
            hap = "".join(seq[site - 1] for site in sites)
            assert hap in haplotype_values, f"{hap=}, {accession=}, {haplotype_values=}"
            record.append(haplotype_values[hap])
        haplotype_records.append(record)
    haplotype_metadata = pd.DataFrame(
        haplotype_records, columns=["accession"] + list(haplotypes)
    )
    if not set(haplotype_metadata["accession"]).issubset(metadata["accession"]):
        raise ValueError(f"{(set(haplotype_metadata['accession']) - set(metadata['accession']))=}")
    metadata = metadata.merge(haplotype_metadata, on="accession", validate="one_to_one")

if mutations_from:
    mutations_from_records = []
    mask_sites = set(snakemake.params.mask_sites)
    assert len(mutations_from) == len(snakemake.input.mutations_from_alignments)
    for name, fasta in zip(mutations_from, snakemake.input.mutations_from_alignments):
        comparator = str(Bio.SeqIO.read(fasta, "fasta").seq).upper()
    for seq in Bio.SeqIO.parse(snakemake.input.alignment_inferred, "fasta"):
        if seq.id.startswith("NODE_"):
            continue  # internal node
        record = [seq.id]
        seq = str(seq.seq).upper()
        for name, fasta in zip(mutations_from, snakemake.input.mutations_from_alignments):
            comparator = str(Bio.SeqIO.read(fasta, "fasta").seq).upper()
            assert len(seq) == len(comparator)
            nmuts = sum(
                x != y
                for (i, (x, y)) in enumerate(zip(seq, comparator))
                if (x != "N") and (y != "N") and (i + 1 not in mask_sites)
            )
            record.append(nmuts)
        mutations_from_records.append(record)
    mutations_from_metadata = pd.DataFrame(
        mutations_from_records, columns=["accession"] + list(mutations_from)
    )
    if not set(mutations_from_metadata["accession"]).issubset(metadata["accession"]):
        raise ValueError(f"{(set(mutations_from_metadata['accession']) - set(metadata['accession']))=}")
    metadata = metadata.merge(mutations_from_metadata, on="accession", validate="one_to_one")

metadata.to_csv(snakemake.output.metadata, index=False)
