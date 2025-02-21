import Bio.SeqIO

import pandas as pd


mask_sites = set(snakemake.params.mask_sites)
max_unmasked_n = snakemake.params.max_unmasked_n

print(f"Reading unmasked alignment from {snakemake.input.alignment_unmasked}")
seqlen = None
dropped_reference_sequence = False
records = []
for seq in Bio.SeqIO.parse(snakemake.input.alignment_unmasked, "fasta"):
    if seq.id == "reference_sequence":
        assert not dropped_reference_sequence
        dropped_reference_sequence = True
        continue
    if seqlen is not None:
        assert len(seq) == seqlen, "all sequences not of same length"
    seqlen = len(seq)
    accession = seq.id
    seq = list(str(seq.seq).upper())
    for site in mask_sites:
        assert 1 <= site <= seqlen
        seq[site - 1] = "N"
    seq = "".join(seq)
    unmasked_n = sum(seq[i] == "N" for i in range(seqlen) if i + 1 not in mask_sites)
    records.append((accession, seq, unmasked_n))

assert dropped_reference_sequence

df = pd.DataFrame(records, columns=["accession", "sequence", "unmasked_n"])
print(f"Read {len(df)} sequences after dropping 'reference_sequence'.")
if any(df["unmasked_n"] > max_unmasked_n):
    raise ValueError(
        f"Some sequences have more than {max_unmasked_n=} unmasked N's:\n"
        + str(
            df
            .sort_values("unmasked_n", ascending=False)
            [["accession", "unmasked_n"]]
            .reset_index()
        )
    )

print(f"Writing sequences to {snakemake.output.alignment_masked}")
with open(snakemake.output.alignment_masked, "w") as f:
    f.write(
        "\n".join(f">{accession}\n{sequence}\n" for (accession, sequence, _) in records)
    )
