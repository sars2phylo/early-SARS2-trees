import Bio.SeqIO

import pandas as pd


print(f"Getting sequences to align into {snakemake.output.to_align}")

seqs = []

if snakemake.wildcards.root != "none":
    root = Bio.SeqIO.read(snakemake.input.root, "fasta")
    root.id = "root"
    seqs.append(root)

reference_sequence = Bio.SeqIO.read(snakemake.input.reference_sequence, "fasta")
reference_sequence.id = "reference_sequence"
seqs.append(reference_sequence)

accession_dates = (
    pd.read_csv(snakemake.input.metadata)
    .assign(
        date=lambda x: pd.to_datetime(x["date"], format="%Y-%m-%d", errors="raise")
    )
    .set_index("accession")
    ["date"]
    .to_dict()
)

assert "root" not in accession_dates
assert "reference_sequence" not in accession_dates

start_date, end_date = snakemake.params.date_range
if start_date is not None:
    start_date = pd.to_datetime(start_date, format="%Y-%m-%d", errors="raise")
if end_date is not None:
    end_date = pd.to_datetime(end_date, format="%Y-%m-%d", errors="raise")

print(f"{start_date=}, {end_date=}")

n = nkept = 0
for seq in Bio.SeqIO.parse(snakemake.input.unaligned_seqs, "fasta"):
    n += 1
    date = pd.to_datetime(accession_dates[seq.id], format="%Y-%m-%d", errors="raise")
    if (start_date is not None) and (date < start_date):
        continue
    if (end_date is not None) and (date > end_date):
        continue
    nkept += 1
    seqs.append(seq)
print(f"Kept {nkept} of {n} sequences as being in date range")

Bio.SeqIO.write(seqs, snakemake.output.to_align, "fasta")
