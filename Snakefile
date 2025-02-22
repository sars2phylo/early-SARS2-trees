"""Snakemake file that runs analysis."""


import re


configfile: "config.yaml"


# process some variables from configuration
sequence_sets = config["sequence_sets"]
roots = config["roots"]
date_ranges = config["date_ranges"]
repo = config["repo"]
auspice_configs = config["auspice_configs"]
mask_sites = [
    site
    for sublist in (
        [item] if isinstance(item, int) else range(item[0], item[1] + 1)
        for item in config["mask_sites"]
    )
    for site in sublist
]

# cannot have underscores in these for Nextstrain Community builds
for key in ["sequence_sets", "roots", "date_ranges", "auspice_configs"]:
    assert all("_" not in val for val in config[key]), f"_ in {config[key]=}"
assert "_" not in repo, repo

wildcard_constraints:
    root="|".join(map(re.escape, roots)),
    date_range="|".join(map(re.escape, date_ranges)),
    sequence_set="|".join(map(re.escape, sequence_sets)),
    auspice_config="|".join(map(re.escape, auspice_configs)),


rule all:
    input:
        expand(
            f"auspice/{repo}_" + "{sequence_set}-{date_range}-{root}-{auspice_config}.json",
            sequence_set=sequence_sets,
            root=roots,
            date_range=date_ranges,
            auspice_config=auspice_configs,
        ),


rule prep_sequences_and_metadata:
    """Concatenate the unaligned sequences and get metadata."""
    input:
       fastas=lambda wc: sequence_sets[wc.sequence_set]["fastas"],
       metadata=lambda wc: sequence_sets[wc.sequence_set]["metadata"],
    output:
        unaligned_seqs="results/{sequence_set}/unaligned.fa",
        metadata="results/{sequence_set}/metadata.csv",
    params:
        metadata_req_cols=config["metadata_req_cols"],
        length_range=lambda wc: sequence_sets[wc.sequence_set]["length_range"],
        drop_accessions=lambda wc: sequence_sets[wc.sequence_set]["drop_accessions"],
    log:
        notebook="results/{sequence_set}/prep_sequences_and_metadata.ipynb",
    notebook:
        "notebooks/prep_sequences_and_metadata.py.ipynb"


rule get_seqs_to_align:
    """Get the set of sequences to align for specific root / date range."""
    input:
        unpack(
            lambda wc: (
                {"root": roots[wc.root]}
                if wc.root != "none"
                else {}
            )
        ),
        reference_sequence=config["reference_sequence"],
        unaligned_seqs="results/{sequence_set}/unaligned.fa",
        metadata="results/{sequence_set}/metadata.csv",
    output:
        to_align="results/{sequence_set}/to_align_root-{root}_{date_range}.fa",
    params:
        date_range=lambda wc: date_ranges[wc.date_range],
    script:
        "scripts/get_seqs_to_align.py"


rule align:
    """Align the sequences."""
    input:
        to_align="results/{sequence_set}/to_align_root-{root}_{date_range}.fa",
    output:
        alignment_unmasked="results/{sequence_set}/root-{root}_{date_range}_unmasked.fa",
    threads: 4
    shell:
        """
        augur align \
            -s {input.to_align} \
            -o {output.alignment_unmasked} \
            --reference-name reference_sequence \
            --fill-gaps \
            --nthreads {threads}
        """


rule mask_alignment:
    """Mask the alignment and drop the reference sequence."""
    input:
        alignment_unmasked="results/{sequence_set}/root-{root}_{date_range}_unmasked.fa",
    output:
        alignment_masked="results/{sequence_set}/root-{root}_{date_range}_masked.fa",
    params:
        mask_sites=mask_sites,
        max_unmasked_n=config["max_unmasked_n"],
    script:
        "scripts/mask_alignment.py"


rule build_tree:
    """Build the phylogenetic tree."""
    input:
        alignment_masked="results/{sequence_set}/root-{root}_{date_range}_masked.fa",
    output:
        exclude_sites="results/{sequence_set}/root-{root}_{date_range}_exclude_sites.txt",
        tree="results/{sequence_set}/root-{root}_{date_range}_initial_tree.newick",
    params:
        exclude_sites="\n".join(map(str, mask_sites)),
        blmin=0.0000001,  # minimum branch length for zero-length branches
        outgroup_args=lambda wc: "" if wc.root == "none" else r"\-o root",
    shell:
        r"""
        echo "{params.exclude_sites}" > {output.exclude_sites}
        augur tree \
            --alignment {input.alignment_masked} \
            --exclude-sites {output.exclude_sites} \
            --method iqtree \
            --substitution-model GTR \
            --tree-builder-args "\-ninit 100 \-n 5 \-me 0.01 \-seed 1 \-blmin {params.blmin} {params.outgroup_args}" \
            --override-default-args \
            --output {output.tree}
        """


rule collapse_and_remove_outgroup:
    """Collapse zero-length branches and remove rooting outgroup if present."""
    input:
        tree="results/{sequence_set}/root-{root}_{date_range}_initial_tree.newick",
    output:
        tree="results/{sequence_set}/root-{root}_{date_range}_collapsed.newick",
    params:
        blmin=rules.build_tree.params.blmin,
        has_outgroup=lambda wc: wc.root != "none",
    script:
        "scripts/collapse_and_remove_outgroup.py"


rule refine_tree:
    """Run the ``augur refine`` command."""
    input:
        tree="results/{sequence_set}/root-{root}_{date_range}_collapsed.newick",
        metadata="results/{sequence_set}/metadata.csv",  # contains more accessions than those in tree
        alignment="results/{sequence_set}/root-{root}_{date_range}_masked.fa",  # still has root and reference
    output:
        tree="results/{sequence_set}/root-{root}_{date_range}_refined.newick",
        branch_lengths="results/{sequence_set}/root-{root}_{date_range}_branch_lengths.json",
    params:
        rooting_params=lambda wc: (
            ""
            if wc.root != "none"
            else r"--keep-root"
        ),
    shell:
        """
        augur refine \
            -a {input.alignment} \
            -t {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns accession \
            --output-tree {output.tree} \
            --output-node-data {output.branch_lengths} \
            {params.rooting_params} \
            --verbosity 2
        """


rule ancestral:
    """Mutations and ancestral sequences."""
    input:
        tree="results/{sequence_set}/root-{root}_{date_range}_refined.newick",
        alignment="results/{sequence_set}/root-{root}_{date_range}_masked.fa",  # still has root and reference
    output:
        nt_muts="results/{sequence_set}/root-{root}_{date_range}_nt-muts.json",
        alignment_inferred="results/{sequence_set}/root-{root}_{date_range}_alignment_inferred.fa",
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.nt_muts} \
            --keep-ambiguous \
            --seed 1
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-sequences {output.alignment_inferred} \
            --infer-ambiguous \
            --seed 1
        """


rule align_for_mutations_from:
    """Align reference to sequences that we use to calculate mutations from."""
    input:
        reference_fasta=config["reference_sequence"],
        mutations_from_fasta=lambda wc: (
            config["metadata_from_alignment"]["mutations_from"][wc.mutations_from]
        ),
    output:
        alignment="results/mutations_from_alignments/{mutations_from}.fa",
    shell:
        """
        augur align \
            -s {input.mutations_from_fasta} \
            -o {output.alignment} \
            --reference-sequence {input.reference_fasta} \
            --remove-reference \
            --fill-gaps
        """


rule metadata_from_alignment:
    """Add to metadata additional fields from alignment w inferred ambiguous identities."""
    input:
        alignment_inferred="results/{sequence_set}/root-{root}_{date_range}_alignment_inferred.fa",
        metadata="results/{sequence_set}/metadata.csv",
        mutations_from_alignments=expand(
            "results/mutations_from_alignments/{mutations_from}.fa",
            mutations_from=config["metadata_from_alignment"]["mutations_from"],
        ),
    output:
        metadata="results/{sequence_set}/root-{root}_{date_range}_metadata.csv",
    params:
        **config["metadata_from_alignment"],
        mask_sites=mask_sites,
    script:
        "scripts/metadata_from_alignment.py"


rule export_tree:
    """Export the tree to JSON for auspice."""
    input:
        tree="results/{sequence_set}/root-{root}_{date_range}_refined.newick",
        nt_muts="results/{sequence_set}/root-{root}_{date_range}_nt-muts.json",
        branch_lengths="results/{sequence_set}/root-{root}_{date_range}_branch_lengths.json",
        metadata="results/{sequence_set}/root-{root}_{date_range}_metadata.csv",
        auspice_config=lambda wc: auspice_configs[wc.auspice_config],
    output:
        tree="results/{sequence_set}/root-{root}_{date_range}_{auspice_config}_no_mutation_branch_labels.json",
    params:
        title=lambda wc: (
            "Early SARS-CoV-2 tree using sequences from "
            + sequence_sets[wc.sequence_set]["name"]
            + f" up to {wc.date_range} rooted on {wc.root}"
        ),
        metadata_columns=" ".join(
            c for c in config["metadata_req_cols"] if c not in {"accession"}
        ),
        color_by_metadata=" ".join(
            ["date", "addtl_annotations"]
            + list(config["metadata_from_alignment"]["haplotypes"])
            + list(config["metadata_from_alignment"]["mutations_from"])
        ),
        addtl_args=" ".join(
            f"--{key} '{val}'" for (key, val) in config["augur_export_args"].items()
        )
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --output {output.tree} \
            --auspice-config {input.auspice_config} \
            --node-data {input.branch_lengths} {input.nt_muts} \
            --metadata {input.metadata} \
            --metadata-id-columns accession \
            --color-by-metadata {params.color_by_metadata} \
            --metadata-columns {params.metadata_columns} \
            --include-root-sequence-inline \
            --panels tree \
            --title "{params.title}" \
            --no-minify-json \
            {params.addtl_args}
        """


rule add_mutation_branch_labels:
    """Add branch labels as mutations."""
    input:
        tree="results/{sequence_set}/root-{root}_{date_range}_no_mutation_branch_labels.json",
    output:
        tree="results/{sequence_set}/root-{root}_{date_range}_{auspice_config}.json",
    script:
        "scripts/add_mutation_branch_labels.py"


rule final_jsons:
    """Copy the final JSONs to a separate subdirectory."""
    input:
        json="results/{sequence_set}/root-{root}_{date_range}_{auspice_config}.json",
    output:
        json=f"auspice/{repo}_" + "{sequence_set}-{date_range}-{root}-{auspice_config}.json",
    shell:
        "cp {input.json} {output.json}"
