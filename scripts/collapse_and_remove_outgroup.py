"""Collapse zero length branches."""


import Bio.Phylo


tree = Bio.Phylo.read(snakemake.input.tree, "newick")

to_prune = []
for node in tree.find_clades(order='postorder'):
    if node.branch_length <= snakemake.params.blmin and node != tree.root and (node.is_terminal() == False):
        to_prune.append(node)
print(f"Collapsing {len(to_prune)} nodes")
for node in to_prune:
    tree.collapse(node)

if snakemake.params.has_outgroup:
    print(f"Removing outgroup")
    assert 1 == sum(clade.name == "root" for clade in tree.get_terminals())
    n_init = len(tree.get_terminals())
    # outgroup should be a single clade off the base clade of tree 
    assert 1 == sum(clade.name == "root" for clade in tree.clade.clades), tree.clade.clades
    tree.root_with_outgroup("root")
    assert tree.root.clades[1].name == "root"
    tree.root = tree.root.clades[0]
    n_final = len(tree.get_terminals())
    assert n_init == n_final + 1
    
assert 0 == sum(clade.name == "root" for clade in tree.get_terminals())

Bio.Phylo.write(tree, snakemake.output.tree, "newick")
