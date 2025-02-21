"""Based on this: https://discussion.nextstrain.org/t/how-to-add-branch-labels-from-node-data/1339/2"""

import json


with open(snakemake.input.tree, "r") as f:
    auspice_json = json.load(f)


def attach_labels(n):
    assert "branch_attrs" in n
    if "labels" not in n["branch_attrs"]:
        n["branch_attrs"]["labels"] = {}
    mutations = n["branch_attrs"]["mutations"]
    if mutations:
        n["branch_attrs"]["labels"]["mutations"] = " ".join(
            m for m in n["branch_attrs"]["mutations"]["nuc"] if not "N" in m
        )
    if "children" in n:
        for c in n["children"]:
            attach_labels(c)

attach_labels(auspice_json["tree"])

if "display_defaults" not in auspice_json["meta"]:
    auspice_json["meta"]["display_defaults"] = {}
auspice_json["meta"]["display_defaults"]["branch_label"] = "mutations"

with open(snakemake.output.tree, "w") as f:
    json.dump(auspice_json, f, indent=2)
