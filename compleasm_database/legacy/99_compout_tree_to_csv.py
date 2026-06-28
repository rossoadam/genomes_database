#!/usr/bin/env python3

import argparse
import csv
import re

import matplotlib.pyplot as plt
from Bio import Phylo


COMP_RE = re.compile(
    r"^>Node\s+(\d+),\s+son of Node\s+(\d+)\s+by branch of length\s+([0-9.eE+-]+)"
)


def parse_compout(compout_path, root_name="ancestral"):
    """
    Parse an nhPhyML .compout file.

    The self-line like:
        >Node 43, son of Node 43 by branch of length ...

    is treated as a root sentinel only and is not kept as a numbered node
    in the final CSV. Instead, anything whose parent is that root sentinel
    is assigned parent = root_name (default: 'ancestral').

    Returns:
        mapping: dict[str, str]
            Internal-node child -> parent relationships, e.g.
            {"node_44": "ancestral", "node_45": "node_44", ...}

        root_id: int
            The nhPhyML root sentinel ID, e.g. 43
    """
    rows = []

    with open(compout_path, "r", encoding="utf-8") as handle:
        for line in handle:
            m = COMP_RE.match(line.strip())
            if m:
                child = int(m.group(1))
                parent = int(m.group(2))
                rows.append((child, parent))

    if not rows:
        raise ValueError("No nhPhyML node lines found in the .compout file.")

    root_id = None
    for child, parent in rows:
        if child == parent:
            root_id = child
            break

    if root_id is None:
        raise ValueError(
            "No root self-line found in .compout "
            "(expected something like '>Node X, son of Node X ...')."
        )

    mapping = {}
    for child, parent in rows:
        if child == parent:
            continue

        child_name = f"node_{child}"
        parent_name = root_name if parent == root_id else f"node_{parent}"
        mapping[child_name] = parent_name

    return mapping, root_id


def assign_internal_ids_from_compout(tree, root_id, root_name="ancestral"):
    """
    Label internal nodes in the tree to match the intended CSV convention:

      - root is named 'ancestral'
      - non-root internal nodes are named:
            node_(root_id + 1), node_(root_id + 2), ...

    Assignment is done in postorder for all non-root internal nodes so that
    deeper nodes receive smaller numbers.

    Example for 43 taxa:
      root_id from compout = 43
      tree root             -> ancestral
      first non-root node   -> node_44
      next                  -> node_45
      ...
    """
    tree.root.name = root_name

    next_id = root_id + 1

    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            continue
        if clade == tree.root:
            continue

        clade.name = f"node_{next_id}"
        next_id += 1


def build_full_mapping(tree, root_name="ancestral"):
    """
    Build child -> parent relationships from the labeled tree.

    Assumes:
      - leaves already have names
      - internal nodes have already been labeled
      - tree.root.name is root_name

    Returns:
        mapping: dict[str, str]
    """
    mapping = {}

    def visit(parent):
        parent_name = parent.name if parent.name else root_name

        for child in parent.clades:
            if not child.name:
                raise ValueError("Encountered unnamed node after labeling step.")
            mapping[child.name] = parent_name
            visit(child)

    visit(tree.root)
    return mapping


def validate_node_counts(tree, comp_internal_map):
    """
    Sanity-check that the number of non-root internal clades in the tree matches
    the number of numbered internal nodes described by the compout mapping.
    """
    non_root_internal_clades = [
        c for c in tree.find_clades(order="postorder")
        if (not c.is_terminal()) and (c != tree.root)
    ]

    comp_nodes = set(comp_internal_map.keys())

    if len(non_root_internal_clades) != len(comp_nodes):
        raise ValueError(
            "Mismatch between tree and compout:\n"
            f"  non-root internal clades in tree: {len(non_root_internal_clades)}\n"
            f"  internal nodes in compout mapping: {len(comp_nodes)}"
        )


def write_two_row_csv(mapping, output_csv):
    """
    Write mapping as:
      row 1 = child nodes
      row 2 = parent nodes
    """
    children = list(mapping.keys())
    parents = [mapping[child] for child in children]

    with open(output_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(children)
        writer.writerow(parents)


def plot_tree_pdf(tree, output_pdf, fig_width=18, fig_height=24, dpi=300):
    """
    Plot the tree left-to-right and save it as a PDF.
    """
    fig = plt.figure(figsize=(fig_width, fig_height), dpi=dpi)
    ax = fig.add_subplot(1, 1, 1)

    def label_func(clade):
        return clade.name if clade.name else ""

    plt.rcParams["font.size"] = 7

    Phylo.draw(
        tree,
        label_func=label_func,
        do_show=False,
        axes=ax,
    )

    ax.set_title("Tree with labeled internal nodes")
    plt.tight_layout()
    plt.savefig(output_pdf, format="pdf", bbox_inches="tight")
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build a two-row child->parent CSV by combining an nhPhyML .compout "
            "file (internal-node relationships) with a rooted Newick tree "
            "(leaf relationships). The root is labeled 'ancestral'."
        )
    )
    parser.add_argument(
        "input_tree",
        help="Path to rooted Newick tree file"
    )
    parser.add_argument(
        "input_compout",
        help="Path to nhPhyML .compout file"
    )
    parser.add_argument(
        "output_csv",
        help="Path to output CSV"
    )
    parser.add_argument(
        "--output-pdf",
        default=None,
        help="Optional path to save a labeled tree PDF"
    )
    parser.add_argument(
        "--root-name",
        default="ancestral",
        help="Name to use for the root node in the CSV and PDF (default: ancestral)"
    )
    parser.add_argument(
        "--fig-width",
        type=float,
        default=18,
        help="Figure width for PDF output (default: 18)"
    )
    parser.add_argument(
        "--fig-height",
        type=float,
        default=24,
        help="Figure height for PDF output (default: 24)"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    comp_internal_map, root_id = parse_compout(
        args.input_compout,
        root_name=args.root_name
    )

    tree = Phylo.read(args.input_tree, "newick")

    assign_internal_ids_from_compout(
        tree,
        root_id=root_id,
        root_name=args.root_name
    )

    validate_node_counts(tree, comp_internal_map)

    full_mapping = build_full_mapping(tree, root_name=args.root_name)

    write_two_row_csv(full_mapping, args.output_csv)

    if args.output_pdf:
        plot_tree_pdf(
            tree,
            args.output_pdf,
            fig_width=args.fig_width,
            fig_height=args.fig_height
        )


if __name__ == "__main__":
    main()
