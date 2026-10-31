"""Independent motif censuses of the bundled food webs, for test/motifs.jl.

igraph counts the connected three- and four-species motifs with its own ESU implementation, and
NetworkX's triad census gives the three-species ones a second time. Both are mapped to the
mfinder IDs used by FoodWebs.jl. Run from the repository root:

    uvx --with igraph --with networkx python -I test/fixtures/make_motif_census.py
"""

import csv
import itertools
import pathlib

import igraph
import networkx

ROOT = pathlib.Path("data")
COLLECTIONS = ["cheddar", "williams_martinez_2000"]


def mfinder_id(edges, k):
    """Smallest adjacency code over relabelings, bit (i*k + j) set when i eats j."""
    best = None
    for p in itertools.permutations(range(k)):
        code = sum(1 << (p[i] * k + p[j]) for i, j in edges)
        best = code if best is None else min(best, code)
    return best


def read_web(directory):
    with open(directory / "nodes.csv", newline="") as f:
        species = [row[0] for row in list(csv.reader(f))[1:]]
    index = {s: n for n, s in enumerate(species)}
    with open(directory / "trophic.links.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    # predator -> prey, without cannibalism
    edges = {(index[r["consumer"]], index[r["resource"]]) for r in rows if r["consumer"] != r["resource"]}
    return len(species), sorted(edges)


def igraph_census(n, edges, k):
    g = igraph.Graph(n=n, edges=edges, directed=True)
    counts = g.motifs_randesu(size=k)
    result = {}
    for cls, c in enumerate(counts):
        if c != c:  # NaN: a disconnected class
            continue
        motif = igraph.Graph.Isoclass(k, cls, directed=True)
        result[mfinder_id(motif.get_edgelist(), k)] = int(c)
    return result


def networkx_census(n, edges):
    g = networkx.DiGraph()
    g.add_nodes_from(range(n))
    g.add_edges_from(edges)
    census = networkx.triadic_census(g)
    # the triad type of each connected three-species motif
    result = {}
    for name, count in census.items():
        for edges3 in itertools.chain.from_iterable(
            itertools.combinations([(i, j) for i in range(3) for j in range(3) if i != j], m) for m in range(2, 7)
        ):
            t = networkx.DiGraph()
            t.add_nodes_from(range(3))
            t.add_edges_from(edges3)
            if networkx.is_weakly_connected(t) and networkx.triad_type(t) == name:
                result[mfinder_id(edges3, 3)] = count
                break
    return result


def main():
    out = pathlib.Path(__file__).parent / "motif_census.csv"
    with open(out, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["web", "size", "motif", "igraph", "networkx"])
        for collection in COLLECTIONS:
            for directory in sorted((ROOT / collection).iterdir()):
                if not (directory / "nodes.csv").exists():
                    continue
                n, edges = read_web(directory)
                web = f"{collection}/{directory.name}"
                triads = networkx_census(n, edges)
                for k in (3, 4):
                    for motif, count in sorted(igraph_census(n, edges, k).items()):
                        nx = triads.get(motif, 0) if k == 3 else ""
                        writer.writerow([web, k, motif, count, nx])
    print(f"wrote {out}")


main()
