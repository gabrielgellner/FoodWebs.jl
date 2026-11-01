"""Independent modularity values of fixed partitions of the bundled food webs, for test/compartments.jl.

NetworkX and igraph each compute the directed modularity (Leicht & Newman) of the web and the
undirected modularity (Newman & Girvan) of the web with every link kept as an undirected edge, so
that two species that eat each other are joined by two edges. Links run from prey to predator and
cannibalistic links are left out. Weighted rows give the link from prey i to predator j the weight
1 + (3i + 7j) mod 5, with 0-based species indices. Run from the repository root:

    uvx --with igraph --with networkx python -I test/fixtures/make_modularity.py
"""

import csv
import pathlib

import igraph
import networkx

ROOT = pathlib.Path("data")
COLLECTIONS = ["cheddar", "williams_martinez_2000"]


def read_web(directory):
    with open(directory / "nodes.csv", newline="") as f:
        species = [row[0] for row in list(csv.reader(f))[1:]]
    index = {s: n for n, s in enumerate(species)}
    with open(directory / "trophic.links.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    # prey -> predator, without cannibalism
    edges = {(index[r["resource"]], index[r["consumer"]]) for r in rows if r["consumer"] != r["resource"]}
    return len(species), sorted(edges)


def partitions(n):
    return {
        "thirds": [3 * i // n for i in range(n)],
        "mod3": [i % 3 for i in range(n)],
    }


def weight(i, j):
    return 1 + (3 * i + 7 * j) % 5


def networkx_q(n, edges, groups, directed, weighted, gamma):
    g = networkx.DiGraph() if directed else networkx.MultiGraph()
    g.add_nodes_from(range(n))
    for i, j in edges:
        g.add_edge(i, j, weight=weight(i, j) if weighted else 1)
    communities = [{v for v in range(n) if groups[v] == c} for c in set(groups)]
    return networkx.community.modularity(g, communities, weight="weight", resolution=gamma)


def igraph_q(n, edges, groups, directed, weighted, gamma):
    g = igraph.Graph(n=n, edges=edges, directed=directed)
    weights = [weight(i, j) for i, j in edges] if weighted else None
    return g.modularity(groups, weights=weights, resolution=gamma, directed=directed)


def main():
    out = pathlib.Path(__file__).parent / "modularity.csv"
    with open(out, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["web", "partition", "directed", "weighted", "gamma", "networkx", "igraph"])
        for collection in COLLECTIONS:
            for directory in sorted((ROOT / collection).iterdir()):
                if not (directory / "nodes.csv").exists():
                    continue
                n, edges = read_web(directory)
                for name, groups in partitions(n).items():
                    for directed in (True, False):
                        for weighted in (False, True):
                            for gamma in (1.0, 0.5):
                                writer.writerow([
                                    f"{collection}/{directory.name}", name, directed, weighted, gamma,
                                    repr(networkx_q(n, edges, groups, directed, weighted, gamma)),
                                    repr(igraph_q(n, edges, groups, directed, weighted, gamma)),
                                ])


if __name__ == "__main__":
    main()
