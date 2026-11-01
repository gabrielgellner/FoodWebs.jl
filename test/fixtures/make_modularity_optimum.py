"""Maximum modularity of the small bundled food webs, for test/compartments.jl.

igraph's cluster_optimal finds the partition of maximum modularity exactly, by integer
programming, for the directed modularity of the web (Leicht & Newman) and the undirected modularity
of the web with every link kept as an undirected edge, so that two species that eat each other
are joined by two edges (summed into one edge of weight 2). Links run from prey to predator and
cannibalistic links are left out. Only the webs with at most 110 links are solved; denser webs
take hours. Run from the repository root:

    uvx --with igraph python -I test/fixtures/make_modularity_optimum.py
"""

import csv
import pathlib

import igraph

ROOT = pathlib.Path("data")
COLLECTIONS = ["cheddar", "williams_martinez_2000"]
MAX_LINKS = 110


def read_web(directory):
    with open(directory / "nodes.csv", newline="") as f:
        species = [row[0] for row in list(csv.reader(f))[1:]]
    index = {s: n for n, s in enumerate(species)}
    with open(directory / "trophic.links.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    # prey -> predator, without cannibalism
    edges = {(index[r["resource"]], index[r["consumer"]]) for r in rows if r["consumer"] != r["resource"]}
    return len(species), sorted(edges)


def optimum(n, edges, directed):
    g = igraph.Graph(n=n, edges=edges, directed=directed)
    g.es["weight"] = 1.0
    if not directed:
        g.simplify(multiple=True, loops=True, combine_edges="sum")
    partition = g.community_optimal_modularity(weights="weight")
    return g.modularity(partition.membership, weights="weight", directed=directed)


def main():
    out = pathlib.Path(__file__).parent / "modularity_optimum.csv"
    with open(out, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["web", "directed", "optimum"])
        for collection in COLLECTIONS:
            for directory in sorted((ROOT / collection).iterdir()):
                if not (directory / "nodes.csv").exists():
                    continue
                n, edges = read_web(directory)
                if len(edges) > MAX_LINKS:
                    continue
                for directed in (True, False):
                    writer.writerow([f"{collection}/{directory.name}", directed, repr(optimum(n, edges, directed))])
                    f.flush()


if __name__ == "__main__":
    main()
