"""The Chesapeake Bay carbon flow network, for the energy channel test in test/compartments.jl.

Reads the trophic models of the R package enaR (Borrett & Lau 2014, GPL-3) from its data file at a
pinned commit, without R, and writes the Chesapeake Bay network (Baird & Ulanowicz 1989) as the
flows between its 36 compartments (mg C m⁻² yr⁻¹, from compartment `from` to compartment `to`) and
its nodes, with whether each is living and its external input. Zhao et al. (2018) use this network
as their worked example. Run from the repository root:

    uvx --with rdata python -I test/fixtures/make_chesapeake_bay.py
"""

import csv
import pathlib
import tempfile
import urllib.request

import rdata

COMMIT = "46cf03e07620a34b34a8336d08f7b17b2b204814"
URL = f"https://github.com/SEELab/enaR/raw/{COMMIT}/data/troModels.rda"
MODEL = "Chesapeake Bay"


def main():
    out = pathlib.Path(__file__).parent
    with tempfile.TemporaryDirectory() as tmp:
        path = pathlib.Path(tmp) / "troModels.rda"
        urllib.request.urlretrieve(URL, path)
        models = rdata.conversion.convert(rdata.parser.parse_file(path), default_encoding="utf8")["troModels"]
    name = next(k for k in models if str(k).strip() == MODEL)
    network = models[name]
    network = network.value if hasattr(network, "value") else network
    with open(out / "chesapeake_bay_nodes.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["index", "name", "living", "input"])
        for i, v in enumerate(network["val"], 1):
            writer.writerow([i, str(v["vertex.names"][0]), bool(v["living"][0]), repr(float(v["input"][0]))])
    with open(out / "chesapeake_bay_flows.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["from", "to", "flow"])
        for e in network["mel"]:
            writer.writerow([int(e["outl"][0]), int(e["inl"][0]), repr(float(e["atl"]["flow"][0]))])


if __name__ == "__main__":
    main()
