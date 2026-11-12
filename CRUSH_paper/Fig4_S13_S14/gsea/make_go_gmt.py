#!/usr/bin/env python3
import gzip
from pathlib import Path

HERE = Path(__file__).resolve().parent
RAW = HERE / "ref" / "go_raw"


def ontology():
    terms, cur = {}, None
    for line in open(RAW / "go-basic.obo"):
        line = line.rstrip("\n")
        if line == "[Term]":
            cur = {"parents": set()}
        elif line.startswith("[") or not line:
            if cur and "id" in cur and cur.get("ns") == "biological_process" and not cur.get("obsolete"):
                terms[cur["id"]] = cur
            cur = None if line.startswith("[") else cur
        elif cur is not None:
            k, _, v = line.partition(": ")
            if k == "id":
                cur["id"] = v
            elif k == "name":
                cur["name"] = v
            elif k == "namespace":
                cur["ns"] = v
            elif k == "is_obsolete" and v == "true":
                cur["obsolete"] = True
            elif k == "is_a":
                cur["parents"].add(v.split()[0])
            elif k == "relationship" and v.startswith("part_of "):
                cur["parents"].add(v.split()[1])
    return terms


def main():
    terms = ontology()
    anc = {}

    def up(t):
        if t not in anc:
            anc[t] = {t}.union(*(up(p) for p in terms[t]["parents"] if p in terms))
        return anc[t]

    genes = {t: set() for t in terms}
    with gzip.open(RAW / "goa_human.gaf.gz", "rt") as f:
        for line in f:
            if line.startswith("!"):
                continue
            x = line.split("\t")
            if x[8] != "P" or "NOT" in x[3] or x[6] == "ND" or x[4] not in terms:
                continue
            for a in up(x[4]):
                genes[a].add(x[2])
    with open(HERE / "ref" / "GO_BP_full.gmt", "w") as out:
        for t, g in sorted(genes.items()):
            if g:
                out.write("\t".join([f"{terms[t]['name']} ({t})", ""] + sorted(g)) + "\n")
    print(len(terms), "BP terms;", sum(1 for g in genes.values() if g), "with genes;",
          "hair follicle development:", len(genes.get("GO:0001942", ())), "genes")


if __name__ == "__main__":
    main()
