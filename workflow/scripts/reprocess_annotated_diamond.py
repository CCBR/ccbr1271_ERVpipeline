import argparse
import os
import sys
import gzip
import time, random
from Bio import Entrez

Entrez.email = "CCBR_Pipeliner@mail.nih.gov"


def checkfile(f):
    if not os.path.exists(f):
        sys.exit("File %s does not exist!" % (f))
    if not os.access(f, os.R_OK):
        sys.exit("File %s is not readable!" % (f))
    return True


def get_taxid(accid):
    tid = "Unknown"
    time.sleep(random.choice(list(range(10, 20))))
    try:
        handle = Entrez.esummary(db="protein", id=accid, retmode="xml")
        records = Entrez.parse(handle)
        x = 0
        for r in records:
            x += 1
            if x != 1:
                break
            tid = str(int(r["TaxId"]))
    except:
        tid = "Unknown"
    return tid


def main():
    parser = argparse.ArgumentParser(
        description="Find missing annotations in pre-annotated DIAMOND output"
    )
    parser.add_argument(
        "-l",
        dest="taxid2lineage",
        required=True,
        help="taxid to lineage lookup gzip file",
    )
    parser.add_argument(
        "-a",
        dest="annotateddiamondoutput",
        required=True,
        help="annotated DIAMOND output TSV file",
    )
    parser.add_argument(
        "-o", dest="outfile", required=True, help="reannotated output TSV file"
    )
    args = parser.parse_args()

    checkfile(args.taxid2lineage)
    taxid2lineage = dict()
    with gzip.open(args.taxid2lineage) as fin:
        for l in fin:
            l = l.decode("utf-8")
            l = l.strip().split(",")
            tid = l.pop(0)
            taxid2lineage[tid] = ",".join(l)
    taxid2lineage["Unknown"] = "Unknown"
    print("Done reading taxid2lineage info!")

    checkfile(args.annotateddiamondoutput)

    of = open(args.outfile, "w")
    with open(args.annotateddiamondoutput) as fin:
        for l in fin:
            l = l.strip().split("\t")
            # Expected columns:
            # 1       qseqid
            # 2       accessionversion
            # 3       pident
            # 4       length
            # 5       mismatch
            # 6       gapopen
            # 7       qstart
            # 8       qend
            # 9       sstart
            # 10      send
            # 11      evalue
            # 12      bitscore
            # 13      taxid
            # 14      title
            # 15      lineage
            tid = l[12]
            lin = l[14]
            if tid == "Unknown" or lin == "Unknown":
                print("Unknown found in line:")
                print("\t".join(l))
                new_tid = get_taxid(l[1])
                try:
                    new_lineage = taxid2lineage[new_tid]
                except:
                    new_lineage = "Unknown"
                l[12] = new_tid
                l[14] = new_lineage
                print("new_taxid=", new_tid)
                print("new_lineage=", new_lineage)
            of.write("%s\n" % ("\t".join(l)))
    of.close()


if __name__ == "__main__":
    main()
