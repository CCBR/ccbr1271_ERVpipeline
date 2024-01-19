import pickle
import argparse
import os
import sys
import gzip


def checkfile(f):
    if not os.path.exists(f):
        sys.exit("File %s does not exist!" % (f))
    if not os.access(f, os.R_OK):
        sys.exit("File %s is not readable!" % (f))
    return True


def main():
    """
    Reads in 3 different lookup tables and adds annotations to DIAMOND tab-delimted output.
    """

    parser = argparse.ArgumentParser(
        description="Annotate DIAMOND tab-delimited output with taxid and lineage"
    )
    parser.add_argument(
        "-d", dest="diamondoutput", required=True, help="DIAMOND output TSV file"
    )
    parser.add_argument(
        "-p", dest="accid2taxid", required=True, help="accession_id to tax_id pkl file"
    )
    parser.add_argument(
        "-l",
        dest="taxid2lineage",
        required=True,
        help="taxid to lineage lookup gzip file",
    )
    parser.add_argument(
        "-t", dest="nrtitle", required=True, help="nr title lookup pkl file"
    )
    parser.add_argument(
        "-o", dest="outfile", required=True, help="annotated output TSV file"
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

    checkfile(args.nrtitle)
    with open(args.nrtitle, "rb") as input_file:
        nrtitle = pickle.load(input_file)
    print("Done reading NR title info!")

    checkfile(args.accid2taxid)
    with open(args.accid2taxid, "rb") as input_file:
        accid2taxid = pickle.load(input_file)
    print("Done reading accid2taxid info!")
    # accid2taxid is dict()

    of = open(args.outfile, "w")
    outheader = [
        "qseqid",
        "accessionversion",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
    ]
    outheader.append("taxid")
    outheader.append("title")
    outheader.append(
        "lineage(superkingdom,phylum,class,order,family,genus,species,biotype,clade,clade1,clade10,clade11,clade12,clade13,clade14,clade15,clade16,clade17,clade18,clade19,clade2,clade3,clade4,clade5,clade6,clade7,clade8,clade9,cohort,forma,forma specialis,forma specialis1,genotype,infraclass,infraorder,isolate,kingdom,morph,no rank,no rank1,no rank2,no rank3,no rank4,no rank5,parvorder,pathogroup,section,series,serogroup,serotype,species group,species subgroup,strain,subclass,subcohort,subfamily,subgenus,subkingdom,suborder,subphylum,subsection,subspecies,subtribe,superclass,superfamily,superorder,superphylum,tribe,varietas)"
    )
    of.write("%s\n" % ("\t".join(outheader)))
    with open(args.diamondoutput, "r") as fin:
        for l in fin:
            l = l.strip().split("\t")
            accid = l[1]
            # print(accid)
            try:
                taxid = accid2taxid[accid]
            except:
                taxid = "Unknown"
            # print(taxid)
            try:
                title = nrtitle[accid]
            except:
                title = "Unknown"
            # print(title)
            try:
                lineage = taxid2lineage[taxid]
            except:
                lineage = "Unknown"
            # print(lineage)
            l.extend([taxid, title, lineage])
            of.write("%s\n" % ("\t".join(l)))
    of.close()


if __name__ == "__main__":
    main()
