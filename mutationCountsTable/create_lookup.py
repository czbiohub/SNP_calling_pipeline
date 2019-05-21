import argparse
from collections import defaultdict
import pickle

import gffutils
import intervaltree


def arg_parse():
        parser = argparse.ArgumentParser()
        parser.add_argument("--ref-gtf",
                            type=argparse.FileType('r'),
                            required=True)
        parser.add_argument("--output",
                            type=argparse.FileType('wb'),
                            required=True)
        return parser.parse_args()


def main():
    args = arg_parse()
    fn = gffutils.example_filename(args.ref_gtf.name)
    db = gffutils.create_db(fn,
                            ":memory:",
                            keep_order=True,
                            disable_infer_genes=True,
                            disable_infer_transcripts=True)

    interval_map = defaultdict(intervaltree.IntervalTree)
    
    for i in db.all_features():
        gene_names = i.attributes.get("gene_name")
        assert len(gene_names) == 1
        gene_name = gene_names[0]
        interval_map[i.chrom][i.start:i.stop+1] = gene_name

    with open(args.output.name, "wb") as fh:
        pickle.dump(interval_map, fh)


if __name__ == "__main__":
    main()
