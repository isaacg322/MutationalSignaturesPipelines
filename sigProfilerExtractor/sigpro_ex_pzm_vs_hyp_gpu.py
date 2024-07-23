#!/usr/bin/env python

print("Entred python...")

import sys, argparse

parser = argparse.ArgumentParser()
parser.add_argument("input", help="input matrix file or folder contating vcf files")
parser.add_argument("-o", "--output", type=str, default="output_folder", help="output folder")
parser.add_argument("-d", "--data_type", type=str, default="matrix")
parser.add_argument("-c", "--context_type", type=str, default="96")
parser.add_argument("-r", "--reference", type=str, default="GRCh38")
parser.add_argument("-n", "--cpus", type=int, default=8)
parser.add_argument("-s", "--min_sigs", type=int, default=2)
parser.add_argument("-x", "--max_sigs", type=int, default=6)
args = parser.parse_args()

print(args)

from SigProfilerExtractor import sigpro as sig

print("Starting sigprofiler analysis...")

if __name__ == '__main__':
    sig.sigProfilerExtractor(input_type=args.data_type,
    input_data=args.input,
    output=args.output,
    reference_genome=args.reference,
    opportunity_genome=args.reference,
    context_type=args.context_type,
    minimum_signatures=args.min_sigs,
    cpu=args.cpus,
    maximum_signatures=args.max_sigs,
    nmf_replicates=100, gpu=True,
    cosmic_version="3.3")
