#!/usr/bin/env python

print("Entred python...")

import sys, argparse

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--samples_path", help="input matrix")
parser.add_argument("-i", "--dnsignatures_path", help="input pyr probability matrix per de novo signature")
parser.add_argument("-db", "--dbsignatures_path", help="input pyr probability matrix per reference signature")
parser.add_argument("-o", "--output_path", type=str, default="output_folder", help="output folder")
parser.add_argument("-d", "--data_type", type=str, default="matrix")
parser.add_argument("-r", "--reference", type=str, default="GRCh38")
args = parser.parse_args()


from SigProfilerAssignment import Analyzer as Analyze

print("Starting decompose fit...")

Analyze.decompose_fit(samples=args.samples_path,
                      output=args.output_path,
                      input_type=args.data_type,
                      signatures=args.dnsignatures_path,
                      signature_database=args.dbsignatures_path,
                      genome_build=args.reference)
