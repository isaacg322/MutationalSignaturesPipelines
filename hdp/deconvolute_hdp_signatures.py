#!/usr/bin/env python

print("Entred python...")

import sys, argparse, os, re

parser = argparse.ArgumentParser()
parser.add_argument("input_folder", help="input folder with hdp results")
parser.add_argument("-e", "--experiment_name", type=str, help="experiment name")
parser.add_argument("-r", "--reference", type=str, default="GRCh38")
parser.add_argument("-t", "--costhreshold", type=float, default=0.8)
parser.add_argument("-s", "--sigthreshold", type=float, default=0.8)
args = parser.parse_args()

print(args)

from SigProfilerAssignment import Analyzer as Analyze

os.chdir(args.input_folder)

signatures_name=args.experiment_name+"_hdp2deconv_cosinesim_he_"+re.sub("\.", "", str(args.costhreshold))+".txt"

print("HDP signatures are: ", signatures_name)

cosmic_name=args.experiment_name+"_cosmic2deconv_cosinesim_he_"+re.sub("\.", "", str(args.costhreshold))+".txt"

print("COSMIC reference to decompose is: ", cosmic_name)

samples_name=args.experiment_name+"_sigpro_matrix_4deconv.txt"

print("Sample matrix is: ", samples_name)

out_name=args.experiment_name+"sigpro_deconvolution_cos_he"+re.sub("\.", "", str(args.costhreshold))

print("Starting decompose fit...")

if __name__ == '__main__':
  Analyze.decompose_fit(samples=samples_name,
  output=out_name,
  input_type="matrix",
  signatures=signatures_name,
  signature_database=cosmic_name,
  genome_build=args.reference,
  verbose=False, make_plots=True,
  new_signature_thresh_hold=args.sigthreshold,
  exome=False)
