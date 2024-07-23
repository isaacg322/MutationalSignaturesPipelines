#!/usr/bin/env bash

source ~/.bash_profile

conda activate sigpro_extract

echo "${LSB_JOBINDEX}"

run_what=${1}

curr_line=$(sed -n "${LSB_JOBINDEX}p" ${MANIFEST})

if [ ${run_what} == "extract" ]; then
  # Expects a column with two lines: full path to output dir and full path to
  # input matrix
  outdir=$(echo "${curr_line}" | cut -f 1)       # A NON existent directory
  input_matrix=$(echo "${curr_line}" | cut -f 2) # A counts matrix (96pyr x Samples)

  echo "Starting sigProfilerExtractor script..."

  python /lustre/scratch125/casm/team294rr/ig4/pzms_gel/scripts/sigpro_ex_pzm_vs_hyp_gpu.py \
    ${input_matrix} \
    --output ${outdir} \
    --cpus ${CPUS} --min_sigs ${MINSIGS} --max_sigs ${MAXSIGS}

elif [ ${run_what} == "denovofit" ]; then

  lookInFolder=$(echo "${curr_line}" | cut -f 1) # Should be same as output folder above
  signature_solution=$(echo "${curr_line}" | cut -f 2)
  outDNFit=${lookInFolder}/deNovo_fit
  dnSigName=SBS96_S${signature_solution}_Signatures
  DNSignatures=${lookInFolder}/SBS96/All_Solutions/SBS96_${signature_solution}_Signatures/Signatures/${dnSigName}.txt
  samples_file=${lookInFolder}/SBS96/Samples.txt

  [ ! -d ${outDNFit} ] && echo "Directory deNovo_fit DOES NOT exists...creating"
  mkdir ${outDNFit}

  outDNFit=${lookInFolder}/deNovo_fit/${dnSigName}

  [ ! -d ${outDNFit} ] && echo "Directory for deNovo_fit (S${signature_solution}) DOES NOT exists...creating"
  mkdir ${outDNFit}

  echo "Starting deNovo fit script..."

  python /lustre/scratch125/casm/team294rr/ig4/pzms_gel/scripts/deNovo_fit.py \
    --samples_path ${samples_file} \
    --dnsignatures_path ${DNSignatures} \
    --output_path ${outDNFit}

elif [ ${run_what} == "decompose" ]; then

  lookInFolder=$(echo "${curr_line}" | cut -f 1)       # Should be same as output folder in extractor mode
  signature_solution=$(echo "${curr_line}" | cut -f 2) # An index for an existing solution
  sigDB=$(echo "${curr_line}" | cut -f 3)              # Path to sigDB
  sigDBID=$(echo "${curr_line}" | cut -f 4)            # To be used in output folder name
  outDecFit=${lookInFolder}/decomposed_fit
  dnSigName=SBS96_S${signature_solution}_Signatures
  DNSignatures=${lookInFolder}/SBS96/All_Solutions/SBS96_${signature_solution}_Signatures/Signatures/${dnSigName}.txt
  samples_file=${lookInFolder}/SBS96/Samples.txt

  [ ! -d ${outDecFit} ] && echo "Directory decompose_fit DOES NOT exists...creating"
  mkdir ${outDecFit}

  outDecFit=${lookInFolder}/decomposed_fit/${dnSigName}

  [ ! -d ${outDecFit} ] && echo "Directory for decompose_fit (S${signature_solution}) DOES NOT exists...creating"
  mkdir ${outDecFit}

  outDecFit=${lookInFolder}/decomposed_fit/${dnSigName}/fit2dbID_${sigDBID}

  [ ! -d ${outDecFit} ] && mkdir ${outDecFit}

  echo "Starting decompose fit script..."

  python /lustre/scratch125/casm/team294rr/ig4/pzms_gel/scripts/decompose_fit.py \
    --samples_path ${samples_file} \
    --dnsignatures_path ${DNSignatures} \
    --output_path ${outDecFit} \
    --dbsignatures_path ${sigDB}

else
  echo "Valid options are: extract, denovofit, decompose"
fi
