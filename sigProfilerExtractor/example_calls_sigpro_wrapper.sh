# Examples for calling sigpro_ex_pzm_vs_hyp_gpu_launcher.sh wrapper to run
# de novo mutational signature extraction and mutational signature decompostion

# Mutational signature extraction
# MANIFEST HAS TO BE A 2 column plain text file. First column is the path to the output directory.
# Final output directory should not exist yet!.
# Second column is the path to mutation counts matrix (sigprofiler format!)
export MANIFEST=/lustre/scratch125/casm/team294rr/ig4/pzms_gel/manifest_all_pzm_combinations.txt
logs_dir=/lustre/scratch125/casm/team294rr/ig4/pzms_gel/job_logs
job_name=pzm_sigpro
total_jobs=$(wc -l ${MANIFEST} | cut -d' ' -f 1)
export CPUS=16
#total_jobs=2

bsub -J "${job_name}[1-${total_jobs}]%2" \
    -o ${logs_dir}/${job_name}_log.%J-%I \
    -e ${logs_dir}/${job_name}_err.%J-%I \
    -q gpu-huge -gpu - \
    -R 'select[mem>=32000] rusage[mem=32000]' \
    -M32000 -env "all" \
    "bash /lustre/scratch125/casm/team294rr/ig4/pzms_gel/scripts/sigpro_ex_pzm_vs_hyp_gpu_launcher.sh extract"

# Mutational signature decomposition
# MANIFEST has to be a 4 column plain text file. First column is
# output directory used for the "extract" step. Second column is the solution index chosen from manual inspection of
# results from the "extract" step. Third column is path to the reference signature database (e.g. cosmic).
# Fourth column is the ID of such database (e.g cosmic_v2), will be used as suffix in outputs from the decompose commmand
export MANIFEST=/lustre/scratch125/casm/team294rr/ig4/pzms_gel/sig3_4_decomposition4Hypermutators_manifest.txt
logs_dir=/lustre/scratch125/casm/team294rr/ig4/pzms_gel/job_logs
job_name=decompose_fit
total_jobs=$(wc -l ${MANIFEST} | cut -d' ' -f 1)
#total_jobs=2

bsub -J "${job_name}[2-${total_jobs}]%10" \
    -o ${logs_dir}/${job_name}_log.%J-%I \
    -e ${logs_dir}/${job_name}_err.%J-%I \
    -q normal \
    -R 'span[hosts=1] select[mem>=16000] rusage[mem=16000]' \
    -M16000 -env "all" \
    "bash /lustre/scratch125/casm/team294rr/ig4/pzms_gel/scripts/sigpro_ex_pzm_vs_hyp_gpu_launcher.sh decompose"
