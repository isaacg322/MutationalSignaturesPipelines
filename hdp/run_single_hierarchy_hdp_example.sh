## Set paths to directories containing scripts and inputs
main_dir=/lustre/scratch125/casm/team294rr/ig4/s126_main/dnms_smokers
pzms_dir=/lustre/scratch125/casm/team294rr/ig4/pzms_gel/re_analysis2024

export input_dir=${pzms_dir}

# Running single chain based de novo signature annotation
export single_chain_rscript=${main_dir}/hdp_try/hdp_single_chain.R 
export exp_name=dnms_vs_pzms025 ## Will be used as prefix in outputs
export lower_threshold=0 # Remove samples containing less than ${lower_threshold} mutations
# Conda environment for hdp (R v4.1)
export hdp_env=/lustre/scratch126/casm/team294rr/rs30/anaconda3/envs/hdp_env

## Running jobs in LSF HPC
## Set stout and stderr names and output paths
job_name=hdp_${exp_name}
logs_dir=${input_dir}/logs
allocate_ram=60000 
total_jobs=20 ## 20 chains - running as array jobs

bsub -J "${job_name}[1-${total_jobs}]%20" \
-o ${logs_dir}/${job_name}_${exp_name}_stdout.%J-%I \
-e ${logs_dir}/${job_name}_${exp_name}_stderr.%J-%I \
-q normal \
-n 2 \
-R "span[hosts=1] select[mem>${allocate_ram}] rusage[mem=${allocate_ram}]" \
-M ${allocate_ram} \

'
source ~/.bash_profile
conda activate ${hdp_env}

nchain=${LSB_JOBINDEX}

Rscript --vanilla ${single_chain_rscript} "${nchain}" \
"${input_dir}" \
"${exp_name}" \
"${lower_threshold}"
'
