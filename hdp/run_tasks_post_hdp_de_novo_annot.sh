## Next describes 3 post de novo signature annotation tasks
## 1) combine HDP chains and get raw components (i.e. de novo sigantures)
## 2) get all cosine similarities (vs COSMIC.v3), extract top n based on cosine sim threshold
## 3) run mutational signature deconvolution (sigprofiler based) based on most similar signatures in (2)

########################################
## Combine chains, get raw components ##
########################################

## Set paths to directories containing scripts and inputs
main_dir=/lustre/scratch125/casm/team294rr/ig4/s126_main/dnms_smokers
pzms_dir=/lustre/scratch125/casm/team294rr/ig4/pzms_gel/re_analysis2024

export combine_chain_rscript=${main_dir}/hdp_try/hdp_combine.v3.R
export input_dir=${pzms_dir}

export exp_name=dnms_vs_pzms025 # Finding files with this prefix/suffix
export is_multi_hierarchy="n" ## y for multihierarchy files
export lower_threshold=0 ## remove samples with less than ${lower_threshold} mutations
export n_chains=20 ## n_chains running as array jobs

# Conda environment for hdp (R v4.1)
export hdp_env=/lustre/scratch126/casm/team294rr/rs30/anaconda3/envs/hdp_env

## Running in LSF HPC
job_name=hdp_combine_${exp_name}
logs_dir=${input_dir}/logs
allocate_ram=20000
total_jobs=20

bsub -J "${job_name}" \
-o ${logs_dir}/${job_name}_stdout.%J \
-e ${logs_dir}/${job_name}_stderr.%J \
-q normal \
-R "span[hosts=1] select[mem>${allocate_ram}] rusage[mem=${allocate_ram}]" \
-M ${allocate_ram} \
'
source ~/.bash_profile
conda activate ${hdp_env}

Rscript --vanilla ${combine_chain_rscript} "${n_chains}" \
"${input_dir}" \
"${exp_name}" \
"${is_multi_hierarchy}" \
"${lower_threshold}"
'

##################################################
## Run cosine similarity, prepare sigfit inputs ##
##################################################

main_dir=/lustre/scratch125/casm/team294rr/ig4/s126_main/dnms_smokers
pzms_dir=/lustre/scratch125/casm/team294rr/ig4/pzms_gel/re_analysis2024

export cosine_similarity_rscript=${main_dir}/hdp_try/hdp_cosine_similarity.R
export input_dir=${pzms_dir}
export exp_name=dnms_vs_pzms025
export is_multi_hierarchy="n"
export cosine_threshold=0.8
export cosmic_location=/lustre/scratch125/casm/team294rr/ig4/resources/cosmic38/COSMIC_v3.3.1_SBS_GRCh38.txt

## Running in LSF HPC
job_name=hdp_cosine_sim_${exp_name}
logs_dir=${input_dir}/logs
allocate_ram=8000

# Conda environment for my personal (R v4.2)
export my_r_env=/lustre/scratch125/casm/team294rr/ig4/software/anaconda3/envs/R4.2

bsub -J "${job_name}" \
-o ${logs_dir}/${job_name}_stdout.%J \
-e ${logs_dir}/${job_name}_stderr.%J \
-q normal \
-R "span[hosts=1] select[mem>${allocate_ram}] rusage[mem=${allocate_ram}]" \
-M ${allocate_ram} \
'
source ~/.bash_profile
conda activate ${my_r_env}

Rscript --vanilla ${cosine_similarity_rscript} "${cosmic_location}" \
"${input_dir}" \
"${exp_name}" \
"${is_multi_hierarchy}" \
"${cosine_threshold}"
'

##############################
## Run sigfit deconvolution ##
##############################

main_dir=/lustre/scratch125/casm/team294rr/ig4/s126_main/dnms_smokers
pzms_dir=/lustre/scratch125/casm/team294rr/ig4/pzms_gel/re_analysis2024

export input_dir=${pzms_dir}

export deconvolution_script=${main_dir}/hdp_try/deconvolute_hdp_signatures.py
export exp_name=dnms_vs_pzms025
export costhreshold="0.7" # Threshold used to select SBS to deconvolute
export multi_hierarchy=0 # Set to 1 for multi-hierarchy runs

## Running in LSF HPC
job_name=sigfit_hdp_${exp_name}
logs_dir=${input_dir}/logs
allocate_ram=16000

export sigpro_extract=/lustre/scratch125/casm/team294rr/ig4/software/anaconda3/envs/sigpro_extract

bsub -J "${job_name}" \
-o ${logs_dir}/${job_name}_stdout.%J-%I \
-e ${logs_dir}/${job_name}_stderr.%J-%I \
-q normal \
-R "span[hosts=1] select[mem>${allocate_ram}] rusage[mem=${allocate_ram}]" \
-M ${allocate_ram} \
'
source ~/.bash_profile
conda activate ${sigpro_extract}

if [ "${multi_hierarchy}" -eq 1 ]
then
  exp_dir=${input_dir}/${exp_name}_hdp_run_multi_hierarchy
else 
  exp_dir=${input_dir}/${exp_name}_hdp_run
fi

if [ -d "${exp_dir}" ]
then
  echo "Entring folder containing HDP output..."
  echo "Starting sigProfiler decomponseFit script..."
  python ${deconvolution_script} ${exp_dir} \
  --experiment_name ${exp_name} \
  --costhreshold ${costhreshold}
else 
  echo "Specified HDP output folder does NOT exist!"
fi
'
