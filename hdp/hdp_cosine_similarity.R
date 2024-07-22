#!/usr/bin/env Rscript

#conda activate /lustre/scratch125/casm/team294rr/ig4/software/anaconda3/envs/R4.2

library(lsa)
library(tidyverse)
library(lattice)
library(gtools)
library(vroom)

args = commandArgs(trailingOnly=TRUE)
cosmic_location=args[1]
#cosmic_location="/lustre/scratch125/casm/team294rr/ig4/resources/cosmic38/COSMIC_v3.3.1_SBS_GRCh38.txt"
input_dir=args[2]
#input_dir="/lustre/scratch125/casm/team294rr/ig4/s126_main/dnms_smokers/hdp_try"
exp_name=args[3]
#exp_name="kenichisSamplesOnly_no_unmatchREF"
is_multi_hierarchy=args[4]
#is_multi_hierarchy="y"
CosThres=args[5]
#CosThres=0.8
CosThres=as.numeric(CosThres)

if(is_multi_hierarchy=="y"){
  out_dir_suffix="_hdp_run_multi_hierarchy"
}else{
  out_dir_suffix="_hdp_run"
}

ref=read.table(cosmic_location ,header=T, stringsAsFactors = F, row.names = 1, sep = '\t')
sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
full_vec = paste0(rep(c("A","C","G","T"),each=4),"[",rep(sub_vec,each=16),"]",rep(c("A","C","G","T"),times=4))
ref=ref[full_vec,]



out_dir = file.path(input_dir, paste0(exp_name, out_dir_suffix))
message("Moving to directory containing hdp output...")
setwd(out_dir)

# Check this for original Mimy's input 
# test=read.table("/lustre/scratch126/casm/team294rr/mp29/Emily_blood/signatures/HDP_w_contaminated/hdp_PD_method_woPD50306/run/hdp_sigs_new.txt")

hdp_sigs=read.table(file.path(paste0("hdp_sigs_", exp_name, ".txt")), sep=" ")
tinuc_sort = c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")

rownames(hdp_sigs) = tinuc_sort
hdp_sigs=hdp_sigs[rownames(ref),] # Reorder to match COSMIC pyr subs order

# Name components
max_components=length(hdp_sigs)-1
colnames(hdp_sigs)= paste0("component", 0:max_components)

cosine_matrix=data.frame(matrix(nrow=ncol(hdp_sigs), ncol=ncol(ref)))
rownames(cosine_matrix)=colnames(hdp_sigs)
colnames(cosine_matrix)=colnames(ref)

for (n in 1:nrow(cosine_matrix)) {
  for (m in 1:ncol(cosine_matrix)) {
    cosine_matrix[n,m] <- cosine(x=hdp_sigs[,rownames(cosine_matrix)[n]],
                                 y=ref[,colnames(cosine_matrix)[m]])
  }
}

## Make a readabe version for the cosine similarities matrix

all_cosine_similarities = as_tibble(cosine_matrix) %>%
  mutate(hdp_sig=rownames(cosine_matrix)) %>%
  pivot_longer(-hdp_sig, names_to = "cosmic_sbs", values_to = "cosine_similarity") %>%
  arrange(hdp_sig, desc(cosine_similarity))

message("Writting cosine similarities report...")
write.table(all_cosine_similarities, 
            file.path(out_dir, paste0(exp_name,"_all_cosine_smilarities.tsv")), 
            quote = F,  row.names = F, sep="\t")

message("Writting cosine similarities plot...")
pdf(file.path(out_dir, paste0(exp_name,"_all_cosine_similarities.pdf")), height=5, 
    width=15)
color.palette = colorRampPalette(c("white", "orange", "purple"))
levelplot(t(cosine_matrix[dim(cosine_matrix)[1]:1,]),col.regions=color.palette, 
          aspect="fill", scales=list(x=list(rot=90)))
dev.off()

message("Writting matrices for signatures with cosine similarities >= ", CosThres)

all_cos_heCosThres=all_cosine_similarities %>%
  filter(cosine_similarity>=CosThres) %>%
  mutate(cosine_flag=paste0("cos>=", CosThres))

excluded_components = setdiff(unique(all_cosine_similarities[["hdp_sig"]]) ,
        unique(all_cos_heCosThres[["hdp_sig"]]))

if(length(excluded_components)>0){
  message("Found exluded components...will recover top 2 signatures for these!")
  message("Excluded components are: ", excluded_components)
  exc_recovered = all_cosine_similarities %>%
  filter(hdp_sig %in% excluded_components) %>%
    arrange(hdp_sig, desc(cosine_similarity)) %>%
    group_by(hdp_sig) %>%
    top_n(n=2) %>% 
    ungroup() %>%
    mutate(cosine_flag=paste0("top2_only"))

  all_cos_heCosThres=rbind(all_cos_heCosThres, exc_recovered) %>% arrange(hdp_sig)

} else{
  message("No components were excluded due to insufficient cosine similarity!")
  all_cos_heCosThres = all_cos_heCosThres %>%
    arrange(hdp_sig)
}

hdp_components = unique(all_cos_heCosThres[["hdp_sig"]] )
names(hdp_components)=paste0("SBS96",LETTERS[c(1:length(hdp_components))])
hdp_components_dict= data.frame(hdp_sig=hdp_components, 
                                sbs96_name=names(hdp_components))

all_cos_heCosThres = all_cos_heCosThres %>%
  left_join(., hdp_components_dict,  by="hdp_sig")

vroom_write(all_cos_heCosThres, 
            file =  paste0(exp_name,"_selected_cosine_smilarities_cosinesim_he_", 
                           gsub(CosThres, pattern="\\.", replacement=""), 
                                      ".tsv"), delim="\t")

## Select COSMIC SBS to deconvolute
cosmic_sbs2deconv_sbsids = mixedsort(unique(c(all_cos_heCosThres[["cosmic_sbs"]], "SBS1", "SBS5")))
cosmic_sbs2deconv=ref[cosmic_sbs2deconv_sbsids] %>%
  rownames_to_column(var = "Type")

write.table(cosmic_sbs2deconv, 
            file.path(out_dir, paste0(exp_name,"_cosmic2deconv_cosinesim_he_", 
                                      gsub(CosThres, pattern="\\.", replacement=""), 
                                      ".txt")), quote = F,  row.names = F, sep="\t")

## Select HDP SBS to deconvolute
hdp_sbs2deconv_ids = mixedsort(unique(all_cos_heCosThres[["hdp_sig"]]))
hdp_sbs2deconv = hdp_sigs[hdp_sbs2deconv_ids] %>%
  rownames_to_column(var = "Type")

colnames(hdp_sbs2deconv) = c("Type", names(hdp_components))

write.table(hdp_sbs2deconv, 
            file.path(out_dir, paste0(exp_name,"_hdp2deconv_cosinesim_he_", 
                                      gsub(CosThres, pattern="\\.", replacement=""), 
                                      ".txt")), quote = F, row.names = F, sep="\t")

message("Preparing sample input for sigprofiler...")

hdp_matrix = vroom(file.path(input_dir, paste0(exp_name, "_hdp_matrix.tsv")),
                 delim= "\t") %>%
  column_to_rownames("Sample") %>% as.matrix() %>% 
  t()


hdp_matrix_out = as_tibble(hdp_matrix[rownames(ref), ]) %>%
  mutate(MutationType=rownames(hdp_matrix))

last_column = length(hdp_matrix_out)
pre_last_column = length(hdp_matrix_out)-1

message("Writing sample matrix in sigprofiler style...")

hdp_matrix_out[,c(last_column,1:pre_last_column)] %>%
  vroom_write(., file.path(out_dir, paste0(exp_name, "_sigpro_matrix_4deconv.txt")))

message("Finished programme...")
