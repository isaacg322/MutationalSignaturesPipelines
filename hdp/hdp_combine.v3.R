#!/usr/bin/env Rscript

#conda activate /lustre/scratch126/casm/team294rr/rs30/anaconda3/envs/hdp_env

options(stringsAsFactors = F)
library(hdp)
library("RColorBrewer")

args = commandArgs(trailingOnly=TRUE)

n_chains = args[1]
#n_chains = 20
input_dir=args[2]
#input_dir="/lustre/scratch125/casm/team294rr/ig4/s126_main/dnms_smokers/hdp_try"
exp_name=args[3]
#exp_name="kenichisOnly_no_unmatchREF_byDonor"
is_multi_hierarchy=args[4]
#is_multi_hierarchy="y"
lower_threshold=args[5]
#lower_threshold=0

if(is_multi_hierarchy=="y"){
  out_dir_suffix="_hdp_run_multi_hierarchy"
  chain_prefix="_multi_hierarchy.Rdata"
}else{
  out_dir_suffix="_hdp_run"
  chain_prefix=".Rdata"
}


current_matrix = file.path(input_dir, paste0(exp_name, "_hdp_matrix.tsv"))
current_keytable = file.path(input_dir, paste0("key_table_", exp_name, ".txt"))

out_dir = file.path(input_dir, paste0(exp_name, out_dir_suffix))
message("Moving to directory containing chains...")
setwd(out_dir)

message("Reading single chains...")

chlist <- vector("list", n_chains)

for (i in 1:n_chains){
  if(is_multi_hierarchy=="y"){
    curr_chain = paste0("hdp_chain_",i, chain_prefix)
  } else{
    curr_chain = paste0("hdp_chain_",i, "_", exp_name, chain_prefix)
  }
  if(file.exists(curr_chain)){
    message("Found chain ", i)
    chlist[[i]] <- readRDS(curr_chain)
  } else{
    stop("Missing chain: ", i, " check your single chain run!")
  }
}
if(any(unlist(lapply(chlist,is.null)))) chlist=chlist[-which(unlist(lapply(chlist,is.null)))]

message("Plotting QC...")
mut_example_multi <- hdp_multi_chain(chlist)
pdf(paste0("QC_plots_chain_", exp_name, ".pdf"))
par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
p1 <- lapply(chains(mut_example_multi), plot_lik, bty="L", start=1000)
p2 <- lapply(chains(mut_example_multi), plot_numcluster, bty="L")
p3 <- lapply(chains(mut_example_multi), plot_data_assigned, bty="L")
dev.off()

message("Saving multi-chain file...")
mut_example_multi <- hdp_extract_components(mut_example_multi) #This step can take a while. If too long, submit R script as job
saveRDS(mut_example_multi,paste0("HDP_multi_chain", exp_name, ".Rdata"))

### Restart from here if necesary
#mut_example_multi= readRDS(paste0("HDP_multi_chain", exp_name, ".Rdata"))

message("Plotting muts attributed...")
pdf(paste0("muts_attributed_", exp_name, ".pdf"))
plot_comp_size(mut_example_multi, bty="L")
dev.off()

trinuc_context <- sapply(strsplit(colnames(mut_count), '\\.'), `[`, 4)
group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                              each=16))
mut_colours=c("dodgerblue","black","red","grey70","olivedrab3","plum2")

#dev.new(width=12,height=4)
#par(mfrow=c(3,4))


for (i in 0:mut_example_multi@numcomp){
  pdf(paste0("hdp_component_",i,"_", exp_name, ".pdf"),width=12,height=4)

  plot_comp_distn(mut_example_multi, cat_names=trinuc_context,
                  grouping=group_factor, col=mut_colours,comp=i,
                  col_nonsig="grey80", show_group_labels=TRUE)
  dev.off()
}

nb.cols <- n_chains
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

pdf(paste0("hierarchy_exposures_", exp_name, ".pdf"))
plot_dp_comp_exposure(mut_example_multi,
                      dpindices=2:4, incl_numdata_plot=FALSE,
                      col=mycolors,
                      incl_nonsig=TRUE, cex.names=0.8,
                      ylab_exp = 'Signature exposure', leg.title = 'Signature')
dev.off()

message("Reading original matrix and key table...")
mutations=read.table(current_matrix, header=T, check.names=FALSE, sep="\t",quote = "", row.names=1)
key_table=read.table(current_keytable, header=T, check.names=FALSE, sep="\t",quote = "")

#If requiring a minimum number of mutations:
message("Removing samples with less than ", lower_threshold , " mutations")
message(nrow(key_table), " original samples...")

sample_remove=rownames(mutations)[rowSums(mutations)<lower_threshold]
mutations=mutations[!rownames(mutations)%in%sample_remove,]
key_table=key_table[!key_table[["Sample"]]%in%sample_remove,]

message(nrow(key_table), " samples remain...")

freq=table(key_table[["Sample"]])

message("Plotting sample exposures...")
pdf(paste0("signature_attribution_", exp_name, ".pdf"),width=10,height=8)
plot_dp_comp_exposure(mut_example_multi,
                      dpindices=(length(freq)+2):length(mut_example_multi@comp_dp_counts),
                      incl_nonsig = T, ylab_exp = 'Signature exposure', leg.title = 'Signature', col=mycolors, incl_numdata_plot=F)
dev.off()


dp_distn <- comp_dp_distn(mut_example_multi)
ndp <- nrow(dp_distn$mean)
ncomp <- ncol(dp_distn$mean)

message("Saving summary tables...")
mean_assignment=as.data.frame(comp_dp_distn(mut_example_multi)$mean)
write.table(mean_assignment,paste0("mean_assignment_hdp_", exp_name,".txt"))
mean_sigs=as.data.frame(t(comp_categ_distn(mut_example_multi)$mean))
write.table(mean_sigs,paste0("hdp_sigs_", exp_name,".txt"))

message("Programme finished!")
