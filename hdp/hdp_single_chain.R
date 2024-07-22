options(stringsAsFactors = F)
library(hdp)

args = commandArgs(trailingOnly=TRUE)

n = args[1]
n = as.numeric(n)
input_dir=args[2]
exp_name=args[3]
lower_threshold=args[4]

current_matrix = file.path(input_dir, paste0(exp_name, "_hdp_matrix.tsv"))
current_keytable = file.path(input_dir, paste0("key_table_", exp_name, ".txt"))

if(file.exists(current_matrix)){
  message("Processing matrix:\n", 
          current_matrix)
  message("n chains is: ", n)
} else{
  stop(current_matrix, " doesn't exist!")
}

out_dir = file.path(input_dir, paste0(exp_name, "_hdp_run"))

if(dir.exists(out_dir)){
  message("Outdir is:\n", 
          out_dir)
} else{
   message("Creating outdir:\n",
   out_dir)
  dir.create(out_dir)
}


mutations=read.table(current_matrix, header=T,check.names =F, sep="\t",
                     quote = "", row.names=1)

tinuc_sort <- c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")
mutations <- mutations[tinuc_sort]

if(file.exists(current_keytable)){
  message("Key table is:\n", 
          current_keytable)  
} else{
  stop(current_keytable, " doesn't exist!")
}


key_table=read.table(current_keytable, header=T, check.names=FALSE, sep="\t",
                     quote = "")

if(length(key_table)==2){
  message("Assuming first column is sample!")
  message("Assuming second column is hierarchy!")
  colnames(key_table) = c("Sample", "category")
} else{
  stop("Key table is NOT a two column file!")
}


#If requiring a minimum number of mutations:
message("Removing samples with less than ", lower_threshold , " mutations")
message(nrow(key_table), " original samples...")

sample_remove=rownames(mutations)[rowSums(mutations)<lower_threshold]
mutations=mutations[!rownames(mutations)%in%sample_remove,]
key_table=key_table[!key_table[["Sample"]]%in%sample_remove,]

message(nrow(key_table), " samples remain...")

key_table[["Type"]] <- factor(key_table[["category"]])

freq=table(key_table[["category"]])

message("Starting hdp chain step 1...")

#with Just type as parent (1 hierarchy)
hdp_mut <- hdp_init(ppindex = c(0, rep(1,length(freq)),rep(2:(length(freq)+1), times=freq)), # index of parental node
                    cpindex = c(1, rep(2,length(freq)),rep(3:(length(freq)+2), times=freq)), # index of the CP to use
                    hh = rep(1, 96), # prior is uniform over 96 categories
                    alphaa = rep(1,length(freq)+2), # shape hyperparameters for 2 CPs
                    alphab = rep(1,length(freq)+2))  # rate hyperparameters for 2 CPs

hdp_mut <- hdp_setdata(hdp_mut,
                       dpindex = (length(freq)+2):numdp(hdp_mut), # index of nodes to add data to
                       mutations)

message("Starting hdp chain step 2...")

hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), initcc=10,seed=n*300)

chain=hdp_posterior(hdp_activated,
                    burnin=30000,
                    n=100,
		    seed=n*1000,
                    space=200,
                    cpiter=3)

out_file = file.path(out_dir, paste0("hdp_chain_", n, "_", exp_name, ".Rdata"))
message("Saving chain to:\n", 
        out_file)

saveRDS(chain,out_file)

message("Finished hdp chain run!")
