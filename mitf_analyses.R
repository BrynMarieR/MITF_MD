#### HELPER FUNCTIONS

# On average, the slow-down is more dramatic for the 
# MITF:DNA:7 simulated complex than for the MITF:DNA:8 
# complex, with an average number of swaps of 66 and 144, 
# respectively, while the number of swaps for the simulated 
# unliganded MITF:DNA is on average 164. 

setwd("/Users/brynr/Desktop/mitf_paper/MITF_revisions")
source("utils.R")


#### LOGIC

data_dirs <- paste0("./distance_files05072025/", 
       list.files("./distance_files05072025"))

#data_dirs <- paste0("./distance_files_noions/",
#                    list.files("./distance_files_noions"))

helix1_data <- data.frame(lapply(data_dirs, function(x) {
  f1 <- read.delim(paste0(x,"/helix1.out"), skip=0,
                   header=FALSE,
                   comment.char="#", sep="")
  as.numeric(f1$V2)
}))

helix2_data <- data.frame(lapply(data_dirs, function(x) {
  f1 <- read.delim(paste0(x,"/helix2.out"), skip=0,
                   header=FALSE,
                   comment.char="#", sep="")
  as.numeric(f1$V2)
}))

#base_names = c("1bzkq_o1_noions", 
#               "1gybv_o2_noions",
#               "mitf_nolig_noions")

base_names = c("c8_o1_nodna", "c8_o1_withdna",
               "c8_o2_nodna", "c8_o2_withdna",
               "c7_o1_nodna", "c7_o1_withdna",
               "c7_o2_nodna", "c7_o2_withdna",
               "nolig_nodna", "nolig_withdna")

colnames(helix1_data) <- paste0(base_names, "_h1")
colnames(helix2_data) <- paste0(base_names, "_h2")

all_helix_data = cbind(helix1_data, helix2_data)
# re-order to interleave columns
all_helix_data <- all_helix_data[rep(1:length(base_names), 
                                     each = 2) + (0:1) * length(base_names)]


my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 50)
truncate_first <- create_truncated_indices(500)
pdf("correlation_heatmap.pdf", width=6, height=6)
pheatmap(cor(apply(all_helix_data[truncate_first,],
             2, rolling_average, 50)),
         cluster_rows=FALSE, cluster_cols=FALSE,
         col=my_palette, breaks=seq(-1, 1, l=50),
         cellwidth = 15, cellheight=15)
dev.off()

result_list = list()
data_list = list()
par(mfrow=c(1,1))

average_helix_data = apply(all_helix_data[truncate_first,], 2, 
                           rolling_average, 50)

for (base_name in base_names) {
  vals_per_run = dim(average_helix_data)[1]/5
  for (run_num in c(1:5)) {
    rows_of_interest = (1+(run_num-1)*vals_per_run):(run_num*vals_per_run)
    
    tmp_plot_df = data.frame(average_helix_data[rows_of_interest,
                                     c(paste0(base_name, "_h1"),
                                       paste0(base_name, "_h2"))])
    colnames(tmp_plot_df) <- c("Helix_1", "Helix_2")
    tmp_plot_df$time <- (c(1:vals_per_run)/10)
    
    title_str = gsub("_withdna", ", with DNA", 
                     gsub("_nodna", ", no DNA", 
                          gsub("_o", ", orientation ", 
                               gsub("nolig", "MITF, no compound",
                                gsub("c", "MITF, compound ", 
                                    base_name)))))
    
    p <- ggplot(tmp_plot_df, 
                aes(x=time)) + 
      geom_line(aes(y=Helix_1), color="red") + 
      geom_line(aes(y=Helix_2), color="blue") + 
      xlab("Time (ns)") +
      ylab("Helix length (Angstroms)") +
      ggtitle(paste0("Helix length trace for ", title_str)) +
      theme_pubr()
    print(p)
    ggsave(paste0(base_name, "_run_", run_num, "_helix_trace.pdf"), width = 8, height = 3)
  }
  
  # count crossings per run
  swap_count = c(0,0,0,0,0)
  for (run_num in c(1:5)) {
    rows_of_interest = (1+(run_num-1)*vals_per_run):(run_num*vals_per_run)
    diffs = average_helix_data[rows_of_interest,paste0(base_name, "_h1")] - 
      average_helix_data[rows_of_interest,paste0(base_name, "_h2")]
    binarized_diffs = diffs / abs(diffs)
    swap_count[run_num] = sum(rolling_difference(binarized_diffs)/2, na.rm=TRUE)
  }

  result_list[[base_name]] = swap_count
}

result_df = sapply(data.frame(result_list), as.numeric)


ggplot_results <- data.frame(
  compound = c(rep("8",4),rep("7",4),"none", "none"),
  DNA = c("absent", "present"),
  orientation = c(rep(c("orient1", "orient1",
                        "orient2", "orient2"),2), "none", "none"),
  run1 = result_df[1,],
  run2 = result_df[2,],
  run3 = result_df[3,],
  run4 = result_df[4,],
  run5 = result_df[5,]
)

ggplot_reshaped <- reshape(ggplot_results, 
                           direction = "long",
                           varying = list(names(ggplot_results)[4:8]),
                           v.names = "swaps",
                           idvar = c("compound", "DNA", "orientation"),
                           timevar = "Run",
                           times = c("run1","run2","run3","run4","run5"))

ggplot_reshaped$compound <- factor(ggplot_reshaped$compound,
                                   levels=c("none","7","8"))
ggplot_reshaped$orientation <- factor(ggplot_reshaped$orientation,
                                   levels=c("none", "orient1", "orient2"))
ggplot_reshaped$DNA <- factor(ggplot_reshaped$DNA,
                                      levels=c("absent","present"))

p <- ggplot(ggplot_reshaped, aes(x=DNA, y=swaps)) + 
  ylim(0, 260) + 
  geom_boxplot(width=0.1) + 
  geom_jitter(shape=16, position=position_jitter(0.04)) + 
  theme_pubr(base_size=20) + 
  geom_pwc(method="wilcox_test",
           label="p.signif",
           label.size=6,
           tip.length=0.02)
print(p)
ggsave("swaps_by_DNA.pdf", width = 5, height = 4)

ggplot_reshaped$compound_and_orientation = paste0(
  ggplot_reshaped$compound, "_", ggplot_reshaped$orientation
)
ggplot_reshaped$compound_and_orientation <- factor(
  ggplot_reshaped$compound_and_orientation, 
  ordered=TRUE,
  levels=c("none_none", "7_orient1", "7_orient2",
           "8_orient1", "8_orient2")
)
p <- ggplot(ggplot_reshaped, aes(x=compound_and_orientation,
                                 y=swaps, fill=orientation)) + 
  xlab("compound and orientation") + 
  geom_boxplot() + 
  scale_fill_manual(values=c("grey", "orange", "forestgreen")) + 
  theme_pubr(base_size=20) + 
  theme(axis.ticks = element_blank()) +
  scale_x_discrete(labels=c("none_none" = "no ligand", 
                            "7_orient1" = "cmpd 7",
                            "7_orient2" = "",
                            "8_orient1" = "cmpd 8",
                            "8_orient2" = "")) +
  geom_pwc(method = "wilcox_test",
           label="p.signif",
           ref.group="none_none",
           p.adjust.method = "fdr",
           tip.length = 0.02, 
           group.by = "x.var",
           hide.ns = TRUE,
           vjust=0.6,
           label.size=6
  )
print(p)

ggsave("swaps_by_compound.pdf", width = 6, height = 4)

dev.off()
