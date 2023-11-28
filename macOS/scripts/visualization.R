suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))
theme_set(theme_bw())

# Define the command-line arguments
parser <- ArgumentParser(description = "HT candidates visualization")
parser$add_argument("--rep", help = "Path to summary.tsv file")
parser$add_argument("--nonrep", help = "Path to summary.tsv file")
parser$add_argument("--output1", help = "Path to output1 file")
parser$add_argument("--output2", help = "Path to output2 file")

# Parse the command-line arguments
args <- parser$parse_args()

# Read the input file1
candidates <- read_tsv(args$rep, show_col_types = FALSE, col_names = c("all", "bp", "x", "y", "z")) %>%
  select(-x, -y, -z) %>%
  separate(all, into = c("cluster","credibility","n"), sep = "-") %>% type_convert() %>%
  mutate(cluster = str_replace(cluster, ".consensus", "")) %>%
  mutate(cluster = str_remove(cluster, "cluster_")) %>%
  mutate(cred_value = case_when(credibility <= 0.5 ~ "low", credibility > 0.5 & credibility <= 1.5 ~ "optimal", credibility > 1.5 ~ "excessive"))

candidates$cred_value <- factor(candidates$cred_value, levels = c("low", "optimal", "excessive"))

# Read the input file2
non_rep <- read_tsv(args$nonrep, show_col_types = FALSE, col_names = c("all", "bp", "x", "y", "z")) %>%
  select(-x, -y, -z) %>%
  separate(all, into = c("cluster","end","credibility"), sep = "-") %>%
  mutate(cluster = paste0(cluster, "-", end)) %>% select(-end) %>% type_convert() %>%
  mutate(cred_value = case_when(credibility <= 0.5 ~ "low", credibility > 0.5 & credibility <= 1.5 ~ "optimal", credibility > 1.5 ~ "excessive"))

non_rep$cred_value <- factor(non_rep$cred_value, levels = c("low", "optimal", "excessive"))

# Visualize
plot_rep <- ggplot(candidates, aes(bp, n, color=cred_value)) +
  geom_point(size=5, alpha=0.5)+
  scale_color_manual(values=c("red", "darkgreen", "purple"))+
  labs(y="number of sequences", x="consensus length (bp)", color="credibility")+
  geom_text(aes(label = cluster), size=3, color="black", fontface="bold")

plot_nonrep <- ggplot(non_rep, aes(bp, credibility, color=cred_value)) +
  geom_point(size=5, alpha=0.5)+
  scale_color_manual(values=c("red", "darkgreen", "purple"))+
  labs(y="credibility", x="consensus length (bp)", color="credibility")+
  geom_text(data = filter(non_rep, cred_value=="optimal" & bp > 5000), aes(label = cluster), size=2, color="black", fontface="bold")

# Save visualization
ggsave(plot_rep, filename=args$output1)
ggsave(plot_nonrep, filename=args$output2)