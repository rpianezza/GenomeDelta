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
  separate(all, into = c("x","cluster","credibility","n"), sep="_") %>%
  mutate(cluster = gsub(".consensus", "", cluster)) %>% select(-x) %>% type_convert()

# Read the input file2
non_rep <- read_tsv(args$nonrep, show_col_types = FALSE, col_names = c("all", "bp", "x", "y", "z")) %>%
  select(-x, -y, -z) %>% separate(all, into = c("cluster","credibility"), sep="_") %>% type_convert()
  
#print(candidates)
#print(non_rep)

# Visualize
plot_rep <- ggplot(candidates, aes(bp, n, color=credibility)) +
  geom_point(size=5, alpha=0.5)+
  scale_colour_gradient2(low = "red", mid = "darkgreen", high = "purple", midpoint = 0, limits = c(-1, 1))+
  labs(y="number of sequences", x="consensus length (bp)", color="coverage bias")+
  geom_text(aes(label = cluster), size=3, color="black", fontface="bold")

plot_nonrep <- ggplot(non_rep, aes(bp, credibility, color=credibility)) +
  geom_point(size=5, alpha=0.5)+
  scale_colour_gradient2(low = "red", mid = "darkgreen", high = "purple", midpoint = 0, limits = c(-1, 1))+
  labs(y="coverage bias", x="consensus length (bp)", color="coverage bias")

# Save visualization
ggsave(plot_rep, filename=args$output1)
ggsave(plot_nonrep, filename=args$output2)