library(tidygraph)
library(ggraph)
library(dplyr)
library(mitolina)

sites <- 100
mu <- rep(1e-3, sites)

set.seed(1)

sim <- sample_mtdna_geneology_varying_size(
  population_sizes_females = rep(5, 10),
  population_sizes_males = c(0L, rep(5, 9)),
  #generations_full = 8,
  #generations_return = 8,
  generations_full = 1,
  generations_return = 1,
  progress = FALSE)
pop <- sim$population

peds <- build_pedigrees(pop, progress = FALSE)

end_gen_peds <- lapply(sim$end_generation_male_individuals, get_pedigree_from_individual)
pedids <- sapply(end_gen_peds, get_pedigree_id)
id_tab <- sort(table(pedids))
ped_id_endgen <- as.integer(names(id_tab)[1]) # smallest

get_founder_mito <- function() {
  sample(c(FALSE, TRUE), sites, replace = TRUE)
}

pedigrees_all_populate_haplotypes_custom_founders(
  pedigrees = peds,
  mutation_rates = mu,
  get_founder_haplotype = get_founder_mito, 
  progress = FALSE)

g <- as_tbl_graph(peds) %>% 
  activate(nodes) %>% 
  filter(ped_id == ped_id_endgen) %>% 
  mutate(num_vars = sapply(haplotype, sum))


##################################################

#data("mtdna_partitions")

#tps <- do.call(rbind, g %>% activate(nodes) %>% pull(haplotype))


library(tidyverse)
d_phylo_rle <- tribble(
  ~id, ~Pos_start, ~Pos_end,
  1L, 1, 100,
  2L, 101, 200, 
  3L, 201, 300, 
  4L, 301, 400, 
  5L, 401, 500, 
  6L, 501, 600) %>% 
  I()
d_phylo_rle

#stopifnot(all.equal(d_phylo_rle %>% pull(id), g %>% activate(nodes) %>% pull(haplotype_id) %>% unique() %>% sort()))
cls <- RColorBrewer::brewer.pal(nrow(d_phylo_rle), name = "Set1")

dir.create("tmp")

for (i in g %>% activate(nodes) %>% pull(haplotype_id) %>% unique()) {
  #i <- 1
  vls <- rep("grey", d_phylo_rle %>% pull(id) %>% length())
  names(vls) <- d_phylo_rle %>% pull(id)
  
  mark <- case_when(
    i == 1L ~ list(c(1, 2, 4)),
    i == 2L ~ list(c(2)),
    i == 3L ~ list(c(2, 3)),
    i == 4L ~ list(c(2, 4)),
    i == 5L ~ list(c(2, 4, 5))
  ) %>% unlist()
  
  for (j in mark) {
    vls[j] <- cls[j]
  }
  
  pmtdna <- ggplot(d_phylo_rle) +
    geom_rect(aes(xmin = Pos_start - 1, xmax = Pos_end + 1, ymin = 1, ymax = 2, group = as.character(id), fill = as.character(id)), size = 10, show.legend = FALSE) +
    #scale_x_continuous(breaks = brks) +
    scale_y_continuous(limits = c(0, 3)) +
    scale_fill_manual(values = vls) +
    coord_polar() + 
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0, "cm"),
          panel.spacing = margin(0, 0, 0, 0, "cm"),
          plot.background = element_rect(fill = "transparent", colour = NA)) +
    labs(x = NULL, y = NULL) +
    theme(line = element_blank(),
          text = element_blank(),
          title = element_blank()) +
    NULL
  #pmtdna
  ggsave(pmtdna, file = paste0("tmp/mtdna-", i, ".png"), width = 1, height = 1, bg = "transparent")
  system(paste0("convert tmp/mtdna-", i, ".png -trim tmp/mtdna-", i, "-crop.png"), wait = TRUE)
}

library(igraph)
library(png)


#https://stackoverflow.com/questions/21423854/png-images-as-vertices-in-r-igraph

hap_id <- g %>% activate(nodes) %>% pull(haplotype_id)

gi <- as.igraph(g)

V(gi)$raster <- vector("list", vcount(gi))
#options(warn=1)
for (i in unique(hap_id)) {
  #print(i)
  img <- readPNG(paste0("tmp/mtdna-", i, "-crop.png"))
  V(gi)$raster[which(hap_id == i)] <- replicate(sum(hap_id == i), img, simplify = FALSE)
}

unlink("tmp", recursive = TRUE, force = TRUE)


lo <- layout_as_tree(graph = gi)

png("paper-fig-simulation.png", width = 3000, height = 1800)
par(mar = c(0, 0, 0, 0))
plot(gi, layout = lo, 
     xlim = range(lo[, 1]), 
     ylim = range(lo[, 2]), 
     vertex.shape = "raster", vertex.label = NA, vertex.size = 50, vertex.size2 = 50,
     edge.arrow.size = 3, 
     rescale = FALSE)
dev.off()

