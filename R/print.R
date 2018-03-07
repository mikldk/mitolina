#' @export
print.mitolina_population_abort <-
  function(x, ...) {
    if (!is(x, "mitolina_population_abort")) stop("x must be a mitolina_population_abort object")
    cat("Operation was cancelled, hence the assembly was not finished.\n")    
    return(invisible(NULL))
  }

#' @export
print.mitolina_population <-
  function(x, ...) {
    if (!is(x, "mitolina_population")) stop("x must be a mitolina_population object")
    
    cat("Population with ", formatC(pop_size(x), big.mark = ","), " individuals\n", sep = "")
    
    return(invisible(NULL))
  }

#' @export
print.mitolina_pedigreelist <-
  function(x, ...) {
    if (!is(x, "mitolina_pedigreelist")) stop("x must be a mitolina_pedigreelist object")
    
    sizes <- unlist(lapply(1L:pedigrees_count(x), function(i) pedigree_size(x[[i]])))
    sizes_str <- ""
    
    max_print <- 6L
    
    if (length(sizes) > 0L) {
      if (length(sizes) <= max_print) {
        sizes_str <- paste0(" (of size ", paste0(sizes, collapse = ", "), ")")
      } else {
        sizes_str <- paste0(" (of size ", paste0(sizes[1L:max_print], collapse = ", "), ", ...)")
      }
    }
    
    cat("List of ", formatC(pedigrees_count(x), big.mark = ","), " pedigrees", sizes_str, "\n", sep = "")
    
    return(invisible(NULL))
  }

stop_invalid_id <- function(id) {
  if (length(id) != 1L || !is.numeric(id) || id <= 0L || round(id) != id) {
    stop("Invalid id: ", id)
  }
}
  
#' @export
`[[.mitolina_pedigreelist` <- function(x, ...) {
  i <- ..1
  stop_invalid_id(i)
  
  #if (length(i) != 1L || !is.integer(i) || i[1L] <= 0L || i > pedigrees_count(x)) {
  if (i > pedigrees_count(x)) {
    stop("Wrong pedigree selected (not that many pedigrees exist)")
  }
  
  p <- get_pedigree(x, i - 1L) # -1 to go to 0-based indexing
  return(p)
}

#' @export
`[[.mitolina_population` <- function(x, ...) {
  pid <- ..1
  stop_invalid_id(pid)
  
  #if (!is.integer(pid)) {
  #  pid <- as.integer(pid)
  #  warning("Converting to integer explicitely (remember L postfix)")
  #}
  
  p <- get_individual(x, pid)
  return(p)
}

#' @export
print.mitolina_pedigree <-
  function(x, ...) {
    if (!is(x, "mitolina_pedigree")) stop("x must be a mitolina_pedigree object")
    
    print_pedigree(x)
    
    return(invisible(NULL))
  }
  


#' @importFrom igraph graph_from_data_frame plot.igraph union layout_as_tree layout.reingold.tilford vcount V
#' @import tibble
#' @importFrom graphics par
#' @importFrom methods is
#' @importFrom utils head
#' @export
pedigree_as_igraph <-
  function(x, ...) {
    if (!is(x, "mitolina_pedigree")) stop("x must be a mitolina_pedigree object")
    
    ginfo <- get_pedigree_as_graph(x)
    g <- igraph::graph_from_data_frame(ginfo$edgelist, directed = TRUE, vertices = ginfo$nodes)
    
    #co <- igraph::layout_nicely(g, dim = 2)
    co <- igraph::layout_as_tree(g, mode = "out")
    attr(g, "layout") <- co

    return(g)
  }

#   
# #' @export  
# #plot_pedigrees <-
# plot.mitolina_pedigreelist <-
#   function(x, ...) {
#     pedigrees <- x
#     peds_gs <- lapply(1L:pedigrees_count(pedigrees), function(i) pedigree_as_igraph(pedigrees[[i]]))
#     
#     big_graph <- do.call(igraph::union, peds_gs)
# 
#     #http://stackoverflow.com/questions/15558218/draw-multiple-discrete-networks-in-r-using-igraph
#     
#     roots <- sapply(lapply(peds_gs, igraph::topological.sort), head, n = 1)
#     coords <- mapply(FUN = igraph::layout.reingold.tilford, peds_gs, root = roots, SIMPLIFY = FALSE)
#     
#     ## Put the graphs side by side, roots on the top
#     width <- sapply(coords, function(x) { r <- range(x[, 1]); r[2] - r[1] })
#     gap <- 0.5
#     shift <- c(0, cumsum(width[-length(width)] + gap))
#     ncoords <- mapply(FUN=function(mat, shift) {
#       mat[,1] <- mat[,1] - min(mat[,1]) + shift
#       mat[,2] <- mat[,2] - max(mat[,2])
#       mat
#     }, coords, shift, SIMPLIFY=FALSE)
#     
#    ## Put together the coordinates for the original graph,
#     ## based on the names of the vertices
#     lay <- matrix(0, ncol = 2, nrow = igraph::vcount(big_graph))
#     for (i in seq_along(peds_gs)) {
#       lay[match(igraph::V(peds_gs[[i]])$name, igraph::V(big_graph)$name),] <- ncoords[[i]]
#     }
#     
#     ## Plot everything
#     old_mar <- par("mar")
#     par(mar = c(0, 0, 0, 0))
#     igraph::plot.igraph(big_graph, layout = lay, ...)
#     par(mar = old_mar)
#         
#     return(invisible(NULL))
#     #eturn(g)
#   }




#' @export
#tidy_graph_tbl <- function(x, ...) {
#  if (!is(x, "mitolina_pedigreelist")) stop("x must be a mitolina_pedigreelist object")
#  
#  ret <- get_pedigrees_tidy(peds)
#  
#  d_edges <- bind_rows(lapply(seq_along(ret$ped_ids), function(i) {
#    tibble(ped_id = ret$ped_ids[[i]], 
#           from = ret$edgelists[[i]][, 1], 
#           to = ret$edgelists[[i]][, 2])
#  }))
#  #d_edges
#  
#  d_indv <- bind_rows(lapply(seq_along(ret$ped_ids), function(i) {
#    tibble(pid = ret$pids[[i]], is_female = ret$is_female[[i]], haplotype = ret$haplotypes[[i]])
#  }))
#  #d_indv
#  
#  d <- d_edges %>% 
#    left_join(d_indv, by = c("from" = "pid")) %>% 
#    rename(from_haplotype = haplotype,
#           from_is_female = is_female) %>% 
#    left_join(d_indv, by = c("to" = "pid")) %>% 
#    rename(to_haplotype = haplotype,
#           to_is_female = is_female)
#  
#  return(d)
#}

#' @importFrom magrittr "%>%"
#' @importFrom dplyr mutate
#' @export
get_nodes_edges <- function(x, ...) {
  if (!is(x, "mitolina_pedigreelist")) stop("x must be a mitolina_pedigreelist object")
  
  ret <- get_pedigrees_tidy(x)
  
  d_edges <- dplyr::bind_rows(lapply(seq_along(ret$ped_ids), function(i) {
    tibble(from = ret$edgelists[[i]][, 1], 
           to = ret$edgelists[[i]][, 2])
  })) %>%
  dplyr::mutate(from = as.character(from),
         to = as.character(to))
  #d_edges
  
  d_indv <- dplyr::bind_rows(lapply(seq_along(ret$ped_ids), function(i) {
    tibble(name = ret$pids[[i]], 
           gens_from_final = ret$generation[[i]], 
           ped_id = ret$ped_ids[[i]], 
           sex = factor(ifelse(ret$is_female[[i]], "Female", "Male"), levels = c("Female", "Male")), 
           haplotype = ret$haplotypes[[i]])
  })) %>%
  dplyr::mutate(name = as.character(name))
  #d_indv
      
  return(list(nodes = d_indv, edges = d_edges))
}


#' @importFrom tidygraph as_tbl_graph tbl_graph
#' @export
as_tbl_graph.mitolina_pedigreelist <- function(x, ...) {
  if (!is(x, "mitolina_pedigreelist")) stop("x must be a mitolina_pedigreelist object")
  
  VE <- get_nodes_edges(x)
  
  g <- tidygraph::tbl_graph(nodes = VE$nodes, edges = VE$edges)
    
  return(g)
}


  
# mark_pids vector: use mark_color (with reuse)
#' @export  
plot.mitolina_pedigree <-
  function(x, ids = TRUE, haplotypes = FALSE, locus_sep = " ", mark_pids = NULL, label_color = "black", node_color = "lightgray", mark_color = "orange", ...) {
  
    if (!is(x, "mitolina_pedigree")) stop("x must be a mitolina_pedigree object")
    
    x_pids <- get_pids_in_pedigree(x)
    x_is_females <- get_is_female_in_pedigree(x)
    
    vertex_shapes <- rep("circle", length(x_pids))
    vertex_shapes[!x_is_females] <- "square"

    vertex_label <- rep("", length(x_pids))
    
    if (ids) {
      vertex_label <- x_pids
    }
    
    if (haplotypes) {      
      haps <- get_haplotypes_in_pedigree(x)
      
      vertex_label <- unlist(lapply(seq_along(haps), function(h_i) {
        h <- haps[[h_i]]
        prefix <- ""
        
        if (ids) {
          prefix <- paste0(x_pids[h_i], ": ")
        }
        
        paste0(strwrap(paste0(prefix, paste0(h, collapse = locus_sep)), width = 15), collapse = "\n")
      }))
    }
    
    vertex_colors <- rep(node_color, length(vertex_label))
    if (!is.null(mark_pids)) {
      #vertex_colors[x_pids %in% mark_pids] <- mark_color
      
      if (length(mark_color) == 1L) {
        #mark_color <- rep(mark_color, length(mark_pids))
        mark_color <- rep(mark_color, length(mark_pids))
      } else if (length(mark_color) != length(mark_pids)) {
        stop("Expected mark_color of length 1 or same length as mark_pids")
      }
      
      for (m_id in seq_along(mark_pids)) {
        vertex_colors[which(x_pids == mark_pids[m_id])] <- mark_color[m_id]
      }
    }
    
    vertex_colors[!x_is_females] <- "orange"
    
    g <- pedigree_as_igraph(x)
    igraph::V(g)$color <- vertex_colors
    
    
    old_mar <- par("mar")
    par(mar = c(0, 0, 0, 0))        
    igraph::plot.igraph(g, 
                        vertex.shape = vertex_shapes,
                        vertex.label = vertex_label, 
                        vertex.label.cex = 0.75, 
                        vertex.label.color = label_color, 
                        layout = igraph::layout_as_tree(graph = g),
                        ...)
    par(mar = old_mar)
    
    return(invisible(NULL))
    #eturn(g)
  }
  


