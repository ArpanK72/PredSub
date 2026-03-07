sourceCpp("process_simplices.cpp")

# ---- Load raw data ----
nverts    <- fread("~/Data/DBLP/coauth-DBLP-full-nverts.txt", header = FALSE)$V1
simplices <- fread("~/Data/DBLP/coauth-DBLP-full-simplices.txt", header = FALSE)$V1
times     <- fread("~/Data/DBLP/coauth-DBLP-full-times.txt", header = FALSE)$V1

# ---- Create edge list ----
dblp_df <- create_edge_list_cpp(nverts, simplices, times)
cat("Edge list created\n")
cat("Total edges:", nrow(dblp_df), "\n")

# ---- Define time periods ----
era1 <- c(2011, 2014)
era2 <- c(2015, 2018)

# ---- Filter and aggregate by era ----
era_1_df <- dblp_df %>%
  filter(year >= era1[1] & year <= era1[2]) %>%
  mutate(p1 = pmin(author_1, author_2), p2 = pmax(author_1, author_2)) %>%
  group_by(p1, p2) %>%
  summarise(weight = n(), .groups = "drop") %>%
  rename(author_1 = p1, author_2 = p2)

era_2_df <- dblp_df %>%
  filter(year >= era2[1] & year <= era2[2]) %>%
  mutate(p1 = pmin(author_1, author_2), p2 = pmax(author_1, author_2)) %>%
  group_by(p1, p2) %>%
  summarise(weight = n(), .groups = "drop") %>%
  rename(author_1 = p1, author_2 = p2)

cat(sprintf("Unique collaborations in Era 1 (2011-2014): %d\n", nrow(era_1_df)))
cat(sprintf("Unique collaborations in Era 2 (2015-2018): %d\n", nrow(era_2_df)))

# ---- Find common node set ----
nodes_era_1    <- union(era_1_df$author_1, era_1_df$author_2)
nodes_era_2    <- union(era_2_df$author_1, era_2_df$author_2)
common_node_set <- intersect(nodes_era_1, nodes_era_2)

cat(sprintf("Total unique authors active in both eras: %d\n", length(common_node_set)))

# ---- Filter to common nodes ----
era_1_filtered_df <- era_1_df %>%
  filter(author_1 %in% common_node_set & author_2 %in% common_node_set)

era_2_filtered_df <- era_2_df %>%
  filter(author_1 %in% common_node_set & author_2 %in% common_node_set)

# ---- Build graphs ----
g1 <- graph_from_data_frame(d = era_1_filtered_df, directed = FALSE, vertices = common_node_set)
g1 <- simplify(g1)

g2 <- graph_from_data_frame(d = era_2_filtered_df, directed = FALSE, vertices = common_node_set)
g2 <- simplify(g2)

# ---- Remove low degree nodes ----
isolated_nodes <- union(
  V(g1)[degree(g1) < 5]$name,
  V(g2)[degree(g2) < 5]$name
)

g1 <- delete_vertices(g1, isolated_nodes)
g2 <- delete_vertices(g2, isolated_nodes)

# ---- Align vertex ordering ----
desired_order <- sort(V(g1)$name)
g1 <- permute(g1, match(desired_order, V(g1)$name))
g2 <- permute(g2, match(desired_order, V(g2)$name))

cat(sprintf("Final network size: %d nodes\n", vcount(g1)))
cat(sprintf("Vertex sets identical: %s\n", setequal(V(g1)$name, V(g2)$name)))

# ---- Extract adjacency matrices ----
A1 <- as_adjacency_matrix(g1, sparse = TRUE)
A2 <- as_adjacency_matrix(g2, sparse = TRUE)

cat(sprintf("A1 dimensions: %d x %d\n", nrow(A1), ncol(A1)))
cat(sprintf("A2 dimensions: %d x %d\n", nrow(A2), ncol(A2)))

cat("Adjacency matrices created\n")
