# ---- Load raw data ----
all_edges <- read.table("Cannes2013_multiplex_edges.txt",
                        col.names = c("layer", "source", "target", "weight"),
                        header = FALSE)
cat("Edge data loaded\n")

# Layers: 1 = Retweet, 2 = Mention, 3 = Reply
# ---- Separate by layer and convert to undirected ----
rt_edges <- all_edges %>%
  filter(layer == 1) %>%
  mutate(p1 = pmin(source, target), p2 = pmax(source, target)) %>%
  group_by(p1, p2) %>%
  summarise(weight = sum(weight), .groups = "drop") %>%
  rename(source = p1, target = p2)

mt_edges <- all_edges %>%
  filter(layer == 2) %>%
  mutate(p1 = pmin(source, target), p2 = pmax(source, target)) %>%
  group_by(p1, p2) %>%
  summarise(weight = sum(weight), .groups = "drop") %>%
  rename(source = p1, target = p2)

cat("Edge data separated by layer\n")

# ---- Find common nodes ----
rt_nodes     <- unique(c(rt_edges$source, rt_edges$target))
mt_nodes     <- unique(c(mt_edges$source, mt_edges$target))
common_nodes <- intersect(rt_nodes, mt_nodes)

cat(sprintf("Nodes common to both Retweet and Mention layers: %d\n", length(common_nodes)))

# ---- Filter to common nodes ----
rt_common_edges <- rt_edges %>%
  filter(source %in% common_nodes & target %in% common_nodes)

mt_common_edges <- mt_edges %>%
  filter(source %in% common_nodes & target %in% common_nodes)

common_nodes_df <- data.frame(id = common_nodes)

# ---- Build graphs ----
g_rt <- graph_from_data_frame(
  d        = rt_common_edges[, c("source", "target")],
  directed = FALSE,
  vertices = common_nodes_df
)
g_rt <- simplify(g_rt)

g_mt <- graph_from_data_frame(
  d        = mt_common_edges[, c("source", "target")],
  directed = FALSE,
  vertices = common_nodes_df
)
g_mt <- simplify(g_mt)

# ---- Remove isolated nodes ----
isolated_nodes <- union(
  V(g_rt)[degree(g_rt) < 1]$name,
  V(g_mt)[degree(g_mt) < 1]$name
)

g_rt <- delete_vertices(g_rt, isolated_nodes)
g_mt <- delete_vertices(g_mt, isolated_nodes)

# ---- Align vertex ordering ----
desired_order <- sort(V(g_rt)$name)
g_rt <- permute(g_rt, match(desired_order, V(g_rt)$name))
g_mt <- permute(g_mt, match(desired_order, V(g_mt)$name))

cat(sprintf("Final network size: %d nodes\n", vcount(g_rt)))
cat(sprintf("Vertex sets identical: %s\n", setequal(V(g_rt)$name, V(g_mt)$name)))

# ---- Extract adjacency matrices ----
A1 <- as_adjacency_matrix(g_rt, sparse = TRUE)
A2 <- as_adjacency_matrix(g_mt, sparse = TRUE)

cat(sprintf("A1 (Retweet) dimensions: %d x %d\n", nrow(A1), ncol(A1)))
cat(sprintf("A2 (Mention) dimensions: %d x %d\n", nrow(A2), ncol(A2)))