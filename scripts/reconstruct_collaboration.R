### SOME DEPENDENCY CHANGE OR UPDATE IN THE CONDA ENVIRONMENT HAS GONE UNNOTICED.
### EXACT REPRODUCIBILITY DOES NOT WORK ANYMORE. TO REPRODUCE REGRESSION RESULTS, USE HARD COPIES OF MODELS.

#workspace
library(igraph)
devtools::install_version("Bergm", version = "5.0.6", repos = "https://cloud.r-project.org") #not available on conda-forge 
library(Bergm)
setwd("/m/cs/scratch/ecanet/malkama5/mediaAccess/replication")

#attributes
a = read.table("attribute.csv", header = TRUE, sep = ";")
a$type[a$type == "corporate_private"] = "corporate"
a$type[a$type == "corporate_public"] = "corporate"
a$type[a$type == "city"] = "governmental"
a$type[a$type == "government_agency"] = "governmental"
a$type[a$type == "government_ministry"] = "governmental"
a$type[a$type == "political_party"] = "political_party"
a$type[a$type == "research_institute"] = "research"
a$type[a$type == "university"] = "research"
a$type[a$type == "trade_union"] = "interest_group"

###
#C0
d = read.table("collaboration0.csv", header = TRUE, sep = ";")
d = subset(d, to %in% a$actor[a$roster0 == 1])
g = simplify(graph_from_data_frame(d, directed = TRUE), remove.loops = TRUE, remove.multiple = TRUE)

#spurious vertices
or = V(g)$name[degree(g, mode = "out") > ((vcount(g) - 1) * (1 / 2))]
td = list()
for (e in E(g)){
    so = ends(g, e)[1, 1]
    if (so %in% or){
        td = append(td, e)}}
g = delete_edges(g, td)

#influence
d = read.table("influence0.csv", header = TRUE, sep = ";")
d = subset(d, to %in% a$actor[a$roster0 == 1])
f = simplify(graph_from_data_frame(d, directed = TRUE), remove.loops = TRUE, remove.multiple = TRUE)

#matrix
d = as_adjacency_matrix(as.undirected(g, mode = "collapse"), sparse = FALSE)
n = V(g)$name[degree(g, mode = "out") == 0]
length(n) / length(a$actor[a$roster0 == 1]) #share of spurious or missing vertices
d[n, n] = NA

#drop spurious or missing vertices
p_n = V(g)$name[degree(g, mode = "out") != 0]
p_d = d[p_n, p_n]
p_nw = network(p_d, directed = FALSE, loops = FALSE, matrix.type = "adjacency")

#network
nw = network(d, directed = FALSE, loops = FALSE, matrix.type = "adjacency")
nw %v% "type" = merge(x = data.frame("vertex" = rownames(as.matrix(nw))), y = data.frame("vertex" = a$actor, "variable" = a$type), by = "vertex", all.x = TRUE, sort = FALSE)$variable
nw %v% "influence" = merge(x = data.frame("vertex" = rownames(as.matrix(nw))), y = data.frame("vertex" = V(f)$name, "variable" = degree(f, mode = "in")), by = "vertex", all.x = TRUE, sort = FALSE)$variable

#pre-imputation descriptives
round(c(network.size(p_nw), network.edgecount(p_nw), network.density(p_nw), mean(sna::degree(p_nw)), sd(sna::degree(p_nw)), sna::gtrans(p_nw, mode = "graph")), 2)
round(c(network.size(nw), network.edgecount(nw), network.density(nw), mean(sna::degree(nw)), sd(sna::degree(nw)), sna::gtrans(nw, mode = "graph")), 2)

#model
fo = nw ~ edges + gwesp(decay = 2.0, fixed = TRUE) + nodecov("influence") + nodefactor("type", levels = "governmental") + nodematch("type", diff = TRUE, levels = c("corporate", "governmental", "interest_group", "research")) + gwdegree(decay = 1.0, fixed = TRUE)
pM = c(-4, rep(0, 8))
pS = diag(4, length(pM))

m = bergmM(fo
, nImp = 500
, gamma = 0.1
, nchains = length(pM) * 4
, burn.in = 500
, main.iters = 3000
, aux.iters = 3000
, prior.mean = pM
, prior.sigma = pS
, missingUpdate = network.naedgecount(nw) * 4
, constraints = ~ bd(minout = min(sna::degree(nw, cmode = "outdegree")))
, seed = 0)

saveRDS(m, "rc_m0.rds")
svg("rc_m0_bgof.svg", width = 30, height = 10)
bgof(m)
dev.off()

#reconstruct
#m = readRDS("rc_m0.rds")
m = readRDS("m0_rc.rds") #hard copy
rc = matrix(0, nrow = nrow(as.matrix(nw)), ncol = ncol(as.matrix(nw)))
dimnames(rc) = list(rownames(as.matrix(nw)), colnames(as.matrix(nw)))
for (i in m$impNets){
    rc = rc + as.matrix(i)}
rc = rc / length(m$impNets)
rc = ifelse(rc > 0.5, 1, 0)
g = graph_from_adjacency_matrix(rc, mode = "upper", diag = FALSE)
V(g)$bonpow = power_centrality(g, exponent = 1, sparse = FALSE) #add bonacich

E(g)$layer = 0
write_graph(g, file = "c0.xml", format = "graphml")

#descriptives
rc = network(rc, directed = FALSE, loops = FALSE, matrix.type = "adjacency")
round(c(network.size(rc), network.edgecount(rc), network.density(rc), mean(sna::degree(rc)), sd(sna::degree(rc)), sna::gtrans(rc, mode = "graph")), 2)

###
#C1
d = read.table("collaboration1.csv", header = TRUE, sep = ";")
d = subset(d, to %in% a$actor[a$roster1 == 1])
g = simplify(graph_from_data_frame(d, directed = TRUE), remove.loops = TRUE, remove.multiple = TRUE)

#spurious vertices
or = V(g)$name[degree(g, mode = "out") > ((vcount(g) - 1) * (1 / 2))]
td = list()
for (e in E(g)){
    so = ends(g, e)[1, 1]
    if (so %in% or){
        td = append(td, e)}}
g = delete_edges(g, td)

#influence
d = read.table("influence1.csv", header = TRUE, sep = ";")
d = subset(d, to %in% a$actor[a$roster1 == 1])
f = simplify(graph_from_data_frame(d, directed = TRUE), remove.loops = TRUE, remove.multiple = TRUE)

#matrix
d = as_adjacency_matrix(as.undirected(g, mode = "collapse"), sparse = FALSE)
n = V(g)$name[degree(g, mode = "out") == 0]
length(n) / length(a$actor[a$roster1 == 1]) #share of spurious or missing vertices
d[n, n] = NA

#drop spurious or missing vertices
p_n = V(g)$name[degree(g, mode = "out") != 0]
p_d = d[p_n, p_n]
p_nw = network(p_d, directed = FALSE, loops = FALSE, matrix.type = "adjacency")

#network
nw = network(d, directed = FALSE, loops = FALSE, matrix.type = "adjacency")
nw %v% "type" = merge(x = data.frame("vertex" = rownames(as.matrix(nw))), y = data.frame("vertex" = a$actor, "variable" = a$type), by = "vertex", all.x = TRUE, sort = FALSE)$variable
nw %v% "influence" = merge(x = data.frame("vertex" = rownames(as.matrix(nw))), y = data.frame("vertex" = V(f)$name, "variable" = degree(f, mode = "in")), by = "vertex", all.x = TRUE, sort = FALSE)$variable

#pre-imputation descriptives
round(c(network.size(p_nw), network.edgecount(p_nw), network.density(p_nw), mean(sna::degree(p_nw)), sd(sna::degree(p_nw)), sna::gtrans(p_nw, mode = "graph")), 2)
round(c(network.size(nw), network.edgecount(nw), network.density(nw), mean(sna::degree(nw)), sd(sna::degree(nw)), sna::gtrans(nw, mode = "graph")), 2)

#model
fo = nw ~ edges + gwesp(decay = 2.0, fixed = TRUE) + nodecov("influence") + nodefactor("type", levels = "governmental") + nodematch("type", diff = TRUE, levels = c("corporate", "governmental", "interest_group", "research")) + gwdegree(decay = 1.0, fixed = TRUE)
pM = c(-4, rep(0, 8))
pS = diag(4, length(pM))

m = bergmM(fo
, nImp = 500
, gamma = 0.1
, nchains = length(pM) * 4,
, burn.in = 500
, main.iters = 3000
, aux.iters = 3000
, prior.mean = pM
, prior.sigma = pS
, missingUpdate = network.naedgecount(nw) * 4
, constraints = ~ bd(minout = min(sna::degree(nw, cmode = "outdegree")))
, seed = 0)

saveRDS(m, "rc_m1.rds")
svg("rc_m1_bgof.svg", width = 30, height = 10)
bgof(m)
dev.off()

#reconstruct
#m = readRDS("rc_m1.rds")
m = readRDS("m1_rc.rds") #hard copy
rc = matrix(0, nrow = nrow(as.matrix(nw)), ncol = ncol(as.matrix(nw)))
dimnames(rc) = list(rownames(as.matrix(nw)), colnames(as.matrix(nw)))
for (i in m$impNets){
    rc = rc + as.matrix(i)}
rc = rc / length(m$impNets)
rc = ifelse(rc > 0.5, 1, 0)
g = graph_from_adjacency_matrix(rc, mode = "upper", diag = FALSE)
V(g)$bonpow = power_centrality(g, exponent = 1, sparse = FALSE) #add bonacich

E(g)$layer = 1
write_graph(g, file = "c1.xml", format = "graphml")

#descriptives
rc = network(rc, directed = FALSE, loops = FALSE, matrix.type = "adjacency")
round(c(network.size(rc), network.edgecount(rc), network.density(rc), mean(sna::degree(rc)), sd(sna::degree(rc)), sna::gtrans(rc, mode = "graph")), 2)
