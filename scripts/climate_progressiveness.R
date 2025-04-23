############
## MIRT 0 ##
############

#workspace
library(mirt)
setwd("/m/cs/scratch/ecanet/malkama5/mediaAccess/replication/")

#data
bf = read.table("rb0.csv", header = TRUE, sep = ";")
rownames(bf) = bf$actor
bf = bf[1:ncol(bf) - 1]

#item
ds = psych::describe(bf)
ds = ds[order(ds$sd, decreasing = TRUE), 1:ncol(ds)]
bf = as.matrix(bf[rownames(ds)[ds$sd > 0.5]])

#mirt
m = mirt(bf, model = 4, itemtype = "graded", method = "MHRM", SE.type = "MHRM", optimizer = "NR1", technical = list(set.seed = 0))
summary(m, rotate = "oblimin", suppress = 0.4)

xx = paste0(c("natcmp", "secene", "tmeffo"), collapse = "|")
bf = bf[1:nrow(bf), colnames(bf)[grepl(xx, colnames(bf))]]
m = mirt(bf, model = 1, itemtype = "graded", method = "MHRM", SE.type = "MHRM", optimizer = "NR1", technical = list(set.seed = 0))
summary(m, rotate = "oblimin", suppress = 0.4)

#M2(m, QMC = FALSE, type = "C2")
fs = fscores(m, method = "EAP", full.scores.SE = TRUE, QMC = FALSE)
empirical_rxx(fs)

#store
bf = data.frame("actor" = rownames(bf), "procli" = scales::rescale(as.data.frame(fs)$F1, to = c(1, 0)))
#write.table(bf, "procli0.csv", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ";") #for exact reproducibility, use hard copy

############
## MIRT 1 ##
############

#workspace
library(mirt)

#data
bf = read.table("rb1.csv", header = TRUE, sep = ";")
rownames(bf) = bf$actor
bf = bf[1:ncol(bf) - 1]

#item
ds = psych::describe(bf)
ds = ds[order(ds$sd, decreasing = TRUE), 1:ncol(ds)]
bf = as.matrix(bf[rownames(ds)[ds$sd > 0.5]])

#mirt
m = mirt(bf, model = 4, itemtype = "graded", method = "MHRM", SE.type = "MHRM", optimizer = "NR1", technical = list(set.seed = 0))
summary(m, rotate = "oblimin", suppress = 0.4)

xx = paste0(c("natcmp", "secene", "tmeffo"), collapse = "|")
bf = bf[1:nrow(bf), colnames(bf)[grepl(xx, colnames(bf))]]
m = mirt(bf, model = 1, itemtype = "graded", method = "MHRM", SE.type = "MHRM", optimizer = "NR1", technical = list(set.seed = 0))
summary(m, rotate = "oblimin", suppress = 0.4)

#M2(m, QMC = FALSE, type = "C2")
fs = fscores(m, method = "EAP", full.scores.SE = TRUE, QMC = FALSE)
empirical_rxx(fs)

#store
bf = data.frame("actor" = rownames(bf), "procli" = scales::rescale(as.data.frame(fs)$F1, to = c(1, 0)))
#write.table(bf, "procli1.csv", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ";") #for exact reproducibility, use hard copy
