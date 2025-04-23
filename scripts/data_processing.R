#########
### 0 ###
#########

#workspace
library(igraph)
setwd("/m/cs/scratch/ecanet/malkama5/mediaAccess/replication/")

#data
dt = read.table("dt.csv", header = TRUE, sep = ";")
dt = subset(dt, period == 0)

#brokerage
g = read.graph("c0.xml", format = "graphml")
a = as_adjacency_matrix(g, attr = NULL, sparse = FALSE)
a = a[dt$actor, dt$actor]
b = sna::brokerage(a, dt$coalit)
dt$coordi = scales::rescale(as.data.frame(b$raw.nli)$w_I, to = c(0, 1)) #coordinator
dt$liaiso = scales::rescale(as.data.frame(b$raw.nli)$b_O, to = c(0, 1)) #liaison

#brokerage filter
g = read.graph("rc0_f.xml", format = "graphml")
a = as_adjacency_matrix(g, attr = NULL, sparse = FALSE)
a = a[dt$actor, dt$actor]
b = sna::brokerage(a, dt$coalit)
dt$coordi_f = scales::rescale(as.data.frame(b$raw.nli)$w_I, to = c(0, 1)) #coordinator
dt$liaiso_f = scales::rescale(as.data.frame(b$raw.nli)$b_O, to = c(0, 1)) #liaison

#X activity
k = paste(scan("keys.csv", character(), encoding = "UTF-8"), collapse = "|")
d = readRDS("twitter2013_16.rds") #not provided for privacy reasons
d = subset(d, type == "retweeted" & !is.na(ref_org) & level != "3" & ref_level != "3" & created_at >= as.Date("2013-01-01") & created_at <= as.Date("2016-12-31"))
d = subset(d, ref_org %in% dt$actor)
d = d[grep(k, d$text, ignore.case = TRUE), 1:ncol(d)]
d$counter = 1
d = aggregate(d[which(colnames(d) == "counter")], list(d$org), sum)
colnames(d) = c("actor", "twact")
dt = merge(x = dt, y = d, by = "actor", sort = FALSE, all.x = TRUE)
dt$twact[is.na(dt$twact)] = 0
dt$twact_log = scales::rescale(log(1 + dt$twact), to = c(0, 1))

#variables
dt$type[dt$type == "climate_sceptic"] = "non_governmental"
dt$type[dt$type == "corporate_public"] = "corporate"
dt$type[dt$type == "corporate_private"] = "corporate"
dt$type[dt$type == "trade_union"] = "interest_group"
dt$type[dt$type == "government_ministry"] = "governmental"
dt$type[dt$type == "government_agency"] = "governmental"
dt$type[dt$type == "city"] = "governmental"
dt$type[dt$type == "research_institute"] = "research"
dt$type[dt$type == "university"] = "research"
dt$period = factor(dt$period)
dt$actor = factor(dt$actor)
dt$type = factor(dt$type)
dt$coalit = factor(dt$coalit)
d0 = dt



#########
### 1 ###
#########

#workspace
library(igraph)

#data
dt = read.table("dt.csv", header = TRUE, sep = ";")
dt = subset(dt, period == 1)

#brokerage
g = read.graph("c1.xml", format = "graphml")
a = as_adjacency_matrix(g, attr = NULL, sparse = FALSE)
a = a[dt$actor, dt$actor]
b = sna::brokerage(a, dt$coalit)
dt$coordi = scales::rescale(as.data.frame(b$raw.nli)$w_I, to = c(0, 1)) #coordinator
dt$liaiso = scales::rescale(as.data.frame(b$raw.nli)$b_O, to = c(0, 1)) #liaison

#brokerage filter
g = read.graph("rc1_f.xml", format = "graphml")
a = as_adjacency_matrix(g, attr = NULL, sparse = FALSE)
a = a[dt$actor, dt$actor]
b = sna::brokerage(a, dt$coalit)
dt$coordi_f = scales::rescale(as.data.frame(b$raw.nli)$w_I, to = c(0, 1)) #coordinator
dt$liaiso_f = scales::rescale(as.data.frame(b$raw.nli)$b_O, to = c(0, 1)) #liaison

#online activity
k = paste(scan("keys.csv", character(), encoding = "UTF-8"), collapse = "|")
d = readRDS("twitter2017_22.rds") #not provided for privacy reasons
d = subset(d, type == "retweeted" & !is.na(ref_org) & level != "3" & ref_level != "3" & created_at >= as.Date("2017-01-01") & created_at <= as.Date("2020-12-31"))
d = subset(d, ref_org %in% dt$actor)
d = d[grep(k, d$text, ignore.case = TRUE), 1:ncol(d)]
d$counter = 1
d = aggregate(d[which(colnames(d) == "counter")], list(d$org), sum)
colnames(d) = c("actor", "twact")
dt = merge(x = dt, y = d, by = "actor", sort = FALSE, all.x = TRUE)
dt$twact[is.na(dt$twact)] = 0
dt$twact_log = scales::rescale(log(1 + dt$twact), to = c(0, 1))

#variables
dt$type[dt$type == "climate_sceptic"] = "non_governmental"
dt$type[dt$type == "corporate_public"] = "corporate"
dt$type[dt$type == "corporate_private"] = "corporate"
dt$type[dt$type == "trade_union"] = "interest_group"
dt$type[dt$type == "government_ministry"] = "governmental"
dt$type[dt$type == "government_agency"] = "governmental"
dt$type[dt$type == "city"] = "governmental"
dt$type[dt$type == "research_institute"] = "research"
dt$type[dt$type == "university"] = "research"
dt$period = factor(dt$period)
dt$actor = factor(dt$actor)
dt$type = factor(dt$type)
dt$coalit = factor(dt$coalit)
d1 = dt



###########
## JOINT ##
###########

#store
dt = rbind(d0, d1)
saveRDS(dt, paste0(path, "jointData.rds"))
write.table(dt, paste0(path, "jointData.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ";")
