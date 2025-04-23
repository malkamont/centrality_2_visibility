##############
## DESCRIBE ##
##############

#workspace
library(ggplot2)
setwd("/m/cs/scratch/ecanet/malkama5/mediaAccess/replication/")

#data
dt = readRDS("jointData.rds")

#continuous-x-categorical
vr = dt[c(which(names(dt) == "vtot"), which(names(dt) == "popula"), which(names(dt) == "closen"), which(names(dt) == "eigenv"), which(names(dt) == "betwee"), which(names(dt) == "coalea"), which(names(dt) == "actot"), which(names(dt) == "procli"), which(names(dt) == "twact"), which(names(dt) == "bonpow"), which(names(dt) == "constr"))]

cat = rbind(aggregate(vr, list(dt$period), mean), aggregate(vr, list(dt$type), mean), aggregate(vr, list(dt$coalit), mean))
cat[2:ncol(cat)] = round(cat[2:ncol(cat)], 2)
nm = c("Category", "Visibility", "Indegree", "Closeness", "Eigenvector", "Betweenness", "Community eigentrust", "Autocorrelation", "Climate progressiveness", "X activity", "Bonacich's power", "Burt's constraint")
colnames(cat) = nm
rownames(cat) = c("Pre-Paris", "Post-Paris", "Corporation", "Governmental", "Interest group", "Non-governmental", "Political party", "Research", "Economy coalition", "Ecology coalition")
con = round(psych::describe(vr), 2)
rownames(con) = nm[c(2:length(nm))]
vr = rbind(cat[2:ncol(cat)], t(con[c(3:4, 8:9, 11:12)]))
print(xtable::xtable(t(vr)), type = "html", file = "descriptive.html")

#visibility
cpie = c("#b92e34", "#1b03a3")
tb = data.frame(table(dt$vtot[dt$period == 0])); colnames(tb) = c("vtot", "count")
t0 = merge(x = data.frame("vtot" = 0:max(dt$vtot), "period" = 0), y = tb, by = "vtot", all.x = TRUE, sort = FALSE)
tb = data.frame(table(dt$vtot[dt$period == 1])); colnames(tb) = c("vtot", "count")
t1 = merge(x = data.frame("vtot" = 0:max(dt$vtot), "period" = 1), y = tb, by = "vtot", all.x = TRUE, sort = FALSE)
tb = rbind(t0, t1); tb$period = factor(tb$period)
tb[is.na(tb)] = 0

p = ggplot(tb, aes(x = vtot, y = count, fill = period, color = period)) +
geom_bar(alpha = 1, width = 0.5, stat = "identity", position = position_dodge(width = 1)) +
scale_fill_manual(values = cpie, labels = c("Pre-Paris", "Post-Paris")) +
scale_color_manual(values = cpie, labels = c("Pre-Paris", "Post-Paris")) +
theme(legend.title = element_blank(), text = element_text(size = 8, color = "black", family = "Times"), axis.text = element_text(color = "black", family = "Times"), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) +
labs(x = "Visibility", y = "Count")
ggsave(plot = p, filename = "visibility.png", width = 5, height = 5, dpi = 1000)

#correlation
vr = dt[c(which(names(dt) == "vtot"), which(names(dt) == "popula"), which(names(dt) == "closen"), which(names(dt) == "eigenv"), which(names(dt) == "betwee"), which(names(dt) == "coalea"), which(names(dt) == "actot"), which(names(dt) == "procli"), which(names(dt) == "twact_log"), which(names(dt) == "bonpow"), which(names(dt) == "constr"))]
names(vr) = c("Visibility", "Indegree [H1]", "Closeness [H2A]", "Eigenvector [H2B]", "Betweenness  [H3]", "Community eigentrust  [H4]", "Autocorrelation", "Climate progressiveness", "X activity [log]", "Bonacich's power", "Burt's constraint")
vr = round(cor(vr, method = "kendall"), 2)
vr[upper.tri(vr, diag = TRUE)] = NA
vr = as.data.frame(as.table(as.matrix(vr)))
colnames(vr) = c("x", "y", "cr")

pcon = ggplot(vr, aes(x, y, fill = cr)) +
geom_tile(color = "black") +
geom_text(aes(label = sprintf("%0.2f", cr)), size = 2, color = "black", family = "Times") +
scale_fill_gradient(low = cpie[2], high = cpie[1], na.value = "black") +
labs(x = NULL, y = NULL) +
theme(text = element_text(color = "black", family = "Times", size = 8), axis.text.x = element_text(color = "black", family = "Times", angle = 30, vjust = 1, hjust = 1), axis.text.y = element_text(color = "black", family = "Times"), legend.position = "none", axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())

eta = function(x, y){
m = mean(x, na.rm = TRUE)
sct = sum((x - m) ^ 2, na.rm = TRUE)
n = table(y)
mk = tapply(x, y, mean)
sce = sum(n * (mk - m) ^ 2)
return(ifelse(sct > 0, sce / sct, 0))}

qlt = list("Period" = dt$period, "Type" = dt$type, "Coalition" = dt$coalit)
qnt = list("Indegree [H1]" = dt$popula, "Closeness [H2A]" = dt$closen, "Eigenvector [H2B]" = dt$eigenv, "Betweenness [H3]" = dt$betwee, "Community eigentrust [H4]" = dt$coalea, "Autocorrelation" = dt$actot, "Climate progressiveness" = dt$procli, "X activity [log]" = dt$twact_log, "Bonacich's power" = dt$bonpow, "Burt's constraint" = dt$constr)

et = round(matrix(mapply(eta, rep(qnt, length(qlt)), rep(qlt, each = length(qnt))), ncol = length(qlt)), 2)
dimnames(et) = list(names(qnt), names(qlt))
et = t(et)
vr = as.data.frame(as.table(as.matrix(et)))
names(vr) = c("x", "y", "et")

pcat = ggplot(vr, aes(x, y, fill = et)) +
geom_tile(color = "black") +
geom_text(aes(label = sprintf("%0.2f", et)), size = 2, color = "black", family = "Times") +
scale_fill_gradient(low = cpie[2], high = cpie[1]) +
labs(x = NULL, y = NULL) +
theme(text = element_text(color = "black", family = "Times", size = 8), axis.text.x = element_text(color = "black", family = "Times", angle = 30, vjust = 1, hjust = 1), axis.text.y = element_text(color = "#ffffff", family = "Times"), legend.position = "none", axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())

ge = ggpubr::ggarrange(pcon, pcat, ncol = 2, nrow = 1, align = "hv") #0.06 medium effect size; 0.14 large effect size
ggsave(plot = ge, filename = "correlation.png", width = 7, height = 3, dpi = 1000)
