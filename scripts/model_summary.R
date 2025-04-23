#workspace
setwd("/m/cs/scratch/ecanet/malkama5/mediaAccess/replication/")

#model summary and evaluation
m0 = readRDS("m0.rds")
m1 = readRDS("m1.rds")
m2 = readRDS("m2.rds")
m3 = readRDS("m3.rds")
m4 = readRDS("m4.rds")
m5 = readRDS("m5.rds")
m6 = readRDS("m6.rds")
m7 = readRDS("m7.rds")
m1i = readRDS("m1i.rds")
m2i = readRDS("m2i.rds")
m3i = readRDS("m3i.rds")
m4i = readRDS("m4i.rds")
m5i = readRDS("m5i.rds")
m6i = readRDS("m6i.rds")
m7i = readRDS("m7i.rds")
m0z = readRDS("m0z.rds")
m0j = readRDS("m0j.rds")
mall = readRDS("mall.rds")

lb = c("m0", "m1", "m2", "m3", "m4", "m5", "mall", "m6", "m7", "m1i", "m2i", "m3i", "m4i", "m5i", "m6i", "m7i", "m0j") #j
ms = list(m0, m1, m2, m3, m4, m5, mall, m6, m7, m1i, m2i, m3i, m4i, m5i, m6i, m7i, m0j)
sjPlot::tab_model(ms, dv.labels = lb, show.intercept = TRUE, digits.rsq = 2, pred.labels = NULL, string.est = "ICC", file = "modelSummary.html")

lm = readRDS("loo.rds") #elpd_diff / se_diff > 5?
bf = bayestestR::bayesfactor_models(m0, m1, m2, m3, m4, m5, mall, m6, m7, m1i, m2i, m3i, m4i, m5i, m6i, m7i, m0z) #z
ft = merge(x = bf, y = lm$diffs, by = "row.names", sort = FALSE, all.x = TRUE)
colnames(ft)[1:2] = c("model", "call")
ft[3:ncol(ft)] = round(ft[3:ncol(ft)], 2)
print(xtable::xtable(ft), type = "html", file = "modelEvaluation.html")

m0 = readRDS("f_m0.rds")
m1 = readRDS("f_m1.rds")
m2 = readRDS("f_m2.rds")
m3 = readRDS("f_m3.rds")
m4 = readRDS("f_m4.rds")
m5 = readRDS("f_m5.rds")
m6 = readRDS("f_m6.rds")
m7 = readRDS("f_m7.rds")
mall = readRDS("f_mall.rds")

lb = c("f_m0", "f_m1", "f_m2", "f_m3", "f_m4", "f_m5", "mall", "f_m6", "f_m7")
ms = list(m0, m1, m2, m3, m4, m5, mall, m6, m7)
sjPlot::tab_model(ms, dv.labels = lb, show.intercept = TRUE, digits.rsq = 2, pred.labels = NULL, string.est = "ICC", file = "modelSummaryFilter.html")

lm = readRDS("f_loo.rds") #elpd_diff / se_diff > 5?
bf = bayestestR::bayesfactor_models(m0, m1, m2, m3, m4, m5, mall, m6, m7)
ft = merge(x = bf, y = lm$diffs, by = "row.names", sort = FALSE, all.x = TRUE)
colnames(ft)[1:2] = c("model", "call")
ft[3:ncol(ft)] = round(ft[3:ncol(ft)], 2)
print(xtable::xtable(ft), type = "html", file = "modelEvaluationFilter.html")

m0 = readRDS("c_m0.rds")
m1 = readRDS("c_m1.rds")
m2 = readRDS("c_m2.rds")
m3 = readRDS("c_m3.rds")
m4 = readRDS("c_m4.rds")
m5 = readRDS("c_m5.rds")
mall = readRDS("c_mall.rds")

lb = c("c_m0", "c_m1", "c_m2", "c_m3", "c_m4", "c_m5", "c_mall")
ms = list(m0, m1, m2, m3, m4, m5, mall)
sjPlot::tab_model(ms, dv.labels = lb, show.intercept = TRUE, digits.rsq = 2, pred.labels = NULL, string.est = "ICC", file = "modelSummaryCheck.html")

lm = readRDS("c_loo.rds") #elpd_diff / se_diff > 5?
bf = bayestestR::bayesfactor_models(m0, m1, m2, m3, m4, m5, mall)
ft = merge(x = bf, y = lm$diffs, by = "row.names", sort = FALSE, all.x = TRUE)
colnames(ft)[1:2] = c("model", "call")
ft[3:ncol(ft)] = round(ft[3:ncol(ft)], 2)
print(xtable::xtable(ft), type = "html", file = "modelEvaluationCheck.html")

#############
## EFFECTS ##
#############

#workspace
library(brms)
library(dplyr)
library(ggplot2)

#summary function
brmsSummary = function(model, formatOptions = list(digits = 2, nsmall = 2), round = 2, seed = 0){
    fe = brms::fixef(model)
    partables_formatted = data.frame("Covariate" = rownames(fe), do.call(format, list(x = round(fe, round), unlist(formatOptions))))
    rownames(partables_formatted) = NULL
    colnames(partables_formatted)[4:5] = c("l-95 CI", "u-95 CI")
    #hypotheses
    testnams = mapply(function(coefname, sign)
    ifelse(sign > 0, paste0(coefname, ">0"), paste0(coefname, "<0")), rownames(fe), sign(fe[, "Estimate"]))
    hypetests = brms::hypothesis(x = model, testnams, seed = seed)
    ER = hypetests$hypothesis$Evid.Ratio
    partables_formatted$pvals = ER / (ER + 1)
    partables_formatted$pvals[is.infinite(ER)] = 1
    partables_formatted$pvals = do.call(format, list(x = round(partables_formatted$pvals, round), formatOptions[[1]]))
    colnames(partables_formatted)[ncol(partables_formatted)] = "B<>0"
    partables_formatted[, "Evi.Ratio"] = do.call(format, list(x = round(ER, round), formatOptions[[1]]))
    partables_formatted$Notable = hypetests$hypothesis$Star
    partables_formatted = cbind(partables_formatted, bayestestR::equivalence_test(model)[5:6])
    colnames(partables_formatted)[9:10] = c("ROPE Perc.", "ROPE Equi.")
    rownames(partables_formatted) = partables_formatted$Covariate
    return(partables_formatted[2:ncol(partables_formatted)])}

#notable
m = readRDS("m1.rds")
n = brmsSummary(m)
print(xtable::xtable(n), type = "html", file = "hypothesesModel1.html")

#quartile
m = readRDS("m1.rds")
cf = conditional_effects(m, effects = "popula")
x = data.frame(
cf$popula$popula,
cf$popula$type,
cf$popula$estimate__,
cf$popula$lower__,
cf$popula$upper__)
x[-2] = round(x[-2], 2)
x = x %>% mutate(qrt = ntile(cf$popula$popula, 4))
x = aggregate(x[c(1, 3:5)], list(x$qrt), mean)
x$type= "corporate"
cen = x

cf = conditional_effects(m, effects = "type")
x = data.frame(
cf$type$popula,
cf$type$type,
cf$type$estimate__,
cf$type$lower__,
cf$type$upper__)
x[-2] = round(x[-2], 2)
x = x[c(2, 1, 3:5)]
x$group = NA
typ = x

colnames(cen) = c("quartile", "popula", "est", "lower", "upper", "type")
cen = cen[c(1, 6, 2:5)]
names(typ) = c("type", "popula", "est", "lower", "upper", "quartile")
typ = typ[c(6, 1:5)]
x = rbind(cen, typ)
rownames(x) = 1:nrow(x)
print(xtable::xtable(x), type = "html", file = "effectsModel1.html")

##################
## PLOT EFFECTS ##
##################

#theme
th = theme(text = element_text(size = 10, color = "black", family = "Times"), axis.text.x = element_text(color = "black", family = "Times", angle = 30, hjust = 1, vjust = 1), axis.text.y = element_text(color = "black", family = "Times"), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), legend.position = "right", legend.title = element_blank(), plot.title = element_text(hjust = 0.0))
ymin = 0
ymax = 5

###########
#effects m1
m = readRDS("m1.rds")
cf = conditional_effects(m)

#popula
p_popula_1 = ggplot(cf$popula, aes(x = popula, y = estimate__)) +
geom_smooth(method = "loess", color = "#000000") +
geom_ribbon(aes(ymin = lower__, ymax = upper__), linewidth = 0.2, linetype = 1, alpha = 0.1, color = "#000000") +
labs(title = "Popular partner [indegree, H1] [m1]", x = NULL, y = "Visibility") +
th + ylim(ymin, 30) #+ theme(axis.text.y = element_text(color = "#b92e34"), axis.title.y = element_text(color = "#b92e34"))

#type
cf$type$type = c("Corporate", "Governmental", "Interest group", "Non-governmental", "Political party", "Research")
p_type_1 = ggplot(cf$type, aes(x = type, y = estimate__)) +
geom_point(size = 4, color = "#000000") +
geom_errorbar(ymin = cf$type$lower__, ymax = cf$type$upper__, linewidth = 1, width = 0.0) +
labs(title = "Type [m1]", x = NULL, y = "Visibility") +
th + ylim(ymin, 30)

#online activity
p_twact_log_1 = ggplot(cf$twact_log, aes(x = twact_log, y = estimate__)) +
geom_smooth(method = "loess", color = "#000000") +
geom_ribbon(aes(ymin = lower__, ymax = upper__), linewidth = 0.2, linetype = 1, alpha = 0.1, color = "#000000") +
labs(title = "Online activity [log] [m1]", x = NULL, y = "Visibility") +
th + ylim(ymin, ymax)

#period
cf$period$period = factor(c("Pre-Paris", "Post-Paris"), levels = c("Pre-Paris", "Post-Paris"))
p_period_1 = ggplot(cf$period, aes(x = period, y = estimate__)) +
geom_point(size = 4, color = "#000000") +
geom_errorbar(ymin = cf$period$lower__, ymax = cf$period$upper__, linewidth = 1, width = 0.0) +
labs(title = "Period [m1]", x = NULL, y = "Visibility") +
th + ylim(ymin, ymax)

#climate progressiveness
p_procli_1 = ggplot(cf$procli, aes(x = procli, y = estimate__)) +
geom_smooth(method = "loess", color = "#000000") +
geom_ribbon(aes(ymin = lower__, ymax = upper__), linewidth = 0.2, linetype = 1, alpha = 0.1, color = "#000000") +
labs(title = "Climate progressiveness [m1]", x = NULL, y = "Visibility") +
th + ylim(ymin, ymax)

#climate progressiveness by period
levels(cf$"procli:period"$period) = c("Pre-Paris", "Post-Paris")
p_procli_period_1 = ggplot(cf$"procli:period", aes(x = procli, y = estimate__)) +
geom_smooth(aes(color = period), method = "loess") +
geom_ribbon(aes(color = period, ymin = lower__, ymax = upper__), linewidth = 0.2, linetype = 1, alpha = 0.1) +
scale_color_manual(values = c("#b92e34", "#1b03a3")) +
labs(title = "Climate progressiveness by period [m1]", x = NULL, y = "Visibility") +
th + ylim(ymin, ymax)

#autocorrelation
p_actot_1 = ggplot(cf$actot, aes(x = actot, y = estimate__)) +
geom_smooth(method = "loess", color = "#000000") +
geom_ribbon(aes(ymin = lower__, ymax = upper__), linewidth = 0.2, linetype = 1, alpha = 0.1, color = "#000000") +
labs(title = "Autocorrelation [m1]", x = NULL, y = "Visibility") +
th + ylim(ymin, ymax)

#empty
p_empty = ggplot(cf$popula, aes(x = popula, y = estimate__)) +
geom_smooth(method = "loess", color = "#ffffff") +
labs(title = NULL, x = NULL, y = NULL) +
theme(axis.text = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.ticks = element_blank())

###########
#effects m2
m = readRDS("m2.rds")
cf = conditional_effects(m)

#closen
p_closen_2 = ggplot(cf$closen, aes(x = closen, y = estimate__)) +
geom_smooth(method = "loess", color = "#000000") +
geom_ribbon(aes(ymin = lower__, ymax = upper__), linewidth = 0.2, linetype = 1, alpha = 0.1, color = "#000000") +
labs(title = "Connector [closeness, H2A] [m2]", x = NULL, y = "Visibility") +
th + ylim(ymin, ymax)

###########
#effects m3
m = readRDS("m3.rds")
cf = conditional_effects(m)

#eigenv
p_eigenv_3 = ggplot(cf$eigenv, aes(x = eigenv, y = estimate__)) +
geom_smooth(method = "loess", color = "#000000") +
geom_ribbon(aes(ymin = lower__, ymax = upper__), linewidth = 0.2, linetype = 1, alpha = 0.1, color = "#000000") +
labs(title = "Authority [eigenvector, H2B] [m3]", x = NULL, y = "Visibility") +
th + ylim(ymin, ymax)

###########
#effects m4
m = readRDS("m4.rds")
cf = conditional_effects(m)

#constr
p_constr_4 = ggplot(cf$constr, aes(x = constr, y = estimate__)) +
geom_smooth(method = "loess", color = "#000000") +
geom_ribbon(aes(ymin = lower__, ymax = upper__), linewidth = 0.2, linetype = 1, alpha = 0.1, color = "#000000") +
labs(title = "Broker [Burt's constraint, H3] [m4]", x = NULL, y = "Visibility") +
th + ylim(ymin, ymax)

###########
#effects m5
m = readRDS("m5.rds")
cf = conditional_effects(m)

#coalea
p_coalea_5 = ggplot(cf$coalea, aes(x = coalea, y = estimate__)) +
geom_smooth(method = "loess", color = "#000000") +
geom_ribbon(aes(ymin = lower__, ymax = upper__), linewidth = 0.2, linetype = 1, alpha = 0.1, color = "#000000") +
labs(title = "Coalition leader [community eigentrust, H4] [m5]", x = NULL, y = "Visibility") +
th + ylim(ymin, ymax)

########
#retheme
th = theme(text = element_text(size = 14, color = "black", family = "Times"), axis.text.x = element_text(color = "black", family = "Times", angle = 30, hjust = 1, vjust = 1), axis.text.y = element_text(color = "black", family = "Times"), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), legend.position = "right", legend.title = element_blank(), plot.title = element_text(hjust = 0.0))
ymin = 0
ymax = 45

############
#effects m1i
m = readRDS("m1i.rds")
cf = conditional_effects(m)

#popula by type
levels(cf$"popula:retype"$retype) = c("Institutional", "Non-institutional")
p_popula_retype_1 = ggplot(cf$"popula:retype", aes(x = popula, y = estimate__)) +
geom_smooth(aes(color = retype), method = "loess") +
geom_ribbon(aes(color = retype, ymin = lower__, ymax = upper__), linewidth = 0.2, linetype = 1, alpha = 0.1) +
scale_color_manual(values = c("#b92e34", "#1b03a3")) +
labs(title = "Popular partner by type [m1i]", x = NULL, y = "Visibility") +
th + ylim(ymin, ymax)

############
#effects m2i
m = readRDS("m2i.rds")
cf = conditional_effects(m)

#closeness by type
levels(cf$"closen:retype"$retype) = c("Institutional", "Non-institutional")
p_closen_retype_2 = ggplot(cf$"closen:retype", aes(x = closen, y = estimate__)) +
geom_smooth(aes(color = retype), method = "loess") +
geom_ribbon(aes(color = retype, ymin = lower__, ymax = upper__), linewidth = 0.2, linetype = 1, alpha = 0.1) +
scale_color_manual(values = c("#b92e34", "#1b03a3")) +
labs(title = "Connector by type [m2i]", x = NULL, y = "Visibility") +
th + ylim(ymin, ymax)

############
#effects m3i
m = readRDS("m3i.rds")
cf = conditional_effects(m)

#eigenvector by type
levels(cf$"eigenv:retype"$retype) = c("Institutional", "Non-institutional")
p_eigenv_retype_3 = ggplot(cf$"eigenv:retype", aes(x = eigenv, y = estimate__)) +
geom_smooth(aes(color = retype), method = "loess") +
geom_ribbon(aes(color = retype, ymin = lower__, ymax = upper__), linewidth = 0.2, linetype = 1, alpha = 0.1) +
scale_color_manual(values = c("#b92e34", "#1b03a3")) +
labs(title = "Authority by type [m3i]", x = NULL, y = "Visibility") +
th + ylim(ymin, ymax)

############
#effects m4i
m = readRDS("m4i.rds")
cf = conditional_effects(m)

#constraint by type
levels(cf$"constr:retype"$retype) = c("Institutional", "Non-institutional")
p_constr_retype_4 = ggplot(cf$"constr:retype", aes(x = constr, y = estimate__)) +
geom_smooth(aes(color = retype), method = "loess") +
geom_ribbon(aes(color = retype, ymin = lower__, ymax = sqrt(upper__)), linewidth = 0.2, linetype = 1, alpha = 0.1) +
scale_color_manual(values = c("#b92e34", "#1b03a3")) +
labs(title = "Broker by type [m4i]", x = NULL, y = "Visibility [square root upper bound]") +
th + ylim(ymin, ymax)

############
#effects m5i
m = readRDS("m5i.rds")
cf = conditional_effects(m)

#coalition leader by type
levels(cf$"coalea:retype"$retype) = c("Institutional", "Non-institutional")
p_coalea_retype_5 = ggplot(cf$"coalea:retype", aes(x = coalea, y = estimate__)) +
geom_smooth(aes(color = retype), method = "loess") +
geom_ribbon(aes(color = retype, ymin = lower__, ymax = upper__), linewidth = 0.2, linetype = 1, alpha = 0.1) +
scale_color_manual(values = c("#b92e34", "#1b03a3")) +
labs(title = "Coalition leader by type [m5i]", x = NULL, y = "Visibility") +
th + ylim(ymin, ymax)

######
#store
ge = ggpubr::ggarrange(p_type_1, p_popula_1, p_closen_2, p_eigenv_3, p_constr_4, p_coalea_5, ncol = 2, nrow = 3, align = "hv")
ggsave(plot = ge, filename = "conditionalEffects.png", width = 10, height = 12, dpi = 300)
ggsave(plot = ge, filename = "conditionalEffects.tiff", width = 10, height = 12, dpi = 500)
ps.options(family = "Times")
ggsave(plot = ge, filename = "conditionalEffects.eps", width = 10, height = 12, device = cairo_ps)

ge = ggpubr::ggarrange(p_popula_retype_1, p_closen_retype_2, p_eigenv_retype_3, p_constr_retype_4, p_coalea_retype_5, ncol = 5, nrow = 1, align = "hv")
ggsave(plot = ge, filename = "conditionalEffectsInteraction.png", width = 30, height = 10, dpi = 300)

ge = ggpubr::ggarrange(p_type_1, p_popula_1, p_twact_log_1, p_period_1, p_procli_1, p_procli_period_1, p_actot_1, p_empty, ncol = 4, nrow = 2, align = "hv")
ggsave(plot = ge, filename = "conditionalEffectsModel1.png", width = 20, height = 10, dpi = 300)

#################
## DIAGNOSTICS ##
#################

m = readRDS("m1.rds")
print(m)
pp = pp_check(m)
ps = pp_check(m, type = "stat_2d")
pbzero = function(y) mean(y == 0)
pz = pp_check(m, type = "stat", stat = "pbzero")

ge = ggpubr::ggarrange(pp, ps, pz, ncol = 3, nrow = 1, align = "hv")
ggsave(plot = ge, filename="posteriorPredictiveModel1.png", width = 12, height = 4, dpi = 300)
