#workspace
setwd("/m/cs/scratch/ecanet/malkama5/mediaAccess/replication/")

#load
dt = readRDS("jointData.rds")

#retype
dt$retype = NA
dt$retype[dt$type == "governmental" | dt$type == "political_party"] = "institutional"
dt$retype[dt$type == "corporate" | dt$type == "research"] = "non_institutional"
dt$retype[dt$type == "non_governmental" | dt$type == "interest_group"] = "non_institutional"
dt$retype = factor(dt$retype)

#brm
library(brms)
dat = dt
iter = 40000
sp = save_pars(all = TRUE)
spr = "no"
ctrl = list(adapt_delta = 0.999, max_treedepth = 20)
fam = negbinomial(link = "log", link_shape = "log")
silent = 2
pg = FALSE
seed = 0
cores = as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE", parallel::detectCores()))
print(paste("Running", cores, "cores"))

#prior
mpr = c(
set_prior("normal(0, 30)", class = "b"), #applies to intercept as well since we treat it as a fixed effect via "0 + Intercept"
set_prior("student_t(3, 0, 3)", class = "sd"),
set_prior("gamma(0.01, 0.001)", class = "shape"))

###############
## ALL EDGES ##
###############

#m0
f = bf(vtot ~ 0 + Intercept + type + twact_log + period + procli + period:procli + actot + (1|actor))
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "m0.rds")

#m0z zero-inflated
f = bf(vtot ~ 0 + Intercept + type + twact_log + period + procli + period:procli + actot + (1|actor))
m = brm(f, family = zero_inflated_negbinomial(link = "log", link_shape = "log"), data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "m0z.rds")

#m0j journalistic content
f = bf(vjou ~ 0 + Intercept + type + twact_log + period + procli + period:procli + acjou + (1|actor))
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "m0j.rds")

#m1
f = bf(vtot ~ 0 + Intercept + type + twact_log + period + procli + period:procli + actot + (1|actor) + popula)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "m1.rds")

f = bf(vtot ~ 0 + Intercept + retype + twact_log + period + procli + period:procli + actot + (1|actor) + popula + retype:popula)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "m1i.rds")

#m2
f = bf(vtot ~ 0 + Intercept + type + twact_log + period + procli + period:procli + actot + (1|actor) + closen)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "m2.rds")

f = bf(vtot ~ 0 + Intercept + retype + twact_log + period + procli + period:procli + actot + (1|actor) + closen + retype:closen)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "m2i.rds")

#m3
f = bf(vtot ~ 0 + Intercept + type + twact_log + period + procli + period:procli + actot + (1|actor) + eigenv)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "m3.rds")

f = bf(vtot ~ 0 + Intercept + retype + twact_log + period + procli + period:procli + actot + (1|actor) + eigenv + retype:eigenv)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "m3i.rds")

#m4
f = bf(vtot ~ 0 + Intercept + type + twact_log + period + procli + period:procli + actot + (1|actor) + constr)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "m4.rds")

f = bf(vtot ~ 0 + Intercept + retype + twact_log + period + procli + period:procli + actot + (1|actor) + constr + retype:constr)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "m4i.rds")

#m5
f = bf(vtot ~ 0 + Intercept + type + twact_log + period + procli + period:procli + actot + (1|actor) + coalea)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "m5.rds")

f = bf(vtot ~ 0 + Intercept + retype + twact_log + period + procli + period:procli + actot + (1|actor) + coalea + retype:coalea)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "m5i.rds")

#m6
f = bf(vtot ~ 0 + Intercept + type + twact_log + period + procli + period:procli + actot + (1|actor) + bonpow)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "m6.rds")

f = bf(vtot ~ 0 + Intercept + retype + twact_log + period + procli + period:procli + actot + (1|actor) + bonpow + retype:bonpow)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "m6i.rds")

#m7
f = bf(vtot ~ 0 + Intercept + type + twact_log + period + procli + period:procli + actot + (1|actor) + betwee)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "m7.rds")

f = bf(vtot ~ 0 + Intercept + retype + twact_log + period + procli + period:procli + actot + (1|actor) + betwee + retype:betwee)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "m7i.rds")

#mall
f = bf(vtot ~ 0 + Intercept + type + twact_log + period + procli + period:procli + actot + (1|actor) + popula + eigenv + constr + coalea)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "mall.rds")

#######################
## SIGNIFICANT EDGES ##
#######################

#m0
f = bf(vtot ~ 0 + Intercept + type + twact_log + period + procli + period:procli + actot_f + (1|actor))
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "f_m0.rds")

#m0z zero-inflated
f = bf(vtot ~ 0 + Intercept + type + twact_log + period + procli + period:procli + actot_f + (1|actor))
m = brm(f, family = zero_inflated_negbinomial(link = "log", link_shape = "log"), data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "f_m0z.rds")

#m0j journalistic content
f = bf(vjou ~ 0 + Intercept + type + twact_log + period + procli + period:procli + acjou_f + (1|actor))
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "f_m0j.rds")

#m1
f = bf(vtot ~ 0 + Intercept + type + twact_log + period + procli + period:procli + actot_f + (1|actor) + popula)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "f_m1.rds")

f = bf(vtot ~ 0 + Intercept + retype + twact_log + period + procli + period:procli + actot_f + (1|actor) + popula + retype:popula)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "f_m1i.rds")

#m2
f = bf(vtot ~ 0 + Intercept + type + twact_log + period + procli + period:procli + actot_f + (1|actor) + closen_f)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "f_m2.rds")

f = bf(vtot ~ 0 + Intercept + retype + twact_log + period + procli + period:procli + actot_f + (1|actor) + closen_f + retype:closen_f)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "f_m2i.rds")

#m3
f = bf(vtot ~ 0 + Intercept + type + twact_log + period + procli + period:procli + actot_f + (1|actor) + eigenv_f)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "f_m3.rds")

f = bf(vtot ~ 0 + Intercept + retype + twact_log + period + procli + period:procli + actot_f + (1|actor) + eigenv_f + retype:eigenv_f)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "f_m3i.rds")

#m4
f = bf(vtot ~ 0 + Intercept + type + twact_log + period + procli + period:procli + actot_f + (1|actor) + constr_f)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "f_m4.rds")

f = bf(vtot ~ 0 + Intercept + retype + twact_log + period + procli + period:procli + actot_f + (1|actor) + constr_f + retype:constr_f)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "f_m4i.rds")

#m5
f = bf(vtot ~ 0 + Intercept + type + twact_log + period + procli + period:procli + actot_f + (1|actor) + coalea_f)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "f_m5.rds")

f = bf(vtot ~ 0 + Intercept + retype + twact_log + period + procli + period:procli + actot_f + (1|actor) + coalea_f + retype:coalea_f)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "f_m5i.rds")

#m6
f = bf(vtot ~ 0 + Intercept + type + twact_log + period + procli + period:procli + actot_f + (1|actor) + bonpow_f)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "f_m6.rds")

f = bf(vtot ~ 0 + Intercept + retype + twact_log + period + procli + period:procli + actot_f + (1|actor) + bonpow_f + retype:bonpow_f)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "f_m6i.rds")

#m7
f = bf(vtot ~ 0 + Intercept + type + twact_log + period + procli + period:procli + actot_f + (1|actor) + betwee_f)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "f_m7.rds")

f = bf(vtot ~ 0 + Intercept + retype + twact_log + period + procli + period:procli + actot_f + (1|actor) + betwee_f + retype:betwee_f)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "f_m7i.rds")

#mall
f = bf(vtot ~ 0 + Intercept + type + twact_log + period + procli + period:procli + actot_f + (1|actor) + popula + eigenv_f + constr_f + coalea_f)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "f_mall.rds")

##########################
## MEDIA STRATEGY CHECK ##
##########################

#reshape
ms = read.table("mstrat.csv", header = TRUE, sep = ";")
dt = merge(x = dt, y = ms, by = c("actor", "period"), sort = FALSE)
dt$mstrat = factor(dt$mstrat)
dat = dt

#m0
f = bf(vtot ~ 0 + Intercept + type + mstrat + period + procli + period:procli + actot + (1|actor))
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "c_m0.rds")

#m1
f = bf(vtot ~ 0 + Intercept + type + mstrat + period + procli + period:procli + actot + (1|actor) + popula)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "c_m1.rds")

#m2
f = bf(vtot ~ 0 + Intercept + type + mstrat + period + procli + period:procli + actot + (1|actor) + closen)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "c_m2.rds")

#m3
f = bf(vtot ~ 0 + Intercept + type + mstrat + period + procli + period:procli + actot + (1|actor) + eigenv)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "c_m3.rds")

#m4
f = bf(vtot ~ 0 + Intercept + type + mstrat + period + procli + period:procli + actot + (1|actor) + constr)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "c_m4.rds")

#m5
f = bf(vtot ~ 0 + Intercept + type + mstrat + period + procli + period:procli + actot + (1|actor) + coalea)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "c_m5.rds")

#mall
f = bf(vtot ~ 0 + Intercept + type + mstrat + period + procli + period:procli + actot + (1|actor) + popula + eigenv + constr + coalea)
m = brm(f, family = fam, data = dat, iter = iter, save_pars = sp, control = ctrl, cores = cores, seed = seed, silent = silent, open_progress = pg, sample_prior = "yes", prior = mpr)
saveRDS(m, "c_mall.rds")

###################
## LEAVE-ONE-OUT ##
###################

#loo
m0 = readRDS("m0.rds")
m0z = readRDS("m0z.rds")
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
mall = readRDS("mall.rds")

lm = loo(m0, m0z, m1, m2, m3, m4, m5, m6, m7, m1i, m2i, m3i, m4i, m5i, m6i, m7i, mall, moment_match = TRUE, cores = cores, reloo = FALSE)
saveRDS(lm, "loo.rds")

#loo filter
m0 = readRDS("f_m0.rds")
m0z = readRDS("f_m0z.rds")
m1 = readRDS("f_m1.rds")
m2 = readRDS("f_m2.rds")
m3 = readRDS("f_m3.rds")
m4 = readRDS("f_m4.rds")
m5 = readRDS("f_m5.rds")
m6 = readRDS("f_m6.rds")
m7 = readRDS("f_m7.rds")
m1i = readRDS("f_m1i.rds")
m2i = readRDS("f_m2i.rds")
m3i = readRDS("f_m3i.rds")
m4i = readRDS("f_m4i.rds")
m5i = readRDS("f_m5i.rds")
m6i = readRDS("f_m6i.rds")
m7i = readRDS("f_m7i.rds")
mall = readRDS("f_mall.rds")

lm = loo(m0, m0z, m1, m2, m3, m4, m5, m6, m7, m1i, m2i, m3i, m4i, m5i, m6i, m7i, mall, moment_match = TRUE, cores = cores, reloo = FALSE)
saveRDS(lm, "f_loo.rds")

#loo media strategy check
m0 = readRDS("c_m0.rds")
m1 = readRDS("c_m1.rds")
m2 = readRDS("c_m2.rds")
m3 = readRDS("c_m3.rds")
m4 = readRDS("c_m4.rds")
m5 = readRDS("c_m5.rds")
mall = readRDS("c_mall.rds")

lm = loo(m0, m1, m2, m3, m4, m5, mall, moment_match = TRUE, cores = cores, reloo = FALSE)
saveRDS(lm, "c_loo.rds")
