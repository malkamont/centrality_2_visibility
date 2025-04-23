#############
## PREPARE ##
#############

#path
import sys
path = "/m/cs/scratch/ecanet/malkama5/mediaAccess/replication/"
sys.path.append(path)

#workspace
from scipy.spatial.distance import pdist, squareform
import graph_tool.all as gt
import matplotlib as mpl
from filtFunc import *
import igraph as ig
import pandas as pd
import numpy as np

#######
## 0 ##
#######

#attributes
at = pd.read_csv(path + "attribute.csv", delimiter = ";", header = 0)

#data
dt = pd.DataFrame({"period": 0, "actor": at[at["roster0"] == 1]["actor"]})
dt = pd.merge(dt, at.loc[:, ["actor", "name", "type"]], on = "actor", how = "left")

#visibility
vs = pd.read_csv(path + "vtot0.csv", delimiter = ";", header = 0)
dt = pd.merge(dt, vs, on = "actor", how = "left").fillna(0).astype({"vjou": "int", "vopi": "int", "vtot": "int"})

#influence
ri = pd.read_csv(path + "influence0.csv", delimiter = ";", header = 0)
ri = ig.Graph.DataFrame(ri, directed = True, use_vids = False).simplify(multiple = True, loops = True)
dt = pd.merge(dt, pd.DataFrame({"actor": ri.vs["name"], "repinf": ri.degree(mode = "in")}), on = "actor", how = "left")
dt["repinf"] = scaler(dt["repinf"])

#popularity
po = pd.read_csv(path + "collaboration0.csv", delimiter = ";", header = 0)
po = ig.Graph.DataFrame(po, directed = True, use_vids = False).simplify(multiple = True, loops = True)
dt = pd.merge(dt, pd.DataFrame({"actor": po.vs["name"], "popula": po.degree(mode = "in")}), on = "actor", how = "left")
dt["popula"] = scaler(dt["popula"])

#belief
co = ig.Graph.Read_GraphML(path + "c0.xml").to_graph_tool(vertex_attributes = {"name": "string"}) #reconstruct
bn = pd.read_csv(path + "beliefN.csv", delimiter = ";", header = 0)
bf = pd.read_csv(path + "belief0.csv", delimiter = ";", header = 0)
bf.index = bf["actor"].tolist()
bf.columns = ["actor", "status"] + bn["label"].tolist()
bf = bf[(bf["status"] != "non-respondent") & (bf["status"] != "invalid")].iloc[:, np.r_[0, 2:8, 12:20, 21, 25:33]]

im = pd.merge(pd.DataFrame({"actor": dt["actor"]}), bf, on = "actor", how = "left")
im.index = im["actor"]
im = im.iloc[:, 1:]
bf.index = bf["actor"]
bf = bf.iloc[:, 1:]
bf.isna().sum().sum() / (bf.shape[0] * bf.shape[1]) #percent blank
for row in im.index:
    for col in im.columns:
        if np.isnan(im.loc[row, col]):
            alterlist = []
            for e in co.edges():
                so, ta = e
                if (co.vp.name[so] == row):
                    alterlist.append(co.vp.name[ta])
                if (co.vp.name[ta] == row):
                    alterlist.append(co.vp.name[so])
            x = im.loc[alterlist, col].dropna().mean()
            if np.isnan(x):
                x = 3
            im.loc[row, col] = x
im = im.round().astype("int")
im["actor"] = im.index
np.savetxt(path + "rb0.csv", im.values, delimiter = ";", fmt = "%s", header = ";".join([str(e) for e in im.columns.tolist()]), comments = "")

#belief homophily
bf = pd.read_csv(path + "procli0.csv", delimiter = ";", header = 0)
dt = pd.merge(dt, bf, on = "actor", how = "left")
ds = pdist(pd.DataFrame({"d0": dt["procli"], "d1": 0}).values, metric = "euclidean")
ds = pd.DataFrame(squareform(ds)).max() - pd.DataFrame(squareform(ds))
ds.index = dt["actor"]
ds.columns = dt["actor"]

co = ig.Graph.Read_GraphML(path + "c0.xml")
co.es["weight"] = 0
for e in co.es():
    so = e.source
    ta = e.target
    e["weight"] = int(np.around(ds.loc[co.vs["name"][so], co.vs["name"][ta]] * 10) + 1)
cf = co.copy() #filter

#autocorrelation
ac = pd.DataFrame(co.get_adjacency())
ac.index = co.vs["name"]
ac.columns = co.vs["name"]
ac = ac.loc[dt["actor"], dt["actor"]]
ac.index = np.arange(len(ac))
ac.columns = np.arange(len(ac))

dt["acjou"] = scaler((ac / ac.sum(axis = 1)).transpose().dot(dt["vjou"]))
dt["acopi"] = scaler((ac / ac.sum(axis = 1)).transpose().dot(dt["vopi"]))
dt["actot"] = scaler((ac / ac.sum(axis = 1)).transpose().dot(dt["vtot"]))

#centrality
dg = scaler(co.degree())
cn = scaler(co.coreness())
bu = scaler(co.constraint(weights = None))
bp = scaler(co.vs["bonpow"])
co = co.to_graph_tool(vertex_attributes = {"name": "string"}, edge_attributes = {"weight": "int"})
cl = scaler(gt.closeness(co, weight = None, norm = True).a)
ev = scaler(gt.eigenvector(co, weight = None)[1].a)
bw = scaler(gt.betweenness(co, weight = None, norm = True)[0].a)
et = scaler(gt.eigentrust(co, trust_map = co.ep.weight, norm = True).a)
dt = pd.merge(dt, pd.DataFrame({"actor": co.vp.name, "degree": dg, "corene": cn, "constr": bu, "bonpow": bp, "closen": cl, "eigenv": ev, "betwee": bw, "eigent": et}), on = "actor", how = "left")

#filter
cf = marginal_likelihood(g = cf)
cf.delete_edges([e.index for e in cf.es if e["p_value"] >= 0.05])
ig.write(cf, filename = path + "rc0_f.xml", format = "graphml")

#filter autocorrelation
ac = pd.DataFrame(cf.get_adjacency())
ac.index = cf.vs["name"]
ac.columns = cf.vs["name"]
ac = ac.loc[dt["actor"], dt["actor"]]
ac.index = np.arange(len(ac))
ac.columns = np.arange(len(ac))

dt["acjou_f"] = scaler((ac / ac.sum(axis = 1)).transpose().dot(dt["vjou"]))
dt["acopi_f"] = scaler((ac / ac.sum(axis = 1)).transpose().dot(dt["vopi"]))
dt["actot_f"] = scaler((ac / ac.sum(axis = 1)).transpose().dot(dt["vtot"]))

#filter centrality
dg = scaler(cf.degree())
cn = scaler(cf.coreness())
bu = scaler(cf.constraint(weights = None))
bp = scaler(cf.vs["bonpow"])
cf = cf.to_graph_tool(vertex_attributes = {"name": "string"}, edge_attributes = {"weight": "int"})
cl = scaler(gt.closeness(cf, weight = None, norm = True).a)
ev = scaler(gt.eigenvector(cf, weight = None)[1].a)
bw = scaler(gt.betweenness(cf, weight = None, norm = True)[0].a)
et = scaler(gt.eigentrust(cf, trust_map = cf.ep.weight, norm = True).a)
dt = pd.merge(dt, pd.DataFrame({"actor": cf.vp.name, "degree_f": dg, "corene_f": cn, "constr_f": bu, "bonpow_f": bp, "closen_f": cl, "eigenv_f": ev, "betwee_f": bw, "eigent_f": et}), on = "actor", how = "left")

#community
gt.seed_rng(1); np.random.seed(1)
s = gt.minimize_blockmodel_dl(cf, state = gt.PPBlockState, multilevel_mcmc_args = dict(B_min = 1, B_max = 2, niter = 10))
p = gt.sfdp_layout(cf, groups = s.b, gamma = 0.01, kappa = 20)
s.draw(pos = p)
values, counts = np.unique(s.b.a, return_counts = True)
value_to_label = {value: label for label, (value, count) in enumerate(sorted(zip(values, counts), key = lambda x: - x[1]))}
s.b.a = np.array([value_to_label[value] for value in s.b.a])
cf.vp.b = s.b
b = co.new_vp("int")
b.a = cf.vp.b.a
co.vp.b = b
dt = pd.merge(dt, pd.DataFrame({"actor": cf.vp.name, "coalit": s.b}), on = "actor", how = "left")

#community centrality
cc = co.copy()
td = []
for v in cc.vertices():
    if cc.vp.b[v] != 0:
        td.append(v)
for v in sorted(td, reverse = True):
    cc.remove_vertex(v)
et = scaler(gt.eigentrust(cc, trust_map = cc.ep.weight, norm = True).a)
cc = cf.copy()
td = []
for v in cc.vertices():
    if cc.vp.b[v] != 0:
        td.append(v)
for v in sorted(td, reverse = True):
    cc.remove_vertex(v)
ef = scaler(gt.eigentrust(cc, trust_map = cc.ep.weight, norm = True).a)
cm = pd.DataFrame({"actor": cc.vp.name, "coalea": et, "coalea_f": ef})

cc = co.copy()
td = []
for v in cc.vertices():
    if cc.vp.b[v] != 1:
        td.append(v)
for v in sorted(td, reverse = True):
    cc.remove_vertex(v)
et = scaler(gt.eigentrust(cc, trust_map = cc.ep.weight, norm = True).a)
cc = cf.copy()
td = []
for v in cc.vertices():
    if cc.vp.b[v] != 1:
        td.append(v)
for v in sorted(td, reverse = True):
    cc.remove_vertex(v)
ef = scaler(gt.eigentrust(cc, trust_map = cc.ep.weight, norm = True).a)
cm = pd.concat([cm, pd.DataFrame({"actor": cc.vp.name, "coalea": et, "coalea_f": ef})], axis = 0, ignore_index = True)

dt = pd.merge(dt, cm, on = "actor", how = "left")
d0 = dt.copy()

#plot
frame = pd.merge(pd.DataFrame({"actor": cf.vp.name}), dt.loc[:, ["actor", "coalit", "coalea"]], on = "actor", how = "left")
color = cf.new_vp("string")
vsize = cf.new_vp("float")
for v in np.arange(cf.num_vertices()):
    if (frame["coalit"][v] == 0):
        color[v] = "#b92e34"
    if (frame["coalit"][v] == 1):
        color[v] = "#1b03a3"
    vsize[v] = frame["coalea"][v]
gt.graph_draw(cf, pos = p, vertex_color = color, vertex_fill_color = vsize, vcmap = mpl.cm.magma, vorder = vsize, vertex_size = gt.prop_to_size(vsize, 8, 22), edge_pen_width = 0.2, vertex_text = None, output = path + "coalea0.svg")



#######
## 1 ##
#######

#attributes
at = pd.read_csv(path + "attribute.csv", delimiter = ";", header = 0)

#data
dt = pd.DataFrame({"period": 1, "actor": at[at["roster1"] == 1]["actor"]})
dt = pd.merge(dt, at.loc[:, ["actor", "name", "type"]], on = "actor", how = "left")

#visibility
vs = pd.read_csv(path + "vtot1.csv", delimiter = ";", header = 0)
dt = pd.merge(dt, vs, on = "actor", how = "left").fillna(0).astype({"vjou": int, "vopi": int, "vtot": int})

#influence
ri = pd.read_csv(path + "influence1.csv", delimiter = ";", header = 0)
ri = ig.Graph.DataFrame(ri, directed = True, use_vids = False).simplify(multiple = True, loops = True)
dt = pd.merge(dt, pd.DataFrame({"actor": ri.vs["name"], "repinf": ri.degree(mode = "in")}), on = "actor", how = "left")
dt["repinf"] = scaler(dt["repinf"])

#popularity
po = pd.read_csv(path + "collaboration1.csv", delimiter = ";", header = 0)
po = ig.Graph.DataFrame(po, directed = True, use_vids = False).simplify(multiple = True, loops = True)
dt = pd.merge(dt, pd.DataFrame({"actor": po.vs["name"], "popula": po.degree(mode = "in")}), on = "actor", how = "left")
dt["popula"] = scaler(dt["popula"])

#belief
co = ig.Graph.Read_GraphML(path + "c1.xml").to_graph_tool(vertex_attributes = {"name": "string"}) #reconstruct
bn = pd.read_csv(path + "beliefN.csv", delimiter = ";", header = 0)

b0 = pd.read_csv(path + "belief0.csv", delimiter = ";", header = 0) #adopt
b0.index = b0["actor"].tolist()
b0.columns = ["actor", "status"] + bn["label"].tolist()
b0 = b0[(b0["status"] != "non-respondent") & (b0["status"] != "invalid")].iloc[:, np.r_[0, 2:8, 12:20, 21, 25:33]]

bf = pd.read_csv(path + "belief1.csv", delimiter = ";", header = 0)
bf.index = bf["actor"].tolist()
bf.columns = ["actor", "status"] + bn["label"].tolist()
bf = bf[(bf["status"] != "non-respondent") & (bf["status"] != "invalid")].iloc[:, np.r_[0, 2:8, 12:20, 21, 25:33]]

for row in bf.index:
    for col in bf.columns[1:]:
        if np.isnan(bf.loc[row, col]):
            if row in b0.index:
                if col in b0.columns:
                    bf.loc[row, col] = b0.loc[row, col]

im = pd.merge(pd.DataFrame({"actor": dt["actor"]}), bf, on = "actor", how = "left")
im.index = im["actor"]
im = im.iloc[:, 1:]
bf.index = bf["actor"]
bf = bf.iloc[:, 1:]
bf.isna().sum().sum() / (bf.shape[0] * bf.shape[1]) #percent blank
for row in im.index:
    for col in im.columns:
        if np.isnan(im.loc[row, col]):
            alterlist = []
            for e in co.edges():
                so, ta = e
                if (co.vp.name[so] == row):
                    alterlist.append(co.vp.name[ta])
                if (co.vp.name[ta] == row):
                    alterlist.append(co.vp.name[so])
            x = im.loc[alterlist, col].dropna().mean()
            if np.isnan(x):
                x = 3
            im.loc[row, col] = x
im = im.round().astype("int")
im["actor"] = im.index
np.savetxt(path + "rb1.csv", im.values, delimiter = ";", fmt = "%s", header = ";".join([str(e) for e in im.columns.tolist()]), comments = "")

#belief homophily
bf = pd.read_csv(path + "procli1.csv", delimiter = ";", header = 0)
dt = pd.merge(dt, bf, on = "actor", how = "left")
ds = pdist(pd.DataFrame({"d0": dt["procli"], "d1": 0}).values, metric = "euclidean")
ds = pd.DataFrame(squareform(ds)).max() - pd.DataFrame(squareform(ds))
ds.index = dt["actor"]
ds.columns = dt["actor"]

co = ig.Graph.Read_GraphML(path + "c1.xml")
co.es["weight"] = 0
for e in co.es():
    so = e.source
    ta = e.target
    e["weight"] = int(np.around(ds.loc[co.vs["name"][so], co.vs["name"][ta]] * 10) + 1)
cf = co.copy() #filter

#autocorrelation
ac = pd.DataFrame(co.get_adjacency())
ac.index = co.vs["name"]
ac.columns = co.vs["name"]
ac = ac.loc[dt["actor"], dt["actor"]]
ac.index = np.arange(len(ac))
ac.columns = np.arange(len(ac))

dt["acjou"] = scaler((ac / ac.sum(axis = 1)).transpose().dot(dt["vjou"]))
dt["acopi"] = scaler((ac / ac.sum(axis = 1)).transpose().dot(dt["vopi"]))
dt["actot"] = scaler((ac / ac.sum(axis = 1)).transpose().dot(dt["vtot"]))

#centrality
dg = scaler(co.degree())
cn = scaler(co.coreness())
bu = scaler(co.constraint(weights = None))
bp = scaler(co.vs["bonpow"])
co = co.to_graph_tool(vertex_attributes = {"name": "string"}, edge_attributes = {"weight": "int"})
cl = scaler(gt.closeness(co, weight = None, norm = True).a)
ev = scaler(gt.eigenvector(co, weight = None)[1].a)
bw = scaler(gt.betweenness(co, weight = None, norm = True)[0].a)
et = scaler(gt.eigentrust(co, trust_map = co.ep.weight, norm = True).a)
dt = pd.merge(dt, pd.DataFrame({"actor": co.vp.name, "degree": dg, "corene": cn, "constr": bu, "bonpow": bp, "closen": cl, "eigenv": ev, "betwee": bw, "eigent": et}), on = "actor", how = "left")

#filter
cf = marginal_likelihood(g = cf)
cf.delete_edges([e.index for e in cf.es if e["p_value"] >= 0.05])
ig.write(cf, filename = path + "rc1_f.xml", format = "graphml")

#filter autocorrelation
ac = pd.DataFrame(cf.get_adjacency())
ac.index = cf.vs["name"]
ac.columns = cf.vs["name"]
ac = ac.loc[dt["actor"], dt["actor"]]
ac.index = np.arange(len(ac))
ac.columns = np.arange(len(ac))

dt["acjou_f"] = scaler((ac / ac.sum(axis = 1)).transpose().dot(dt["vjou"]))
dt["acopi_f"] = scaler((ac / ac.sum(axis = 1)).transpose().dot(dt["vopi"]))
dt["actot_f"] = scaler((ac / ac.sum(axis = 1)).transpose().dot(dt["vtot"]))

#filter centrality
dg = scaler(cf.degree())
cn = scaler(cf.coreness())
bu = scaler(cf.constraint(weights = None))
bp = scaler(cf.vs["bonpow"])
cf = cf.to_graph_tool(vertex_attributes = {"name": "string"}, edge_attributes = {"weight": "int"})
cl = scaler(gt.closeness(cf, weight = None, norm = True).a)
ev = scaler(gt.eigenvector(cf, weight = None)[1].a)
bw = scaler(gt.betweenness(cf, weight = None, norm = True)[0].a)
et = scaler(gt.eigentrust(cf, trust_map = cf.ep.weight, norm = True).a)
dt = pd.merge(dt, pd.DataFrame({"actor": cf.vp.name, "degree_f": dg, "corene_f": cn, "constr_f": bu, "bonpow_f": bp, "closen_f": cl, "eigenv_f": ev, "betwee_f": bw, "eigent_f": et}), on = "actor", how = "left")

#community
gt.seed_rng(1); np.random.seed(1)
s = gt.minimize_blockmodel_dl(cf, state = gt.PPBlockState, multilevel_mcmc_args = dict(B_min = 1, B_max = 2, niter = 10))
p = gt.sfdp_layout(cf, groups = s.b, gamma = 0.01, kappa = 20)
s.draw(pos = p)
values, counts = np.unique(s.b.a, return_counts = True)
value_to_label = {value: label for label, (value, count) in enumerate(sorted(zip(values, counts), key = lambda x: - x[1]))}
s.b.a = np.array([value_to_label[value] for value in s.b.a])
cf.vp.b = s.b
b = co.new_vp("int")
b.a = cf.vp.b.a
co.vp.b = b
dt = pd.merge(dt, pd.DataFrame({"actor": cf.vp.name, "coalit": s.b}), on = "actor", how = "left")

#community centrality
cc = co.copy()
td = []
for v in cc.vertices():
    if cc.vp.b[v] != 0:
        td.append(v)
for v in sorted(td, reverse = True):
    cc.remove_vertex(v)
et = scaler(gt.eigentrust(cc, trust_map = cc.ep.weight, norm = True).a)
cc = cf.copy()
td = []
for v in cc.vertices():
    if cc.vp.b[v] != 0:
        td.append(v)
for v in sorted(td, reverse = True):
    cc.remove_vertex(v)
ef = scaler(gt.eigentrust(cc, trust_map = cc.ep.weight, norm = True).a)
cm = pd.DataFrame({"actor": cc.vp.name, "coalea": et, "coalea_f": ef})

cc = co.copy()
td = []
for v in cc.vertices():
    if cc.vp.b[v] != 1:
        td.append(v)
for v in sorted(td, reverse = True):
    cc.remove_vertex(v)
et = scaler(gt.eigentrust(cc, trust_map = cc.ep.weight, norm = True).a)
cc = cf.copy()
td = []
for v in cc.vertices():
    if cc.vp.b[v] != 1:
        td.append(v)
for v in sorted(td, reverse = True):
    cc.remove_vertex(v)
ef = scaler(gt.eigentrust(cc, trust_map = cc.ep.weight, norm = True).a)
cm = pd.concat([cm, pd.DataFrame({"actor": cc.vp.name, "coalea": et, "coalea_f": ef})], axis = 0, ignore_index = True)

dt = pd.merge(dt, cm, on = "actor", how = "left")
d1 = dt.copy()

#plot
frame = pd.merge(pd.DataFrame({"actor": cf.vp.name}), dt.loc[:, ["actor", "coalit", "coalea"]], on = "actor", how = "left")
color = cf.new_vp("string")
vsize = cf.new_vp("float")
for v in np.arange(cf.num_vertices()):
    if (frame["coalit"][v] == 0):
        color[v] = "#b92e34"
    if (frame["coalit"][v] == 1):
        color[v] = "#1b03a3"
    vsize[v] = frame["coalea"][v]
gt.graph_draw(cf, pos = p, vertex_color = color, vertex_fill_color = vsize, vcmap = mpl.cm.magma, vorder = vsize, vertex_size = gt.prop_to_size(vsize, 8, 22), edge_pen_width = 0.2, vertex_text = None, output = path + "coalea1.svg")



###########
## STORE ##
###########

#save
dt = pd.concat([d0, d1], axis = 0, ignore_index = True)
np.savetxt(path + "dt.csv", dt.values, delimiter = ";", fmt = "%s", header = ";".join([str(e) for e in dt.columns.tolist()]), comments = "")
