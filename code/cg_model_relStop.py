import numpy as np
import pandas as pd
from gurobipy import *
import time
#test comment

def DFS(adj):
    color = {vert : 'white' for vert in adj.keys()}
    cc = []
    for vert in adj.keys():
        if color[vert] == 'white':
            comp = []
            cc.append(process(comp, adj, vert, color))
    return cc
def process(comp, adj, vert, color):
    color[vert] = 'black'
    comp.append(vert)
    for v in adj[vert]:
        if color[v] == 'white':
            comp = process(comp, adj, v, color)
    return comp

def callback(model, where):
    if where == GRB.Callback.MIP:
        runtime = model.cbGet(GRB.Callback.RUNTIME)
        objbst = model.cbGet(GRB.Callback.MIP_OBJBST)
        #objbnd = model.cbGet(GRB.Callback.MIP_OBJBND)
        
        if runtime > 60 and objbst < 0:
            model.terminate()
    elif where == GRB.Callback.MIPSOL:
        a_sol = model.cbGetSolution(model._a)
        inNtwk = []
        ids = []
        for i in range(model._n):
            if a_sol[i] >= 0.5:
                inNtwk.append(model._raw_in.iloc[i].name)
                ids.append(i)
        adj = {dmu: set() for dmu in inNtwk}
        for dmu in inNtwk:
            for dmu1 in inNtwk:
                if model._distr_int[model._distr_int["BUSINESS_UNIT"] == dmu][model._distr_int["SHIPTO_NO"] == dmu1]["TE"].sum() > 0:
                    adj[dmu].add(dmu1)
                    adj[dmu1].add(dmu)
        cc = DFS(adj)
        if len(cc) > 1:
            model.cbLazy(quicksum(model._a[i] for i in ids) <= len(ids) - 1)

class DEAcg:
    def __init__(self, raw_in, raw_out, groups, distr_int):

        self.obj_stats = []
        self.red_csts = []
        self.master_sols = []

        self.raw_in = raw_in - 2 + 1 #handling undesirable
        self.raw_out = raw_out - 1 + 2 #handling undesirable
        self.distr_int = distr_int
        self.n = raw_in.shape[0]

        X_lst = []
        for dmu in list(set(groups["group0"])):
            row_dict = {}
            row_dict["dmu"] = dmu
            row_dict["tot_flow"] = raw_in.loc[groups[groups["group0"] == dmu].index]["tot_flow"].sum()
            row_dict["veh_TE"] = raw_in.loc[groups[groups["group0"] == dmu].index]["veh_TE"].sum()
            #--------------undesirable outputs------------------
            row_dict["not_kept"] = raw_out.loc[groups[groups["group0"] == dmu].index]["not_kept"].sum()
            X_lst.append(row_dict)
        X = pd.DataFrame(X_lst).set_index("dmu")

        Y_lst = []
        for dmu in list(set(groups["group0"])):
            row_dict = {}
            row_dict["dmu"] = dmu
            row_dict["help"] = raw_out.loc[groups[groups["group0"] == dmu].index]["help"].sum()
            #--------------undesirable inputs------------------
            row_dict["ncust"] = raw_in.loc[groups[groups["group0"] == dmu].index]["ncust"].sum()
            row_dict["cdiff"] = np.dot(raw_in.loc[groups[groups["group0"] == dmu].index]["ncust"].sum(), raw_in.loc[groups[groups["group0"] == dmu].index]["cdiff"].sum())/row_dict["ncust"]
            Y_lst.append(row_dict)
        Y = pd.DataFrame(Y_lst).set_index("dmu")

        #binary constants defining subntwks in master
        ai_lst = []
        for k in X.index:
            t = [0 for x in range(self.n)]
            i = 0
            for dmu in raw_in.index:
                if groups.loc[dmu, "group0"] == k:
                    t[i] = 1
                i += 1
            ai_lst.append(t)
        ai = pd.DataFrame(ai_lst, columns=raw_in.index) 

        self.X = X
        self.Y = Y
        self.ai = ai

        self.n_in = raw_in.shape[1] - 2 + 1 #handling undesirable
        self.n_out = raw_out.shape[1] - 1 + 2 #handling undesirable
        self.K = self.X.shape[0]
        self.c = [] #DEA scores
        self.old_c = []

        self.master = Model("Partition")
        self.z = [] #master primal variables
        self.dmus = [] #master dual variables

    def setupDeaInit(self, e=0.001):
        #inputs and outputs transformations
        self.X_cur = self.X.copy()
        self.Y_cur = self.Y.copy()

        #normalization
        self.X_cur = (self.X_cur - self.X_cur.min())/(self.X_cur.max() - self.X_cur.min())  
        self.Y_cur = (self.Y_cur - self.Y_cur.min())/(self.Y_cur.max() - self.Y_cur.min())

        self.deaInit = Model("DEA")
        self.deaInit.setParam("OutputFlag", 0)
        self.nu_init = self.deaInit.addVars(self.n_in, lb=e, vtype=GRB.CONTINUOUS, name="nuIn")
        self.mu_init = self.deaInit.addVars(self.n_out, lb=e, vtype=GRB.CONTINUOUS, name="muIn")
        self.deaInit.update()
        self.deaInit.addConstrs(np.dot(self.mu_init.values(), self.Y_cur.loc[j]) - np.dot(self.nu_init.values(), self.X_cur.loc[j]) <= 0 for j in list(set(self.X_cur.index)))

    def updDEAscores(self, e=0.001):

        self.old_c.append(self.c)
        self.c = []
        #inputs and outputs transformations
        self.X_cur = self.X.copy()
        self.Y_cur = self.Y.copy()

        #normalization
        self.X_cur = (self.X_cur - self.X_cur.min())/(self.X_cur.max() - self.X_cur.min())  
        self.Y_cur = (self.Y_cur - self.Y_cur.min())/(self.Y_cur.max() - self.Y_cur.min())

        self.deaInit = Model("DEA")
        self.deaInit.setParam("OutputFlag", 0)
        self.nu_init = self.deaInit.addVars(self.n_in, lb=e, vtype=GRB.CONTINUOUS, name="nuIn")
        self.mu_init = self.deaInit.addVars(self.n_out, lb=e, vtype=GRB.CONTINUOUS, name="muIn")
        self.deaInit.update()
        self.deaInit.addConstrs(np.dot(self.mu_init.values(), self.Y_cur.loc[k]) - np.dot(self.nu_init.values(), self.X_cur.loc[k]) <= 0 for k in list(set(self.X_cur.index)))
        
        for k in self.X.index:
            tmpConstr = self.deaInit.addConstr(np.dot(self.nu_init.values(), self.X_cur.loc[k]) == 1)
            self.deaInit.setObjective(np.dot(self.mu_init.values(), self.Y_cur.loc[k]), GRB.MAXIMIZE)
            self.deaInit.optimize()
            self.c.append(self.deaInit.objVal)
            self.deaInit.remove(tmpConstr)
            self.deaInit.reset(0)

    def add_row(self):
        
        #add column (binary params)
        self.ai.loc[self.K] = [self.a[i].x for i in range(self.n)]

        #add input row
        self.X.loc[self.K, "tot_flow"] = self.raw_in.iloc[[i for i in range(self.n) if self.a[i].x > 0.5]]["tot_flow"].sum()
        self.X.loc[self.K, "veh_TE"] = self.raw_in.iloc[[i for i in range(self.n) if self.a[i].x > 0.5]]["veh_TE"].sum()
        #-------------undesirable outputs-----------------
        self.X.loc[self.K, "not_kept"] = self.raw_out.iloc[[i for i in range(self.n) if self.a[i].x > 0.5]]["not_kept"].sum()

        #add output rpw
        self.Y.loc[self.K, "help"] = self.raw_out.iloc[[i for i in range(self.n) if self.a[i].x > 0.5]]["help"].sum()
        #-------------undesirable inputs-----------------
        self.Y.loc[self.K, "ncust"] = self.raw_in.iloc[[i for i in range(self.n) if self.a[i].x > 0.5]]["ncust"].sum()
        self.Y.loc[self.K, "cdiff"] = np.dot(self.raw_in.iloc[[i for i in range(self.n) if self.a[i].x > 0.5]]["ncust"], self.raw_in.iloc[[i for i in range(self.n) if self.a[i].x > 0.5]]["cdiff"])/self.Y.loc[self.K, "ncust"] 

        self.K += 1

    def getXY(self):
        return self.X, self.Y

    def DFS(self, adj):
        color = {vert : 'white' for vert in adj.keys()}
        cc = []
        for vert in adj.keys():
            if color[vert] == 'white':
                comp = []
                cc.append(self.process(comp, adj, vert, color))
        return cc

    def process(self, comp, adj, vert, color):
        color[vert] = 'black'
        comp.append(vert)
        for v in adj[vert]:
            if color[v] == 'white':
                comp = self.process(comp, adj, v, color)
        return comp
    
    def setupSub(self, e=0.001):

        self.subproblem = Model("DEAgen")
        self.subproblem.setParam("OutputFlag", 0)
        self.subproblem.setParam("NonConvex", 2)
        self.subproblem.setParam("BestBdStop", 0)

        #setup subproblem
        self.a = self.subproblem.addVars(self.n, lb=0, ub=1, vtype=GRB.BINARY, name="a")
        self.w = self.subproblem.addVars(self.n, vtype=GRB.CONTINUOUS, name="w1") #cdiff linearization
        self.nu = self.subproblem.addVars(self.n_in, lb=e, vtype=GRB.CONTINUOUS, name="nu")
        self.mu = self.subproblem.addVars(self.n_out, lb=e, vtype=GRB.CONTINUOUS, name="mu")
        self.th = self.subproblem.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name="theta") #dea dual obj
        self.ld = self.subproblem.addVars(self.K, vtype=GRB.CONTINUOUS, name="l") #dea dual var
        self.alpha = self.subproblem.addVars(self.n_in, vtype=GRB.CONTINUOUS, name="alpha")
        self.beta = self.subproblem.addVars(self.n_out, vtype=GRB.CONTINUOUS, name="beta")
        self.f1 = self.subproblem.addVar(vtype=GRB.CONTINUOUS, name="f1") #denotes raw cdiff (linearization applies to raw cdiff)
        self.f2 = self.subproblem.addVar(vtype=GRB.CONTINUOUS, name="f2") #denotes transformed cdiff
        self.f3 = self.subproblem.addVar(vtype=GRB.CONTINUOUS, name="f3") #denotes normalized cdiff
        self.f4 = self.subproblem.addVar(vtype=GRB.CONTINUOUS, name="f4") #denotes normalized cdiff
        self.f5 = self.subproblem.addVar(vtype=GRB.CONTINUOUS, name="f5") #denotes normalized cdiff
        self.f6 = self.subproblem.addVar(vtype=GRB.CONTINUOUS, name="f6") #denotes normalized cdiff
        self.switch1 = self.subproblem.addVar(vtype=GRB.BINARY, name="switch1")
        self.switch2 = self.subproblem.addVar(vtype=GRB.BINARY, name="switch2")
        self.jswitch1 = self.subproblem.addVars(self.n - 1, vtype=GRB.BINARY, name="jswitch1")
        self.jswitch2 = self.subproblem.addVars(self.n - 1, vtype=GRB.BINARY, name="jswitch2")
        self.neqVars = self.subproblem.addVars(self.K, self.n, vtype=GRB.BINARY, name="neqVar") #switchers for banning same subnetworks
        self.curScore = self.subproblem.addVar(name="curScore")
        #tmp vars to check is inputs and outputs match in sub and updDea
        self.tmp_in0 = self.subproblem.addVar(lb=-GRB.INFINITY,name = "tmp_in0")
        self.tmp_in1 = self.subproblem.addVar(lb=-GRB.INFINITY,name = "tmp_in1")
        self.tmp_in2 = self.subproblem.addVar(lb=-GRB.INFINITY,name = "tmp_in2")
        # self.tmp_in3 = self.subproblem.addVar(lb=-GRB.INFINITY,name = "tmp_in3")
        self.tmp_out0 = self.subproblem.addVar(lb=-GRB.INFINITY,name = "tmp_out0")
        self.tmp_out1 = self.subproblem.addVar(lb=-GRB.INFINITY,name = "tmp_out1")
        self.tmp_out2 = self.subproblem.addVar(lb=-GRB.INFINITY,name = "tmp_out2")

        self.subproblem.update()

        #i/o transformations
        #inputs and outputs transformations
        X = self.X.copy()
        Y = self.Y.copy()
        '''#unwanted inputs/outputs
        self.max_ncust = 2000#max(X["ncust"])
        self.min_ncust = min(X["ncust"])
        self.max_cdiff = max(X["cdiff"])
        self.min_cdiff = min(X["cdiff"])
        X["ncust"] = -X["ncust"] + (self.max_ncust + self.min_ncust)
        X["cdiff"] = -X["cdiff"] + (self.max_cdiff + self.min_cdiff)
        self.min_nk = min(0, min(Y["not_kept"]))
        Y["not_kept"] = Y["not_kept"] - self.min_nk'''

        #linear expressions for current subnetwork inputs and outputs
        self.ncust = np.dot(self.raw_in["ncust"], self.a.values())
        self.not_kept = np.dot(self.raw_out["not_kept"], self.a.values())

        self.x0 = [np.dot(self.raw_in["tot_flow"], self.a.values()), np.dot(self.raw_in["veh_TE"], self.a.values()), self.not_kept]
        self.y0 = [np.dot(self.raw_out["help"], self.a.values()), self.ncust, 0]
        self.x0 = (self.x0 - X.min())/(X.max() - X.min()) 
        self.y0[2] = self.f1 
        self.y0 = (self.y0 - Y.min())/(Y.max() - Y.min())

        #normalization
        ymin = Y.min()
        ymax = Y.max()
        X = (X - X.min())/(X.max() - X.min())  
        Y = (Y - Y.min())/(Y.max() - Y.min())
    
        #nonnegative flow
        for i in range(self.n):
            tot_out = self.distr_int[self.distr_int["BUSINESS_UNIT"] == self.raw_in.iloc[i].name]["TE"].sum()
            tot_in = self.distr_int[self.distr_int["SHIPTO_NO"] == self.raw_in.iloc[i].name]["TE"].sum()
            wts1 = []
            wts2 = []
            for j in range(self.n):
                te = self.distr_int[self.distr_int["BUSINESS_UNIT"] == self.raw_in.iloc[i].name][self.distr_int["SHIPTO_NO"] == self.raw_in.iloc[j].name]["TE"].sum()# + self.distr_int[self.distr_int["BUSINESS_UNIT"] == raw_in.iloc[j].name][self.distr_int["SHIPTO_NO"] == raw_in.iloc[i].name]["TE"].sum()
                wts1.append(te) #t_ij
                te =  + self.distr_int[self.distr_int["BUSINESS_UNIT"] == self.raw_in.iloc[j].name][self.distr_int["SHIPTO_NO"] == self.raw_in.iloc[i].name]["TE"].sum()
                wts2.append(te) #t_ji
            for j in range(i):
                tot_outj = self.distr_int[self.distr_int["BUSINESS_UNIT"] == self.raw_in.iloc[j].name]["TE"].sum()
                tot_inj = self.distr_int[self.distr_int["SHIPTO_NO"] == self.raw_in.iloc[j].name]["TE"].sum()
                self.subproblem.addConstr(quicksum(wts1[j] * self.a[i] * self.a[j] for j in range(self.n)) >= 0.02 * (tot_inj + 0.01) * self.a[i] - 100000 * (1 - self.jswitch1[j]))
                self.subproblem.addConstr(quicksum(wts2[j] * self.a[i] * self.a[j] for j in range(self.n)) >= 0.02 * (tot_outj + 0.01) * self.a[i] - 100000 * (1 - self.jswitch2[j]))
            for j in range(i + 1, self.n):
                tot_outj = self.distr_int[self.distr_int["BUSINESS_UNIT"] == self.raw_in.iloc[j].name]["TE"].sum()
                tot_inj = self.distr_int[self.distr_int["SHIPTO_NO"] == self.raw_in.iloc[j].name]["TE"].sum()
                self.subproblem.addConstr(quicksum(wts1[j] * self.a[i] * self.a[j] for j in range(self.n)) >= 0.02 * (tot_inj + 0.01)* self.a[i] - 100000 * (1 - self.jswitch1[j - 1]))
                self.subproblem.addConstr(quicksum(wts2[j] * self.a[i] * self.a[j] for j in range(self.n)) >= 0.02 * (tot_outj + 0.01) * self.a[i] - 100000 * (1 - self.jswitch2[j - 1]))
            self.subproblem.addConstr(wts1[j] * self.a[i] * self.a[j] >= 0.02 * (tot_out + 0.01) * self.a[i] - 100000 * (1 - self.switch1)) #epsilon
            self.subproblem.addConstr(wts2[j] * self.a[i] * self.a[j] >= 0.02 * (tot_in + 0.01) * self.a[i] - 100000 * (1 - self.switch2))#'''
        self.subproblem.addConstr(self.switch1 + self.switch2 + self.jswitch1.sum() + self.jswitch2.sum() >= 1)

        #ban same subnetworks
        self.sameSubConstrs1 = self.subproblem.addConstrs(self.a[i] + self.ai.iloc[k, i] >= self.neqVars[k, i] for i in range(self.n) for k in range(self.K))
        self.sameSubConstrs2 = self.subproblem.addConstrs(self.a[i] + self.ai.iloc[k, i] <= 2 - self.neqVars[k, i] for i in range(self.n) for k in range(self.K))
        self.sameSubConstrs3 = self.subproblem.addConstrs(self.neqVars.sum(k, '*') >= 1 for k in range(self.K))

        self.subproblem.update()
        self.n_keep = len(self.subproblem.getConstrs())
        self.n_keepQ = len(self.subproblem.getQConstrs())
    
        #strong duality
        self.sdc = self.subproblem.addConstr(np.dot(self.mu.values(), self.y0) == self.th - e * (sum(self.alpha.values()) + sum(self.beta.values()))) #switch to linearize
        #primal feasibility
        self.pc1 = self.subproblem.addConstrs(np.dot(self.mu.values(), Y.loc[k]) - np.dot(self.nu.values(), X.loc[k]) <= 0 for k in list(set(X.index)))
        self.pc2 = self.subproblem.addConstr(np.dot(self.nu.values(), self.x0) == 1) #switch to linearize
        #dual feasibility
        self.dc1 = self.subproblem.addConstrs(np.dot(X.iloc[:, l], self.ld.values()) <= self.th * self.x0[l] - self.alpha[l] for l in range(self.n_in))
        self.dc2 = self.subproblem.addConstrs(np.dot(Y.iloc[:, r], self.ld.values()) >= self.y0[r] + self.beta[r] for r in range(self.n_out))
        #cdiff linearization
        self.cdiffc1 = self.subproblem.addConstr(self.f1 * self.ncust == sum([self.raw_in["ncust"].iloc[i] * self.raw_in["cdiff"].iloc[i] * self.a[i] for i in range(self.n)]))
        #transformation
        self.subproblem.addConstr(self.f3 == (self.f1 - ymin["cdiff"])/(ymax["cdiff"] - ymin["cdiff"]))
        self.csc = self.subproblem.addConstr(self.curScore == np.dot(self.mu.values(), self.y0))

        self.subproblem._a = self.a
        self.subproblem._n = self.n
        self.subproblem._distr_int = self.distr_int
        self.subproblem._raw_in = self.raw_in
        self.subproblem.Params.lazyConstraints = 1
    
    def updSub(self, e=0.01, niter=0):
        
        self.ld[self.K - 1] = self.subproblem.addVar(vtype=GRB.CONTINUOUS, name="l[%d]" % (self.K - 1))
        for i in range(self.n):
            self.neqVars[self.K - 1, i] = self.subproblem.addVar(vtype=GRB.BINARY, name="neqVar[%d,%d]" % (self.K - 1, i))

        self.subproblem.remove(self.subproblem.getConstrs()[self.n_keep:])
        self.subproblem.remove(self.subproblem.getQConstrs()[self.n_keepQ:])

        self.subproblem.update()

        #i/o transformations
        #inputs and outputs transformations
        X = self.X.copy()
        Y = self.Y.copy()

        #linear expressions for current subnetwork inputs and outputs
        self.ncust = np.dot(self.raw_in["ncust"], self.a.values())
        self.not_kept = np.dot(self.raw_out["not_kept"], self.a.values())

        self.x0 = [np.dot(self.raw_in["tot_flow"], self.a.values()), np.dot(self.raw_in["veh_TE"], self.a.values()), self.not_kept]
        self.y0 = [np.dot(self.raw_out["help"], self.a.values()), self.ncust, 0]
        self.x0 = (self.x0 - X.min())/(X.max() - X.min()) 
        self.y0[2] = self.f1 
        self.y0 = (self.y0 - Y.min())/(Y.max() - Y.min())

        #normalization
        ymin = Y.min()
        ymax = Y.max()
        X = (X - X.min())/(X.max() - X.min())  
        Y = (Y - Y.min())/(Y.max() - Y.min())
    
        #ban same subnetworks
        self.subproblem.addConstrs(self.a[i] + self.ai.iloc[self.K - 1, i] >= self.neqVars[self.K - 1, i] for i in range(self.n))
        self.subproblem.addConstr(self.neqVars.sum(self.K - 1, '*') >= 1)
        self.subproblem.addConstrs(self.a[i] + self.ai.iloc[self.K - 1, i] <= 2 - self.neqVars[self.K - 1, i] for i in range(self.n))

        self.subproblem.update()
        self.n_keep = len(self.subproblem.getConstrs())
    
        #strong duality
        self.sdc = self.subproblem.addConstr(np.dot(self.mu.values(), self.y0) == self.th - e * (sum(self.alpha.values()) + sum(self.beta.values()))) #switch to linearize
        #primal feasibility
        self.pc1 = self.subproblem.addConstrs(np.dot(self.mu.values(), Y.loc[k]) - np.dot(self.nu.values(), X.loc[k]) <= 0 for k in list(set(X.index)))
        self.pc2 = self.subproblem.addConstr(np.dot(self.nu.values(), self.x0) == 1) #switch to linearize
        #dual feasibility
        self.dc1 = self.subproblem.addConstrs(np.dot(X.iloc[:, l], self.ld.values()) <= self.th * self.x0[l] - self.alpha[l] for l in range(self.n_in))
        self.dc2 = self.subproblem.addConstrs(np.dot(Y.iloc[:, r], self.ld.values()) >= self.y0[r] + self.beta[r] for r in range(self.n_out))
        #cdiff linearization
        self.cdiffc1 = self.subproblem.addConstr(self.f1 * self.ncust == sum([self.raw_in["ncust"].iloc[i] * self.raw_in["cdiff"].iloc[i] * self.a[i] for i in range(self.n)]))
        #transformation
        self.subproblem.addConstr(self.f3 == (self.f1 - ymin["cdiff"])/(ymax["cdiff"] - ymin["cdiff"]))
        self.csc = self.subproblem.addConstr(self.curScore == np.dot(self.mu.values(), self.y0))


    def setup(self, M=10000):

        self.M = M #replace by adequate value
        
        self.master.setParam("OutputFlag", 0)
        #setup master problem
        for k in range(self.K):
            self.z.append(self.master.addVar(lb=0, ub=1, vtype=GRB.CONTINUOUS, name = "z%d" % k))
        self.updDEAscores()

        self.master.setObjective(np.dot(self.c, self.z))# - quicksum(np.max(0, quicksum(self.ai.iloc[k] * self.z[k])) for k in range(self.K)))
        self.master.update()

        for i in self.raw_in.index:
            self.dmus.append(self.master.addConstr(np.dot(self.ai.loc[:, i], self.z) >= 1))

    def solveLP(self, maxiter=200):
        self.last10 = []
        self.master.reset()
        self.master.optimize()
        
        self.obj_stats.append(self.master.objVal)
        self.master_sols.append([self.z[i].x for i in range(len(self.z))])

        self.setupSub()
        self.subTime = []
        self.subproblem.setObjective(np.dot(self.mu.values(), self.y0) - np.dot(self.a.values(), [self.dmus[i].pi for i in range(self.n)]), GRB.MINIMIZE)
        self.sObj1 = self.subproblem.addVar(lb=-GRB.INFINITY)
        self.sObj2 = self.subproblem.addVar(lb=-GRB.INFINITY)
        self.subproblem.update()
        self.subproblem.addConstr(self.sObj1 == np.dot(self.mu.values(), self.y0))
        self.subproblem.addConstr(self.sObj2 == - np.dot(self.a.values(), [self.dmus[i].pi for i in range(self.n)]))
        start = time.time()
        glob_start = time.time()
        self.subproblem.optimize(callback)
        stop = time.time()
        self.subTime.append(stop - start)
        self.red_csts.append(self.subproblem.objVal)
        self.add_row()
        self.last10.append(self.master.objVal)

        niter = 0
        while self.subproblem.objVal < 0:
            self.z.append(self.master.addVar(lb=0, ub=1, vtype=GRB.CONTINUOUS, name="z_%g" % len(self.z)))
            self.updDEAscores()
            for i in range(self.n):
                self.master.chgCoeff(self.dmus[i], self.z[-1], self.ai.iloc[-1, i])
            self.master.setObjective(np.dot(self.c, self.z))# + 0.001 * sum(self.z))
            self.master.update()
            self.master.optimize()
            self.obj_stats.append(self.master.objVal)
            self.master_sols.append([self.z[i].x for i in range(len(self.z))])

            #with linearization:
            #self.subproblem.setObjective(np.dot(self.raw_out["not_kept"], self.ma0.values()) + np.dot(self.raw_out["help"], self.ma1.values()) - np.dot(self.a.values(), [self.dmus[i].pi for i in range(self.n)]))
            #base case
            self.updSub(niter=niter * 0)
            self.subproblem.setObjective(np.dot(self.mu.values(), self.y0) - np.dot(self.a.values(), [self.dmus[i].pi for i in range(self.n)]), GRB.MINIMIZE)
            self.subproblem.update()
            start = time.time()
            self.subproblem.optimize(callback)
            stop = time.time()
            self.subTime.append(stop - start)
            self.red_csts.append(self.subproblem.objVal)
            self.add_row()

            if len(self.last10) >= 10:
                del self.last10[0]
                if (self.last10[0] - self.last10[-1])/self.last10[-1] <= 0.05:
                    break
            self.last10.append(self.master.objVal)
            
            niter += 1
            print("Iteration %d, time elapsed %g" % (niter, time.time() - glob_start))
            if niter == maxiter:
                break
        
    def solve_h(self, lpIter=500):

        self.solveLP(maxiter=lpIter)
        lpObj = self.master.objVal
        
        for zVar in self.z:
            zVar.vtype = GRB.BINARY

        self.master.optimize()
        self.ISol = [self.z[i].x for i in range(len(self.z))]
        mipObj = self.master.objVal
        print(mipObj)
        print(lpObj)
        gap = (mipObj - lpObj)/mipObj * 100
        print(gap)