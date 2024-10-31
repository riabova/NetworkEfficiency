import numpy as np
import pandas as pd
import time
# from matplotlib import pyplot as plt
from gurobipy import *
# import folium
# from folium.features import DivIcon
from cg_model import DEAcg
import warnings
warnings.simplefilter("ignore", UserWarning)

def getDMUs(ccn):
    return set(branches[branches["cc_num"] == ccn].index)

if __name__ == "__main__":

    df_int = pd.read_csv("..\data\distr_int.csv", dtype={"BUSINESS_UNIT": "int64", "SHIPTO_NO": "int64"})
    df_ext = pd.read_csv("..\data\distr_ext.csv", dtype={"BUSINESS_UNIT": "int64", "SHIPTO_NO": "int64"})

    branches = pd.read_csv("..\data\\branches_cc_cur80_sz.csv", dtype={"cc_num": "int64"}).set_index("dmu")

    cc_sizes = []
    for cc in list(set(branches["cc_num"])):
        cc_sizes.append(len(branches[branches["cc_num"] == cc]))

    cc_nums = list(set(branches["cc_num"]))

    for cc in [16]:
        raw_in = pd.read_csv("..\data\\raw_in.csv").set_index("dmu")
        raw_out = pd.read_csv("..\data\\raw_out.csv").set_index("dmu")
        groups = pd.read_csv("..\data\groups_clean.csv").set_index("dmu")
        groups["group0"] = groups.index
        dmus = getDMUs(cc)
        raw_in = raw_in.loc[list(dmus.intersection(raw_in.index))]
        raw_out = raw_out.loc[list(dmus.intersection(raw_in.index))]
        groups = groups.loc[list(dmus.intersection(raw_in.index))]
        
        dea = DEAcg(raw_in, raw_out, groups, df_int)
        dea.setup()
        dea.solve_h(lpIter=50)
        
        with open("..\out\scores_reg%d.txt" % cc, "w") as f:
            for score in dea.c:
                f.write(str(np.round(score, 2)) + "\n")
        with open("..\out\master_obj_reg%d.txt" % cc, "w") as f:
            for objval in [np.dot(dea.c[:len(dea.master_sols[i])], dea.master_sols[i]) for i in range(len(dea.master_sols))]:
                f.write(str(np.round(objval, 2)) + "\n")
        with open("..\out\\red_csts_reg%d.txt" % cc, "w") as f:
            for rc in dea.red_csts:
                f.write(str(np.round(rc, 2)) + "\n")
        with open("..\out\sub_time_reg%d.txt" % cc, "w") as f:
            for t in dea.subTime:
                f.write(str(np.round(t, 2)) + "\n")
        dea.ai.to_csv("..\out\\ai_reg%d.csv" % cc)
        dea.X.to_csv("..\out\X_reg%d.csv" % cc)
        dea.Y.to_csv("..\out\Y_reg%d.csv" % cc)
        with open("..\out\lp_sol_reg%d.txt" % cc, "w") as f:
            for x in dea.master_sols[-1]:
                f.write(str(x) + "\n")
        with open("..\out\ip_sol_reg%d.txt" % cc, "w") as f:
            for z in dea.z:
                f.write(str(int(z.x)) + "\n")