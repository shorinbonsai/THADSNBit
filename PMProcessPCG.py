import sys
import os
import csv
import math
import re
import string
from operator import itemgetter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

wordCheck = "-fitness\n"
finame = "best.lint"
outp = "./Output/PMpcgTest/"
# os.mkdir(outp)
samps = 30
lower_better = True
precision = 6
col_width = 6 + precision


def getFits(dir_path: str, ascending: bool):
    data = []
    with open(dir_path) as f:
        lines = f.readlines()
        for line in lines:
            spl = line.split(' ')
            if wordCheck in spl:
                data.append(float(spl[0]))
        pass
    data.sort()  # Ascending
    if not ascending:
        data.reverse()
        pass
    return data


def writeStat(data: [], out):
    mean = float(np.mean(data))
    mean = round(mean, precision)
    std = float(np.std(data, ddof=0))
    std = round(std, precision)  # Population standard deviation
    diff = 1.96 * std / math.sqrt(30)  # 95% CI
    diff = round(diff, precision)
    if lower_better:
        maxima = float(min(data))
        pass
    else:
        maxima = float(max(data))
        pass
    maxima = round(maxima, precision)
    out.write(str(mean).ljust(col_width))
    out.write(str(std).ljust(col_width))
    out.write(u'\u00B1' + str(diff).ljust(col_width))
    out.write(str(maxima).ljust(col_width))
    return mean, maxima

#ED
def toFile(data: [], fname: string, exp: int, stats: [], inf: str):
    out = open(outp + fname + str(exp) + ".dat", "w")
    out.write(inf + '\n')
    out.write("Mean: " + str(stats[0]) + '\n')
    out.write("Best: " + str(stats[1]) + '\n')
    for d in data:
        out.write(str(d) + "\n")
        pass
    out.close()
    pass

#Profiles
# def toFile(data: [], fname: string, exp: int):
#     out = open(outp + fname + str(exp) + ".dat", "w")
#     for d in data:
#         out.write(str(d) + "\n")
#         pass
#     out.close()
#     pass


def main():
    alphas = [30]
    edges = [3]
    # densities = ["ps1","ps2","ps3","ps4","ps5","ps6","ps7","ps8"]
    densities = ["ps1","ps2","ps3","ps4","ps5","ps6","ps7","ps8","ps9","ps10","ps11","ps12","ps13","ps14","ps15","ps16","ps17","ps18","ps19","ps20","ps21","ps22","ps23","ps24","ps25","ps26","ps27","ps28","ps29"]
    profiles = [1,2,3,4,5,6,7,8,9]

    # data = [[[[] for _ in range(8)] for _ in range(len(edges))] for _ in range(len(alphas))]
    # print(data)
    data2 = [[[[[] for _ in range(len(densities))] for _ in range(len(edges))] for _ in range(len(profiles))] for _ in range(len(alphas))]
    tables = []
    tab_root = outp + "table_"
    # for pridx, pr_dat in enumerate(alphas):
    #     out_name = tab_root + "P" + str(pr_dat) + ".dat"
    #     tables.append(open(out_name, "w"))
        
    for dirpath, dirnames, files in os.walk('.'):
        for file_name in files:
            if file_name.endswith(".lint"):
                direc = os.path.join(dirpath, file_name)
                exper = []
                exper = dirpath.split('/')
                if exper[-4][-1] =='t' and exper[-4][0] == 'P':
                    if exper[-3][0] == '3':
                        alp = 0
                    if exper[-3][0] == '4':
                        alp = 1
                    edg = 0
                    if exper[-4][2] == '1':
                        edg = 0
                    if exper[-4][2] == '3':
                        edg = 1
                    if exper[-4][2] == '5':
                        edg = 2
                    if "P1" in exper[-1]:
                        pr = 0
                    if "P2" in exper[-1]:
                        pr = 1
                    if "P3" in exper[-1]:
                        pr = 2
                    if "P4" in exper[-1]:
                        pr = 3
                    if "P5" in exper[-1]:
                        pr = 4
                    if "P6" in exper[-1]:
                        pr = 5
                    if "P7" in exper[-1]:
                        pr = 6
                    if "P8" in exper[-1]:
                        pr = 7
                    if "P9" in exper[-1]:
                        pr = 8
                    ps = int(exper[-2][-1])-1
                    data2[alp][pr][edg][ps] = getFits(direc, lower_better)
                    # print(data2)

#Profiles
    col_ws = [7, 5, 6, 4, 5]
    # means = [[] for _ in range(len(profiles))]
    # bests = [[] for _ in range(len(profiles))]
    means = []
    bests = []
    exp = 1
    data_1d = []
    # for pridx, pr_dat in enumerate(data):
    #     for alidx in range(len(alphas)):
    #         col_idx = 0
    #         tables[pridx].write(base_lbls[bidx].ljust(col_ws[col_idx]))
    
    out = open(outp + "PM_pcg.dat", "w", encoding='utf-16')
    out.write("EXP".ljust(col_ws[0]))
    out.write("A".ljust(col_ws[1]))
    out.write("Prof".ljust(col_ws[2]))
    out.write("E".ljust(col_ws[3]))
    out.write("PS".ljust(col_ws[4]))
    out.write("Mean".ljust(col_width))
    out.write("SD".ljust(col_width))
    out.write("95%CI".ljust(col_width))
    out.write("Best".ljust(col_width))
    out.write('\n')
    for alindx, al_dat in enumerate(data2):
        col_idx = 1
        al_info = str(str(alphas[alindx]) + "A").ljust(col_ws[col_idx])
        for pridx, pr_dat in enumerate(al_dat):
            col_idx = 2
            pr_info = al_info + "P" + str(profiles[pridx]).ljust(col_ws[col_idx])
            for edidx, ed_dat in enumerate(pr_dat):
                col_idx = 3
                ed_info = pr_info + str(str(edges[edidx]) + "E").ljust(col_ws[col_idx])
                for psidx, dat in enumerate(ed_dat):
                    col_idx = 4
                    if len(dat) < 30:
                        print(pr_info + str(densities[psidx]))
                   
                    assert len(dat) == samps
                    # raise ValueError("ValueError exception thrown")
                    data_1d.append(dat)
                    all_info = ed_info + str(densities[psidx]).ljust(col_ws[col_idx])
                    out.write(str("EXP" + str(exp)).ljust(col_ws[0]))
                    out.write(all_info)

                    vals = writeStat(dat, out)
                    out.write('\n')
                    toFile(dat, "EXP", exp, vals, all_info)
                    exp += 1
                    means.append([vals[0], [alindx,pridx, edidx, psidx]])
                    bests.append([vals[1], [alindx,pridx, edidx, psidx]])
                    pass
                pass
            pass
        pass
    out.close()




    means.sort(key=itemgetter(0))
    bests.sort(key=itemgetter(0))
                    
    # out = open(outp + "EDpcgbest.dat", "w")
    # # for idx in range(len(profs)):
    # # out.write("Profile " + str(profs[idx]) + "\n")
    # out.write("Best Mean: ")
    # out.write(str(means[0][0]) + "\n")
    # out.write("Alpha: " + str(alphas[means[0][1][0]]) + "\n")
    # out.write("Edges: " + str(edges[means[0][1][1]]) + "\n")
    # out.write("PS: " + str(densities[means[0][1][2]]) + "\n")
    # out.write("Best Fitness: ")
    # out.write(str(bests[0][0]) + "\n")
    # out.write("Alpha: " + str(alphas[bests[0][1][0]]) + "\n")
    # out.write("Edges: " + str(edges[bests[0][1][1]]) + "\n")
    # out.write("PS: " + str(densities[bests[0][1][2]]) + "\n")
    # out.write("\n")
    # out.close()


    all_data30 = [[] for _ in range(len(profiles))]
    x_labels30 = [[] for _ in range(len(profiles))]
    all_data40 = [[] for _ in range(len(profiles))]
    x_labels40 = [[] for _ in range(len(profiles))]
    for pridx, pr in enumerate(profiles):
        for psidx, ps in enumerate(densities):
            for edidx, edg in enumerate(edges):
                all_data30[pridx].append(data2[0][pridx][edidx][psidx])
                x_labels30[pridx].append("ps=" + ps + " E=" + str(edg))
    for pridx, pr in enumerate(profiles):
        for psidx, ps in enumerate(densities):
            for edidx, edg in enumerate(edges):
                all_data40[pridx].append(data2[1][pridx][edidx][psidx])
                x_labels40[pridx].append("ps=" + ps + " E=" + str(edg))

    fig_root = outp + "boxplot_30Alph_P"
    fig_titles = ["Profile 1 - 30alpha PCG", "Profile 2 - 30alpha PCG", "Profile 3 - 30alpha PCG", \
         "Profile 4 - 30alpha PCG", "Profile 5 - 30alpha PCG", "Profile 6 - 30alpha PCG", "Profile 7 - 30alpha PCG", \
         "Profile 8 - 30alpha PCG", "Profile 9 - 30alpha PCG"]

    fig_titles2 = ["Profile 1 - 40alpha PCG", "Profile 2 - 40alpha PCG", "Profile 3 - 40alpha PCG", \
         "Profile 4 - 40alpha PCG", "Profile 5 - 40alpha PCG", "Profile 6 - 40alpha PCG", "Profile 7 - 40alpha PCG", \
         "Profile 8 - 40alpha PCG", "Profile 9 - 40alpha PCG"]

    fig_titles3 = ["Profile 1 - PCG", "Profile 2 - PCG", "Profile 3 - PCG", \
         "Profile 4 - PCG", "Profile 5 - PCG", "Profile 6 - PCG", "Profile 7 - PCG", \
         "Profile 8 - PCG", "Profile 9 - PCG"]
    
    plt.rc('xtick', labelsize=4)
    plt.rc('ytick', labelsize=6)
    figs = [plt.figure() for _ in range(len(profiles))]
    for f in figs:
        f.set_dpi(600)
        f.set_figheight(3)

    plts = [figs[idx].add_subplot(111) for idx in range(len(profiles))]
    for idx in range(len(all_data30)):
        plts[idx].boxplot(all_data30[idx])
        plts[idx].set_xticklabels(x_labels30[idx], rotation=75)

        figs[idx].suptitle(fig_titles[idx], fontsize=12)
        plts[idx].set_xlabel("Parameter Setting (E=Starting Edges, ps=Density Profiles)", fontsize=10)
        plts[idx].set_ylabel("PM - Fitness", fontsize=10)
        
        figs[idx].subplots_adjust(left=.08, bottom=.1, right=.98, top=.91, wspace=0, hspace=0)
        figs[idx].tight_layout()
        figs[idx].savefig(fig_root + str(profiles[idx]) + ".png")

        ################
    fig_root2 = outp + "boxplot_40Alph_P"
    plt.rc('xtick', labelsize=4)
    plt.rc('ytick', labelsize=6)
    figs = [plt.figure() for _ in range(len(profiles))]
    for f in figs:
        f.set_dpi(600)
        f.set_figheight(3)

    plts = [figs[idx].add_subplot(111) for idx in range(len(profiles))]
    for idx in range(len(all_data40)):
        plts[idx].boxplot(all_data40[idx])
        plts[idx].set_xticklabels(x_labels40[idx], rotation=75)

        figs[idx].suptitle(fig_titles2[idx], fontsize=12)
        plts[idx].set_xlabel("Parameter Setting (E=Starting Edges, ps=Density Profiles)", fontsize=10)
        plts[idx].set_ylabel("PM - Fitness", fontsize=10)
        
        figs[idx].subplots_adjust(left=.08, bottom=.1, right=.98, top=.91, wspace=0, hspace=0)
        figs[idx].tight_layout()
        figs[idx].savefig(fig_root2 + str(profiles[idx]) + ".png")

    # x_labels = []
    # for stidx, st in enumerate(alphas):
    #     for poidx, p in enumerate(edges):
    #         for midx, m in enumerate(densities):
    #             x_labels.append(" E=" + str(p) + " PS=" + str(m))
    #             pass
    #         pass
    #     pass

    # fig_root = outp + "boxplotEDpcg"
    # plt.rc('xtick', labelsize=4)
    # plt.rc('ytick', labelsize=6)
    # f = plt.figure()
    # f.set_dpi(600)
    # f.set_figheight(5)
    # plot = f.add_subplot(111)

    # plot.boxplot(data_1d[:24])
    # plot.set_xticklabels(x_labels, rotation=75)
    # f.suptitle("ED PCG 30 Alpha", fontsize=10)
    # plot.set_xlabel("Parameter Setting (E=Edges, PS=densities)", fontsize=8)
    # plot.set_ylabel("Distribution of Fitness", fontsize=8)
    # f.subplots_adjust(left=.08, bottom=.1, right=.98, top=.91, wspace=0, hspace=0)
    # f.savefig(fig_root + "30.png")



    # fig_root = outp + "boxplotEDpcg"
    # plt.rc('xtick', labelsize=4)
    # plt.rc('ytick', labelsize=6)
    # f = plt.figure()
    # f.set_dpi(600)
    # f.set_figheight(5)
    # plot = f.add_subplot(111)

    # plot.boxplot(data_1d[24:])
    # plot.set_xticklabels(x_labels, rotation=75)
    # f.suptitle("ED PCG 40 Alpha", fontsize=10)
    # plot.set_xlabel("Parameter Setting (E=Edges, PS=densities)", fontsize=8)
    # plot.set_ylabel("Distribution of Fitness", fontsize=8)
    # f.subplots_adjust(left=.08, bottom=.1, right=.98, top=.91, wspace=0, hspace=0)
    # f.savefig(fig_root + "40.png")
    return 0

        
if '__main__' == __name__:
    main()