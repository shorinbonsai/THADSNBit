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
outp = "../Output/EDpcgTest/"
# os.mkdir(outp)
samps = 30
lower_better = False
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

def init_fit(dir_path: str):
    with open(dir_path) as f:
        lines = f.readlines()
        line = lines[0].strip()
        first = line.split(' ')
        result = first[-1]
        return result



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
    alphas = [30, 40]
    edges = [1, 3, 5]
    # densities = ["ps1","ps2","ps3","ps4","ps5","ps6","ps7","ps8"]
    densities = ["ps1","ps2","ps3","ps4","ps5","ps6","ps7","ps8","ps9","ps10","ps11","ps12","ps13","ps14","ps15","ps16","ps17","ps18","ps19","ps20","ps21","ps22","ps23","ps24","ps25","ps26","ps27","ps28","ps29"]
    profiles = ["pr1", "pr2", "pr3", "pr4", "pr5", "pr6", "pr7", "pr8", "pr9"]
    initfile = "initGraph.dat"

    data = [[] for _ in range(len(densities))]
    print(data)
    data2 = [[[[[] for _ in range(8)] for _ in range(len(edges))] for _ in range(9)] for _ in range(len(alphas))]
    init_data = []
    
    #ED
    # for dirpath, dirnames, files in os.walk('.'):
    #     for file_name in files:
    #         if file_name.endswith(".lint"):
    #             direc = os.path.join(dirpath, file_name)
    #             input_file = open(direc,"r")
    #             exper = []
    #             exper = dirpath.split('/')
    #             if exper[-4][-1] =='r' and exper[-4][0] == 'E':
    #                 if exper[-3][0] == '3':
    #                     alp = 0
    #                 if exper[-3][0] == '4':
    #                     alp = 1
    #                 if exper[-4][2] == '1':
    #                     edg = 0
    #                 if exper[-4][2] == '3':
    #                     edg = 1
    #                 if exper[-4][2] == '5':
    #                     edg = 2
    #                 ps = int(exper[-2][-1])-1
    #                 data[alp][edg][ps] = getFits(direc, lower_better)
                        # print(data)
                        # print(direc)
                        # raise ValueError("ValueError exception thrown")
    for dirpath, dirnames, files in os.walk('.'):
        for file_name in files:
            if file_name.endswith(".lint"):
                direc = os.path.join(dirpath, file_name)
                # input_file = open(direc,"r")
                exper = []
                exper = dirpath.split('/')
                ps = int(exper[-2][2:])-1
                direc2 = os.path.join(dirpath, initfile)
                init_data.append(float(init_fit(direc2)))
                data[ps] = getFits(direc, lower_better)

    # for dirpath, dirnames, files in os.walk('.'):
    #     for file_name in files:
    #         if file_name.endswith(".lint"):
    #             direc = os.path.join(dirpath, file_name)
    #             input_file = open(direc,"r")
    #             exper = []
    #             exper = dirpath.split('/')
    #             if exper[-4][-1] =='r' and exper[-4][0] == 'P':
    #                 if exper[-3][0] == '3':
    #                     alp = 0
    #                 if exper[-3][0] == '4':
    #                     alp = 1
    #                 if exper[-4][2] == '1':
    #                     edg = 0
    #                 if exper[-4][2] == '3':
    #                     edg = 1
    #                 if exper[-4][2] == '5':
    #                     edg = 2
    #                 if "P1" in exper[-1]:
    #                     pr = 0
    #                 if "P2" in exper[-1]:
    #                     pr = 1
    #                 if "P3" in exper[-1]:
    #                     pr = 2
    #                 if "P4" in exper[-1]:
    #                     pr = 3
    #                 if "P5" in exper[-1]:
    #                     pr = 4
    #                 if "P6" in exper[-1]:
    #                     pr = 5
    #                 if "P7" in exper[-1]:
    #                     pr = 6
    #                 if "P8" in exper[-1]:
    #                     pr = 7
    #                 if "P9" in exper[-1]:
    #                     pr = 8
    #                 ps = int(exper[-2][-1])-1
    #                 data2[alp][pr][edg][ps] = getFits(direc, lower_better)
                    # print(data2)

#Profiles
    # col_ws = [6, 5, 6, 4, 5]
    # means = []
    # bests = []
    # exp = 1
    # data_1d = []
    # out = open(outp + "PM_ring.dat", "w", encoding='utf-16')
    # out.write("EXP".ljust(col_ws[0]))
    # out.write("A".ljust(col_ws[1]))
    # out.write("E".ljust(col_ws[2]))
    # out.write("PS".ljust(col_ws[3]))
    # out.write("Prof".ljust(col_ws[4]))
    # out.write("Mean".ljust(col_width))
    # out.write("SD".ljust(col_width))
    # out.write("95%CI".ljust(col_width))
    # out.write("Best".ljust(col_width))
    # out.write('\n')
    # for alindx, al_dat in enumerate(data2):
    #     col_idx = 1
    #     al_info = str(str(alphas[alindx]) + "A").ljust(col_ws[col_idx])
    #     for pridx, pr_dat in enumerate(al_dat):
    #         col_idx = 3
    #         pr_info = al_info + "P" + str(profiles[pridx]).ljust(col_ws[col_idx])
    #         for edidx, ed_dat in enumerate(pr_dat):
    #             col_idx = 2
    #             ed_info = al_info + str(str(edges[edidx]) + "E").ljust(col_ws[col_idx])
    #             for psidx, dat in enumerate(pr_dat):
    #                 col_idx = 4
    #                 if len(dat) < 30:
    #                     print(pr_info + str(densities[psidx]))

                    

    #                 assert len(dat) == samps
    #                 data_1d.append(dat)
    #                 all_info = pr_info + str(densities[psidx]).ljust(col_ws[col_idx])
    #                 out.write(str("EXP" + str(exp)).ljust(col_ws[0]))
    #                 out.write(all_info)
    #                 vals = writeStat(dat, out)
    #                 out.write('\n')
    #                 toFile(dat, "EXP", exp, vals, all_info)
    #                 exp += 1
    #                 means.append([vals[0], [alindx, edidx, psidx]])
    #                 bests.append([vals[1], [alindx, edidx, psidx]])
    #                 pass
    #             pass
    #         pass
    #     pass
    # out.close()


#ED
    # col_ws = [6, 5, 6, 4]
    # means = []
    # bests = []
    # exp = 1
    # data_1d = []
    # out = open(outp + "ED_pcgTest.dat", "w", encoding='utf-16')
    # out.write("EXP".ljust(col_ws[0]))
    # # out.write("A".ljust(col_ws[1]))
    # # out.write("E".ljust(col_ws[2]))
    # out.write("PS".ljust(col_ws[3]))
    # out.write("Mean".ljust(col_width))
    # out.write("SD".ljust(col_width))
    # out.write("95%CI".ljust(col_width))
    # out.write("Best".ljust(col_width))
    # out.write('\n')
    # for alindx, al_dat in enumerate(data):
    #     col_idx = 1
    #     al_info = str(str(alphas[alindx]) + "A").ljust(col_ws[col_idx])
    #     for edidx, ed_dat in enumerate(al_dat):
    #         col_idx = 2
    #         ed_info = al_info + str(str(edges[edidx]) + "E").ljust(col_ws[col_idx])
    #         for psidx, dat in enumerate(ed_dat):
    #             col_idx = 3
    #             assert len(dat) == samps
    #             data_1d.append(dat)
    #             all_info = ed_info + str(densities[psidx]).ljust(col_ws[col_idx])
    #             out.write(str("EXP" + str(exp)).ljust(col_ws[0]))
    #             out.write(all_info)
    #             vals = writeStat(dat, out)
    #             out.write('\n')
    #             toFile(dat, "EXP", exp, vals, all_info)
    #             exp += 1
    #             means.append([vals[0], [alindx, edidx, psidx]])
    #             bests.append([vals[1], [alindx, edidx, psidx]])
    #             pass
    #         pass
    #     pass
    # out.close()
    # print(data)
    col_ws = [6, 5, 6, 4]
    means = []
    bests = []
    exp = 1
    data_1d = []
    out = open(outp + "ED_pcgTest.dat", "w", encoding='utf-16')
    out.write("EXP".ljust(col_ws[0]))
    # out.write("A".ljust(col_ws[1]))
    # out.write("E".ljust(col_ws[2]))
    out.write("PS".ljust(col_ws[3]))
    out.write("Mean".ljust(col_width))
    out.write("SD".ljust(col_width))
    out.write("95%CI".ljust(col_width))
    out.write("Best".ljust(col_width))
    out.write('\n')
    for psidx, dat in enumerate(data):
        col_idx = 1
        all_info = str(densities[psidx]).ljust(col_ws[col_idx])
        out.write(str("EXP" + str(exp)).ljust(col_ws[0]))
        out.write(all_info)
        if len(dat) < 30:
            print(str(densities[psidx]) + str(len(dat))+ "   "+ str(dat))
        assert len(dat) == samps
        data_1d.append(dat)
        vals = writeStat(dat, out)
        out.write('\n')
        toFile(dat, "EXP", exp, vals, all_info)
        exp += 1
        means.append([vals[0], [psidx]])
        bests.append([vals[1], [psidx]])
        pass
    out.close()

    means.sort(key=itemgetter(0), reverse=True)
    bests.sort(key=itemgetter(0), reverse=True)
                    
    out = open(outp + "EDpcgbesttest.dat", "w")
    # for idx in range(len(profs)):
    # out.write("Profile " + str(profs[idx]) + "\n")
    out.write("Top 8 Best Mean: " + "\n")
    for i in range(8):
        out.write("Mean  " + str(means[i][0]) )
        print("testing" + str(means[i][1][0]))
        out.write("   PS: " + str(densities[means[i][1][0]]) + "\n")
    out.write("Top 8 Best Fitness: " + "\n")
    for i in range(8):
        out.write("Fitness  " + str(bests[i][0]) )
        out.write("   PS: " + str(densities[bests[i][1][0]]) + "\n")
    out.write("\n")
    out.close()

    x_labels = []
    for psidx, ps in enumerate(densities):
        x_labels.append("PS=" + str(ps))
    fig_root = outp + "boxplotEDpcg"
    # x_labels = []
    # for stidx, st in enumerate(alphas):
    #     for poidx, p in enumerate(edges):
    #         for midx, m in enumerate(densities):
    #             x_labels.append(" E=" + str(p) + " PS=" + str(m))
    #             pass
    #         pass
    #     pass

    # fig_root = outp + "boxplotEDpcg"
    plt.rc('xtick', labelsize=4)
    plt.rc('ytick', labelsize=6)
    f = plt.figure()
    f.set_dpi(600)
    f.set_figheight(5)
    plot = f.add_subplot(111)

    plot.boxplot(data_1d)
    plt.axhline(y=init_data[0], color='r', linestyle='-')
    plot.set_xticklabels(x_labels, rotation=75)
    f.suptitle("ED PCG Densities test", fontsize=10)
    plot.set_xlabel("Parameter Setting (PS=densities)", fontsize=8)
    plot.set_ylabel("Distribution of Fitness", fontsize=8)
    f.subplots_adjust(left=.08, bottom=.1, right=.98, top=.91, wspace=0, hspace=0)
    f.savefig(fig_root + ".png")



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
    # return 0

        
if '__main__' == __name__:
    main()