import math
import re
import string
from operator import itemgetter

import matplotlib.pyplot as plt
import numpy as np

"""
Processes folders of output in the form "Output - [Profile Number] w [States], [Population Size], [Mutations]" to 
create boxplots, tables, fitness files, and information about the best parameter settings.
"""

inp = "C:/Users/micha/OneDrive - University of Guelph/Coding Projects/BitSprayer/Output/"
outp = "./Output/"
finame = "best.lint"
samps = 30
lower_better = False
precision = 6
col_width = 6 + precision


def getFits(dir_path: str, ascending: bool):
    data = []
    with open(dir_path + finame) as f:
        lines = f.readlines()
        for line in lines:
            if line.__contains__(str("fitness")):
                data.append(float(re.findall("\d+\.\d+", line)[0]))
                pass
            pass
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


def main():
    states = [12, 16]
    pops = ["024", "036", "048"]
    muts = [1, 2]

    # Collect the data
    f_root = inp + "Output - "

    # Create data[profile][states][popsize][mutations] = 30 fitness vals
    # tables = []
    # tab_root = outp + "table_"
    data = [[] for _ in range(len(states))]
    for stidx, s_dat in enumerate(states):
        f_st = f_root + str(s_dat) + "S, "
        data[stidx] = [[] for _ in range(len(pops))]
        for poidx, po_dat in enumerate(pops):
            f_pop = f_st + str(po_dat) + "P, "
            data[stidx][poidx] = [[] for _ in range(len(muts))]
            for midx, m in enumerate(muts):
                folder = f_pop + str(m) + "M/"
                data[stidx][poidx][midx] = getFits(folder, lower_better)
                pass
            pass
        pass

    # Process the data and make table
    col_ws = [6, 5, 6, 4]
    means = []
    bests = []
    exp = 1
    data_1d = []
    out = open(outp + "EXP Summary.dat", "w", encoding='utf-16')
    out.write("EXP".ljust(col_ws[0]))
    out.write("S".ljust(col_ws[1]))
    out.write("P".ljust(col_ws[2]))
    out.write("M".ljust(col_ws[3]))
    out.write("Mean".ljust(col_width))
    out.write("SD".ljust(col_width))
    out.write("95%CI".ljust(col_width))
    out.write("Best".ljust(col_width))
    out.write('\n')
    for stidx, s_dat in enumerate(data):
        col_idx = 1
        st_info = str(str(states[stidx]) + "S").ljust(col_ws[col_idx])
        for poidx, po_dat in enumerate(s_dat):
            col_idx = 2
            po_info = st_info + str(str(pops[poidx]) + "P").ljust(col_ws[col_idx])
            for midx, dat in enumerate(po_dat):
                col_idx = 3
                assert len(dat) == samps
                data_1d.append(dat)
                all_info = po_info + str(str(muts[midx]) + "M").ljust(col_ws[col_idx])
                out.write(str("EXP" + str(exp)).ljust(col_ws[0]))
                out.write(all_info)
                vals = writeStat(dat, out)
                out.write('\n')
                toFile(dat, "EXP", exp, vals, all_info)
                exp += 1
                means.append([vals[0], [stidx, poidx, midx]])
                bests.append([vals[1], [stidx, poidx, midx]])
                pass
            pass
        pass
    out.close()

    means.sort(key=itemgetter(0))
    bests.sort(key=itemgetter(0))

    out = open(outp + "best.dat", "w")
    # for idx in range(len(profs)):
    # out.write("Profile " + str(profs[idx]) + "\n")
    out.write("Best Mean: ")
    out.write(str(means[0][0]) + "\n")
    out.write("States: " + str(states[means[0][1][0]]) + "\n")
    out.write("Population: " + str(pops[means[0][1][1]]) + "\n")
    out.write("Mutations: " + str(muts[means[0][1][2]]) + "\n")
    out.write("Best Fitness: ")
    out.write(str(bests[0][0]) + "\n")
    out.write("States: " + str(states[bests[0][1][0]]) + "\n")
    out.write("Population: " + str(pops[bests[0][1][1]]) + "\n")
    out.write("Mutations: " + str(muts[bests[0][1][2]]) + "\n")
    out.write("\n")
    out.close()

    x_labels = []
    for stidx, st in enumerate(states):
        for poidx, p in enumerate(pops):
            for midx, m in enumerate(muts):
                x_labels.append("S=" + str(st) + " P=" + str(p) + " M=" + str(m))
                pass
            pass
        pass

    fig_root = outp + "boxplot"
    plt.rc('xtick', labelsize=6)
    plt.rc('ytick', labelsize=6)

    f = plt.figure()
    f.set_dpi(500)
    f.set_figheight(3)

    plot = f.add_subplot(111)
    plot.boxplot(data_1d)
    plot.set_xticklabels(x_labels, rotation=90, va='bottom')

    f.suptitle("Dungeon", fontsize=12)
    plot.set_xlabel("Parameter Setting (S=states, P=population, M=max. mutations)", fontsize=10)
    plot.set_ylabel("Distribution of Fitness", fontsize=10)
    # figs[idx].tight_layout()
    f.subplots_adjust(left=.08, bottom=.1, right=.98, top=.91, wspace=0, hspace=0)
    f.savefig(fig_root + ".png")
    return 0


main()
