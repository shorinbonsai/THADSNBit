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

inp = "./Input/"
outp = "./Output/"
finame = "best.lint"
samps = 30
lower_better = True
precision = 6
col_width = 4 + precision


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
    out.write(str(diff).ljust(col_width))
    out.write(str(maxima).ljust(col_width))
    return mean, maxima


def toFile(data: [], fname: string, exp: int):
    out = open(outp + fname + str(exp) + ".dat", "w")
    for d in data:
        out.write(str(d) + "\n")
        pass
    out.close()
    pass


def main():
    profs = [1, 7]
    states = [8, 12, 16]
    pops = [10, 100, 1000]
    muts = [1, 2, 3, 4]

    # Collect the data
    f_root = inp + "Output - "
    base_dir_strs = ["Ring Baseline Profile ", "PLC Baseline Profile "]
    base_data = [[] for _ in range(len(profs))]
    base_lbls = ["Ring", "PLC"]
    for pridx, pr_dat in enumerate(profs):
        for bstr in base_dir_strs:
            base_data[pridx].append(getFits(f_root + bstr + str(pr_dat) + "/", True))
            pass
        pass

    # Create data[profile][states][popsize][mutations] = 30 fitness vals
    data = [[] for _ in range(len(profs))]
    tables = []
    tab_root = outp + "table_"
    for pridx, pr_dat in enumerate(profs):
        f_prof = f_root + "P" + str(pr_dat) + " w "
        out_name = tab_root + "P" + str(pr_dat) + ".dat"
        tables.append(open(out_name, "w"))
        data[pridx] = [[] for _ in range(len(states))]
        for stidx, s_dat in enumerate(states):
            f_st = f_prof + str(s_dat) + "S, "
            data[pridx][stidx] = [[] for _ in range(len(pops))]
            for poidx, po_dat in enumerate(pops):
                f_pop = f_st + str(po_dat) + "P, "
                data[pridx][stidx][poidx] = [[] for _ in range(len(muts))]
                for midx, m in enumerate(muts):
                    folder = f_pop + str(m) + "M/"
                    data[pridx][stidx][poidx][midx] = getFits(folder, lower_better)
                    pass
                pass
            pass
        pass

    # Process the data and make tables
    col_ws = [6, 5, 7, 4]
    means = [[] for _ in range(len(profs))]
    bests = [[] for _ in range(len(profs))]
    exp = 1
    for pridx, pr_dat in enumerate(data):
        for bidx in range(len(base_dir_strs)):
            col_idx = 0
            tables[pridx].write(base_lbls[bidx].ljust(col_ws[col_idx]))
            col_idx += 1
            for _ in range(1, 4):
                tables[pridx].write("NA".ljust(col_ws[col_idx]))
                col_idx += 1
                pass
            writeStat(base_data[pridx][bidx], tables[pridx])
            tables[pridx].write("\n")
            toFile(base_data[pridx][bidx], "BASE" + str(profs[pridx]) + "_", bidx - 1)
            pass
        col_idx = 0
        pr_info = str("P" + str(profs[pridx])).ljust(col_ws[col_idx])
        for stidx, s_dat in enumerate(pr_dat):
            col_idx = 1
            st_info = pr_info + str(str(states[stidx]) + "S").ljust(col_ws[col_idx])
            for poidx, po_dat in enumerate(s_dat):
                col_idx = 2
                po_info = st_info + str(str(pops[poidx]) + "P").ljust(col_ws[col_idx])
                for midx, dat in enumerate(po_dat):
                    col_idx = 3
                    all_info = po_info + str(str(muts[midx]) + "M").ljust(col_ws[col_idx])
                    tables[pridx].write(all_info)
                    vals = writeStat(dat, tables[pridx])
                    assert len(dat) == samps
                    toFile(dat, "EXP", exp)
                    exp += 1
                    means[pridx].append([vals[0], [pridx, stidx, poidx, midx]])
                    bests[pridx].append([vals[1], [pridx, stidx, poidx, midx]])
                    tables[pridx].write("\n")
                    pass
                pass
            pass
        pass

    for f in tables:
        f.close()
        pass

    for li in means:
        li.sort(key=itemgetter(0))
        pass
    for li in bests:
        li.sort(key=itemgetter(0))
        pass

    out = open(outp + "best.dat", "w")
    for idx in range(len(profs)):
        out.write("Profile " + str(profs[idx]) + "\n")
        out.write("Best Mean: ")
        out.write(str(means[idx][0][0]) + "\n")
        out.write("States: " + str(states[means[idx][0][1][1]]) + "\n")
        out.write("Population: " + str(pops[means[idx][0][1][2]]) + "\n")
        out.write("Mutations: " + str(muts[means[idx][0][1][3]]) + "\n")
        out.write("Best Fitness: ")
        out.write(str(bests[idx][0][0]) + "\n")
        out.write("States: " + str(states[bests[idx][0][1][1]]) + "\n")
        out.write("Population: " + str(pops[bests[idx][0][1][2]]) + "\n")
        out.write("Mutations: " + str(muts[bests[idx][0][1][3]]) + "\n")
        out.write("\n")
        pass
    out.close()

    all_data = [[] for _ in range(len(profs))]
    x_labels = [[] for _ in range(len(profs))]
    for pridx, pr in enumerate(profs):
        for bidx in range(len(base_lbls)):
            all_data[pridx].append(base_data[pridx][bidx])
            x_labels[pridx].append(base_lbls[bidx] + " Baseline")
            pass
        for stidx, st in enumerate(states):
            for poidx, p in enumerate(pops):
                for midx, m in enumerate(muts):
                    all_data[pridx].append(data[pridx][stidx][poidx][midx])
                    x_labels[pridx].append("S=" + str(st) + " P=" + str(p) + " M=" + str(m))
                    pass
                pass
            pass
        pass

    fig_root = outp + "boxplot_P"
    fig_titles = ["Unimodal Profile", "Bimodal Profile"]
    plt.rc('xtick', labelsize=6)
    plt.rc('ytick', labelsize=6)

    figs = [plt.figure() for _ in range(len(profs))]
    for f in figs:
        f.set_dpi(500)
        f.set_figheight(3)

    plts = [figs[idx].add_subplot(111) for idx in range(len(profs))]

    for idx in range(len(all_data)):
        plts[idx].boxplot(all_data[idx])
        plts[idx].set_xticklabels(x_labels[idx], rotation=90, va='bottom')
        if idx == 0:
            for xl in plts[idx].xaxis.get_majorticklabels()[:2]:
                xl.set_va('top')
                xl.set_y(1)
            for xl in plts[idx].xaxis.get_majorticklabels()[2:]:
                xl.set_y(0.05)
        else:
            for xl in plts[idx].xaxis.get_majorticklabels()[:1]:
                xl.set_va('top')
                xl.set_y(1)
            for xl in plts[idx].xaxis.get_majorticklabels()[1:]:
                xl.set_y(0.05)

        figs[idx].suptitle(fig_titles[idx], fontsize=12)
        plts[idx].set_xlabel("Parameter Setting (S=states, P=population, M=max. mutations)", fontsize=10)
        plts[idx].set_ylabel("Distribution of RMS Error", fontsize=10)
        # figs[idx].tight_layout()
        figs[idx].subplots_adjust(left=.08, bottom=.1, right=.98, top=.91, wspace=0, hspace=0)
        figs[idx].savefig(fig_root + str(profs[idx]) + ".png")
        pass
    pass


main()
