import statistics
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import csv
import math
import re
import string
from operator import itemgetter

import matplotlib
matplotlib.use('Agg')
# from pyvis.network import Network

densities = ["ps1", "ps2", "ps3", "ps4", "ps5", "ps6", "ps7", "ps8", "ps9", "ps10", "ps11", "ps12", "ps13", "ps14",
             "ps15", "ps16", "ps17", "ps18", "ps19", "ps20", "ps21", "ps22", "ps23", "ps24", "ps25", "ps26", "ps27", "ps28", "ps29"]


def get_fits(dir_path: str):
    # data = [[] for _ in range(30)]
    RI = []
    fit = []
    means = []
    with open(dir_path, encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            spl = line.split(' ')
            spl = list(filter(None, spl))
            RI.append(int(spl[0]))
            # fit.append(float(spl[-1]))  # Epidemic Duration
            fit.append(float(spl[-3]))  # Profile matching
            means.append(float(spl[1]))
    return RI, fit, means


def main():
    count = 0
    pathy = ""
    dens = ""
    overall_x = [x for x in range(101)]
    overall_y = [[] for _ in range(101)]
    overall_z = [[] for _ in range(101)]
    for dirpath, dirnames, files in os.walk('.'):
        for file_name in files:
            if file_name.startswith("ru") and file_name.endswith("dat"):
                if count == 0:
                    pathy = dirpath.split('/')
                    dens = pathy[-2]
                filey = file_name.split(".")
                pathy = dirpath.split('/')
                print(pathy[-2] + " " + dens)
                if pathy[-2] != dens:
                    plt.title(dens + " - summary")
                    plt.ylabel("Fitness")
                    plt.xlabel("RI")
                    meany = []
                    meanz = []
                    for x in range(len(overall_y)):
                        meany.append(statistics.mean(overall_y[x]))
                        meanz.append(statistics.mean(overall_z[x]))
                    plt.plot(overall_x, meany, meanz)
                    plt.grid(b=True, which='major',
                             color='#666666', linestyle='-')
                    plt.minorticks_on()
                    plt.grid(b=True, which='minor', color='#999999',
                             linestyle='-', alpha=0.2)
                    plt.savefig(dens + "overall.png")
                    plt.close()
                    dens = pathy[-2]
                    overall_y = [[] for _ in range(101)]
                    overall_z = [[] for _ in range(101)]
                # raise ValueError("ValueError exception thrown")
                direc = os.path.join(dirpath, file_name)
                vals = get_fits(direc)
                x = vals[0]
                y = vals[1]
                z = vals[2]
                for i in range(len(y)):
                    overall_y[i].append(y[i])
                    overall_z[i].append(z[i])

                title_name = pathy[-2] + " - " + filey[0]
                plt.title(title_name)
                plt.ylabel("Fitness")
                plt.xlabel("RI")
                plt.plot(x, y, z)
                plt.grid(b=True, which='major', color='#666666', linestyle='-')
                plt.minorticks_on()
                plt.grid(b=True, which='minor', color='#999999',
                         linestyle='-', alpha=0.2)
                plt.savefig(pathy[-2]+filey[0]+".png")
                plt.close()
                count += 1
            # plt.title(pathy[1] + " - summary")
            # plt.ylabel("Fitness")
            # plt.xlabel("RI")
            # meany = []
            # meanz = []
            # # print(overall_y)
            # # raise ValueError("ValueError exception thrown")
            # for x in range(len(overall_y)):
            #     meany.append(statistics.mean(overall_y[x]))
            #     meanz.append(statistics.mean(overall_z[x]))
            # print(meanz)
            # # raise ValueError("ValueError exception thrown")
            # plt.plot(overall_x,meany, meanz)

            # plt.savefig(pathy[-2] +"overall.png")
            # plt.close()
    print("##################################################################################")
    plt.title(dens + " - summary")
    plt.ylabel("Fitness")
    plt.xlabel("RI")
    meany = []
    meanz = []
    # print(overall_y)
    # raise ValueError("ValueError exception thrown")
    for x in range(len(overall_y)):
        meany.append(statistics.mean(overall_y[x]))
        meanz.append(statistics.mean(overall_z[x]))
        # raise ValueError("ValueError exception thrown")
    plt.plot(overall_x, meany, meanz)
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.savefig(dens + "overall.png")
    plt.close()


if '__main__' == __name__:
    main()
