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
import statistics
# from pyvis.network import Network


def get_fits(dir_path: str):
    # data = [[] for _ in range(30)]
    RI = []
    fit = []
    means = []
    with open(dir_path,encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            spl = line.split(' ')
            spl = list(filter(None, spl))
            # raise ValueError("ValueError exception thrown")
            RI.append(int(spl[0]))
            fit.append(float(spl[-1]))
            means.append(float(spl[1]))
    return RI, fit, means

def main():
    count = 0
    overall_x = [ x for x in range(101)]
    overall_y = [[] for _ in range(101)]
    overall_z = [[] for _ in range(101)]
    for dirpath, dirnames, files in os.walk('.'):
        for file_name in files:
            if file_name.startswith("ru") and file_name.endswith("dat"):
                direc = os.path.join(dirpath, file_name)
                vals = get_fits(direc)
                x = vals[0]
                y = vals[1]
                z = vals[2]
                for i in range(len(y)):
                    overall_y[i].append(y[i])
                    overall_z[i].append(z[i])

                filey = file_name.split(".")
                plt.title(filey[0])
                plt.ylabel("Fitness")
                plt.xlabel("RI")
                plt.plot(x,y,z)
                plt.grid(b=True, which='major', color='#666666', linestyle='-')
                plt.minorticks_on()
                plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
                
                # plt.savefig("aug8test.png")
                plt.savefig(filey[0]+".png")
                plt.close()
                # raise ValueError("ValueError exception thrown")
                # plt.show()
                # print(y)
                count +=1
    plt.title("summary")
    plt.ylabel("Fitness")
    plt.xlabel("RI")
    meany = []
    meanz = []
    # print(overall_y)
    # raise ValueError("ValueError exception thrown")
    for x in range(len(overall_y)):
        meany.append(statistics.mean(overall_y[x]))
        meanz.append(statistics.mean(overall_z[x]))
    print(meanz)
    # raise ValueError("ValueError exception thrown")
    plt.plot(overall_x,meany, meanz)
                
    plt.savefig("overall.png")
    plt.close()
    print(count)




if '__main__' == __name__:
    main()