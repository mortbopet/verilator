import argparse
import os
import matplotlib
import matplotlib.pyplot as plt
import subprocess
import re
import sys
import numpy as np
from datetime import datetime
from collections import defaultdict
import progressbar


def doPlot(filepath):
    runs = 0
    values = []
    tests = defaultdict(lambda: {})
    splitter = re.compile(r"(s:)?(\d+)=([0-9]*[.]?[0-9]+)")

    # Parse
    with open(filepath) as f:
        for line in f.readlines():
            res = re.match(splitter, line)
            if res == None:
                print("Could not parse test output")
                sys.exit(1)
            isSpec = res[1] != None
            if isSpec not in tests[res[2]]:
                tests[res[2]][isSpec] = []
            tests[res[2]][isSpec].append(float(res[3]))

    # Average
    for thread, inner in tests.items():
        for isSpeculative, vals in inner.items():
            runs = len(vals)
            tests[thread][isSpeculative] = sum(vals)/len(vals)


    labels = list(tests.keys())
    x = np.arange(len(labels))  # the label locations
    fig, ax = plt.subplots()

    w = 0.8
    width = w / len(labels)
    i = 0

    nonSpecValues = []
    specValues = []
    for threads in labels:
        specValues.append(tests[threads][True])
        nonSpecValues.append(tests[threads][False])

    ax.bar(x - width/2, nonSpecValues, width, label="Normal")
    ax.bar(x + width/2, specValues, width, label="Speculative")


    ax.set_ylabel('Time [s]')
    ax.set_xlabel('Threads [#]')
    ax.set_title(f'Synthetic Test (avg. over {runs} runs)')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    fig.tight_layout()

    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Verilator benchmark runner')
    parser.add_argument('--f', type=str, required=True,
                        help='Spec test result file')

    args = parser.parse_args()
    doPlot(args.f)