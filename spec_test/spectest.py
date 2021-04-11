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

RESULT_DIR = os.path.join(os.curdir, "spec_results")
VERILATOR_THREAD_VAR = "VERILATOR_THREAD_CMD"
VERILATOR_ROOT= os.environ["VERILATOR_ROOT"]

VERILATOR_BIN = os.path.join(VERILATOR_ROOT, "bin")
VERILATOR_EXEC = os.path.join(VERILATOR_BIN, "verilator_bin_dbg")
OUT_DIR = os.path.join(VERILATOR_BIN, "obj_dir")

def ensureDir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)


def doPlot():
    runs = 0
    values = []
    tests = defaultdict(lambda: {})
    splitter = re.compile(r"(s:)?(\d+)=([0-9]*[.]?[0-9]+)")

    # Parse
    for testname in os.listdir(RESULT_DIR):
        with open(os.path.join(RESULT_DIR, testname)) as f:
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


def runTest(vl_file, cpp_file, nthreads, speculate):
    # Compilation
    pargs = [VERILATOR_EXEC, "--Mdir", OUT_DIR, "--no-skip-identical", "-o", "out", "--cc", f"--threads", f"{nthreads}", vl_file, "--exe", "--build", cpp_file]
    if speculate:
        pargs.append("--speculate")
    process = subprocess.Popen(pargs, stdout=subprocess.PIPE)
    (output, err) = process.communicate()
    exit_code = process.wait()
    if exit_code != 0:
        print("Error: Could not verilate file")
    
    # Test execution
    pargs = [os.path.join(OUT_DIR, "out")]
    process = subprocess.Popen(pargs, stdout=subprocess.PIPE)
    (output, err) = process.communicate()
    exit_code = process.wait()

    return parseOutput(output)


def parseOutput(output):
    """ Parses the cycles/second result determined by the benchmark"""
    output = output.decode().strip()
    return float(output)

if __name__ == "__main__":
    if VERILATOR_ROOT == "":
        print("Environment variable 'VERILATOR_ROOT' not set")
        sys.exit(1)

    parser = argparse.ArgumentParser(description='Verilator benchmark runner')
    parser.add_argument(
        '--N', help='Number of test runs for each config', type=int, required=True)
    parser.add_argument('--threads', type=str, default="1",
                        help='Comma-separated list of # of threads to test')
    parser.add_argument('--prefix', type=str, default="new",
                        help='Prefix for test results (to separate runs)')
    parser.add_argument('--fvl', type=str, required=True,
                        help='verilog test file')
    parser.add_argument('--fcpp', type=str, default="new", required=True,
                        help='cpp test file')
    args = parser.parse_args()

    ensureDir(RESULT_DIR)

    now = datetime.now() # current date and time
    timestamp = now.strftime("%Y%d%m%H%M%S")

    threadsToTest = [int(t) for t in args.threads.split(',')]
    nRuns = len(threadsToTest)*args.N*2

    # Run tests
    i = 1
    with progressbar.ProgressBar(max_value=nRuns, redirect_stdout=True) as bar:
        with open(os.path.join(RESULT_DIR, f"{args.prefix}_{timestamp}.txt"), "w") as resfile:
            for n_threads in threadsToTest:
                results = []
                for i in range(0, args.N):
                    for doSpec in [True, False]:
                        print(f"Test: T:{n_threads} N:{i+1} S:{doSpec}")
                        results.append(runTest(args.fvl, args.fcpp, n_threads, doSpec))
                        resfile.write((f"s:" if doSpec else "") + f"{n_threads}={sum(results)/len(results)}\n")
                        bar.update(i)
                        i += 1
    
    # Plot
    doPlot()
