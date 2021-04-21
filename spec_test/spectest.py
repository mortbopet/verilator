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

from plotspectest import *

RESULT_DIR = os.path.join(os.curdir, "spec_results")
RESFILEPATH = ""
VERILATOR_THREAD_VAR = "VERILATOR_THREAD_CMD"
VERILATOR_ROOT= os.environ["VERILATOR_ROOT"]

VERILATOR_BIN = os.path.join(VERILATOR_ROOT, "bin")
VERILATOR_EXEC = os.path.join(VERILATOR_BIN, "verilator_bin_dbg")
OUT_DIR = os.path.join(VERILATOR_BIN, "obj_dir")

def ensureDir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)

def runTest(vl_file, cpp_file, nthreads, speculate, additionalArgs):
    # Compilation
    pargs = [VERILATOR_EXEC,
        "--Mdir", OUT_DIR,
        "--no-skip-identical", "-o", "out",
        "--cc", f"--threads", f"{nthreads}"]
    if(additionalArgs != None):
        pargs += additionalArgs.split(" ")

    pargs += [vl_file, "--exe", "--build", cpp_file]
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
    parser.add_argument('--additional', type=str, default=None, required=False,
                        help='Additonal arguments to pass to the verilator command line call')
    args = parser.parse_args()

    ensureDir(RESULT_DIR)

    now = datetime.now() # current date and time
    timestamp = now.strftime("%Y%d%m%H%M%S")

    threadsToTest = [int(t) for t in args.threads.split(',')]
    nRuns = len(threadsToTest)*args.N*2

    # Run tests
    run = 0
    RESFILEPATH = os.path.join(RESULT_DIR, f"{args.prefix}_{timestamp}.txt")
    with progressbar.ProgressBar(max_value=nRuns, redirect_stdout=True) as bar:
        with open(RESFILEPATH, "w") as resfile:
            for n_threads in threadsToTest:
                results = []
                for i in range(0, args.N):
                    for doSpec in [True, False]:
                        bar.update(run)
                        run += 1
                        print(f"Test: T:{n_threads} N:{i+1}/{args.N} S:{doSpec}")
                        results.append(runTest(args.fvl, args.fcpp, n_threads, doSpec, args.additional))
                        resfile.write((f"s:" if doSpec else "") + f"{n_threads}={sum(results)/len(results)}\n")
    
    # Plot
    doPlot(RESFILEPATH)
