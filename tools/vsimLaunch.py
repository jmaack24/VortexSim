#!/usr/bin/python

# Script for doing multiple vsim runs--builds initial data too

import argparse
import numpy as np
import os
import subprocess
import time

import pvparse as pvp
import stat_equil as se

# Set locations for various things
home = os.path.expanduser("~")
VORTEX = "/Vortex/"
SCRIPTS = home + VORTEX + "tools/"
RESULTS = home + VORTEX + "output/"
LOG = home + VORTEX + "log/"
DATA = home + VORTEX + "data/"
BIN = home + VORTEX + "bin/"
MOVIE = home + VORTEX + "movies/"

# Generally useful values
COV = 0.40
NUMPVS = 1000
TIME = 500.0
STEP = 1e-2
dim = 2

##### Functions #####

def buildInput(domain, nruns, npvs, seed, distarr, directory, filename,
               start=0):
    status = 0

    if not os.path.exists(directory):
        os.makedirs(directory)
    fname = directory + filename

    # if distarr[0] == "energy":
    #     init_script = "generateICs.py"
    #     distarr = distarr[1:]
    # else:
    #     init_script = "genInit.py"
    if domain[0] == se.Domain.TWO_LAYER:
        genscript = "genTwoLayerICs.py"
        domain = ["--defrad", domain[1]]
        circ = float(2.0/npvs)
    else:
        genscript = "generateICs.py"
        domain = ["--domain"] + domain
        circ = float(1.0/npvs)

    pargs = ["python", SCRIPTS + genscript,
             "-n", str(npvs),
             "-v", str(circ),
             "-r", str(nruns),
             "--dist"]
    pargs += distarr
    pargs += domain
    pargs += ["--seed", str(seed),
              "--plot", "--start", str(start),
              fname]
    retcode = subprocess.call(pargs)
    if retcode != 0:
        print pargs
        print  genscript, "failed -- aborting"
        status = 1
    return status

def runVsim(domain, nruns, nprocs, infile, logfile, outfile, orate=1.0,
            tfinal=10.0, step=1e-2, moviefile=None, mrate=-1.0, start=0):
    complete = start
    count = complete
    procs = []
    error = False
    found = False

    if not os.path.exists(BIN + "vsim"):
        print "Executable vsim not found in", BIN
        return -1

    directory = os.path.dirname(outfile)
    if not os.path.exists(directory):
        os.makedirs(directory)

    if moviefile != None:
        directory = os.path.dirname(moviefile)
        if not os.path.exists(directory):
            os.makedirs(directory)

    while complete < nruns:
        # Check if we can launch a process
        N = len(procs)
        if (N < nprocs) and (complete + N < nruns):
            # Launch process
            print "Launching vsim run", count
            inf = infile % count
            outf = outfile % count
            logf = logfile % count
            pargs = [BIN + "vsim", "-d"]
            pargs += domain
            pargs += ["-init", inf,
                      "-t", str(tfinal),
                      "-step", str(step),
                      "-log", logf,
                      # "-c", str(1.0/step),
                      "-o", str(orate/step), outf]
            if (moviefile != None):
                movf = moviefile % count
                pargs += ["-m", str(mrate/step), movf]

            proc = subprocess.Popen(pargs)
            procs.append(proc)
            count += 1
            # Otherwise check if there is a process to clean up
        else:
            for k in range(N):
                retcode = procs[k].poll()
                if retcode != None:
                    found = True
                    # Process finished -- check for error
                    if retcode != 0:
                        print "vsim failed -- aborting"
                        error = True
                    else:
                        # No error -- remove the process the list and continue
                        procs.pop(k)
                        # print "len(procs) = ", len(procs)
                        break

        if (error == True):
            # Error out
            for k in range(len(procs)):
                procs[k].terminate()
            break
        elif (found == True):
            complete += 1
            found = False
        else:
            # No processes exited so we sleep
            time.sleep(1.0)

    return int(error)

def processOutput(dom, nruns, tfinal, step, root, orate=1.0):
    nentries = int(tfinal/step * step/orate) + 1
    root = RESULTS + root + "/" + root
    pargs = ["python", SCRIPTS + "buildAvgOutputPlot.py",
             "-n", str(nentries),
             "-r", str(nruns),
             "-d", dom,
             root]

    retcode = subprocess.call(pargs)
    return retcode

def createMovie(nfiles, pos_file, start=0):
    retcode = 0
    for k in range(start, nfiles):
        mf = pos_file % k + ".pos"
        pargs = ["python", SCRIPTS + "buildPositionPlots.py",
                 "--ylim", str(4.0),
                 "--xlim", str(4.0),
                 "--movie", mf]
        retcode += subprocess.call(pargs)
    return retcode

#def runOptimalClosure(tfinal, step, nfiles, npvs, dom, ftemplate, obsnum):
#    nentries = int(tfinal/step * step/1.0) + 1
#    pargs = ["python", SCRIPTS + "buildClosurePlots.py",
#             "-n", str(nentries),
#             "-r", str(nfiles),
#             "-d"]
#    pargs += dom
#    pargs += ["--npvs", str(npvs),
#              "--obs", str(obsnum),
#              # "-m", # "--no-positions",
#             ftemplate]
#    retcode = subprocess.call(pargs)
#    return retcode

def renameFiles(old_nruns, new_nruns, root_name):
    pargs = ["python", SCRIPTS + "renumberResults.py",
             "-d", DATA, RESULTS, "-o", str(old_nruns),
             "-n", str(new_nruns), root_name]
    retcode = subprocess.call(pargs)
    return retcode

##### Script #####

# Create command line argument parser
parser = argparse.ArgumentParser()

# Define accepted arguments
parser.add_argument("root", help="root name for files produced by these runs")
parser.add_argument("-d", "--domain", default=None, nargs="+",
                    help="domain and interaction type for PV simulation")
parser.add_argument("--dist", help="distribution for initial condition",
                    default=None, nargs="+")
parser.add_argument("-r","--nruns", help="total number of runs",
                    type=int, default=10)
parser.add_argument("-j","--nprocs", type=int, default=1,
                    help="number of vsim instances running at one time")
parser.add_argument("-s", "--seed", type=int, default=None,
                    help="seed to initialize random number generator")
parser.add_argument("-v", "--npvs", type=int, default=None,
                    help="number of point vortices in simulations")
parser.add_argument("-m", "--movie", default=-1.0, type=float,
                    help="generate movie of run with given time between frames")
parser.add_argument("-p", "--positions", default=1.0, type=float,
                    help="save positions data from vsim run")
parser.add_argument("-o", "--outrate", default=1.0, type=float,
                    help="time interval between observable outputs for DNS")
parser.add_argument("-a", "--append", action="store_true",
                    help="perform additional runs to go with the given sets")
parser.add_argument("--input", action="store_true",
                    help="run the input building part of the script")
parser.add_argument("--no-input", action="store_true",
                    help="omit running the input building part of the script")
parser.add_argument("--vsim", action="store_true",
                    help="run the vsim part of the script")
parser.add_argument("--output", action="store_true",
                    help="run the output processing part of the script")
#parser.add_argument("--opt-close", action="store_true",
#                    help="run the optimal closure processing")
#parser.add_argument("--no-opt-close", action="store_true",
#                    help="omit running the optimal closure part of the script")
parser.add_argument("--cov", help="specify value for covariance",
                    type=float,default=None)
parser.add_argument("--time", help="specify end time of simulation",
                    type=float, default=None)
parser.add_argument("--step", help="time step size for DNS",
                    type=float, default=None)
#parser.add_argument("--obs", type=int, default=0,
#                    help="observable group number")

# Parse arguments
args = parser.parse_args()

nruns = args.nruns
nprocs = args.nprocs
if args.seed == None:
    seed = NUMPVS
else:
    seed = args.seed

if args.domain == None:
    domain = ["plane"]
else:
    domain = args.domain

if args.cov != None:
    COV = args.cov
if args.time != None:
    TIME = args.time
if args.step != None:
    STEP = args.step
if args.npvs != None:
    NUMPVS = args.npvs

if (args.dist == None):
    darr = ["norm", "0.0", "0.0", "1.0", str(COV), "1.0"]
else:
    # Off load argument processing on genInit script
    darr = args.dist

# Used as file name base throughout this script
root = args.root + "_%0" + str(len(str(nruns))) + "d"

# Figure out which parts of the script to run
runall = not(args.input or args.vsim or args.output)# or args.opt_close)
run_input = (args.input or runall) and not args.no_input
run_vsim = args.vsim or runall
run_output = args.output or runall
#run_opt_close = (args.opt_close or runall) and not args.no_opt_close

# Check if we are appending to an existing set of files
start_num = 0
if (args.append):
    print "Renumbering old files..."
    import glob
    old_nruns = len(glob.glob(RESULTS + args.root + "/*.out"))
    start_num = old_nruns
    # Rename everything with the rename script
    status = renameFiles(old_nruns, nruns, args.root)
    if (status != 0):
        print "Renaming old run files -- returned", status
        quit()
    else:
        print "Done!"

# Build the initial data files
if (run_input):
    dirstr = DATA + args.root + "/"
    print "Building initial data sets..."
    status = buildInput(domain, nruns, NUMPVS, seed, darr, dirstr, args.root,
                        start=start_num)
    if (status != 0):
        print "Error building input -- returned", status
        quit()
    else:
        print "Done!!!"

# Do the vsim runs
if (run_vsim):
    print "Performing runs...."
    infile = DATA + args.root + "/" + root + ".pv"
    output = RESULTS + args.root + "/" + root + ".out"
    logdir = LOG + args.root + "/"
    if not os.path.exists(logdir):
        os.makedirs(logdir)
    log = logdir + root
    if (args.movie > 0.0):
        args.positions = args.movie
    if (args.positions > 0.0):
        mfile = MOVIE + args.root + "/" + root + ".pos"
    else:
        mfile = None
    status = runVsim(domain, nruns, nprocs, infile, log, output,
                     orate=args.outrate, tfinal=TIME, step=STEP,
                     moviefile=mfile, mrate=args.positions, start=start_num)
    if (status != 0):
        print "Error conducting vsim runs -- returned", status
        quit()
    else:
        print "Runs complete!!!"

# Process the output
if (run_output):
    print "Processing output..."
    status = processOutput(domain[0], nruns, TIME, STEP, args.root,
                           args.outrate)
    if (status != 0):
        print "Error processing output -- returned", status
        quit()
    else:
        print "Done!!!"

if (args.movie > 0.0 and run_output):
    print "Creating movie..."
    mfile = MOVIE + args.root + "/" + root
    status = createMovie(nruns, mfile, start=start_num)
    if (status != 0):
        print "Error creating movie -- returned", status
        quit()
    else:
        print "Done!!"

## Run optimal closure stuff
#if (run_opt_close):
#    print "Running optimal closure..."
#    status = runOptimalClosure(TIME, STEP, nruns, NUMPVS, domain, args.root,
#                               args.obs)
#    if (status != 0):
#        print "Error running optimal closure -- returned", status
#        quit()
#    else:
#        print "Done!!!"
