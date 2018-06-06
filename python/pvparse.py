
import math
import numpy as np
import sys

NUM_ENTRIES = "Number of Entries"
TIME = "Simulation Time"
CENTER = "Center of Vorticity"
MOMENT = "Moment of Inertia"
HAMIL = "Hamiltonian"
COV = "Covariance Matrix"
ACORR = "Autocorrelation"
SECOND = "Second Moment"
THIRD = "Third Moment"
FOURTH = "Fourth Moment"
MSQ = "Eq Error"
# TwoLayerMoments.cpp denotes second moments from layer 1 (marked in the
# C++ code and input generation code with layer == -1) as "Second Moment 2" and
# the second moments from layer 2 (marked in the C++ and input generation code
# with layer == 1) as "Second Moment 1".  I have only myself to blame...
# Anyway, fix it here so that the naming and associations are the same from
# outside the black box
SECOND1 = "Second Moment 2"
SECOND2 = "Second Moment 1"

CENTER1 = "Center 1"
CENTER2 = "Center 2"

# Maps that convert from input generation and C++ code layers to person layers
LAYER_MAP = {-1:1, 1:2}
INV_LAYER_MAP = {1:-1, 2:1}

def reportProgress(processed, total, width):
    # Check progess and report
    percent = float(processed)/total
    nprint = int(math.floor(percent*width))
    nblank = width - nprint
    pstr = "\r["
    pstr += ":" * nprint
    pstr += " " * nblank
    percent *= 100
    pstr += "]{:>6.2f}% complete".format(percent)
    if percent >= 100.0:
        pstr += "\n"
    sys.stdout.write(pstr)
    sys.stdout.flush()
    return percent


def findNumEntries(data):
    searching = True
    while (searching):
        line = data.readline()
        if(line):
            first, sep, last = line.partition("=")
            if (first == "Number of Entries"):
                searching = False
                numEntries = int(last)
        else:
            print "Reached end of file without finding number of entries"
            numEntries = -1
            searching = False
    return numEntries

def getGeneratingDistribution(datafile):
    data = open(datafile, mode='r')

    searching = True
    retval = None
    while (searching):
        line = data.readline()
        if (line):
            first, sep, last = line.partition(":")
            if (first == "# dist"):
                retval = last.strip().split(" ")
                searching = False
        else:
            print "Reached end of file without finding distribution"
            searching = False
    return retval

def readPositionEntry(data):

    status = 0
    processing = True
    num = dim = time = vorticity = position = 0

    while (processing and status == 0):
        line = data.readline()

        if (line == ""):
            print "File ended before entry completed"
            status = -1

        # print "Parsing the line '", line, "'"

        first, sep, last = line.partition("=")
        # Check for start of a block
        if (first == "Simulation Time"):
            # Get time
            time = float(last)

            # Get number of vortices
            line = data.readline()
            first, sep, last = line.partition("=")
            if (first != "Number"):
                print "Unrecognized token '", first, "' -- expected","'Number'"
                status = 1

            num = int(last)

            # Get dimension--that is the size of the position vector
            line = data.readline()
            first, sep, last = line.partition("=")
            if (first != "Dimension"):
                print "Unrecognized token '",first,"' -- expected","'Dimension'"
                status = 2

            dim = int(last)

            layer = np.zeros((num,), dtype=np.int)
            vorticity = np.zeros((num,))
            position = np.zeros((num, dim))

            # Now get vorticity and positions for each vortex
            for i in range(num):
                line = data.readline()
                first, sep, last = line.partition(" ")
                # Layer -- skips to vorticity for backward compatibility if
                # string does not match
                lstr, sep, lay = first.partition("=")
                if (lstr == "Layer"):
                    layer[i] = int(lay)
                    first, sep, last = last.partition(" ")
                # Vorticity
                vstr, sep, vort = first.partition("=")
                if (vstr != "Vort"):
                    print "Unrecognized token '",vstr,"' -- expected ","'Vort'"
                    status = 3

                v = float(vort)
                vorticity[i] = v
                # Position
                pstr, sep, pvect = last.partition("=")
                if (pstr != "Pos"):
                    print "Unrecognized token '",pstr,"' -- expected ","'Pos'"
                    status = 4

                pvect=pvect.lstrip("(").rstrip(")\n")
                for j in range(dim):
                    first, sep, remain = pvect.partition(",")
                    # print first
                    # print sep
                    # print remain
                    coord = float(first)
                    position[i,j] = coord
                    pvect = remain

                # Check for error--if remain is not empty there were more
                # values than expected
                if (remain):
                    print "Mismatch between number of stated dimensions and number found"
                    status = 5

            # Check for error -- read next line and strip all whitespace
            line = data.readline().lstrip(" ")
            if (line != '\n'):
                print "Found more vortices than expected"
                status = 6
            else:
                processing = False

    return status, num, dim, time, layer, vorticity, position

def readPositionFile(filename, number, dimension, silent=False):
    data = open(filename, mode='r')

    status = 0
    times = positions = 0
    WIDTH = 80

    numEntries = findNumEntries(data)
    if (numEntries < 1):
        print "Found number of entries is", numEntries
        status = -1
    else:
        times = np.zeros((numEntries,))
        layer = np.zeros((number,))
        vorticity = np.zeros((number,))
        positions = np.zeros((numEntries, number, dimension))
        k = 0
        while (k < numEntries and status == 0):
            if (not silent):
                reportProgress(k, numEntries, WIDTH)

            (status,num,dim,times[k],layer,vorticity,positions[k]
            )=readPositionEntry(data)
            if (status != 0):
                print "**ERROR::readPositionEntry returned", status
            elif (num != number):
                print "Expected", number, "vortices and found", num
                status = 1
            elif (dim != dimension):
                print "Expected dimension", dimension, "and found", dim
                status = 2
            else:
                k += 1

        if (k != numEntries):
            print "Error processing file"
            status = 3
        elif (not silent):
            reportProgress(numEntries, numEntries, WIDTH)

    data.close()

    return status, times, layer, vorticity, positions

def readOutputFile(filename, dim):
    sts = 0
    data = open(filename, mode='r')

    # Find number of entries in the file
    numEntries = findNumEntries(data);
    if (numEntries <= 0):
        sts = -1
        numEntries = 1

    # Build storage spots for each quantity
    time = np.zeros((numEntries,))
    center = np.zeros((numEntries, dim))
    moment = np.zeros((numEntries,))
    hamil = np.zeros((numEntries,))
    cov = np.zeros((numEntries, dim, dim))
    acorr = np.zeros((numEntries, dim))
    msq = np.zeros((numEntries, 4))

    # THESE ASSUME DIMENSION 2!!!
    center1 = np.zeros((numEntries, 2))
    center2 = np.zeros((numEntries, 2))
    second = np.zeros((numEntries, 3))
    second1 = np.zeros((numEntries, 3))
    second2 = np.zeros((numEntries, 3))
    third = np.zeros((numEntries, 4))
    fourth = np.zeros((numEntries, 5))

    # Build return object
    retobj = {TIME:time}

    # Read the file
    if (sts != 0):
        print "Error finding number of entries. Unable to proceed"
    else:
        count = 0
        line = data.readline()
        while(line):
            first, sep, last = line.partition(":")
            if (first == TIME):
                time[count] = float(last)
                line = data.readline()
                while(line != "\n"):
                    line = line.strip(" ").rstrip(" \n")
                    first, sep, last = line.partition(":")

                    if (first == CENTER):

                        if not CENTER in retobj:
                            retobj[CENTER] = center

                        for l in range(dim):
                            last = last.strip(" ()").rstrip(" ()\n")
                            entry, sep, rem = last.partition(",")
                            center[count, l] = float(entry)
                            last = rem

                        if (last):
                            print "Extra values in", CENTER, "--", last
                            sts = 1
                    elif (first == MOMENT):
                        if not MOMENT in retobj:
                            retobj[MOMENT] = moment

                        last = last.strip(" ").rstrip(" \n")
                        moment[count] = float(last)

                    elif (first == HAMIL):
                        if not HAMIL in retobj:
                            retobj[HAMIL] = hamil

                        last = last.strip(" ").rstrip(" \n")
                        hamil[count] = float(last)

                    elif (first == COV):
                        if not COV in retobj:
                            retobj[COV] = cov

                        for i in range(dim):
                            for j in range(i,dim):
                                last = last.strip(" []")
                                entry, sep, rem = last.partition(",")
                                entry = entry.strip("[]").rstrip("[]")
                                cov[count, i, j] = float(entry)
                                last = rem

                        if (last):
                            print "Extra values in", COV, "--", last
                            sts = 1

                    elif (first == ACORR):
                        if not ACORR in retobj:
                            retobj[ACORR] = acorr

                        for i in range(dim):
                            last = last.strip(" ()").rstrip(" ()\n")
                            entry, sep, rem = last.partition(",")
                            acorr[count, i] = float(entry)
                            last = rem

                    elif (first == SECOND):
                        if not SECOND in retobj:
                            retobj[SECOND] = second

                        if (dim != 2):
                            print SECOND, "parsed only for dimension 2"
                        else:
                            last = last.strip(" [").rstrip("]")
                            xx, sep, rem = last.partition(",")
                            xy, sep, yy = rem.partition(",")
                            second[count, 0] = float(xx)
                            second[count, 1] = float(xy)
                            second[count, 2] = float(yy)

                    elif (first == THIRD):
                        if not THIRD in retobj:
                            retobj[THIRD] = third

                        if (dim != 2):
                            print THIRD, "parsed only for dimension 2"
                        else:
                            last = last.strip(" [").rstrip("]")
                            xxx, sep, rem = last.partition(",")
                            xxy, sep, rem = rem.partition(",")
                            xyy, sep, yyy = rem.partition(",")
                            third[count, 0] = float(xxx)
                            third[count, 1] = float(xxy)
                            third[count, 2] = float(xyy)
                            third[count, 3] = float(yyy)

                    elif(first == FOURTH):
                        if not FOURTH in retobj:
                            retobj[FOURTH] = fourth

                        if (dim != 2):
                            print FOURTH, "parsed only for dimension 2"
                        else:
                            last = last.strip(" [").rstrip("]")
                            xxxx, sep, rem = last.partition(",")
                            xxxy, sep, rem = rem.partition(",")
                            xxyy, sep, rem = rem.partition(",")
                            xyyy, sep, yyyy = rem.partition(",")
                            fourth[count, 0] = str(xxxx)
                            fourth[count, 1] = str(xxxy)
                            fourth[count, 2] = str(xxyy)
                            fourth[count, 3] = str(xyyy)
                            fourth[count, 4] = str(yyyy)

                    elif(first == MSQ):
                        if not MSQ in retobj:
                            retobj[MSQ] = msq

                        string = last.strip(" ").rstrip(" \n")
                        mse, sep, rem = string.partition(",")
                        mxe, sep, rem = rem.partition(",")
                        mag, sep, emg = rem.partition(",")
                        msq[count,0] = float(mse)
                        msq[count,1] = float(mxe)
                        msq[count,2] = float(mag)
                        msq[count,3] = float(emg)

                    elif (first == SECOND1):
                        if not SECOND1 in retobj:
                            retobj[SECOND1] = second1

                        if (dim != 2):
                            print SECOND1, "parsed only for dimension 2"
                        else:
                            last = last.strip(" [").rstrip("]")
                            xx, sep, rem = last.partition(",")
                            xy, sep, yy = rem.partition(",")
                            second1[count, 0] = float(xx)
                            second1[count, 1] = float(xy)
                            second1[count, 2] = float(yy)
                    elif (first == SECOND2):
                        if not SECOND2 in retobj:
                            retobj[SECOND2] = second2

                        if (dim != 2):
                            print SECOND2, "parsed only for dimension 2"
                        else:
                            last = last.strip(" [").rstrip("]")
                            xx, sep, rem = last.partition(",")
                            xy, sep, yy = rem.partition(",")
                            second2[count, 0] = float(xx)
                            second2[count, 1] = float(xy)
                            second2[count, 2] = float(yy)
                    elif (first == CENTER1):
                        if not CENTER1 in retobj:
                            retobj[CENTER1] = center1
                        if (dim != 2):
                            print CENTER1, "parsed only for dimension 2"
                        else:
                            last = last.strip(" (").rstrip(")")
                            x, sep, y = last.partition(",")
                            center1[count, 0] = float(x)
                            center1[count, 1] = float(y)
                    elif (first == CENTER2):
                        if not CENTER2 in retobj:
                            retobj[CENTER2] = center2
                        if (dim != 2):
                            print CENTER2, "parsed only for dimension 2"
                        else:
                            last = last.strip(" (").rstrip(")")
                            x, sep, y = last.partition(",")
                            center2[count, 0] = float(x)
                            center2[count, 1] = float(y)
                    else:
                        print "Unrecognized quantity -- ", first
                        sts = 1

                    line = data.readline()
                count += 1
            line = data.readline()

    return sts,retobj

def readOutputFileBatch(nruns, nentries, dim, ftemplate, compute_mean=False):

    #nruns=60

    status = 0
    time = np.zeros((nruns, nentries))
    center = np.zeros((nruns, nentries, dim))
    moment = np.zeros((nruns, nentries))
    hamil = np.zeros((nruns, nentries))
    cov = np.zeros((nruns, nentries, dim, dim))
    acorr = np.zeros((nruns, nentries, dim))
    msq = np.zeros((nruns, nentries, 4))

    # ASSUMES TWO DIMENSIONS!!!!
    center1 = np.zeros((nruns, nentries, 2))
    center2 = np.zeros((nruns, nentries, 2))
    second = np.zeros((nruns, nentries, 3))
    second1 = np.zeros((nruns, nentries, 3))
    second2 = np.zeros((nruns, nentries, 3))
    third = np.zeros((nruns, nentries, 4))
    fourth = np.zeros((nruns, nentries, 5))

    retdict = {}

    for k in range(nruns):
        outfile = ftemplate % k
        sts, rdict = readOutputFile(outfile, dim)
        if sts != 0:
            print "readOutputFile returned", sts, "for file", outfile
        for key in iter(rdict):
            if key not in retdict:
                if key == TIME:
                    retdict[TIME] = time
                elif key == CENTER:
                    retdict[CENTER] = center
                elif key == MOMENT:
                    retdict[MOMENT] = moment
                elif key == HAMIL:
                    retdict[HAMIL] = hamil
                elif key == COV:
                    retdict[COV] = cov
                elif key == ACORR:
                    retdict[ACORR] = acorr
                elif key == SECOND:
                    retdict[SECOND] = second
                elif key == THIRD:
                    retdict[THIRD] = third
                elif key == FOURTH:
                    retdict[FOURTH] = fourth
                elif key == MSQ:
                    retdict[MSQ] = msq
                elif key == SECOND1:
                    retdict[SECOND1] = second1
                elif key == SECOND2:
                    retdict[SECOND2] = second2
                elif key == CENTER1:
                    retdict[CENTER1] = center1
                elif key == CENTER2:
                    retdict[CENTER2] = center2
                else:
                    print "readOutputFileBatch: unrecognized quantity --", key
                    status = 1
                    break

            retdict[key][k] = rdict[key]

    if (compute_mean):
        for key in iter(retdict):
            retdict[key] = np.mean(retdict[key], 0)

    return status, retdict

def writeInitFile(filename, append, header, vorticity, position, layer=None):
    if (append == True):
        out = open(filename, mode="a")
    else:
        out = open(filename + ".pv", mode='w')

    out.write(header)

    (num, dim) = position.shape

    for k in range(num):
        output = "{:.15e}".format(vorticity[k])
        for l in range(dim):
            output += "    {:.15e}".format(position[k,l])
        if isinstance(layer, np.ndarray):
            output += "    {:d}".format(layer[k])
        output += "\n"
        out.write(output)

    out.flush()
    out.close()

def readInitFile(filename, num, dim):
    data = open(filename, mode="r")

    vort = np.zeros((num,))
    pos = np.zeros((num, dim))

    count = 0
    status = 0
    line = data.readline()
    while(line):
        # Vorticity
        line = line.strip(" \t").rstrip(" \t\n")
        first, sep, last = line.partition(" ")
        vort[count] = float(first)
        last = last.strip(" \t")
        # Position
        for k in range(dim):
            first, sep, last = last.partition(" ")
            pos[count, k] = float(first)
            last = last.strip(" \t")

        count += 1
        line = data.readline()

    if (count != num):
        print "Error: Expected", num, "vortices but found", count
        status = -1

    return status, vort, pos

def parseBceFile(filename, nentries):
    data = open(filename, mode='r')
    time = np.zeros(nentries)
    ebc = np.zeros(nentries)

    count = 0
    status = 0
    line = data.readline()
    while(line):
        # Initial line processing
        line = line.strip(" \t").rstrip(" \t\n")
        first, sep, last = line.partition(" ")

        # Time
        grbg, sep, ti = first.partition("=")
        time[count] = float(ti)

        # Baroclinic Energy
        grbg, sep, bce = last.partition("=")
        ebc[count] = float(bce)

        # Loop stuff
        count += 1
        line = data.readline()

    if (count != nentries):
        print "Error: Expected", nentries, "entries but found", count
        status = -1

    return status, time, ebc

def parseApeFile(filename, nentries):
    data = open(filename, mode='r')
    time = np.zeros(nentries)
    ape = np.zeros(nentries)

    count = 0
    status = 0
    line = data.readline()
    while(line):
        # Initial line processing
        line = line.strip(" \t").rstrip(" \t\n")
        first, sep, last = line.partition(" ")

        # Time
        grbg, sep, ti = first.partition("=")
        time[count] = float(ti)

        # Baroclinic Energy
        grbg, sep, pe = last.partition("=")
        ape[count] = float(pe)

        # Loop stuff
        count += 1
        line = data.readline()

    if (count != nentries):
        print "Error: Expected", nentries, "entries but found", count
        status = -1

    return status, time, ape
