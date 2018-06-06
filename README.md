# VortexSim
Direct numerical simulator for 2D point vortices

All source code is located in the src directory.  After compilation the executable 'vsim' will be located in the directory ../bin (relative to the src directory).  You will need to create this directory. 

The doc directory contains a Doxyfile for generating doxygen documentation of the code.

The tools directory contains python scripts for generating input, analyzing output files and (most importantly) running vsim.  This last task is done with the vsimLaunch script (as well as the input generation and output analysis tasks).  You will find the option descriptions on this script far superior to those found in the vsim executable (this was the primary way I interacted with vsim). At the top of vsimLaunch is the VORTEX variable which will need to be modified to the appropriate subdirectory where you have placed the project.

The python directory contains the pvparse module for interacting with vsim input and output files.  This module is used extensively by various scripts.  It may also be used for loading data into a python session.
