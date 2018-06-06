/*******************************************************************************
 * \class VortexSim
 *
 * \brief Class for running point vortex fluid simulations.  Glues together
 * many smaller classes.
 *
 * Author: Jonathan Maack
 *
 ******************************************************************************/

#ifndef VORTEX_SIM
#define VORTEX_SIM

#include "Constants.h"
#include "FlowComputer.h"
#include "Integrator.h"
#include "Logger.h"
#include "OutputHQ.h"
#include "Storage.h"
#include "Timer.h"

class ofstream;

class VortexSim
{
public:
    // -------- Functions -------- //
    /** Constructor **/
    VortexSim();
    /** Destructor **/
    virtual ~VortexSim();

    /**
     * Sets up the simulation run. Pass command line parameters directly.
     * @param argc Number of arguments
     * @param argv Array of arguments
     **/
    bool setup(int argc, char **argv);

    /**
     * Run simulation
     **/
    int run();

private:
    // -------- Attributes -------- //

    /** Logging object to handle all screen/file output **/
    Logger *log;
    /** Storage object to keep track of all the point vortex info **/
    Storage *storage;
    /** Storage object for holding the initial positions **/
    Storage *initial;
    /** Object for computing flows for the integrator **/
    FlowComputer *flowcomp;
    /** Integrator which computes time advances **/
    Integrator *integrator;
    /** Object responsible for outputting everything **/
    OutputHQ *outhq;
    /** Object for recording timing information about the run **/
    Timer *timer;
    /** File for saving position information **/
    std::ofstream *movie;

    unsigned batchSize;

    struct CmdLn
    {
        bool computes;
        bool positions;
        bool warmstart;
        int mode; // -1 = unspecified, 0 = time, 1 = steps
        unsigned number;
        unsigned dimension;
        unsigned numSteps;
        unsigned outputRate;
        unsigned positionRate;
        unsigned seed;
        unsigned precision;
        double timeFinal;
        double stepSize;
        double outputTime;
        double radius;
        double rotation;
        Constants::Domain domain;
        Constants::Integrate integrate;
        std::string initialFile;
        std::string logFile;
        std::string outputFile;
        std::string movieFile;
    } cmdln;

    // -------- Functions -------- //

    /**
     * Routine to parse the command line arguments
     * @param argc Number of arguments
     * @param argv Array of arguments
     * @return True on success
     **/
    bool parseCommandLine(int argc, char **argv);

    bool checkArgs(int argc, char **argv, int numArgs, unsigned current);

    /**
     * Run the simulation.
     * @param steps Number of timesteps to take
     * @return true if clean run
     **/
    bool runNumSteps(unsigned steps);

    /**
     * Summarize setup info into log
     **/
    void summarizeSetup();

    /**
     * Summarize run info to the log and screen
     **/
    void summarizeRun();

    /**
     * Setup movie making output file
     **/
    bool initializeMovieFile();

    /**
     * Output positions to movie file
     **/
    void savePositions();

    /**
     * synch up input simulation time with timestep numbers
     **/
    void synchTimeAndStep();

    void reportProgress(unsigned width, double prog) const;
};

#endif
