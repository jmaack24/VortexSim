/*******************************************************************************
 * \class VortexSim
 *
 * \brief Class for running point vortex fluid simulations.  Glues together
 * many smaller classes.
 *
 * Author: Jonathan Maack
 *
 ******************************************************************************/

#include "VortexSim.h"

// Flow computers
#include "FlowComputer.h"
#include "ScreenedFlow.h"
//#include "SphereFlow.h"
#include "TwoLayerFlow.h"
#include "UnboundedFlow.h"

// Time integrators
#include "Integrator.h"
#include "RuKu4.h"
//#include "RuKu45.h"
#include "AdMoulPC.h"

// Utility classes
#include "Logger.h"
#include "OutputHQ.h"
#include "Storage.h"

#include <math.h>
#include <sstream>
#include <string.h>

VortexSim::VortexSim():
    integrator(0),
    flowcomp(0),
    storage(0),
    initial(0),
    log(0),
    outhq(0),
    movie(0)
{
    log = new Logger();
    timer = new Timer(10000);
    //log->writeMsg(Logger::INFO, "VortexSim construction complete");
}

VortexSim::~VortexSim()
{
    delete outhq;
    outhq = 0;

    delete integrator;
    integrator = 0;

    delete flowcomp;
    flowcomp = 0;

    delete movie;
    movie = 0;

    delete storage;
    storage = 0;

    delete initial;
    initial = 0;

    delete timer;
    timer = 0;

    // Make sure log is the last to go--other things may use it
    delete log;
    log = 0;
}

bool VortexSim::setup(int argc, char **argv)
{
    bool good = parseCommandLine(argc, argv);
    if (!good)
    {
        std::cerr << "Error parsing command line. Aborting..." << std::endl;
    }
    else if ((cmdln.number == 0 && cmdln.initialFile.empty()) ||
             cmdln.dimension == 0)
    {
        std::cerr << "Required command line argument not specified--"
                  << std::endl;
        if (cmdln.number == 0 && cmdln.initialFile.empty())
            std::cerr << "must specify number of vortices (-v or -vor) "
                      << "or an initialization file (-f or -init)"
                      << std::endl;
        if (cmdln.dimension == 0)
            std::cerr << "must specify dimension (-d)" << std::endl;
        good = false;
    }
    else
    {
        // Setup logger
        if (!log->setupLogFile(cmdln.logFile))
        {
            std::cerr << "Error setting up log file '" << cmdln.logFile
                      << "'" << std::endl;
            good = false;
        }

        // Setup storage/vortices
        storage = new Storage(log, cmdln.number, cmdln.dimension);
        if (!cmdln.initialFile.empty())
        {
            if (!storage->initializeFromFile(cmdln.warmstart,
                                             cmdln.initialFile,
                                             cmdln.domain))
            {
                std::cerr << "Error initializing from file '"
                          << cmdln.initialFile
                          << "'" << std::endl;
                good = false;
            }
        }
        else if (!storage->initializeRandom(cmdln.seed))
        {
            std::cerr << "Error setting up vortices with random data"
                      << std::endl;
            good = false;
        }

        // Store copy of initial data
        initial = new Storage(*storage);

        // Setup movie file
        if (!initializeMovieFile())
        {
            std::cerr << "Error setting up movie file" << std::endl;
        }

        // Setup flow computer stuff
        if (cmdln.domain == Constants::PLANE)
        {
            flowcomp = new UnboundedFlow(log, timer);
        }
        else if (cmdln.domain == Constants::SCREENED)
        {
            flowcomp = new ScreenedFlow(log, timer, cmdln.radius);
        }
        else if (cmdln.domain == Constants::TWO_LAYER)
        {
            flowcomp = new TwoLayerFlow(log, timer, cmdln.radius);
        }
        // else if (cmdln.domain == Constants::SPHERE)
        // {
        //     std::cerr << "****WARN: Sphere vortex simulation "
        //               << "not fully implemented" << std::endl;
        //     flowcomp = new SphereFlow(log, timer);
        // }
        else
        {
            std::cerr << "Unrecognized flow type -- " << cmdln.domain
                      << std::endl;
            good = false;
        }

        // Setup time integration stuff
        if (cmdln.integrate == Constants::RK4)
        {
            integrator = new RuKu4(log, timer, flowcomp,
                                   storage, cmdln.stepSize);
        }
        else if (cmdln.integrate == Constants::RK45)
        {
            integrator = new RuKu45(log, timer, flowcomp,
                                    storage, cmdln.stepSize);
        }
        else if (cmdln.integrate == Constants::ADMOUL)
        {
            integrator = new AdMoulPC(log, timer, flowcomp,
                                      storage, cmdln.stepSize);
        }
        else
        {
            std::cerr << "Unrecognized integrator -- " << cmdln.integrate
                      << std::endl;
            good = false;
        }

        // Make sure simulation time and number of steps say the same thing
        synchTimeAndStep();

        // Resolve when things are outputted
        if (cmdln.outputRate == 0)
        {
            cmdln.outputRate = cmdln.numSteps;
        }
        else if (cmdln.outputRate > cmdln.numSteps)
        {
            cmdln.outputRate = cmdln.numSteps;
            std::stringstream ss;
            ss << "Output rate exceeds number of steps -- reset "
               << "output rate to " << cmdln.outputRate << std::endl;
            if (log->writeMsg(Logger::WARN, ss.str()))
                std::cerr << "****WARN: " << ss.str();
        }

        // Setup outputting stuff
        outhq = new OutputHQ(log, timer);
        unsigned nentries = cmdln.numSteps / cmdln.outputRate + 1;
        if (!outhq->setup(cmdln.outputFile, nentries, storage, cmdln.domain,
                cmdln.radius))
        {
            std::cerr << "Error setting up file '"
                      << cmdln.outputFile
                      << "' for output" << std::endl;
            good = false;
        }

        if (movie != 0)
        {
            nentries = cmdln.numSteps / cmdln.positionRate + 1;
            (*movie) << "Number of Entries=" << nentries
                     << std::endl << std::endl;
        }

        if (good)
            summarizeSetup();

        timer->stamp(Timer::SETUP);
    }

    return good;
}

int VortexSim::run()
{
    int retval = 0;

    unsigned steps = batchSize;

    double prog;

    std::cout << std::endl << "Running..." << std::endl;

    // Output initial values
    if (cmdln.computes)
    {
        double time = integrator->stepNum() * integrator->stepSize();
        if (!outhq->processOutput(storage, time))
        {
            std::cerr << "Error computing end values" << std::endl;
            retval = 1;
        }
    }

    // Output initial positions
    if (cmdln.positions)
        savePositions();

    // Setup progress bar
    reportProgress(0, 0.0);

    // Run to end time
    while (integrator->stepNum() < cmdln.numSteps)
    {
        // Make sure we do not run further than needed
        if (cmdln.numSteps - integrator->stepNum() < steps)
            steps = cmdln.numSteps - integrator->stepNum();

        // Run so many steps
        if (!runNumSteps(steps))
        {
            std::cerr << "Error running simulation. Run aborted." << std::endl;
            retval = 2;
            break;
        }

        // Do some house keeping
        timer->condense();

        // Report on progress
        prog = double(integrator->stepNum())/cmdln.numSteps;
        reportProgress(0, prog);
    }

    std::cout << std::endl << "Run complete!!" << std::endl << std::endl;

    summarizeRun();

    return retval;
}

bool VortexSim::checkArgs(int argc, char **argv,
                          int numArgs, unsigned current)
{
    bool retval = true;

    if (argc <= current + numArgs)
    {
        retval = false;
    }
    else
    {
        for (unsigned add = 1; add <= numArgs; ++add)
        {
            if (argv[current + add][0] == '-')
            {
                retval = false;
                break;
            }
        }
    }

    return retval;
}

bool VortexSim::parseCommandLine(int argc, char **argv)
{
    /** @todo implement command line options
        -h help
        -d domain instead of dimension
        -i integration options
        -c options of things to compute and output
        -seed input seed for random initial data
     **/

    // Options are all found below in the giant if else if thingy.
    // They should be alphabetical by short option (if it exists)

    bool ok = true;

    // initialize cmdln arguments
    cmdln.computes = false;
    cmdln.positions = false;
    cmdln.warmstart = false;
    cmdln.mode = -1;
    cmdln.number = 0;
    cmdln.dimension = 0;
    cmdln.numSteps = 0;
    cmdln.outputRate = 0;
    cmdln.positionRate = 0;
    cmdln.seed = 83476234;
    cmdln.precision = 6;
    cmdln.timeFinal = 0.0;
    cmdln.stepSize = 1e-3;
    cmdln.outputTime = 0.25;
    cmdln.radius = 0.0;
    cmdln.rotation = 0.0;
    cmdln.initialFile.clear();
    cmdln.logFile.clear();
    cmdln.outputFile.clear();
    cmdln.movieFile.clear();
    cmdln.domain = Constants::PLANE;
    cmdln.integrate = Constants::RK4;

    unsigned arg = 1;
    char *in;
    while (arg < argc)
    {
        in = argv[arg];
        if (strcmp(in, "-h") == 0 || strcmp(in, "-help") == 0)
        {
            std::string hstr("<-d domain [domain args]> [-f <init file>] [-i <integrator> [integrator args]] [-l <log_file>] [-m <rate> <file>] [-n <num steps>][-o <output rate> <output file>] [-p <precision>] [-s <step size>] [-t <end time>] [-v <num PVs>] [-w]");
            std::cerr << argv[0] << " " << hstr << std::endl;
            ok = false;
            break;
        }
        // else if (strcmp(in, "-c") == 0 || strcmp(in, "-comp") == 0)
        // {
        //     /** @todo add compute quantity name options here **/
        //     ++arg;
        //     in = argv[arg];
        //     if (isdigit(*in))
        //     {
        //         unsigned rate = strtol(in, NULL, 0);
        //         if (cmdln.outputRate == 0 || rate < cmdln.outputRate)
        //             cmdln.outputRate = rate;
        //         if (rate == 0)
        //             cmdln.computes = false;
        //     }
        //     else
        //     {
        //         std::cerr << "Unrecognized argument after -c option: "
        //                   << in << std::endl;
        //         ok = false;
        //         break;
        //     }
        // }
        else if (strcmp(in, "-d") == 0 || strcmp(in, "-dom") == 0)
        {
            if (!checkArgs(argc, argv, 1, arg))
            {
                std::cerr << "-d option requires domain or dimension "
                          << "argument" << std::endl;
                ok = false;
                break;
            }
            ++arg;
            in = argv[arg];
            if (isdigit(*in))
            {
                cmdln.dimension = strtol(in, NULL, 0);
                if (cmdln.dimension == 2)
                    cmdln.domain = Constants::PLANE;
                else if (cmdln.dimension == 3)
                    cmdln.domain = Constants::SPHERE;
                else
                {
                    std::cerr << "Dimension " << cmdln.dimension
                              << " is not currently supported" << std::endl;
                    ok = false;
                    break;

                }
            }
            // else if (strcmp(in, "sph") == 0 ||
            //          strcmp(in, "sphere") == 0 ||
            //          strcmp(in, "S2") == 0)
            // {
            //     cmdln.domain = Constants::SPHERE;
            //     cmdln.dimension = 3;
            // }
            else if (strcmp(in, "plane") == 0 ||
                     strcmp(in, "unbounded") == 0 ||
                     strcmp(in, "R2") == 0)
            {
                cmdln.domain = Constants::PLANE;
                cmdln.dimension = 2;
            }
            else if (strcmp(in, "screened") == 0)
            {
                cmdln.domain = Constants::SCREENED;
                cmdln.dimension = 2;
                if (!checkArgs(argc, argv, 1, arg))
                {
                    std::cerr << "-d screened <deformation radius> "
                              << "::requires argument giving deformation radius"
                              << std::endl;
                    ok = false;
                    break;
                }
                ++arg;
                in = argv[arg];
                cmdln.radius = strtod(in, NULL);
                if (cmdln.radius == 0.0)
                {
                    std::cerr << "Deformation radius given with '-d screened' "
                              << "not understood or non-positive."
                              << std::endl;
                    ok = false;
                    break;
                }
            }
            // else if (strcmp(in, "rotating sphere") == 0 ||
            //          strcmp(in, "RS2") == 0)
            // {
            //     cmdln.domain = Constants::ROTATING_SPHERE;
            //     cmdln.dimension = 3;
            //     if (!checkArgs(argc, argv, 1, arg))
            //     {
            //         std::cerr << "Rotating sphere requires rotation "
            //                   << "rate argument" << std::endl;
            //         ok = false;
            //         break;
            //     }
            //     ++arg;
            //     in = argv[arg];
            //     cmdln.rotation = strtol(in, NULL, 0);
            // }
            else if (strcmp(in, "two-layer") == 0)
            {
                cmdln.domain = Constants::TWO_LAYER;
                cmdln.dimension = 2;
                if (!checkArgs(argc, argv, 1, arg))
                {
                    std::cerr << "-d twolayer <deformation radius> "
                              << "::requires argument giving deformation radius"
                              << std::endl;
                    ok = false;
                    break;
                }
                ++arg;
                in = argv[arg];
                cmdln.radius = strtod(in, NULL);
                if (cmdln.radius == 0.0)
                {
                    std::cerr << "Deformation radius given with '-d twolayer' "
                              << "not understood or non-positive."
                              << std::endl;
                    ok = false;
                    break;
                }
            }
            else
            {
                std::cerr << "Unrecognized domain -- " << in << std::endl
                          << "Currently unbounded (default) and screened "
                          << "are supported" << std::endl;
                ok = false;
                break;
            }
        }
        else if (strcmp(in, "-f") == 0 || strcmp(in, "-init") == 0)
        {
            if(!checkArgs(argc, argv, 1, arg))
            {
                std::cerr << "-f/-init option requires file name" << std::endl;
                ok = false;
                break;
            }
            ++arg;
            cmdln.initialFile = argv[arg];
        }
        else if (strcmp(in, "-i") == 0 || strcmp(in,"-int") == 0)
        {
            if(!checkArgs(argc, argv, 1, arg))
            {
                std::cerr << "-i/-int option requires integrator name"
                          << std::endl;
                ok = false;
                break;
            }

            ++arg;
            in = argv[arg];
            if (strcmp(in, "rk4") == 0 || strcmp(in, "RK4") == 0 ||
                strcmp(in, "std") == 0)
            {
                cmdln.integrate = Constants::RK4;
            }
            else if (strcmp(in, "rk45") == 0 || strcmp(in, "RK45") == 0 ||
                     strcmp(in, "adapt") == 0)
            {
                cmdln.integrate = Constants::RK45;
            }
            else if (strcmp(in, "pc") == 0 || strcmp(in, "PC") == 0 ||
                     strcmp(in, "adams") == 0)
            {
                cmdln.integrate = Constants::ADMOUL;
            }
            else
            {
                std::cerr << "Unrecognized integrator option -- '"
                          << in << "'" << std::endl;
            }
        }
        else if (strcmp(in, "-l") == 0 || strcmp(in, "-log") == 0)
        {
            if (!checkArgs(argc, argv, 1, arg))
            {
                std::cerr << "-l/-log option requires file name" << std::endl;
                ok = false;
                break;
            }
            ++arg;
            cmdln.logFile = argv[arg];
        }
        else if (strcmp(in, "-m") == 0 || strcmp(in, "-movie") == 0)
        {
            if (!checkArgs(argc, argv, 2, arg))
            {
                std::cerr << "-m/-movie requires 2 arguments -- "
                          << "<rate> <file>" << std::endl;
                ok = false;
                break;
            }
            cmdln.positions = true;
            // Rate
            ++arg;
            in = argv[arg];
            unsigned rate = strtol(in, NULL, 0);
            if (cmdln.positionRate == 0 || rate < cmdln.positionRate)
                cmdln.positionRate = rate;
            // File
            ++arg;
            cmdln.movieFile = argv[arg];
        }
        else if (strcmp(in, "-n") == 0 || strcmp(in, "-num") == 0)
        {
            if (!checkArgs(argc, argv, 1, arg))
            {
                std::cerr << "-n option requires number of steps argument"
                          << std::endl;
                ok = false;
                break;
            }
            ++arg;
            cmdln.numSteps = strtol(argv[arg], NULL, 0);
            cmdln.mode = 1;
        }
        else if (strcmp(in, "-o") == 0 || strcmp(in, "-out") == 0)
        {
            if (!checkArgs(argc, argv, 2, arg))
            {
                std::cerr << "-o/-out option requires output rate and file name"
                          << std::endl;
                ok = false;
                break;
            }
            cmdln.computes = true;

            // Rate
            ++arg;
            in = argv[arg];
            unsigned rate = strtol(in, NULL, 0);
            if (cmdln.outputRate == 0 || rate < cmdln.outputRate)
                cmdln.outputRate = rate;
            // File
            ++arg;
            cmdln.outputFile = argv[arg];
        }
        else if (strcmp(in, "-p") == 0 || strcmp(in, "-prec") == 0)
        {
            if (!checkArgs(argc, argv, 1, arg))
            {
                std::cerr << "-p/-prec option requires precision argument"
                          << std::endl;
                ok = false;
                break;
            }
            ++arg;
            cmdln.precision = strtol(argv[arg], NULL, 0);
            std::cout.precision(cmdln.precision);
            std::cerr.precision(cmdln.precision);
        }
        else if (strcmp(in, "-s") == 0 || strcmp(in, "-step") == 0)
        {
            if (!checkArgs(argc, argv, 1, arg))
            {
                std::cerr << "-s/-step requires step size argument"
                          << std::endl;
                ok = false;
                break;
            }
            ++arg;
            cmdln.stepSize = strtod(argv[arg], NULL);
        }
        else if (strcmp(in, "-t") == 0 || strcmp(in, "-time") == 0)
        {
            if (!checkArgs(argc, argv, 1, arg))
            {
                std::cerr << "-t/-time requires end time argument"
                          << std::endl;
                ok = false;
                break;
            }
            ++arg;
            cmdln.timeFinal = strtod(argv[arg], NULL);
            cmdln.mode = 0;
        }
        else if (strcmp(in, "-v") == 0 || strcmp(in, "-vor") == 0)
        {
            if (!checkArgs(argc, argv, 1, arg))
            {
                std::cerr << "-v/-vor option requires number of vortices "
                          << "argument" << std::endl;
                ok = false;
                break;
            }
            ++arg;
            in = argv[arg];
            cmdln.number = strtol(in, NULL, 0);
        }
        else if (strcmp(in, "-w") == 0 || strcmp(in, "-warmstart") == 0)
        {
            if (!checkArgs(argc, argv, 1, arg))
            {
                std::cerr << "-w/-warmstart option requires position "
                          << "file from previous run" << std::endl;
                ok = false;
                break;
            }
            ++arg;
            cmdln.initialFile = argv[arg];
            cmdln.warmstart = true;
        }
        else
        {
            std::cerr << "Unknown input argument: " << in << std::endl;
            ok = false;
            break;
        }
        ++arg;
    }

    return ok;
}

bool VortexSim::runNumSteps(unsigned steps)
{
    bool ok = true;
    unsigned count = 0;
    double time;

    while (count < steps)
    {
        if (integrator->step() == Integrator::FAIL)
        {
            log->writeMsg(Logger::ERROR, "Timestep failed. Run aborted.");
            ok = false;
            break;
        }

        // Output quantities of interest to file
        if (cmdln.computes && integrator->stepNum() % cmdln.outputRate == 0)
        {
            time = integrator->stepNum() * integrator->stepSize();
            if (!outhq->processOutput(storage, time))
            {
                log->writeMsg(Logger::ERROR,
                              "Error computing output. Run aborted.");
                ok = false;
                break;
            }
        }

        // Output positions to file
        if (cmdln.positions && integrator->stepNum() % cmdln.positionRate == 0)
            savePositions();

        ++count;
    }
    return ok;
}

void VortexSim::summarizeSetup()
{
    std::stringstream sum;
    sum.precision(cmdln.precision);
    // Summarize my info -- really storage...
    sum << "Number Vortices: " << storage->size() << std::endl
        << "Domain: ";
    if (cmdln.domain == Constants::PLANE)
    {
        sum << "PLANE";
    }
    else if (cmdln.domain == Constants::SCREENED)
    {
        sum << "SCREENED" << std::endl
            << "Deformation Radius: " << cmdln.radius;
    }
    else if (cmdln.domain == Constants::SPHERE)
    {
        sum << "SPHERE";
    }
    else if (cmdln.domain == Constants::ROTATING_SPHERE)
    {
        sum << "ROTATING SPHERE";
    }
    else if (cmdln.domain == Constants::TWO_LAYER)
    {
        sum << "TWO LAYER" << std::endl
            << "Deformation Radius: " << cmdln.radius;
    }
    else
    {
        sum << "Unrecognized domain: " << cmdln.domain;
    }
    sum << std::endl;
    sum << "Total Vorticity: " << storage->totalvort() << std::endl
        << "Initial Condition: "
        << (cmdln.initialFile.empty() ? "random" : cmdln.initialFile)
        << std::endl
        << "End Time: " << cmdln.timeFinal << std::endl
        << "Stepsize: " << integrator->stepSize() << std::endl
        << "Number of Steps: " << cmdln.numSteps << std::endl
        << "Housekeeping Interval: " << batchSize << std::endl;

    // Summarize flow computer info
    sum << "Flow Type: " << flowcomp->name() << std::endl;

    // Summarize integrator info
    sum << "Integrator: " << integrator->name() << std::endl;

    // Summarize output info
    outhq->summarize(sum);
    sum << "Output Rate: " << cmdln.outputRate << std::endl;

    // Summarize movie output info
    if (!cmdln.movieFile.empty())
    {
        sum << "Movie File: " << cmdln.movieFile << std::endl
            << "Output Rate: " << cmdln.positionRate << std::endl;
    }

    if(log->writeMsg(Logger::INFO, sum.str()))
        std::cout << sum.str();
}

bool VortexSim::initializeMovieFile()
{
    bool ok = true;

    if (cmdln.positions)
    {
        movie = new std::ofstream(cmdln.movieFile);
        if (!movie->is_open())
        {
            ok = false;
            std::stringstream msg;
            msg << "Storage--Failure to open movie file '"
                << cmdln.movieFile << "'";
            log->writeMsg(Logger::ERROR, msg.str());
        }
        movie->precision(15);
    }

    return ok;
}

void VortexSim::savePositions()
{
    double time = integrator->stepNum() * integrator->stepSize();
    (*movie) << "Simulation Time=" << time << std::endl
             << *storage << std::endl;
    timer->stamp(Timer::OUT);
}

void VortexSim::synchTimeAndStep()
{
    if (cmdln.mode == 0)
    {
        cmdln.numSteps = (unsigned)(round(cmdln.timeFinal / cmdln.stepSize));
    }
    else if (cmdln.mode == 1)
    {
        cmdln.timeFinal = cmdln.numSteps * cmdln.stepSize;
    }
    else
    {
        // Default -- go to simulation time 10.0
        cmdln.timeFinal = 10.0;
        cmdln.numSteps = (unsigned)(round(cmdln.timeFinal / cmdln.stepSize));
    }

    // Determine interval between housekeeping tasks
    unsigned steps = cmdln.numSteps/10;
    batchSize = (steps < 1000) ? steps : 1000;
}

void VortexSim::summarizeRun()
{
    std::stringstream sum;
    sum.precision(cmdln.precision);

    // Print initial quantities
    outhq->processOutput(initial, 0.0, sum);

    // Print final quantities
    double time = integrator->stepNum() * integrator->stepSize();
    outhq->processOutput(storage, time, sum);

    // Print timing information to the log and screen
    timer->printTimingInfo(sum);
    if(log->writeMsg(Logger::INFO, sum.str()))
        std::cout << sum.str() << std::endl;
}

void VortexSim::reportProgress(unsigned width, double prog) const
{
    // unsigned pos = int(prog * width);
    // std::cout << "[";
    // for (unsigned k = 0; k < width; ++k)
    // {
    //     if (k < pos)
    //         std::cout << ":";
    //     else
    //         std::cout << " ";
    // }
    // std::cout.precision(3);
    // std::cout << "]" << prog * 100.0 << "%\r";
    // std::cout.precision(cmdln.precision);
    // std::cout.flush();

    char buff[30];
    sprintf(buff, "\rProgress: %4.1f%%", prog * 100.0);
    std::cout << buff;
    std::cout.flush();
}
