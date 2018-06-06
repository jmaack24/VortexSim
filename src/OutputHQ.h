/*******************************************************************************
 * \class OutputHQ
 *
 * \brief Class for knowing and managing all outputs
 *
 * Author: Jonathan Maack
 *
 ******************************************************************************/

#ifndef OUTPUT_HQ
#define OUTPUT_HQ

#include "Output.h"

#include "Constants.h"

#include <fstream>
#include <string>
#include <vector>

class Logger;
class Timer;

class OutputHQ
{
public:

    /** Constructor -- intentionally unimplemented **/
    OutputHQ(){}
    /**
     * Constructor
     * @param ln logger to use
     **/
    OutputHQ(Logger *ln, Timer *tm);
    /** Destructor **/
    virtual ~OutputHQ();

    /**
     * In which all outputs are computed and output
     * @param store collection of vortices to compute things for
     * @param time simulation time of the vortices in given storage object
     **/
    bool processOutput(const Storage *store, double time);

    /**
     * Computes outputs and writes them to then given ostream instead of the
     * internal one.
     * @param store collection of vortices to compute things for
     * @param time simulation time of the vortices in given storage object
     * @param out ostream to output things to
     **/
    bool processOutput(const Storage *store, double time, std::ostream &os);

    /**
     * Setup the output objects
     **/
    bool setup(const std::string &file, unsigned numEntries, const Storage *st,
               Constants::Domain dom, double defrad);

    /**
     * Provide summary of what will be printed to output file
     **/
    void summarize(std::ostream &out) const;

private:
    /** Object for holding all the objects **/
    std::vector<Output *> outputs;

    /** Logger **/
    Logger *log;
    /** Timer **/
    Timer *timer;

    /** Output file **/
    std::string file;
    /**
     * Output stream that will be used -- set to fout if printing to file
     * and to std::cout otherwise
     **/
    std::ostream *out;
    /** Output stream to file **/
    std::ofstream *fout;
};

#endif
