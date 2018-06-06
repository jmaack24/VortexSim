/*******************************************************************************
 * \class OutputHQ
 *
 * \brief Class for knowing and managing all outputs
 *
 * Author: Jonathan Maack
 * 
 ******************************************************************************/

#include "OutputHQ.h"

#include "Logger.h"
#include "Output.h"
#include "Timer.h"

// Computable classes
//#include "Autocorrelation.h"
#include "CenterVort.h"
//#include "Covariance.h"
//#include "EquilibriumError.h"
#include "Hamiltonian.h"
#include "MomentInertia.h"
#include "Moments.h"
#include "ScreenedHamiltonian.h"
#include "TwoLayerCenter.h"
#include "TwoLayerHamiltonian.h"
#include "TwoLayerMoments.h"

#include <fstream>
#include <iostream>
#include <sstream>

OutputHQ::OutputHQ(Logger *ln, Timer *tm):
    log(ln), timer(tm), outputs(0), fout(0)
{
    out = &std::cout;
    //log->writeMsg(Logger::INFO, "OutputHQ construction complete");
}

OutputHQ::~OutputHQ()
{
    unsigned num = outputs.size();
    Output *op;
    for (unsigned k = 0; k < num; ++k)
    {
        op = outputs[k];
        delete op;
    }
    outputs.clear();

    if (fout != 0)
    {
        fout->close();
        delete fout;
        fout = 0;
    }
    out = 0;
    timer = 0;
    log = 0;
}

bool OutputHQ::processOutput(const Storage *store, double time)
{
    bool good = processOutput(store, time, *out);

    return good;
}

bool OutputHQ::processOutput(const Storage *store, double time,
                             std::ostream &os)
{
    bool good = true;

    Output *op;
    unsigned num = outputs.size();
    os << "Simulation Time: " << time << std::endl;
    for (unsigned k = 0; k < num; ++k)
    {
        // "Indent" all outputs
        os << "  ";
        op = outputs[k];
        good = op->computeOutput(store, os);
        if (!good)
        {
            std::stringstream err;
            err << "Error computing " << op->name();
            log->writeMsg(Logger::ERROR, err.str());
            break;
        }
    }
    os << std::endl;

    timer->stamp(Timer::OUT);

    return good;
}

bool OutputHQ::setup(const std::string &file, unsigned numEntries,
                     const Storage *st, Constants::Domain dom, double defrad)
{
    bool good = true;

    if (!file.empty())
    {
        this->file = file;
        fout = new std::ofstream(file.c_str());
        if (!fout->is_open())
        {
            std::stringstream ss;
            ss << "Failed to open output file '" << file << "'";
            log->writeMsg(Logger::ERROR, ss.str());
            good = false;
        }
        else
        {
            // Set output object to file
            out = fout;
        }
    }

    if (good)
    {
        // Build all the outputs
        outputs.push_back(new CenterVort(log));
        outputs.push_back(new MomentInertia(log));
        if (dom == Constants::PLANE)
        {
            outputs.push_back(new Hamiltonian(log));
            outputs.push_back(new Moments(log));
        }
        else if (dom == Constants::SCREENED)
        {
            outputs.push_back(new ScreenedHamiltonian(log, defrad));
            outputs.push_back(new Moments(log));
        }
        else if (dom == Constants::TWO_LAYER)
        {
            outputs.push_back(new TwoLayerHamiltonian(log, defrad));
            outputs.push_back(new TwoLayerCenter(log));
            outputs.push_back(new TwoLayerMoments(log));
        }
        else
        {
            log->writeMsg(Logger::ERROR, "Unimplemented domain type");
            good = false;
        }
        //outputs.push_back(new Covariance(log));
        //outputs.push_back(new Autocorrelation(st, log));
        //outputs.push_back(new EquilibriumError(log, st));

        // Max out the precision on printing
        out->precision(15);

        // Save number of entries
        if (fout != 0)
            (*fout) << "Number of Entries=" << numEntries << std::endl
                    << std::endl;
    }

    return good;
}

void OutputHQ::summarize(std::ostream &out) const
{
    out << "Outputting to: '"
        << (file.empty() ? "terminal" : file) << "'"
        << std::endl
        << "Outputting: ";
    for (unsigned k = 0; k < outputs.size(); ++k)
    {
        out << outputs[k]->name();
        if (k < outputs.size() - 1)
        {
            out << ", ";
        }
    }
    out << std::endl;
}
