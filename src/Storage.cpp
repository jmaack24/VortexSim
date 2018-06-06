/**
 * \class Storage
 *
 * \brief Stores and manages access to point vortices.
 *
 * Author: Jonathan Maack
 *
 * Description: Responsible for storing the point vortex objects of the
 * simulation.  Also controls access to said vortices.
 **/

#include "Storage.h"

#include "PointVortex.h"

#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <random>
#include <vector>

Storage::Storage(const Storage &store)
{
    log = store.log;
    num = store.num;
    dim = store.dim;
    totalVorticity = store.totalVorticity;

    vortices.reserve(num);
    std::vector<PointVortex*>::const_iterator iter;
    for(iter = store.vortices.begin(); iter != store.vortices.end(); ++iter)
    {
        vortices.push_back(new PointVortex(**iter));
    }
}

Storage::Storage(Logger *ln, unsigned N, unsigned d):
    log(ln), num(N), dim(d), vortices(0), totalVorticity(0.0)
{
    vortices.reserve(num);
    for (unsigned k = 0; k < num; ++k)
    {
        vortices.push_back(new PointVortex(dim));
    }

    //log->writeMsg(Logger::INFO, "Storage construction complete");
}

Storage::~Storage()
{
    clear();
    log = 0;
}

void Storage::copyPVPos(const Storage *cp)
{
    std::vector<PointVortex*>::const_iterator citer;
    std::vector<PointVortex*>::iterator uiter = vortices.begin();
    for(citer = cp->vortices.begin(); citer != cp->vortices.end(); ++citer)
    {
        (*uiter)->copyPos(*citer);
        ++uiter;
    }
}

bool Storage::initializeFromInputFile(std::ifstream &in, Constants::Domain dom)
{
    bool ok = true;

    char buffer[1024];
    char *status;
    char *token;
    unsigned count = 0;
    int layer;
    double value;
    PointVortex *pv;

    bool go = true;

    in.getline(buffer, 1023);
    while (go && !in.eof())
    {
        // Check to see if this is a point vortex entry or a comment
        token = strtok(buffer, " \t");
        if (token == 0)
        {
            std::stringstream ss;
            ss << "Incorrect number of arguments to initialize vortex "
               << "-- expected " << dim + 1
               << " arguments and saw 0";
            log->writeMsg(Logger::ERROR, ss.str());
            go = false;
            ok = false;
        }
        else if (buffer[0] != '#')
        {
            // Create Point Vortex object
            pv = new PointVortex(dim);
            vortices.push_back(pv);

            // Get vorticity
            value = strtod(token, &status);
            if (status == token)
            {
                std::stringstream ss;
                ss << "No numeric value found -- saw '"
                   << token << "'";
                log->writeMsg(Logger::ERROR, ss.str());
                go = false;
                ok = false;
                continue;
            }
            else
            {
                totalVorticity += value;
                pv->setVort(value);
            }

            // Get position
            for (unsigned k = 0; k < dim; ++k)
            {
                token = strtok(NULL, " \t");
                if (token == 0)
                {
                    std::stringstream ss;
                    ss << "Incorrect number of arguments to initialize vortex "
                       << "-- expected " << dim + 1
                       << " arguments and saw " << k + 1;
                    log->writeMsg(Logger::ERROR, ss.str());
                    ok = false; // return bad flag
                    go = false; // make sure we exit the while loop
                    break; // exit for loop
                }
                value = strtod(token, &status);
                if (status == token)
                {
                    std::stringstream ss;
                    ss << "No numeric value found -- saw '"
                       << token << "'";
                    log->writeMsg(Logger::ERROR, ss.str());
                    ok = false; // return bad flag
                    go = false; // exit while loop
                    break; // exit for loop
                }
                else
                {
                    pv->setPos(k, value);
                }
            }

            // Get layer if applicable
            if (dom == Constants::TWO_LAYER)
            {
                token = strtok(NULL, " \t");
                layer = strtol(token, &status, 0);
                if (status == token)
                {
                    std::stringstream ss;
                    ss << "Failed to find numberic entry for layer -- saw '"
                       << token << "'";
                    log->writeMsg(Logger::ERROR, ss.str());
                    ok = false;
                }
                else
                {
                    pv->setLayer(layer);
                }
            }
            ++count;
        }
        in.getline(buffer, 1023);
    }

    num = vortices.size();

    return ok;
}

bool Storage::initializeFromPositionFile(std::ifstream &in,
                                         Constants::Domain dom)
{
    bool ok = true;

    char buffer[1024];
    char *status;
    char *token;
    unsigned count = 0;
    unsigned numEntries = 0;
    unsigned number;
    unsigned dimension;
    int layer;
    double value;
    PointVortex *pv;

    // Find the "Number of entries" part of the file
    bool go = true;

    while (go && !in.eof())
    {
        in.getline(buffer, 1023);
        token = strtok(buffer, "=");
        if (strcmp(token, "Number of Entries") == 0)
        {
            token = strtok(NULL, " \t\n");
            numEntries = strtol(token, &status, 0);
            if (status == token)
            {
                std::stringstream ss;
                ss << "No numeric value found for 'Number of Entries' -- saw '"
                   << token << "'";
                log->writeMsg(Logger::ERROR, ss.str());
                ok = false;
            }
            go = false;
        }
    }

    if (ok && numEntries == 0)
    {
        std::stringstream ss;
        ss << "'Number of Entries' not found in file -- unable to initialize";
        log->writeMsg(Logger::ERROR, ss.str());
    }
    else if (ok)
    {
        go = true;
    }

    // Find the last entry to read in
    while (go && !in.eof())
    {
        in.getline(buffer, 1023);
        token = strtok(buffer, "=");
        if (token != 0 && strcmp(token, "Simulation Time") == 0)
            ++count;
        if (count == numEntries)
            go = false;
    }

    // Check if we found the entry we were looking for or if we hit the end
    if (go)
    {
        std::stringstream ss;
        ss << "Expected " << numEntries << " and only found " << count
           << " -- unable to initialize";
        log->writeMsg(Logger::ERROR, ss.str());
        ok = false;
        go = false;
    }

    // Get number of vortices
    in.getline(buffer, 1023);
    token = strtok(buffer, "=");
    if (strcmp(token, "Number") == 0)
    {
        token = strtok(NULL, " \t\n");
        number = strtol(token, &status, 0);
        if (status == token)
        {
            std::stringstream ss;
            ss << "Number of vortices is not numeric -- saw '"
               << token << "'";
            log->writeMsg(Logger::ERROR, ss.str());
            ok = false;
        }
    }
    else
    {
        std::stringstream ss;
        ss << "Failed to find 'Number' entry giving the number of vortices"
           << " -- saw '" << token << "'";
        log->writeMsg(Logger::ERROR, ss.str());
        ok = false;
    }

    // Get dimension
    in.getline(buffer, 1023);
    token = strtok(buffer, "=");
    if (strcmp(token, "Dimension") == 0)
    {
        token = strtok(NULL, " \t\n");
        dimension = strtol(token, &status, 0);
        if (status == token)
        {
            std::stringstream ss;
            ss << "Failed to find numberic entry for dimension -- saw '"
               << token << "'";
            log->writeMsg(Logger::ERROR, ss.str());
            ok = false;
        }
        else if (dimension != dim)
        {
            std::stringstream ss;
            ss << "Dimension of file does not match command line dimension. "
               << "Expected " << dim << " saw " << dimension;
            log->writeMsg(Logger::ERROR, ss.str());
            ok = false;
        }
    }
    else
    {
        std::stringstream ss;
        ss << "Failed to find 'Dimension' entry giving the dimension "
           << "of the simulation -- saw '" << token
           << "'";
        log->writeMsg(Logger::ERROR, ss.str());
        ok = false;
    }

    go = ok;
    count = 0;
    in.getline(buffer, 1023);
    while (go && !in.eof())
    {
        // Create Point Vortex object
        pv = new PointVortex(dim);
        vortices.push_back(pv);

        // Get vorticity
        token = strtok(buffer, "=");
        if (strcmp(token, "Vort") == 0)
        {
            token = strtok(NULL, " ");
            value = strtod(token, &status);
            if (status == token)
            {
                std::stringstream ss;
                ss << "No numeric vorticity value found for vortex " << count
                   << " -- saw '" << token << "'";
                log->writeMsg(Logger::ERROR, ss.str());
                ok = false;
                break;
            }
            totalVorticity += value;
            pv->setVort(value);
        }
        else
        {
            std::stringstream ss;
            ss << "No vorticity value found for vortex " << count
               << " -- saw '" << token << "'";
            log->writeMsg(Logger::ERROR, ss.str());
            ok = false;
            break;
        }

        // Get position -- hacked a bit to get around strtok issues
        // Need the leading '(' removed but I can't break on it and '=' so
        // just separate by the the parenthesis and include the equal sign
        // in the strcmp call...
        token = strtok(NULL, "(");
        if (strcmp(token, "Pos=") == 0)
        {
            for (unsigned k = 0; k < dim; ++k)
            {
                token = strtok(NULL, ",)");
                value = strtod(token, &status);
                if (status == token)
                {
                    std::stringstream ss;
                    ss << "No numeric entry found for coordinate " << k
                       << " of vortex " << count << " -- saw '"
                       << token << "'";
                    log->writeMsg(Logger::ERROR, ss.str());
                    ok = false;
                    go = false; // end the while loop
                    break; // break out of for loop
                }
                pv->setPos(k, value);
            }
        }
        else
        {
            std::stringstream ss;
            ss << "No position value found for vortex " << count
               << " -- saw '" << token << "'";
            log->writeMsg(Logger::ERROR, ss.str());
            ok = false;
            break;
        }

        if (dom == Constants::TWO_LAYER)
        {
            token = strtok(NULL, "=");
            if (strcmp(token, "Layer") == 0)
            {
                token = strtok(NULL, " \t\n");
                layer = strtol(token, &status, 0);
                if (status == token)
                {
                    std::stringstream ss;
                    ss << "Failed to find numberic entry for layer -- saw '"
                       << token << "'";
                    log->writeMsg(Logger::ERROR, ss.str());
                    ok = false;
                }
                else
                {
                    pv->setLayer(layer);
                }
            }
        }

        // Setup for next pass
        in.getline(buffer, 1023);
        ++count;
        if (count == number)
            go = false;
    }

    if (ok && number != vortices.size())
    {
        std::stringstream ss;
        ss << "Expected " << num << " vortices and found " << vortices.size()
           << " -- unable to initialize" << std::endl;
        log->writeMsg(Logger::ERROR, ss.str());
    }

    num = vortices.size();

    return ok;
}

bool Storage::initializeFromFile(bool warmstart, const std::string &file,
                                 Constants::Domain dom)
{
    bool ok = true;

    std::ifstream in(file);

    // Clear anything that came from the normal initialization
    clear();

    if(!in.is_open())
    {
        std::stringstream ss;
        ss << "Failed to open file '" << file << "' -- aborting.";
        log->writeMsg(Logger::ERROR, ss.str());
        ok = false;
    }
    else if(!in.good())
    {
        std::stringstream ss;
        ss << "Error with input file '" << file << "' -- aborting.";
        log->writeMsg(Logger::ERROR, ss.str());
        ok = false;
    }
    else if (warmstart)
    {
        ok = initializeFromPositionFile(in, dom);
    }
    else
    {
        ok = initializeFromInputFile(in, dom);
    }

    return ok;
}

bool Storage::initializeRandom(unsigned seed)
{
    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> unif(-1.0, 1.0);

    PointVortex *pv;
    for (unsigned k = 0; k < num; ++k)
    {
        pv = vortices[k];
        pv->setVort(1.0/num);
        for (unsigned l = 0; l < dim; ++l)
        {
            pv->setPos(l, unif(generator));
        }
    }

    totalVorticity = 1.0;

    return true;
}

PointVortex *Storage::retrieve(unsigned index) const
{
    PointVortex *retval = 0;
    if (index < num)
        retval = vortices[index];
    else
        log->writeMsg(Logger::WARN, "Storage--attempt to access out of range");
    return retval;
}

std::ostream& operator<<(std::ostream &out, const Storage& st)
{
    out << "Number=" << st.num << std::endl
        << "Dimension=" << st.dim << std::endl;
    for (unsigned k = 0; k < st.vortices.size(); ++k)
    {
        out << *(st.vortices[k]) << std::endl;
    }
    return out;
}

void Storage::clear()
{
    for (unsigned k = 0; k < num; ++k)
    {
        delete vortices[k];
    }
    vortices.clear();
    num = 0;
}
