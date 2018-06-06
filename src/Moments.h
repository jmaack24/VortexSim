/*******************************************************************************
 * \class Moments
 *
 * \brief Calculate and output moments of the current PV configuration --
 * ONLY DOES DIMENSION 2!!!!
 *
 * Author: Jonathan Maack
 *
 ******************************************************************************/

#ifndef MOMENTS
#define MOMENTS

#include "Output.h"

#include <string>

class Logger;

class Moments: public Output
{
public:
    /** Default constructor -- intentionally unimplmented **/
    Moments() {}
    /** Constructor **/
    Moments(Logger *ln);
    /** Destructor **/
    virtual ~Moments();

    bool computeOutput(const Storage *st, std::ostream &out);

    virtual inline const char* name() const
    {
        return Name.c_str();
    }
private:
    std::string Name;

    Logger *log;

    double *second;
    double *third;
    double *fourth;

    unsigned dim;

    void buildMomentStorage(unsigned dimension);
    void cleanupStorage();
};

#endif
