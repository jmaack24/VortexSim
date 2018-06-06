/*******************************************************************************
 * \class MomentInertia
 *
 * \brief Output object to compute the moment of inertia
 *
 * Author: Jonathan Maack
 *
 ******************************************************************************/

#ifndef MOMENT_INERTIA
#define MOMENT_INERTIA

#include "Output.h"

#include <string>

class Logger;

class MomentInertia: public Output
{
public:
    /** Constructor **/
    MomentInertia() {}
    /** Constructor **/
    MomentInertia(Logger *ln);
    /** Destructor **/
    virtual ~MomentInertia();

    virtual bool computeOutput(const Storage *st, std::ostream &out);

    inline virtual const char *name () const
    {
        return Name.c_str();
    }

private:
    const std::string Name;

    Logger *log;
};

#endif
