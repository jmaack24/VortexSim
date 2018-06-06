/*******************************************************************************
 * \class ScreenedHamiltonian
 *
 * \brief Computes and outputs the Hamiltonian of the vortex system with a
 * given deformation radius
 *
 * Author: Jonathan Maack
 *
 ******************************************************************************/

#ifndef SCREENED_HAMILTONIAN
#define SCREENED_HAMILTONIAN

#include "Output.h"

#include <string>

class Logger;

class ScreenedHamiltonian: public Output
{
public:
    /** Constructor **/
    ScreenedHamiltonian() {}
    /** Constructor **/
    ScreenedHamiltonian(Logger *ln, double defrad);
    /** Destructor **/
    virtual ~ScreenedHamiltonian();

    virtual bool computeOutput(const Storage *st, std::ostream &out);

    inline virtual const char *name () const
    {
        return Name.c_str();
    }

private:
    const std::string Name;

    double deffreq;

    Logger *theLog;
};

#endif
