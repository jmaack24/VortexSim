/*******************************************************************************
 * \class Hamiltonian
 *
 * \brief Computes and outputs the Hamiltonian of the vortex system
 *
 * Author: Jonathan Maack
 *
 ******************************************************************************/

#ifndef HAMILTONIAN
#define HAMILTONIAN

#include "Output.h"

#include <string>

class Logger;

class Hamiltonian: public Output
{
public:
    /** Constructor **/
    Hamiltonian() {}
    /** Constructor **/
    Hamiltonian(Logger *ln);
    /** Destructor **/
    virtual ~Hamiltonian();

    virtual bool computeOutput(const Storage *st, std::ostream &out);

    inline virtual const char *name () const
    {
        return Name.c_str();
    }

private:
    const std::string Name;

    Logger *theLog;
};

#endif
