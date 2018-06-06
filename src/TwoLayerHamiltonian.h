/*******************************************************************************
 * \class TwoLayerHamiltonian
 *
 * \brief Computes and outputs the Hamiltonian of the vortex system for a two
 * layer quasigeostrophic model
 *
 * Author: Jonathan Maack
 *
 ******************************************************************************/

#ifndef TWO_LAYER_HAMILTONIAN
#define TWO_LAYER_HAMILTONIAN

#include "Output.h"

#include <string>

class Logger;

class TwoLayerHamiltonian: public Output
{
public:
    /** Constructor **/
    TwoLayerHamiltonian() {}
    /** Constructor **/
    TwoLayerHamiltonian(Logger *ln, double defrad);
    /** Destructor **/
    virtual ~TwoLayerHamiltonian();

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
