/*******************************************************************************
 * \class TwoLayerHamiltonian
 *
 * \brief Computes and outputs the Hamiltonian of the vortex system for the
 * two layer quasigeostrophic model
 *
 * Author: Jonathan Maack
 *
 ******************************************************************************/

#include "TwoLayerHamiltonian.h"

#include "Constants.h"
#include "Logger.h"
#include "PointVortex.h"
#include "Storage.h"

#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include <sstream>

TwoLayerHamiltonian::TwoLayerHamiltonian(Logger *ln, double defrad):
    Output(),
    theLog(ln),
    Name("Hamiltonian")
{
    deffreq = sqrt(2.0) / defrad;
}

TwoLayerHamiltonian::~TwoLayerHamiltonian()
{
    theLog = 0;
}

bool TwoLayerHamiltonian::computeOutput(const Storage *st, std::ostream &out)
{
    bool ok = true;

    // Compute
    unsigned N = st->size();
    unsigned dim = st->dimension();

    const PointVortex *pv1;
    const PointVortex *pv2;
    double r, rsq, v1, v2, dx;
    double H = 0.0;
    int sign;

    for (unsigned i = 0; i < N; ++i)
    {
        pv1 = st->retrieve(i);
        v1 = pv1->getVort();
        for (unsigned j = i+1; j < N; ++j)
        {
            pv2 = st->retrieve(j);
            v2 = pv2->getVort();
            rsq = 0.0;
            sign = pv1->getLayer()*pv2->getLayer();
            for (unsigned l = 0; l < dim; ++l)
            {
                dx = pv1->getPos(l) - pv2->getPos(l);
                rsq += dx * dx;
            }
            r = sqrt(rsq);
            H += 0.5*v1*v2 * (log(r) - sign * gsl_sf_bessel_K0(deffreq * r));
        }
    }

    H *= -1.0/Constants::TWO_PI;

    // Output
    out << "Hamiltonian: " << H << std::endl;

    return ok;
}
