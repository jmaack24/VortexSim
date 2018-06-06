/*******************************************************************************
 * \class Hamiltonian
 *
 * \brief Computes and outputs the Hamiltonian of the vortex system
 *
 * Author: Jonathan Maack
 *
 ******************************************************************************/

#include "Hamiltonian.h"

#include "Constants.h"
#include "Logger.h"
#include "PointVortex.h"
#include "Storage.h"

#include <cmath>

Hamiltonian::Hamiltonian(Logger *ln):
    Output(),
    theLog(ln),
    Name("Hamiltonian")
{
}

Hamiltonian::~Hamiltonian()
{
    theLog = 0;
}

bool Hamiltonian::computeOutput(const Storage *st, std::ostream &out)
{
    bool ok = true;

    // Compute
    unsigned N = st->size();
    unsigned dim = st->dimension();

    const PointVortex *pv1;
    const PointVortex *pv2;
    double r, rsq, v1, v2, dx;
    double H = 0.0;

    for (unsigned i = 0; i < N; ++i)
    {
        pv1 = st->retrieve(i);
        v1 = pv1->getVort();
        for (unsigned j = i+1; j < N; ++j)
        {
            pv2 = st->retrieve(j);
            v2 = pv2->getVort();
            rsq = 0.0;
            for (unsigned l = 0; l < dim; ++l)
            {
                dx = pv1->getPos(l) - pv2->getPos(l);
                rsq += dx * dx;
            }
            r = sqrt(rsq);
            H += v1 * v2 * log(r);
        }
    }

    H *= -1.0/Constants::TWO_PI;

    // Output
    out << "Hamiltonian: " << H << std::endl;

    return ok;
}
