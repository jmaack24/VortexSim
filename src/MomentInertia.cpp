/*******************************************************************************
 * \class MomentInertia
 *
 * \brief Output object to compute the moment of inertia
 *
 * Author: Jonathan Maack
 *
 ******************************************************************************/

#include "MomentInertia.h"

#include "Logger.h"
#include "PointVortex.h"
#include "Storage.h"

MomentInertia::MomentInertia(Logger *ln):
    Output(),
    log(ln),
    Name("Moment of Inertia")
{
}

MomentInertia::~MomentInertia()
{
    log = 0;
}

bool MomentInertia::computeOutput(const Storage *st, std::ostream &out)
{
    bool good = true;

    unsigned k;
    unsigned l;
    const PointVortex *pv;
    double vort, xl, magsq;

    unsigned N = st->size();
    unsigned dim = st->dimension();
    double MI = 0.0;

    // Compute
    for (k = 0; k < N; ++k)
    {
        pv = st->retrieve(k);
        vort = pv->getVort();
        magsq = 0.0;
        for (l = 0; l < dim; ++l)
        {
            xl = pv->getPos(l);
            magsq += xl * xl;
        }
        MI += vort * magsq;
    }

    // Output the moment of inertia
    out << "Moment of Inertia: " << MI << std::endl;

    return good;
}
