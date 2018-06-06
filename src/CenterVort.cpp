/*******************************************************************************
 * \class CenterVort
 *
 * \brief Calculate and output center of vorticity
 *
 * Author: Jonathan Maack
 *
 ******************************************************************************/

#include "CenterVort.h"

#include "PointVortex.h"
#include "Storage.h"

CenterVort::CenterVort(Logger *ln):
    log(ln), center(0), dim(0), Name("Center of Vorticity")
{}

CenterVort::~CenterVort()
{
    if (center != 0)
    {
        delete center;
        center = 0;
    }
    log = 0;
}

bool CenterVort::computeOutput(const Storage *st, std::ostream &out)
{
    bool good = true;

    unsigned l;
    unsigned k;
    PointVortex *pv;
    double vort;

    unsigned N = st->size();
    double totalVort = 0.0;

    // Check if we need more space for calculation
    buildCenterSpace(st->dimension());

    // Clear calc space
    for (l = 0; l < dim; ++l)
    {
        center[l] = 0.0;
    }

    // Compute
    for (k = 0; k < N; ++k)
    {
        pv = st->retrieve(k);
        vort = pv->getVort();
        totalVort += vort;
        for (l = 0; l < dim; ++l)
        {
            center[l] += vort*pv->getPos(l);
        }
    }

    // Output the center of vorticity
    out << "Center of Vorticity: (";
    for (l = 0; l < dim; ++l)
    {
        out << center[l]/totalVort;
        if (l < dim - 1)
        {
            out << ", ";
        }
    }
    out << ")" << std::endl;

    return good;
}

void CenterVort::buildCenterSpace(unsigned dimension)
{
    if (dim < dimension)
    {
        delete center;

        center = new double[dimension];
        dim = dimension;
    }
}
