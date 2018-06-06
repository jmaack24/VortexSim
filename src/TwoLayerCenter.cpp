/*******************************************************************************
 * \class TwoLayerMoments
 *
 * \brief Calculate and output moments of the current PV configuration --
 * ONLY DOES DIMENSION 2!!!!
 *
 * Author: Jonathan Maack
 *
 ******************************************************************************/

#include "TwoLayerCenter.h"

#include "PointVortex.h"
#include "Storage.h"

TwoLayerCenter::TwoLayerCenter(Logger *ln):
    Output(),
    Name("Two-Layer Center"),
    center(0),
    dim(0)
{
    buildCenterStorage(2);
}

TwoLayerCenter::~TwoLayerCenter()
{
    cleanupStorage();
}

bool TwoLayerCenter::computeOutput(const Storage *st, std::ostream &out)
{
    // Compute moments
    bool status = true;

    if (dim != st->dimension())
    {
        status = false;
    }
    else
    {
        unsigned l;
        unsigned k;
        unsigned idx;
        unsigned N2 = 0;
        PointVortex *pv;

        unsigned N = st->size();

        // Clear calc space
        for (k = 0; k < 2; ++k)
        {
            for (l = 0; l < dim; ++l)
            {
                center[k][l] = 0.0;
            }
        }

        // Compute
        for (k = 0; k < N; ++k)
        {
            pv = st->retrieve(k);
            idx = ((pv->getLayer() == -1) ? 0 : 1);
            N2 += idx;
            for (l = 0; l < dim; ++l)
            {
                center[idx][l] += pv->getPos(l);
            }
        }

        for (l = 0; l < dim; ++l)
        {
            center[0][l] /= (N-N2);
            center[1][l] /= N2;
        }

        // Output the center of vorticity
        out << "Center 1: (" << center[0][0] << ", "
            << center[0][1] << ")" << std::endl
            << "  Center 2: (" << center[1][0] << ", "
            << center[1][1] << ")" << std::endl;
    }

    return status;
}

void TwoLayerCenter::buildCenterStorage(unsigned dimension)
{
    dim = dimension;
    center = new double *[2];
    for (int i = 0; i < 2; ++i)
        center[i] = new double[2];
}

void TwoLayerCenter::cleanupStorage()
{
    for (int i = 0; i < 2; ++i)
        delete [] center[i];
    delete [] center;
    center = 0;
}
