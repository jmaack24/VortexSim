/*******************************************************************************
 * \class TwoLayerMoments
 *
 * \brief Calculate and output moments of the current PV configuration --
 * ONLY DOES DIMENSION 2!!!!
 *
 * Author: Jonathan Maack
 *
 ******************************************************************************/

#include "TwoLayerMoments.h"

#include "PointVortex.h"
#include "Storage.h"

TwoLayerMoments::TwoLayerMoments(Logger *ln):
    Output(),
    Name("Moments"),
    second(0),
    dim(0)
{
    buildMomentStorage(2);
}

TwoLayerMoments::~TwoLayerMoments()
{
    cleanupStorage();
}

bool TwoLayerMoments::computeOutput(const Storage *st, std::ostream &out)
{
    // Compute moments
    bool status = true;

    if (dim != st->dimension())
    {
        status = false;
    }
    else
    {
        unsigned N = st->size();
        unsigned i, j;
        // Clear calculation buffers
        for (j = 0; j < 2; ++j)
            for (i = 0; i < 3; ++i)
                second[j][i] = 0.0;
        PointVortex *pv;
        double x, y;
        int l;
        // Calculate moments
        j = 0;
        for (i = 0; i < N; ++i)
        {
            pv = st->retrieve(i);
            x = pv->getPos(0);
            y = pv->getPos(1);
            l = pv->getLayer();
            if (l == 1)
            {
                ++j;
                second[0][0] += x * x / N;
                second[0][1] += x * y / N;
                second[0][2] += y * y / N;
            }
            else if (l == -1)
            {
                second[1][0] += x * x / N;
                second[1][1] += x * y / N;
                second[1][2] += y * y / N;
            }
        }
        for (i = 0; i < 3; ++i)
        {
            second[0][i] *= double(N)/double(j);
            second[1][i] *= double(N)/double(N - j);
        }

        // Print stuff
        out << "Second Moment 1: [" << second[0][0] << ", "
            << second[0][1] << ", " << second[0][2] << "]" << std::endl;
        out << "  Second Moment 2: [" << second[1][0] << ", "
            << second[1][1] << ", " << second[1][2] << "]" << std::endl;
    }

    return status;
}

void TwoLayerMoments::buildMomentStorage(unsigned dimension)
{
    dim = 2;
    second = new double *[2];
    for (int i = 0; i < 2; ++i)
        second[i] = new double[3];
}

void TwoLayerMoments::cleanupStorage()
{
    for (int i = 0; i < 2; ++i)
        delete [] second[i];
    delete [] second;
    second = 0;
}
