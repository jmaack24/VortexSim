/*******************************************************************************
 * \class Moments
 *
 * \brief Calculate and output moments of the current PV configuration --
 * ONLY DOES DIMENSION 2!!!!
 *
 * Author: Jonathan Maack
 *
 ******************************************************************************/

#include "Moments.h"

#include "PointVortex.h"
#include "Storage.h"

Moments::Moments(Logger *ln):
    Output(),
    Name("Moments"),
    second(0), third(0), fourth(0),
    dim(0)
{
    buildMomentStorage(2);
}

Moments::~Moments()
{
    cleanupStorage();
}

bool Moments::computeOutput(const Storage *st, std::ostream &out)
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
        for (i = 0; i < 3; ++i)
            second[i] = 0.0;
        for (i = 0; i < 4; ++i)
            third[i] = 0.0;
        for (i = 0; i < 5; ++i)
            fourth[i] = 0.0;

        PointVortex *pv;
        double x, y;
        // Calculate moments
        for (i = 0; i < N; ++i)
        {
            pv = st->retrieve(i);
            x = pv->getPos(0);
            y = pv->getPos(1);
            second[0] += x * x / N;
            second[1] += x * y / N;
            second[2] += y * y / N;

            third[0] += x * x * x / N;
            third[1] += x * x * y / N;
            third[2] += x * y * y / N;
            third[3] += y * y * y / N;

            fourth[0] += x * x * x * x / N;
            fourth[1] += x * x * x * y / N;
            fourth[2] += x * x * y * y / N;
            fourth[3] += x * y * y * y / N;
            fourth[4] += y * y * y * y / N;
        }

        // Print stuff
        out << "Second Moment: [" << second[0] << ", "
            << second[1] << ", " << second[2] << "]" << std::endl;
        out << "  Third Moment: [" << third[0] << ", " << third[1] << ", "
            << third[2] << ", " << third[3] << "]" << std::endl;
        out << "  Fourth Moment: [" << fourth[0] << ", " << fourth[1] << ", "
            << fourth[2] << ", " << fourth[3] << ", " << fourth[4]
            << "]" << std::endl;
    }

    return status;
}

void Moments::buildMomentStorage(unsigned dimension)
{
    dim = 2;
    second = new double[3];
    third = new double[4];
    fourth = new double[5];

    // if (dim < dimension)
    // {
    //     cleanupStorage();
    //     dim = dimension;

    //     second = new double *[dim];
    //     third = new double **[dim];
    //     fourth = new double ***[dim];
    //     for (unsigned i = 0; i < dim; ++i)
    //     {
    //         second[i] = new double[dim];
    //         third[i] = new double *[dim];
    //         fourth[i] = new double **[dim];
    //         for (unsigned j = 0; j < dim; ++j)
    //         {
    //             third[i][j] = new double[dim];
    //             fourth[i][j] = new double *[dim];
    //             for (unsigned k = 0; k < dim; ++k)
    //             {
    //                 fourth[i][j][k] = new double[dim];
    //             }
    //         }
    //     }
    // }
}

void Moments::cleanupStorage()
{
    // for (unsigned i = 0; i < dim; ++i)
    // {
    //     for (unsigned j = 0; j < dim; ++j)
    //     {
    //         for (unsigned k = 0; k < dim; ++k)
    //         {
    //             if (fourth != 0)
    //                 delete [] fourth[i][j][k];
    //         }
    //         if (fourth != 0)
    //             delete [] fourth[i][j];
    //         if (third != 0)
    //             delete [] third[i][j];
    //     }
    //     if (fourth != 0)
    //         delete [] fourth[i];
    //     if (third != 0)
    //         delete [] third[i];
    //     if (second != 0)
    //         delete [] second[i];
    // }

    delete [] fourth;
    delete [] third;
    delete [] second;
    fourth = 0;
    third = 0;
    second = 0;
}
