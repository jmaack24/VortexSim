/*******************************************************************************
 * \class PointVortex
 *
 * \brief Stores information about a point vortex
 *
 * Author: Jonathan Maack
 *
 * Storage class for the point vortex simulation which stores all info needed
 * by the simulation.
 ******************************************************************************/

#include "PointVortex.h"

#include <cstring>
#include <iostream>

PointVortex::PointVortex(unsigned dim):
    position(0), vorticity(0.0), dimension(dim), layer(0)
{
    position = new double[dimension];
}

PointVortex::PointVortex(const PointVortex &pv)
{
    this->layer = pv.layer;
    this->dimension = pv.dimension;
    this->vorticity = pv.vorticity;
    this->position = new double[this->dimension];
    std::memcpy(this->position, pv.position, this->dimension * sizeof(double));
}

PointVortex::~PointVortex()
{
    if (position != 0)
    {
        delete [] position;
        position = 0;
    }
}

std::ostream& operator<<(std::ostream &out, const PointVortex &pv)
{
    out << "Layer=" << pv.layer << " Vort=" << pv.vorticity << " Pos=(";
    for (unsigned k = 0; k < pv.dimension; ++k)
    {
        out << pv.position[k];
        if (k < pv.dimension - 1)
            out << ",";
    }
    out << ")";
    return out;
}
