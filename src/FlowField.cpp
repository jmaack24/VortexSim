/*******************************************************************************
 * \class FlowField
 *
 * \brief Object containing the flow field at the location of the vortices
 *
 * Author: Jonathan Maack
 *
 * Description: Contains the flow field at each vortex location.  The
 * FlowComputer class returns this object from the computeFlow method.
 * The order of storage of the flow matches the order of the vortices in the
 * Storage object used to compute the field.
 ******************************************************************************/

#include "FlowField.h"

#include <cstring>
#include <iostream>

FlowField::FlowField():
    dim(0), num(0), flow(0)
{}

FlowField::~FlowField()
{
    if (flow != 0)
    {
        for(unsigned k = 0; k < num; ++k)
        {
            if (flow[k] != 0)
            {
                delete [] flow[k];
                flow[k] = 0;
            }
        }
        delete flow;
        flow = 0;
    }
}

bool FlowField::initialize(unsigned d, unsigned n)
{
    dim = d;
    num = n;
    flow = new double*[num];
    for (unsigned k = 0; k < num; ++k)
    {
        flow[k] = new double[dim];
    }

    return true;
}

void FlowField::setFlow(unsigned k, double *u)
{
    std::memcpy(flow[k], u, dim*sizeof(double));
}

void FlowField::clearFlow()
{
    for(unsigned i = 0; i < num; ++i)
    {
        for(unsigned j = 0; j < dim; ++j)
        {
            flow[i][j] = 0.0;
        }
    }
}

void FlowField::addFlow(unsigned k, double *u)
{
    for(unsigned l = 0; l < dim; ++l)
        flow[k][l] += u[l];
}

void FlowField::copy(const FlowField *cp)
{
    for (unsigned i = 0; i < num; ++i)
    {
        memcpy(flow[i], cp->flow[i], sizeof(double)*dim);
        // for (unsigned j = 0; j < dim; ++j)
        // {
        //     flow[i][j] = cp->flow[i][j]
        // }
    }
}

std::ostream& operator<<(std::ostream &out, const FlowField &ff)
{
    for (unsigned k = 0; k < ff.num; ++k)
    {
        out << "(";
        for (unsigned l = 0; l < ff.dim; ++l)
        {
            out << ff.flow[k][l];
            if (l < ff.dim - 1)
                out << ", ";
        }
        out << ")" << std::endl;
    }
    return out;
}
