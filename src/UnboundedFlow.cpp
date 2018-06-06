/*******************************************************************************
 * \class UnboundedFlow
 *
 * \brief Computes the flow field for point vortices in an unbounded plane.
 *
 * Author: Jonathan Maack
 *
 * Description: Computes the flow field on the given point vortices for in
 * the R^2 plane with no boundaries.
 ******************************************************************************/

#include "UnboundedFlow.h"

#include "Constants.h"
#include "FlowField.h"
#include "Logger.h"
#include "PointVortex.h"
#include "Storage.h"
#include "Timer.h"

#include <sstream>

UnboundedFlow::UnboundedFlow(Logger *natural, Timer *tm) :
    FlowComputer(),
    log(natural),
    timer(tm),
    Name("Unbounded 2D Flow")
{
    //log->writeMsg(Logger::INFO, "UnboundedFlow construction complete");
}

UnboundedFlow::~UnboundedFlow()
{}

void UnboundedFlow::computeFlow(Storage *store, FlowField *ret)
{
    unsigned numVorts = store->size();
    unsigned numSpots = ret->getNum();
    if (numVorts != numSpots)
    {
        std::stringstream ss;
        ss << "ERR: Have " << numVorts << " vortices and " << numSpots
           << " spots to store flow.";
        log->writeMsg(Logger::ERROR, ss.str());
    }
    else
    {
        PointVortex *pv1;
        PointVortex *pv2;
        double dx, dy, rsq;
        double coeff;
        double u1[2];
        double u2[2];

        // clear flow storage
        ret->clearFlow();

        for(unsigned i = 0; i < numVorts; ++i)
        {
            pv1 = store->retrieve(i);
            for (unsigned j=i+1; j < numVorts; ++j)
            {
                pv2 = store->retrieve(j);
                dx = pv1->getPos(0) - pv2->getPos(0);
                dy = pv1->getPos(1) - pv2->getPos(1);
                rsq = dx*dx + dy*dy;
                coeff = 1.0/(Constants::TWO_PI * rsq);
                u1[0] = -1.0*coeff * pv2->getVort()*dy;
                u1[1] = coeff * pv2->getVort()*dx;
                u2[0] = coeff * pv1->getVort()*dy;
                u2[1] = -1.0*coeff * pv1->getVort()*dx;
                ret->addFlow(i, u1);
                ret->addFlow(j, u2);
            }
        }
    }

    timer->stamp(Timer::FLOW);
}
