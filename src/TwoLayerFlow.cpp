
#include "TwoLayerFlow.h"

#include "Constants.h"
#include "FlowField.h"
#include "Logger.h"
#include "PointVortex.h"
#include "Storage.h"
#include "Timer.h"

#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include <sstream>

TwoLayerFlow::TwoLayerFlow(Logger *natural, Timer *tm, double defrad):
    FlowComputer(),
    log(natural),
    timer(tm),
    Name("Two Layer Flow"),
    deformationRadius(defrad)
{
    deformationFreq = sqrt(2.0) / deformationRadius;
}

TwoLayerFlow::~TwoLayerFlow()
{}

void TwoLayerFlow::computeFlow(Storage *store, FlowField *ret)
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
        double dx, dy, r, rsq;
        double coeff;
        double u1[2];
        double u2[2];
        int sign;

        // clear flow storage
        ret->clearFlow();

        for(unsigned i = 0; i < numVorts; ++i)
        {
            pv1 = store->retrieve(i);
            for (unsigned j=i+1; j < numVorts; ++j)
            {
                pv2 = store->retrieve(j);
                sign = pv1->getLayer() * pv2->getLayer();
                dx = pv1->getPos(0) - pv2->getPos(0);
                dy = pv1->getPos(1) - pv2->getPos(1);
                rsq = dx*dx + dy*dy;
                r = sqrt(rsq);
                coeff = 0.5 / (Constants::TWO_PI * rsq);
                coeff *= (1.0 + sign * deformationFreq * r
                          * gsl_sf_bessel_K1(deformationFreq * r));
                u1[0] = -1.0 * coeff * pv2->getVort() * dy;
                u1[1] = coeff * pv2->getVort() * dx;
                u2[0] = coeff * pv1->getVort() * dy;
                u2[1] = -1.0 * coeff * pv1->getVort() * dx;
                ret->addFlow(i, u1);
                ret->addFlow(j, u2);
            }
        }
    }

    timer->stamp(Timer::FLOW);
}
