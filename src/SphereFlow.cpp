/*******************************************************************************
 * \class SphereFlow
 *
 * \brief Computes the flow field for point vortices on the unit sphere --
 * explicitly assumes 3D
 *
 * Author: Jonathan Maack
 *
 ******************************************************************************/

#include "SphereFlow.h"

#include "Constants.h"
#include "Logger.h"
#include "Timer.h"

SphereFlow::SphereFlow(Logger *ln, Timer *tm):
    FlowComputer(),
    log(ln),
    timer(tm)
{}

SphereFlow::~SphereFlow()
{
    timer = 0;
    log = 0;
}

void SphereFlow::computeFlow(Storage *str, FlowField *ret)
{
    ret->clearFlow();

    unsigned num = str->size();

    unsigned i;
    unsigned j;
    double temp;
    double coefficient;
    double x1, x2, x3, y1, y2, y3;
    double flow1[3];
    double flow2[3];
    PointVortex *pv1;
    PointVortex *pv2;

    for (i = 0; i < num; ++i)
    {
        pv1 = str->retrieve(i);

        x1 = pv1->getPos(0);
        x2 = pv1->getPos(1);
        x3 = pv1->getPos(2);

        for (j = i + 1; j < num; ++j)
        {
            pv2 = str->retrieve(j);

            y1 = pv2->getPos(0);
            y2 = pv2->getPos(1);
            y3 = pv2->getPos(2);

            // Compute common multiplier between the two vortices
            temp = x1*y1 + x2*y2 + x3*y3; // dot product
            coefficient = 1.0/(Constants::FOUR_PI*(1 - temp));

            // Compute flow
            temp = y2*x3 - y3*x2; // cross product
            flow1[0] = pv2->getVort() * temp;
            flow2[0] = pv1->getVort() * temp;

            temp = y3*x1 - y1*x3;
            flow1[1] = pv2->getVort() * temp;
            flow2[1] = pv1->getVort() * temp;

            temp = y1*x2 - y2*x1;
            flow1[2] = pv2->getVort() * temp;
            flow2[2] = pv1->getVort() * temp;

            // Add flows
            ret->addFlow(i, flow1);
            ret->addFlow(j, flow2);
        }
    }

    timer->stamp(Timer::FLOW);
}
