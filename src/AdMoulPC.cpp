/**
 * \class AdMoulPC
 *
 * \brief Implementation of the adams moulton predictor corrector algorithm
 *
 * Author: Jonathan Maack
 **/

#include "AdMoulPC.h"

#include "Logger.h"
#include "FlowComputer.h"
#include "FlowField.h"
#include "RuKu4.h"
#include "Storage.h"
#include "Timer.h"

#include <cmath>

AdMoulPC::AdMoulPC(Logger *ln, Timer *tm, FlowComputer *fc, Storage *st,
                   double stepSize):
    Integrator(stepSize), log(ln), timer(tm), flowcalc(fc),
    store(st), initialized(false),
    Name("Adams-Moulton Predictor-Corrector")
{
    num = st->size();
    dim = st->dimension();

    lstore = new Storage(*st);

    flows.clear();
    for (unsigned i = 0; i < 4; ++i)
    {
        flows.push_back(new FlowField);
        flows[i]->initialize(dim, num);
    }

    initStepper = new RuKu4(ln, tm, fc, st, stepSize);
}

AdMoulPC::~AdMoulPC()
{
    for(unsigned i = 0; i < flows.size(); ++i)
    {
        delete flows[i];
    }
    flows.clear();

    flowcalc = 0;
    store = 0;
    timer = 0;
    log = 0;
}

int AdMoulPC::step(FlowField *ff)
{
    if (!initialized)
    {
        unsigned sn = initStepper->stepNum();
        initStepper->step(flows[3 - sn]);
        if (sn == 3)
        {
            initialized = true;
            delete initStepper;
            initStepper = 0;
        }
    }
    else
    {
        // Adams-Bashforth predictor step
        unsigned i, j;
        double pos, err;
        PointVortex *pv;

        for (i = 0; i < num; ++i)
        {
            pv = lstore->retrieve(i);
            for (j = 0; j < dim; ++j)
            {
                pos = pv->getPos(j);
                pos += dt/24.0 * (
                    55.0*flows[0]->getFlow(i,j) - 59.0*flows[1]->getFlow(i,j)
                    + 37.0*flows[2]->getFlow(i,j) - 9.0*flows[3]->getFlow(i,j));
                pv->setPos(j, pos);
            }
        }

        // Update storage locations
        FlowField *temp = flows[3];
        flows[3] = flows[2];
        flows[2] = flows[1];
        flows[1] = flows[0];
        flows[0] = temp;

        // Adams-Moulton corrector step
        bool go = true;
        while (go)
        {
            go = false;
            timer->stamp(Timer::INTEGRATE);
            // Compute current flow
            flowcalc->computeFlow(lstore, flows[0]);
            for (i = 0; i < num; ++i)
            {
                pv = lstore->retrieve(i);
                for (j = 0; j < dim; ++j)
                {
                    pos = store->retrieve(i)->getPos(j);
                    pos += dt/24.0 * (9.0*flows[0]->getFlow(i,j)
                                      + 19.0*flows[1]->getFlow(i,j)
                                      - 5.0*flows[2]->getFlow(i,j)
                                      + 1.0*flows[3]->getFlow(i,j));
                    err = pv->getPos(j);
                    pv->setPos(j, pos);
                    err -= pos;
                    go = (go || fabs(err) > 1e-8);
                }
            }
        }

        store->copyPVPos(lstore);
    }
    // Update PV positions for the rest of the program
    ++stepnum;
    timer->stamp(Timer::INTEGRATE);

    return Integrator::PASS;
}
