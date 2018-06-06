/**
 * \class RuKu4
 *
 * \brief Implementation of the standard Runge Kutta 4th order algorithm
 *
 * Author: Jonathan Maack
 **/

#include "RuKu4.h"

#include "FlowComputer.h"
#include "Logger.h"
#include "Timer.h"

RuKu4::RuKu4(Logger *ln, Timer *tm, FlowComputer *fc,
             Storage *st, double stepSize):
    Integrator(stepSize),
    log(ln),
    timer(tm),
    flowcalc(fc),
    store(st),
    positions(0),
    flow(0),
    Name("Runge-Kutta 4")
{
    num = st->size();
    dim = st->dimension();

    flow = new FlowField();
    flow->initialize(dim, num);

    positions = new Storage(*st);

    K = new double **[4];
    for (unsigned i = 0; i < 4; ++i)
    {
        K[i] = new double *[num];
        for (unsigned j = 0; j < num; ++j)
        {
            K[i][j] = new double[dim];
        }
    }

    //log->writeMsg(Logger::INFO, "RuKu4 Construction complete");
}

RuKu4::~RuKu4()
{
    unsigned i;
    unsigned j;

    delete flow;
    flow = 0;

    delete positions;
    positions = 0;

    if (K != 0)
    {
        for(i = 0; i < 4; ++i)
        {
            if (K[i] != 0)
            {
                for (j = 0; j < num; ++j)
                {
                    delete [] K[i][j];
                    K[i][j] = 0;
                }
                delete [] K[i];
                K[i] = 0;
            }
        }
        delete [] K;
        K = 0;
    }

    flowcalc = 0;
    store = 0;
    timer = 0;
    log = 0;
}

int RuKu4::step(FlowField *ff)
{
    double x0, x1;
    double *localFlow;
    PointVortex *pv0, *pv1;
    unsigned i,j;

    flowcalc->computeFlow(store, flow);

    if (ff != 0)
    {
        ff->copy(flow);
    }

    // Compute k1
    for (i = 0; i < num; ++i)
    {
        pv0 = store->retrieve(i);
        pv1 = positions->retrieve(i);
        localFlow = flow->getFlow(i);
        for (j = 0; j < dim; ++j)
        {
            x0 = pv0->getPos(j);
            K[0][i][j] = dt * localFlow[j];
            x1 = x0 + 0.5 * K[0][i][j];
            pv1->setPos(j,x1);
        }
    }

    timer->stamp(Timer::INTEGRATE);

    flowcalc->computeFlow(positions, flow);

    // Compute k2
    for (i = 0; i < num; ++i)
    {
        pv0 = store->retrieve(i);
        pv1 = positions->retrieve(i);
        localFlow = flow->getFlow(i);
        for (j = 0; j < dim; ++j)
        {
            x0 = pv0->getPos(j);
            K[1][i][j] = dt * localFlow[j];
            x1 = x0 + 0.5*K[1][i][j];
            pv1->setPos(j, x1);
        }
    }

    timer->stamp(Timer::INTEGRATE);

    flowcalc->computeFlow(positions, flow);

    // Compute k3
    for (i = 0; i < num; ++i)
    {
        pv0 = store->retrieve(i);
        pv1 = positions->retrieve(i);
        localFlow = flow->getFlow(i);
        for (j = 0; j < dim; ++j)
        {
            x0 = pv0->getPos(j);
            K[2][i][j] = dt * localFlow[j];
            x1 = x0 + K[2][i][j];
            pv1->setPos(j, x1);
        }
    }

    timer->stamp(Timer::INTEGRATE);

    flowcalc->computeFlow(positions, flow);

    // Compute k4
    for (i = 0; i < num; ++i)
    {
        pv0 = store->retrieve(i);
        pv1 = pv0; // Want to update actual position
        localFlow = flow->getFlow(i);
        for (j = 0; j < dim; ++j)
        {
            x0 = pv0->getPos(j);
            K[3][i][j] = dt * localFlow[j];
            // Compute actual position update
            x1 = x0 + 1.0/6.0 * (K[0][i][j] + 2.0 * K[1][i][j]
                                 + 2.0 * K[2][i][j] + K[3][i][j]);
            pv1->setPos(j, x1);
        }
    }

    ++stepnum;

    timer->stamp(Timer::INTEGRATE);

    return Integrator::PASS;
}
