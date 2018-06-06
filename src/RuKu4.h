/**
 * \class RuKu4
 *
 * \brief Implementation of the standard Runge-Kutta 4th order algorithm
 *
 * Author: Jonathan Maack
 **/

#ifndef RUKU4
#define RUKU4

#include "Integrator.h"

#include <string>

class FlowComputer;
class FlowField;
class Logger;
class Storage;
class Timer;

class RuKu4: public Integrator
{
public:
    /** Constructor -- intentionally unimplemeted **/
    RuKu4() {}

    /** Copy Constructor -- intentionally unimplemented **/
    RuKu4(const RuKu4 &rk) {}

    /** Constructor
     * @param ln logger to use for this simulation
     * @param tm timer to use for tracking timing info
     * @param fc flow computer to use for computing right-hand side of ODE
     * @param st storage to use for storing current PV positions
     * @param stepSize size of steps integrator should use
     **/
    RuKu4(Logger *ln, Timer *tm, FlowComputer *fc,
          Storage *st, double stepSize);

    /** Destructor **/
    virtual ~RuKu4();

    virtual int step(FlowField *ff=0);

    virtual inline const char *name() const
    {
        return Name.c_str();
    }

private:
    const std::string Name;

    /** Logger for registering any messages/warnings/errors **/
    Logger *log;
    /** Flow computer used for calculating flow--the right hand side **/
    FlowComputer *flowcalc;
    /** Actual positions of vortices of simulation **/
    Storage *store;
    /** Object for calculating running time stuff **/
    Timer *timer;

    /** For storing intermediate positions used by rk4 algorithm **/
    Storage *positions;
    /** For storing intermediate flows used by rk4 algorithm **/
    FlowField *flow;

    unsigned num;
    unsigned dim;
    double ***K;
};

#endif
