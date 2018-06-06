/*******************************************************************************
 * \class SphereFlow
 *
 * \brief Computes the flow field for point vortices on the unit sphere --
 * explicitly assumes 3D
 *
 * Author: Jonathan Maack
 *
 ******************************************************************************/

#ifndef SPHERE_FLOW
#define SPHERE_FLOW

#include "FlowComputer.h"

#include <string>

class Logger;
class Timer;

class SphereFlow: public FlowComputer
{
public:
    /** Constructor -- intentionally unimplemented **/
    SphereFlow() {}
    /** Copy Constructor -- intentionally unimplemented **/
    SphereFlow(const SphereFlow &sf) {}
    /** Constructor **/
    SphereFlow(Logger *ln, Timer *tm);
    /** Destructor **/
    virtual ~SphereFlow();

    virtual void computeFlow(Storage *str, FlowField *ret);

    inline virtual const char *name() const
    {
        return Name.c_str();
    }

private:
    std::string Name;

    Logger *log;
    Timer *timer;
};

#endif
