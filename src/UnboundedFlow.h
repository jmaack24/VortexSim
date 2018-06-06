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

#ifndef UNBOUNDED_FLOW
#define UNBOUNDED_FLOW

#include "FlowComputer.h"

#include "FlowField.h"
#include "Storage.h"

class Logger;
class Timer;

class UnboundedFlow: public FlowComputer
{
public:
    /** Constructor -- intentionally unimplemented **/
    UnboundedFlow() {}
    /** Copy Constructor -- intentionally unimplemented **/
    UnboundedFlow(const UnboundedFlow &uf) {}
    /** Constructor **/
    UnboundedFlow(Logger *natural, Timer *tm);
    /** Destructor **/
    virtual ~UnboundedFlow();

    virtual void computeFlow(Storage *str, FlowField *ret);

    virtual inline const char *name() const
    {
        return Name.c_str();
    }

private:
    const std::string Name;

    Logger *log;
    Timer *timer;
};

#endif
