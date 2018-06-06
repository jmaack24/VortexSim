/*******************************************************************************
 * \class Deformation Flow
 *
 * \brief Computes the screened flow field for point vortices in the unbounded
 * plane with.
 *
 * Author: Jonathan Maack
 *
 * Description: Computes the flow field on the given point vortices in
 * the R^2 plane with no boundaries. The interactions are screened.
 ******************************************************************************/

#ifndef SCREENED_FLOW
#define SCREENED_FLOW

#include "FlowComputer.h"

#include "FlowField.h"
#include "Storage.h"

class Logger;
class Timer;

class ScreenedFlow: public FlowComputer
{
public:
    /** Constructor -- intentionally unimplemented **/
    ScreenedFlow() {}
    /** Copy Constructor -- intentionally unimplemented **/
    ScreenedFlow(const ScreenedFlow &sf) {}
    /** Constructor **/
    ScreenedFlow(Logger *natural, Timer *tm, double defrad);
    /** Destructor **/
    virtual ~ScreenedFlow();

    virtual void computeFlow(Storage *str, FlowField *ret);

    virtual inline const char *name() const
    {
        return Name.c_str();
    }

private:
    const std::string Name;

    double deformationRadius;
    double deformationFreq;

    Logger *log;
    Timer *timer;
};

#endif
