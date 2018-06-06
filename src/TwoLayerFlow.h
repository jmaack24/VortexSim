/*******************************************************************************
 * \class TwoLayerFlow
 *
 * \brief Computes the flow field for point vortices in the two layer
 * quasigeostrophic model
 *
 * Author: Jonathan Maack
 *
 * Description: Computes the flow field on the given point vortices for in
 * the two layer quasigeostrophic twolayer model
 ******************************************************************************/

#ifndef TWO_LAYER_FLOW
#define TWO_LAYER_FLOW

#include "FlowComputer.h"

#include "FlowField.h"
#include "Storage.h"

class Logger;
class Timer;

class TwoLayerFlow: public FlowComputer
{
public:
    /** Constructor -- intentionally unimplemented **/
    TwoLayerFlow();
    /** Copy Constructor -- intentionally unimplemented **/
    TwoLayerFlow(const TwoLayerFlow &tlf);
    /** Constructor **/
    TwoLayerFlow(Logger *natural, Timer *tm, double defrad);
    /** Destructor **/
    virtual ~TwoLayerFlow();

    virtual void computeFlow(Storage *storage, FlowField *ret);

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
