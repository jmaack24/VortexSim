/*******************************************************************************
 * \class FlowComputer
 *
 * \brief Computes the flow field for the given configuration of point vorices.
 *
 * Author: Jonathan Maack
 *
 * Description: Parent class for all different flow field computations.
 * Any child class must implement the computeFlow method.
 ******************************************************************************/

#ifndef FLOW_COMPUTER
#define FLOW_COMPUTER

#include "FlowField.h"
#include "Storage.h"

class FlowComputer
{
public:
    /** Constructor **/
    FlowComputer() {}
    /** Copy Constructor **/
    FlowComputer(const FlowComputer &fc) {}
    /** Destructor **/
    virtual ~FlowComputer() {}

    /**
     * Computes flow at the point vortices in 'str' due to the vortices
     * in 'str' and returns them in the FlowField object
     * @param str object containing point vortices to use for computation
     * @param ret return object which contains the flow
     **/
    virtual void computeFlow(Storage *str, FlowField *ret) = 0;

    /**
     * Provide the name of this flow computer
     * @return objects name
     **/
    virtual const char *name() const = 0;

protected:
private:
};

#endif
