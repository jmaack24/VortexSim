/**
 * \class Integrator
 *
 * \brief Parent class for time integrator function.
 *
 * Author: Jonathan Maack
 *
 * Description: Parent class for time integration of point vortex systems.
 * All time integrators must inherit from this class and implement the step
 * function.
 **/

#ifndef INTEGRATOR
#define INTEGRATOR

class FlowField;
class Storage;

class Integrator
{
public:
    enum{FAIL, PASS};

    /** Constructor--Intentionally unimplemented. **/
    Integrator() {}
    /** Constructor
     * @param stepSize (initial?) timestep size
     **/
    inline Integrator(double stepSize)
    {
        dt = stepSize;
        stepnum = 0;
    }
    /** Destructor **/
    virtual ~Integrator() {}

    /**
    * Takes a single time step forward
    **/
    virtual int step(FlowField *ff=0) = 0;

    /**
     * Gets the name of the integrator object
     **/
    virtual const char *name() const = 0;

    /**
     * Accessor for time step size
     * @return timestep size
     **/
    inline double stepSize() const
    {
        return dt;
    }

    inline unsigned long stepNum() const
    {
        return stepnum;
    }

protected:
    /** Time step amount **/
    double dt;

    /** Number of steps taken **/
    unsigned long stepnum;
};

#endif
