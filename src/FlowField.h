/*******************************************************************************
 * \class FlowField
 *
 * \brief Object containing the flow field at the location of the vortices
 *
 * Author: Jonathan Maack
 *
 * Description: Contains the flow field at each vortex location.  The
 * FlowComputer class returns this object from the computeFlow method.
 * The order of storage of the flow matches the order of the vortices in the
 * Storage object used to compute the field.
 ******************************************************************************/

#ifndef FLOW_FIELD
#define FLOW_FIELD

#include <iostream>

class FlowField
{
public:
    /** Constructor **/
    FlowField();
    /** Copy Constructor -- intentionally unimplemented **/
    FlowField(const FlowField &ff) {}
    /** Destructor **/
    virtual ~FlowField();

    /**
     * Initialize object for storage
     * @param d number of dimensions of flow stored
     * @param n number of voritces
     * @return True if successful **/
    bool initialize(unsigned d, unsigned n);

    /** Get number of dimensions **/
    inline unsigned getDim() const
    {
        return dim;
    }
    /** Get number of vortices **/
    inline unsigned getNum() const
    {
        return num;
    }
    /** Get flow on 'k'th vortex **/
    inline double *const getFlow(unsigned k) const
    {
        return flow[k];
    }
    /** Get flow on 'k'th vortex in 'l'th dimension **/
    inline double getFlow(unsigned k, unsigned l) const
    {
        return flow[k][l];
    }
    /** Get entire flow vector **/
    inline double *const*const getFlow() const
    {
        return flow;
    }

    /**
     * Set flow on 'k'th vortex.  It is assumed that 'k' < 'num'.  No effort
     * is made to ensure this.  The input flow vector is assumed to be of
     * length 'dim'.  It is up to the user to ensure this.  Note the contents
     * of 'u' are copied into the flow vector.
     **/
    void setFlow(unsigned k, double *u);

    /**
     * Set all flow values to zero
     **/
    void clearFlow();

    /**
     * Add to flow on 'k'th vortex.  It is assumed that 'k' < 'num'.  No effort
     * is made to ensure this.  The input flow vector is assumed to be of
     * length 'dim'.  It is up to the user to ensure this.
     **/
    void addFlow(unsigned k, double *u);

    /** Copy the flow information from the given FlowField **/
    void copy(const FlowField *cp);

    friend std::ostream & operator<<(std::ostream &out, const FlowField &ff);

private:
    /** dimension of flow field **/
    unsigned dim;
    /** number of vortices **/
    unsigned num;
    /** flow field **/
    double **flow;
};

#endif
