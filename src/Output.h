/*******************************************************************************
 * \class Output
 *
 * \brief Defines the interface for output calculations and objects
 *
 * Author: Jonathan Maack
 *
 ******************************************************************************/

#ifndef OUTPUT
#define OUTPUT

#include <iostream>

class Storage;

class Output
{
public:
    /** Constructor **/
    Output(){}
    /** Destructor **/
    virtual ~Output(){}

    /**
     * Compute and output whatever this thing does
     * @param st collection of vortices to compute stuff for
     * @param out output stream to use for output
     * @return true if successful
     **/
    virtual bool computeOutput(const Storage *st, std::ostream &out) = 0;

    /**
     * Get the name of the output object--mostly useful for printing things
     * @return char string of the name
     **/
    virtual inline const char *name() const = 0;

private:
};

#endif
