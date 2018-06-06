/*******************************************************************************
 * \class CenterVort
 *
 * \brief Calculate and output center of vorticity
 *
 * Author: Jonathan Maack
 *
 ******************************************************************************/

#ifndef CENTER_VORT
#define CENTER_VORT

#include "Output.h"

class Logger;

class CenterVort: public Output
{
public:
    /** Constructor -- intentionally unimplemented **/
    CenterVort() {}
    /** Constructor **/
    CenterVort(Logger *ln);
    /** Destructor **/
    virtual ~CenterVort();

    bool computeOutput(const Storage *st, std::ostream &out);

    virtual inline const char *name() const
    {
        return Name.c_str();
    }

private:
    /** Name of this output object **/
    const std::string Name;

    /** Logger for logging problems **/
    Logger *log;
    /** Size of storage buffer **/
    unsigned dim;
    /** Center of vorticity storage -- stored so not always creating
        and deleting memory
    **/
    double *center;

    /** Helper function that setups center of vorticity storage buffer **/
    void buildCenterSpace (unsigned dimension);
};

#endif
