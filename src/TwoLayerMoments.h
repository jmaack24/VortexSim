/*******************************************************************************
 * \class TwoLayerMoments
 *
 * \brief Calculate and output moments of the current PV configuration --
 * ONLY DOES DIMENSION 2!!!!
 *
 * Author: Jonathan Maack
 *
 ******************************************************************************/

#ifndef TWO_LAYER_MOMENTS
#define TWO_LAYER_MOMENTS

#include "Output.h"

#include <string>

class Logger;

class TwoLayerMoments: public Output
{
public:
    /** Default constructor -- intentionally unimplmented **/
    TwoLayerMoments() {}
    /** Constructor **/
    TwoLayerMoments(Logger *ln);
    /** Destructor **/
    virtual ~TwoLayerMoments();

    bool computeOutput(const Storage *st, std::ostream &out);

    virtual inline const char* name() const
    {
        return Name.c_str();
    }
private:
    std::string Name;

    Logger *log;

    double **second;

    unsigned dim;

    void buildMomentStorage(unsigned dimension);
    void cleanupStorage();
};

#endif
