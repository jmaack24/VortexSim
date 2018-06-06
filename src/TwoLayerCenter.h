/*******************************************************************************
 * \class TwoLayerCenter
 *
 * \brief Calculate and output center of vorticity of the current PV
 * configuration -- ONLY DOES DIMENSION 2!!!!
 *
 * Author: Jonathan Maack
 *
 ******************************************************************************/

#ifndef TWO_LAYER_CENTER
#define TWO_LAYER_CENTER

#include "Output.h"

#include <string>

class Logger;

class TwoLayerCenter: public Output
{
public:
    /** Default constructor -- intentionally unimplmented **/
    TwoLayerCenter() {}
    /** Constructor **/
    TwoLayerCenter(Logger *ln);
    /** Destructor **/
    virtual ~TwoLayerCenter();

    bool computeOutput(const Storage *st, std::ostream &out);

    virtual inline const char* name() const
    {
        return Name.c_str();
    }
private:
    std::string Name;

    Logger *log;

    double **center;

    unsigned dim;

    void buildCenterStorage(unsigned dimension);
    void cleanupStorage();
};

#endif
