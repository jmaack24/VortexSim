/**
 * \class AdMoulPC
 *
 * \brief Implementation of the adams moulton predictor corrector algorithm
 *
 * Author: Jonathan Maack
 **/

#ifndef AD_MOUL_PC
#define AD_MOUL_PC

#include "Integrator.h"

#include <string>
#include <vector>

class FlowComputer;
class FlowField;
class Logger;
class RuKu4;
class Storage;
class Timer;

class AdMoulPC: public Integrator
{
public:
    /** Default Constructor -- intentionally unimplemented **/
    AdMoulPC() {}
    /** Copy Constructor -- intentionally unimplemented **/
    AdMoulPC(const AdMoulPC &ampc) {}
    /** Constructor
     * @param ln logger to use for this simulation
     * @param tm timer to use for tracking timing info
     * @param fc flow computer to use for computing right-hand side of ODE
     * @param st storage to use for storing current PV positions
     * @param stepSize size of steps integrator should use
     **/
    AdMoulPC(Logger *ln, Timer *tm, FlowComputer *fc, Storage *st,
             double stepSize);
    /** Destructor **/
    virtual ~AdMoulPC();

    int step(FlowField *ff = 0);

    virtual inline const char *name() const
    {
        return Name.c_str();
    }
private:
    const std::string Name;

    /** Logger for registering any messages/warnings/errors **/
    Logger *log;
    /** Flow computer used for calculating flow--the right hand side **/
    FlowComputer *flowcalc;
    /** Actual positions of vortices of simulation **/
    Storage *store;
    /** Object for calculating running time stuff **/
    Timer *timer;

    /** Flow field storage **/
    std::vector<FlowField *> flows;
    /** Local position storage for iterative process **/
    Storage *lstore;

    unsigned num;
    unsigned dim;

    bool initialized;
    RuKu4 *initStepper;
};

#endif
