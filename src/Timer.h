/*******************************************************************************
 * \class Timer
 *
 * \brief Time various aspects of the numerical simulation
 *
 * Author: Jonathan Maack
 *
 ******************************************************************************/

#ifndef TIMER
#define TIMER

#include <iostream>
#include <sys/time.h>
#include <vector>

class Timer
{
public:
    enum stype {SETUP, INTEGRATE, FLOW, OUT, OTHER, NTYPES};

    /** Constructor -- intentionally unimplemented **/
    Timer() {}
    /** Copy Constructor -- intentionally unimplemented **/
    Timer(const Timer &tm) {}
    /** Constructor **/
    Timer(unsigned size);
    /** Destructor **/
    ~Timer();

    /** Record timing information **/
    void stamp(stype type);

    /**
     * Compile info without printing it. Used to keep memory usage
     * from growing too large
     **/
    void condense();

    /** Retrieve timing information **/
    void printTimingInfo(std::ostream &out);

private:
    struct timeInfo
    {
        stype type;
        timeval tm;
    };

    std::vector<timeInfo *> timestamps;
    unsigned next;
    unsigned capacity;
    unsigned max;

    double totals[NTYPES];
    timeval start;
    timeval stop;
};

#endif
