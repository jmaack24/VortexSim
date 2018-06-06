/*******************************************************************************
 * \class Timer
 *
 * \brief Time various aspects of the numerical simulation
 *
 * Author: Jonathan Maack
 *
 ******************************************************************************/

#include "Timer.h"

#include <iostream>

Timer::Timer(unsigned size):
    timestamps(0),
    capacity(size),
    max(0)
{
    // Save start time of simulation
    gettimeofday(&start, NULL);

    // Clear the totals
    for(unsigned count = 0; count < NTYPES; ++count)
    {
        totals[count] = 0.0;
    }

    // Pre allocate the needed entries
    timestamps.reserve(capacity);
    while (timestamps.size() < capacity)
    {
        timestamps.push_back(new timeInfo());
    }

    timeInfo *info = timestamps[0];
    info->type = OTHER;
    info->tm = start;

    next = 1;
}

Timer::~Timer()
{
    for(unsigned k = 0; k < timestamps.size(); ++k)
    {
        delete timestamps[k];
    }
    timestamps.clear();
}

void Timer::stamp(Timer::stype type)
{
    timeInfo *info = timestamps[next];
    info->type = type;
    gettimeofday(&(info->tm), NULL);
    ++next;

    if (next >= capacity)
        condense();
}

void Timer::condense()
{
    double runTime;
    int usec;
    unsigned sec;

    unsigned count = 1;

    timeInfo *last;
    timeInfo *current = timestamps[0];

    // std::cerr << "Timer::next = " << next << std::endl;

    while (count < next)
    {
        last = current;
        current = timestamps[count];

        sec = current->tm.tv_sec - last->tm.tv_sec;
        usec = current->tm.tv_usec - last->tm.tv_usec;

        totals[current->type] += sec;
        totals[current->type] += usec * 1e-6;;
        runTime += sec;
        runTime += usec * 1e-6;

        // if (sec > 0 || usec > 1000000)
        // {
        //     std::cerr << "sec = " << sec << std::endl
        //               << "usec = " << usec << std::endl
        //               << "current->type = " << current->type << std::endl
        //               << "current->tm.tv_sec = "
        //               << current->tm.tv_sec << std::endl
        //               << "current->tm.tv_usec = "
        //               << current->tm.tv_usec << std::endl
        //               << "last->tm.tv_sec = "
        //               << last->tm.tv_sec << std::endl
        //               << "last->tm.tv_usec = "
        //               << last->tm.tv_usec << std::endl;
        // }

        ++count;
    }

    if (next > max)
        max = next;

    // Reset to the beginning
    timestamps[0]->type = current->type;
    timestamps[0]->tm = current->tm;
    next = 1;

    // Save how long this took
    stamp(OTHER);
    stop = timestamps[1]->tm;
}

void Timer::printTimingInfo(std::ostream &out)
{
    // Calculate total time spent in different categories
    condense();

    // Print out info
    unsigned count = 0;
    char buffer[100];
    double rtime = stop.tv_sec - start.tv_sec;
    rtime += (stop.tv_usec - start.tv_usec)*1e-6;
    out << "          Time(sec)    Percentage" << std::endl;
    while (count < NTYPES)
    {
        if (count == SETUP)
            out << "SETUP     ";
        else if (count == INTEGRATE)
            out << "INTEGRATE ";
        else if (count == FLOW)
            out << "FLOW      ";
        else if (count == OUT)
            out << "OUTPUT    ";
        else if (count == OTHER)
            out << "OTHER     ";
        else
            out << "UNKNOWN   ";
        sprintf(buffer, "%8.5g%8.2f",
                totals[count], totals[count]/rtime * 100);
        out << buffer << std::endl;
        ++count;
    }
    out << "Total Runtime: " << rtime << " sec" << std::endl;
    out << "Max Count: " << max << std::endl;
}
