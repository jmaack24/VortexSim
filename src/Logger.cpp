
#include "Logger.h"

#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>

Logger::Logger()
{
    out = &std::cout;
    log = 0;
}

Logger::~Logger()
{
    if (log != 0)
    {
        log->close();
        delete log;
        log = 0;
    }

    out = 0;
}

void Logger::setPrecision(unsigned prec)
{
    out->precision(prec);
}

bool Logger::writeMsg(int code, const std::string &msg)
{
    if (code == INFO)
        (*out) << msg << std::endl;
    else if (code == WARN)
        (*out) << "****WARN: " << msg << std::endl;
    else
        (*out) << "****ERROR: " << msg << std::endl;

    return log != 0;
}

bool Logger::setupLogFile(const std::string &file)
{
    std::stringstream ss;
    if (file.empty())
        ss << "vsim";
    else
        ss << file;

    time_t currentTime;
    tm *localTime;
    time(&currentTime);
    localTime = localtime(&currentTime);
    ss << "_" << localTime->tm_mon + 1 << "_" << localTime->tm_mday
       << "_" << localTime->tm_hour << localTime->tm_min << localTime->tm_sec
       << ".log";

    // Append time stamp to log?
    log = new std::ofstream(ss.str().c_str());
    out = log;
    return out->good();
}
