/**
 * \class Logger
 * \brief Class to handle all the logging and messaging needs.
 *
 * The Logger class is meant to handle all the logging and outputing
 * needs for the vortex simulation.  The writeMsg function accepts a code
 * signifying whether the message is just info or something more serious
 * like a warning or error.
 **/

#ifndef LOGGER
#define LOGGER

#include <fstream>
#include <iostream>
#include <string>

class Logger
{
public:
    enum{INFO, WARN, ERROR};
    enum{FILE, SCREEN, BOTH};

    /**
     * Constructor
     **/
    Logger();
    /**
     * Destructor
     **/
    virtual ~Logger();

    /**
     * Set the precision of the output log
     **/
    void setPrecision(unsigned prec);

    /**
     * Writes messages to the log file (or cerr if no file)
     * @param code Code which signifies message type--currently unused
     * @param msg The message to output
     * @return true when printing to a log file
     **/
    bool writeMsg(int code, const std::string &msg);

    /**
     * Opens the given file and uses it as the log file
     * @param file The name of the log file
     * @return False if log file setup fails
     **/
    bool setupLogFile(const std::string &file);

private:
    /** the output stream used for writing messages **/
    std::ostream *out;
    /** file output stream for the log file--non null only if file is used  **/
    std::ofstream *log;
};

#endif
