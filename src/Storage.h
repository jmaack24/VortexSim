/**
 * \class Storage
 *
 * \brief Stores and manages access to point vortices.
 *
 * Author: Jonathan Maack
 *
 * Description: Responsible for storing the point vortex objects of the
 * simulation.  Also controls access to said vortices.
 **/

#ifndef STORAGE
#define STORAGE

#include "Constants.h"
#include "Logger.h"
#include "PointVortex.h"

#include <vector>

class Storage
{
public:
    /** Constructor--Intentionally Unimplemented **/
    Storage(){}
    /** Copy Constructor **/
    Storage(const Storage &store);
    /** Constructor
     * @param ln logger to use for this storage class
     * @param N number of point vortices to create
     * @param d number of dimensions in simulation
     **/
    Storage(Logger *ln, unsigned N, unsigned d);
    /** Destructor **/
    virtual ~Storage();

    /**
     * Initialize all point vortex information from given file
     **/
    bool initializeFromFile(bool warmstart, const std::string &file,
                            Constants::Domain dom);

    /**
     * Initialize point vortices to random locations
     **/
    bool initializeRandom(unsigned seed);

    /**
     * Copy position information in *cp for all point vortices.  Makes no
     * attempt to verify that the two Storage object have same number of
     * point vortices in storage.
     **/
    void copyPVPos(const Storage *cp);

    /**
     * Gets number of vortices stored
     * @return number of vortices stored
     **/
    inline unsigned size() const
    {
        return num;
    }

    /**
     * Gets the dimension of the simulation this storage is set up for
     * @return dimension of simulation
     **/
    inline unsigned dimension() const
    {
        return dim;
    }

    /**
     * Get the total vorticity of stored point vortices
     * @return total vorticity of stored point vortices
     **/
    inline double totalvort() const
    {
        return totalVorticity;
    }

    /**
     * Retrieves the vortex with 'index' ID number
     * @param index vortex number to retrieve
     * @return pointer to the vortex, null if index is not valid
     **/
    PointVortex *retrieve(unsigned index) const;

    /** Define output operator **/
    friend std::ostream& operator<<(std::ostream &out, const Storage& st);

private:
    /** Pointer to the logger to use to report errors **/
    Logger *log;

    /** Stores all the point vortices **/
    std::vector<PointVortex*> vortices;

    /** Number of point vortices stored **/
    unsigned num;
    /** Dimension of the simulation **/
    unsigned dim;
    /** Total vorticity of the stored point vortices **/
    double totalVorticity;

    void clear();

    bool initializeFromInputFile(std::ifstream &in, Constants::Domain dom);
    bool initializeFromPositionFile(std::ifstream &in, Constants::Domain dom);
};

#endif
