/*******************************************************************************
 * \class PointVortex
 *
 * \brief Stores information about a point vortex
 *
 * Author: Jonathan Maack
 *
 * Storage class for the point vortex simulation which stores all info needed
 * by the simulation.
 ******************************************************************************/

#ifndef POINTVORTEX
#define POINTVORTEX

#include <iostream>
#include <string.h>

class PointVortex
{
public:
    /** Constructor -- intentionally unimplemented **/
    PointVortex() {}
    /** Constructor **/
    PointVortex(unsigned dim);
    /** Copy Constructor **/
    PointVortex(const PointVortex &pv);
    /** Destructor **/
    virtual ~PointVortex();

    inline void copy(const PointVortex *cp)
    {
        this->dimension = cp->dimension;
        memcpy(this->position, cp->position, sizeof(double)*dimension);
        this->vorticity = cp->vorticity;
    }

    inline void copyPos(const PointVortex *cp)
    {
        memcpy(this->position, cp->position, sizeof(double)*dimension);
    }

    /**
     * Accessor for layer
     * @return gives layer of the point vortex
     **/
    inline int getLayer() const
    {
        return layer;
    }

    /**
     * Accessor for layer
     * @return gives layer of the point vortex
     **/
    inline void setLayer(int lay)
    {
        this->layer = lay;
    }

    /**
     * Accessor for position
     * @return gives pointer to this vortex's position array
     **/
    inline const double* const getPos() const
    {
        return position;
    }

    /**
     * Accessor for position
     * @param dim dimension of coordinate to get
     **/
    inline double getPos(unsigned dim) const
    {
        return position[dim];
    }
    /**
     * Accessor for the vorticity
     * @return vorticity of the point vortex
     **/
    inline double getVort() const
    {
        return vorticity;
    }

    /**
     * Accessor for updating position
     **/
    inline void setPos(unsigned dim, double coord)
    {
        this->position[dim] = coord;
    }
    /**
     * Accessor for updating the vorticity
     * @param vort vorticity
     **/
    inline void setVort(double vort)
    {
        this->vorticity = vort;
    }

    friend std::ostream& operator<<(std::ostream &out, const PointVortex &pv);

private:
    /** Number of coordinates **/
    unsigned dimension;
    /** Layer number -- not always used **/
    int layer;
    /** Store position **/
    double *position;
    /** Stores vorticity **/
    double vorticity;
};

#endif
