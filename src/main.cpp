/**********************************************************************
 * /brief Main file wrapper for point vortex simulator class
 *
 * Author: Jonathan Maack
 *
 * Description: Loose wrapper for VortexSim class.  Allows for an
 * executable.
 *********************************************************************/

#include <iostream>

#include "VortexSim.h"

int main(int argc, char** argv)
{
    int retval = 0;
    VortexSim sim;
    bool success = sim.setup(argc, argv);
    if (!success)
    {
        std::cerr << "Setup failed." << std::endl;
        retval = -1;
    }
    else
    {
        sim.run();
        std::cout << "Done!!!!" << std::endl;
    }

    return retval;
}
