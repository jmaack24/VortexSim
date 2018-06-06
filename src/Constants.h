#ifndef CONSTANTS
#define CONSTANTS

class Constants
{
public:
    enum Domain {PLANE, SCREENED, TWO_LAYER, SPHERE, ROTATING_SPHERE,
                 NUM_DOMAINS};
    enum Integrate {RK4, RK45, ADMOUL, SYM, NUM_INTEG};

    static const double PI;
    static const double TWO_PI;
    static const double FOUR_PI;
};

#endif
