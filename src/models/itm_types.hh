#ifndef ITM_TYPES_HH
#define ITM_TYPES_HH

#include <array>

struct prop_type {
    double aref;
    double dist; //effective earth curvature
    std::array<double, 2> hg; // height off ground of TX, RX in meters
    std::array<double, 2> rch;
    double wn;
    double dh;
    double dhd;
    double ens;
    double encc;
    double cch;
    double cd;
    double gme;
    double zgndreal;
    double zgndimag;
    std::array<double, 2> he;
    std::array<double, 2> dl; // horizon distances
    std::array<double, 2> the; // horizon takeoff angle
    double tiw;
    double ght;
    double ghr;
    double rph;
    double hht;
    double hhr;
    double tgh;
    double tsgh;
    double thera;
    double thenr;
    int rpl;
    int kwx;
    int mdp;
    int ptx; // transmit polarity: 0 = h, 1 = v, 2 = circular (TODO: should this be an enum?)
    int los;
};

struct propv_type {
    double sgc;
    int lvar;
    int mdvar;
    int klim;
};

struct propa_type {
    double dlsa;
    double dx;
    double ael;
    double ak1;
    double ak2;
    double aed;
    double emd;
    double aes;
    double ems;
    std::array<double, 2> dls;
    double dla;
    double tha;
};

#endif
