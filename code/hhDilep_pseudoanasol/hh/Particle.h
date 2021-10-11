#ifndef PARTICLE_H
#define PARTICLE_H

#include <stdio.h>
#include <math.h>

#include "TLorentzVector.h"



//class Particle{
class Particle : public TLorentzVector {
  public:

    Particle();
    Particle(double);
    Particle(double, double, double, double);


    double   px, py, pz, E, m, pAbs, pt, p[4];

    void     p4(double, double, double, double);
    void     setPxPyPzM(double, double, double, double);
    void     setPxPyPzE(double, double, double, double);

    void     boost(Particle);
    void     twoBodyDecay(Particle&, Particle&);
    void     Print();
};

#endif
