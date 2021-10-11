#include "Particle.h"

Particle::Particle(){
  px = py = pz = E = m = pAbs = pt = 0.0;
  p[0] = p[1] = p[2] = p[3] = 0.0;

}

Particle::Particle(double mass){
  m = mass;
  px = py = pz = E = pAbs = pt = 0.0;
  p[0] = p[1] = p[2] = p[3] = 0.0;
}

/* Delete //////////////////////////////////////////////////////////////////////// */
Particle::Particle(double tempx, double tempy, double tempz, double tempmass)
{
   px = tempx;
   py = tempy;
   pz = tempz;
   m  = tempmass;

   E  = sqrt(px*px + py*py + pz*pz + m*m);

   p[0] = px; p[1] = py; p[2] = pz; p[3] = E;
   pt      = sqrt(tempx*tempx + tempy*tempy);
   pAbs    = sqrt(tempx*tempx + tempy*tempy + tempz*tempz);
   SetPxPyPzE(px,py,pz,E);
}
void Particle::p4(double momx, double momy, double momz, double energy){
  px = p[0] = momx;
  py = p[1] = momy;
  pz = p[2] = momz;
  E  = p[3] = energy;
  pt      = sqrt(momx*momx + momy*momy);
  pAbs    = sqrt(momx*momx + momy*momy + momz*momz);
}
/*////////////////////////////////////////////////////////////////////////////// */


void Particle::setPxPyPzE(double tempx, double tempy, double tempz, double tempE)
{

   px = tempx;
   py = tempy;
   pz = tempz;
   E  = tempE;

   m  =  sqrt(E*E - (px*px + py*py + pz*pz) );

   p[0] = px; p[1] = py; p[2] = pz; p[3] = E;
   pt      = sqrt(px*px + py*py);
   pAbs    = sqrt(px*px + py*py + pz*pz);

}

void Particle::setPxPyPzM(double tempx, double tempy, double tempz, double tempmass)
{
   px = tempx;
   py = tempy;
   pz = tempz;
   m  = tempmass;

   E  = sqrt(px*px + py*py + pz*pz + m*m);

   p[0] = px; p[1] = py; p[2] = pz; p[3] = E;
   pt      = sqrt(tempx*tempx + tempy*tempy);
   pAbs    = sqrt(tempx*tempx + tempy*tempy + tempz*tempz);
}

void Particle::Print(){
    printf("\n %lf %lf %lf %lf %lf", px, py, pz, E, m);
}


void Particle::boost(Particle parent){

  double betax = ((-1)*parent.px) / parent.E;
  double betay = ((-1)*parent.py) / parent.E;
  double betaz = ((-1)*parent.pz) / parent.E;
  double beta2 = betax*betax + betay*betay + betaz*betaz;
  double gamma = 1.0/sqrt(1.0-beta2);
  double dot   = betax*px + betay*py + betaz*pz;
  double prod  = gamma*( gamma*dot/(1.0+gamma) + E );

  double pX = px + betax*prod;
  double pY = py + betay*prod;
  double pZ = pz + betaz*prod;
  double e  = gamma*(E + dot);

  p4( pX, pY, pZ, e );
}





