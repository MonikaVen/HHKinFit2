#include "HHFitObjectEConstM.h"
#include <cassert>
#include <iostream>
#include <sstream>
#include <math.h>
#include "exceptions/HHEnergyRangeException.h"
#include "exceptions/HHEnergyConstraintException.h"

HHKinFit2::HHFitObjectEConstM::HHFitObjectEConstM(HHLorentzVector const& initial4vector)
  :HHFitObjectE(initial4vector){

}


HHKinFit2::HHLorentzVector
HHKinFit2::HHFitObjectEConstM::constrainEtoMinv(double m, HHLorentzVector const& pset) const{
  HHLorentzVector pmod = getInitial4Vector();
  HHLorentzVector combined = pset + pmod;

  double Mc = m;
  double M1c = pset.M();
  double M2c = pmod.M();
//  double M = combined.M();

  double C = 0.5*(pow(Mc,2)-pow(M1c,2)-pow(M2c,2));

//  int loopCount = 0;
//  while(fabs(M-Mc) > 0.000001){
//    loopCount++;
    double P1x = pset.Px();
    double P1y = pset.Py();
    double P1z = pset.Pz();
    double P1  = pset.P();
    double E1  = pset.E();

    double P2x = pmod.Px();
    double P2y = pmod.Py();
    double P2z = pmod.Pz();
    double P2  = pmod.P();

    double cosa = (P1x*P2x+P1y*P2y+P1z*P2z)/(P1*P2);
    double E2new = -1;

    if(cosa==0){
      E2new=C/E1;
    }
    else{
      double cp=C/(cosa*P1);
      double dp=E1/(cosa*P1);
      double a=pow(dp,2)-1;
      double b=-2*dp*cp;
      double c=pow(cp,2)+pow(M2c,2);

      if (cosa>0) E2new = (1./(2*a))*(-b+sqrt(pow(b,2)-4*a*c));
      if (cosa<0) E2new = (1./(2*a))*(-b-sqrt(pow(b,2)-4*a*c));
    }

    if (isnan(E2new)||isinf(E2new)||E2new<0){
      std::stringstream msg;
      msg << "problem in constraining energy to inv. mass: minv="<<Mc<< " E(set)="<<E2new << " cos(alpha)="<<cosa<<" P1="<<P1<<" P2="<< P2;
      throw(HHEnergyConstraintException(msg.str()));
    }

    return(this->changeE(E2new));

//    double P2new = sqrt(pow(E2new,2)-pow(M2c,2));
//    double Pt2new = P2new * sin(2.*atan(exp(-pmod.Eta())));
//    pmod.SetPtEtaPhiE(Pt2new,pmod.Eta(),pmod.Phi(),E2new);
//    combined = pmod + pset;
//    M = combined.M();
//  }
//  std::cout << "needed " << loopCount << " iterations to constrain E on invariant mass." << std::endl;
//  return(pmod);
}

HHKinFit2::HHLorentzVector
HHKinFit2::HHFitObjectEConstM::changeE(double E) const{
  HHLorentzVector temp = this->getFit4Vector();

  if((E<this->getLowerFitLimitE())||(E>this->getUpperFitLimitE())){
    std::stringstream msg;
    msg << "HHFitObjectEConstM: target energy is out of limits "<<"E(set)="<<E<<" "
	<<"E(limits)=["<<this->getLowerFitLimitE()<<","<< this->getUpperFitLimitE() << "]";
    throw(HHEnergyRangeException(msg.str()));
  }

  temp.SetEkeepM(E);

  return(temp);
}

void
HHKinFit2::HHFitObjectEConstM::print() const{
  std::cout << "---" << std::endl;
  std::cout <<  "energy component fit object with constant mass:" << std::endl;
  this->printInitial4Vector();
  this->printFit4Vector();
  this->printLimits();
//  this->printCovMatrix();
}
