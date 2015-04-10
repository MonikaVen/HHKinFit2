#include "HHTauTauEventGenerator.h"
#include "TRandom.h"
#include <cmath>
#include <iostream>
#include "TVector3.h"
HHTauTauEventGenerator::HHTauTauEventGenerator()
  : m_tau1(),
    m_tau2() {

// Generate lorentzvectors for the two tau in the rest frame of the higgs
  
//Tau 1
  double pi(3.14159265359);
  double mhiggs(125.7);
  double mtau1(1.77682);
  double mtau2(1.77682);
  TRandom randomnumber;
  double cthtau1(randomnumber.Uniform(-1,1)); 
  double etatau1(log(tan(acos(cthtau1)/2)));
  double phitau1(randomnumber.Uniform(0,2*pi));
  double etau1((mhiggs*mhiggs-mtau2*mtau2 +mtau1*mtau1)/(2*mhiggs));
  m_tau1.SetEEtaPhiM(etau1, etatau1, phitau1, mtau1);
  //Tau 2
  double etau2((mhiggs*mhiggs-mtau1*mtau1 +mtau2*mtau2)/(2*mhiggs));
  double etatau2(-etatau1);
  double phitau2(phitau1+pi);
  m_tau2.SetEEtaPhiM(etau2, etatau2, phitau2, mtau2);

  //generate lorentzvectors for isr-jet
  double ptisr(fabs(randomnumber.Gaus(0,30))); // absolut Value of Pt
  double phiisr(randomnumber.Uniform(0,2*pi));
  double etaisr(randomnumber.Uniform(-5,5));
  double misr(0); //mass
  m_isr.SetPtEtaPhiM(ptisr,etaisr,phiisr,misr);
  
  //generate lorentzvector for higgs with lorentzvector from isr-jet
  double cthhiggs(randomnumber.Uniform(-1,1));
  double etahiggs(log(tan(acos(cthhiggs)/2)));
  double pthiggs(ptisr);
  double phihiggs(phiisr+pi);
  m_higgs.SetPtEtaPhiM(pthiggs,etahiggs,phihiggs,mhiggs);
  
  //boosting the two tau vectors
  double bx(m_higgs.Px()/m_higgs.E());
  double by(m_higgs.Py()/m_higgs.E());
  double bz(m_higgs.Pz()/m_higgs.E());
    HHLorentzVector boost1(m_tau1);
  HHLorentzVector boost2(m_tau2);
  boost1.Boost(bx,by,bz);
  m_tau1boosted=boost1; // why  m_tau1boosted(boost1); dont work?
  boost2.Boost(bx,by,bz);
  m_tau2boosted=boost2;
}

HHLorentzVector HHTauTauEventGenerator::getTau1(){
  return(m_tau1);
}

HHLorentzVector HHTauTauEventGenerator::getTau2(){
  return(m_tau2);
}

HHLorentzVector  HHTauTauEventGenerator::getTau1boosted(){
  return(m_tau1boosted);
}

HHLorentzVector  HHTauTauEventGenerator::getTau2boosted(){
  return(m_tau2boosted);
}


