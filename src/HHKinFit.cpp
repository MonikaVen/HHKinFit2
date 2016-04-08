#include "HHKinFit.h"
#include "HHFitConstraint4Vector.h"
#include "HHFitConstraintEHardM.h"
#include "HHFitConstraint.h"
#include "HHFitConstraintLikelihood.h"
#include "HHFitObjectEConstM.h"
#include "HHFitObjectE.h"
#include "HHFitObject.h"
#include "HHFitObjectComposite.h"
#include "PSMath.h"
#include <iomanip>
#include "TAxis.h"
#include <iostream>
#include "exceptions/HHEnergyRangeException.h"

HHKinFit2::HHKinFit::HHKinFit()
: m_fitobjects(std::vector<HHFitObjectE*>()),
  m_constraints(std::vector<HHFitConstraint*>()),
  m_chi2(99999),
  m_convergence(0),
  m_printlevel(0),
  m_maxloops(500){
}

///todo: compare with old kinfit method
void 
HHKinFit2::HHKinFit::fit(){
  
  //  ----------  for PSfit -----
  const int np = m_fitobjects.size();
  double a[np];
  double astart[np];
  double alimit[np][2];
  double aprec[np];
  double daN[np];
  double h[np];
  double chi2iter[1], aMemory[np][5], g[np], H[np * np], Hinv[np * np];
  bool noNewtonShifts = false;
  
  double chi2before;
  double abefore[np];
  double daNbefore[np];
  double hbefore[np];
  double chi2iterbefore[1], aMemorybefore[np][5], gbefore[np], Hbefore[np * np], Hinvbefore[np * np];


  int iter = 0;             //  number of iterations
  int method = 1;           //  initial fit method, see PSfit()
  int mode = 1;             //  mode =1 for start of a new fit by PSfit()
  //   int icallNewton = -1;     //  init start of Newton Method
  //   int iloop = 0;            // counter for falls to fit function
  
  int iterbefore = 0;             //  number of iterations
  int methodbefore = 1;           //  initial fit method, see PSfit()
  int modebefore = 1;             //  mode =1 for start of a new fit by PSfit()
  
  bool exceptionFound = false;
  for (unsigned int i=0; i<m_fitobjects.size(); i++){
    //// fill initial tau fit parameters
    //astart[i] = m_fitobjects[i]->getE();    // energy of first tau
    //aprec[i]  = 0.1;                        //0.1                 // precision for fit
    //// fill initial step width
    //h[i] = 0.2*(m_fitobjects[i]->getUpperFitLimitE()-m_fitobjects[i]->getLowerFitLimitE());   //
    //daN[i] = 1.0;   //0.0                 // initial search direction

    astart[i]    = m_fitobjects[i]->getInitStart();
    aprec[i]     = m_fitobjects[i]->getInitPrecision();
    h[i]         = m_fitobjects[i]->getInitStepWidth();
    daN[i]       = m_fitobjects[i]->getInitDirection();
    alimit[i][0] = 1.00001* m_fitobjects[i]->getLowerFitLimitE();
    alimit[i][1] = 0.99999* m_fitobjects[i]->getUpperFitLimitE();
    
    // tau: check initial values against fit range
    if (astart[i] - h[i] < alimit[i][0]) {
      astart[i] = alimit[i][0] + h[i];
    }
    else if (astart[i] + h[i] > alimit[i][1]) {
      astart[i] = alimit[i][1] - h[i];
    }

    a[i] = astart[i];
 
    aMemory[i][0] = -999.0;
    aMemory[i][1] = -995.0;
    aMemory[i][2] = -990.0;
    aMemory[i][3] = -985.0;
    aMemory[i][3] = -980.0;
  }
  
  for (int iloop = 0; iloop < m_maxloops * 10 && iter < m_maxloops; iloop++) { // FIT loop
    for (unsigned int i=0; i<m_fitobjects.size();i++){
      try{
        m_fitobjects[i]->changeEandSave(a[i]);
      }
      catch(HHKinFit2::HHEnergyRangeException const& e){
        std::cout << e.what() << std::endl;
        throw(e);
      }
    }
    m_chi2=this->getChi2();
    //	    std::cout << iloop << " a[0]: " << a[0] << " chi2: " << std::fixed << std::setprecision(8) << chi2 << std::endl;
    //	    m_fitobjects[0]->print();
    
    
    chi2before = m_chi2;
    chi2iterbefore[0]=chi2iter[0];
    iterbefore = iter;
    methodbefore = method;
    modebefore = mode;
    for (unsigned int i=0; i<m_fitobjects.size();i++){
      abefore[i]=a[i];
      daNbefore[i]=daN[i];
      hbefore[i]=h[i];
      aMemorybefore[i][0]=aMemory[i][0];
      aMemorybefore[i][1]=aMemory[i][1];
      aMemorybefore[i][2]=aMemory[i][2];
      aMemorybefore[i][3]=aMemory[i][3];
      aMemorybefore[i][4]=aMemory[i][4];
      gbefore[i]=g[i];
      Hbefore[i]=H[i];
      Hinvbefore[i]=Hinv[i];
    }




    
    if (m_convergence != 0) break;
    m_convergence = PSMath::PSfit(iloop, iter, method, mode, noNewtonShifts, m_printlevel,
                                  np, a, astart, alimit, aprec,
                                  daN, h, aMemory, m_chi2, chi2iter, g, H,
                                  Hinv);
  }
  // ------ end of FIT loop
  
  if(m_convergence != 0 && m_convergence != 5){
    if(a[0] < (alimit[0][0] + 2*aprec[0]) ){
      if(m_convergence == 3)
        m_convergence = 5;
      else{
        //	        if (logLevel>1) std::cout << "Convergence at lower tau limit!" << std::endl;
        m_convergence = 4;
      }
    }
    if(a[0] > (alimit[0][1] - 2*aprec[0]) ){
      if(m_convergence == 3)
        m_convergence = 5;
      else{
        //		if (logLevel>1)
        //		  std::cout << "Convergence at upper tau limit!" << std::endl;
        m_convergence = 4;
      }
    }
  }
  //	  if (m_logLevel>1)
  //	    std::cout << "Convergence is " << m_convergence << std::endl;
  
}

double
HHKinFit2::HHKinFit::getChi2() const{
  double chi2=0;
  for(std::vector<HHFitConstraint*>::const_iterator it = m_constraints.begin();it != m_constraints.end(); ++it){
    //std::cout << (*it)->getChi2() << std::endl;
    chi2 += (*it)->getChi2();
  }
  //std::cout << chi2 << std::endl;
  //std::cout << "----------------------------------------------------------------------------------------------"<<std::endl;
  return(chi2);
}

void
HHKinFit2::HHKinFit::printChi2() const{
  double chi2=0;
  for(std::vector<HHFitConstraint*>::const_iterator it = m_constraints.begin();it != m_constraints.end(); ++it){
    std::cout << (*it)->getChi2() << std::endl;
    chi2 += (*it)->getChi2();
    (*it)->printChi2();
  }
  std::cout << chi2 << std::endl;
  std::cout << "----------------------------------------------------------------------------------------------"<<std::endl;
}

double
HHKinFit2::HHKinFit::getL() const{
  double L=1;
  for(std::vector<HHFitConstraint*>::const_iterator it = m_constraints.begin();it != m_constraints.end(); ++it){
    //std::cout << (*it)->getChi2() << std::endl;
    L *= (*it)->getLikelihood();
  }
  //std::cout << chi2 << std::endl;
  //std::cout << "----------------------------------------------------------------------------------------------"<<std::endl;
  return(L);
}


int
HHKinFit2::HHKinFit::getConvergence() const{
  return(m_convergence);
}

std::vector<HHKinFit2::HHFitObjectE*>
HHKinFit2::HHKinFit::getListOfFitObjects() const{
  return(m_fitobjects);
}

std::vector<HHKinFit2::HHFitConstraint*>
HHKinFit2::HHKinFit::getListOfConstraints() const{
  return(m_constraints);
}

void
HHKinFit2::HHKinFit::addFitObjectE(HHFitObjectE* fitobject){
  m_fitobjects.push_back(fitobject);
}

void
HHKinFit2::HHKinFit::addConstraint(HHFitConstraint* constraint){
  m_constraints.push_back(constraint);
}

TGraph*
HHKinFit2::HHKinFit::getChi2Function(int steps){
  int npoints(m_fitobjects.size()*steps);
  TGraph* gr = new TGraph(npoints);
  gr->SetName("chi2function");
  double stepsize((m_fitobjects[0]->getUpperFitLimitE() - m_fitobjects[0]->getLowerFitLimitE())/steps);
  for (unsigned int i=0; i<npoints; i++){
    double e = 1.00001*m_fitobjects[0]->getLowerFitLimitE()+ i*stepsize;
    m_fitobjects[0]->changeEandSave(e);
    double chi2(this->getChi2());
    //    std::cout << i << " " << e << " " << chi2 << std::endl;
    gr->SetPoint(i,e,chi2);
  }
  gr->SetMinimum(0);
  gr->GetXaxis()->SetTitle("E_{1}[GeV]");
  gr->GetYaxis()->SetTitle("#chi^{2}");
  return(gr);
}


TGraph*
HHKinFit2::HHKinFit::getLFunction(int steps){
  int npoints(m_fitobjects.size()*steps);
  TGraph* gr = new TGraph(npoints);
  gr->SetName("Lfunction");
  double stepsize((m_fitobjects[0]->getUpperFitLimitE() - m_fitobjects[0]->getLowerFitLimitE())/steps);
  for (unsigned int i=0; i<npoints; i++){
    double e = 1.00001*m_fitobjects[0]->getLowerFitLimitE()+ i*stepsize;
    m_fitobjects[0]->changeEandSave(e);
    double L(this->getL());
    //    std::cout << i << " " << e << " " << chi2 << std::endl;
    gr->SetPoint(i,e,L);
  }
  gr->SetMinimum(0);
  gr->GetXaxis()->SetTitle("E_{1}[GeV]");
  gr->GetYaxis()->SetTitle("L");
  return(gr);
}
