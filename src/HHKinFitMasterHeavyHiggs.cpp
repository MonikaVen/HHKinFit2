#include "HHKinFitMasterHeavyHiggs.h"
#include "HHKinFit.h"
#include "TMatrixD.h"
#include "TRandom3.h"
#include "HHFitObjectEConstM.h"
#include "HHFitObjectEConstBeta.h"
#include "HHFitObjectMET.h"
#include "HHFitObjectComposite.h"
#include "HHFitConstraintLikelihood.h"
#include "HHFitConstraint.h"
#include "HHFitConstraintEHardM.h"
#include "HHFitObjectEConstM.h"
#include "HHFitObjectEConstBeta.h"
#include "HHFitConstraint4Vector.h"

#include "exceptions/HHEnergyRangeException.h"
#include "exceptions/HHLimitSettingException.h"
#include "exceptions/HHCovarianceMatrixException.h"
#include <TMath.h>
#include <cmath>
#include <cstdlib>
#include <iterator>

void HHKinFit2::HHKinFitMasterHeavyHiggs::fit()
{
  HHFitObjectE* tau1Fit = new HHFitObjectEConstM(m_tauvis1);
  HHFitObjectE* tau2Fit = new HHFitObjectEConstM(m_tauvis2);

  HHLorentzVector bJet1min = m_bjet1;
  if (m_bjet1.E() - 5.0*m_sigma_bjet1 <= 0) bJet1min.SetEkeepBeta(10);
  else bJet1min.SetEkeepBeta(m_bjet1.E() - 5.0*m_sigma_bjet1);
  
  HHLorentzVector bJet2min = m_bjet2;
  if (m_bjet2.E() - 5.0*m_sigma_bjet2 <= 0) bJet2min.SetEkeepBeta(10);
  else bJet2min.SetEkeepBeta(m_bjet2.E() - 5.0*m_sigma_bjet2);

  HHFitObjectE* b1Fit = new HHFitObjectEConstBeta(m_bjet1);
  HHFitObjectE* b2Fit = new HHFitObjectEConstBeta(m_bjet2);
  
  //prepare MET object
  HHFitObjectMET* metFit = new HHFitObjectMET(m_MET);

  //TODO: Why cant I hand over the matrix here?
  metFit->setCovMatrix(m_MET_COV[0][0], m_MET_COV[1][1], m_MET_COV[1][0]);
  
  //prepare composite object: Higgs
  HHFitObject* higgs  = new HHFitObjectComposite(tau1Fit, tau2Fit, 
						 b1Fit, b2Fit,
						 metFit);
  HHFitObject* higgs1  = new HHFitObjectComposite(tau1Fit, tau2Fit);  
  HHFitObject* higgs2  = new HHFitObjectComposite(b1Fit, b2Fit);
  
  for(unsigned int i = 0; i < m_hypos.size(); ++i)
  { 
    int mh1 = m_hypos[i].first;
    int mh2 = m_hypos[i].second;
  
    try{
      tau1Fit->setFitLimitsE(tau1Fit->getInitial4Vector(),mh1,tau2Fit->getInitial4Vector());
      tau2Fit->setFitLimitsE(tau2Fit->getInitial4Vector(),mh1,tau1Fit->getInitial4Vector());
    }
    catch(HHLimitSettingException const& e){
      std::cout << "Exception while setting tau limits:" << std::endl;
      std::cout << e.what() << std::endl;
      std::cout << "Tau energies are not compatible with invariant mass constraint." << std::endl;

      m_map_chi2[m_hypos[i]] = -pow(10,10);
      m_map_prob[m_hypos[i]] = -pow(10,10);
      //m_bestHypo = HHFitHypothesisHeavyHiggs(-pow(10,10),-pow(10,10));
      m_bestHypo = HHFitHypothesisHeavyHiggs( -1, -1);
      m_chi2_best = -pow(10,10);
      m_mH_best = -pow(10,10);
      continue;
    }
  
    
    try{
      b1Fit->setFitLimitsE(bJet1min, mh2, bJet2min);
      b2Fit->setFitLimitsE(bJet2min, mh2, bJet1min);
    }
    catch(HHLimitSettingException const& e){
      std::cout << "Exception while setting b-jet limits" << std::endl;
      std::cout << e.what() << std::endl;
      std::cout << "Bjet energies are not compatible within 5 sigma with invariant mass constraint." << std::endl;
      m_map_chi2[m_hypos[i]] = -pow(10,10);
      m_map_prob[m_hypos[i]] = -pow(10,10);
      //m_bestHypo = HHFitHypothesisHeavyHiggs(-pow(10,10),-pow(10,10));
      m_bestHypo = HHFitHypothesisHeavyHiggs( -1, -1);
      m_chi2_best = -pow(10,10);
      m_mH_best = -pow(10,10);
      continue;
    }
    b1Fit->setCovMatrix(pow(m_sigma_bjet1,2));
    b2Fit->setCovMatrix(pow(m_sigma_bjet2,2));
    

    //prepare constraints
    HHFitConstraint* c_invmh1 = new HHFitConstraintEHardM(tau1Fit, tau2Fit, mh1);
    HHFitConstraint* c_invmh2 = new HHFitConstraintEHardM(b1Fit, b2Fit, mh2);
    
    HHFitConstraint* c_b1 = new HHFitConstraint4Vector(b1Fit, false, false, 
						       false, true);
    HHFitConstraint* c_b2 = new HHFitConstraint4Vector(b2Fit, false, false, 
						       false, true);
    HHFitConstraint* c_balance = new HHFitConstraint4Vector(higgs, true, true, 
							    false, false);

    //fit
    HHKinFit2::HHKinFit* fitObject = new HHKinFit2::HHKinFit();

    tau1Fit->setInitStart((tau1Fit->getUpperFitLimitE()+tau1Fit->getLowerFitLimitE())/2);
    tau1Fit->setInitPrecision(0.1);
    tau1Fit->setInitStepWidth(0.1*(tau1Fit->getUpperFitLimitE() - tau1Fit->getLowerFitLimitE()));
    tau1Fit->setInitDirection(1.0);

    b1Fit->setInitStart(b1Fit->getInitial4Vector().E());
    b1Fit->setInitPrecision(0.002*b1Fit->getInitial4Vector().E());
    b1Fit->setInitStepWidth(0.5*m_sigma_bjet1);
    b1Fit->setInitDirection(1.0);

    fitObject->addFitObjectE(tau1Fit);
    fitObject->addFitObjectE(b1Fit);

    fitObject->addConstraint(c_invmh1);
    fitObject->addConstraint(c_invmh2);
    fitObject->addConstraint(c_b1);
    fitObject->addConstraint(c_b2);
    fitObject->addConstraint(c_balance);

    fitObject->fit();

    m_map_convergence[m_hypos[i]] = fitObject->getConvergence();    
    
    double chi2 = fitObject->getChi2();
    m_map_chi2[m_hypos[i]] = chi2;
    m_map_prob[m_hypos[i]] = TMath::Prob(chi2, 2);

    if(chi2 < m_chi2_best)
    {
      m_bestHypo = m_hypos[i];
      m_chi2_best = chi2;
      m_mH_best = higgs->getFit4Vector().M();
    }
    
    m_map_mH[m_hypos[i]] = higgs->getFit4Vector().M();
    m_map_chi2BJet1[m_hypos[i]] = c_b1->getChi2();
    m_map_chi2BJet2[m_hypos[i]] = c_b2->getChi2();
    m_map_chi2Balance[m_hypos[i]] = c_balance->getChi2();
    
    TLorentzVector fittedTau1 =  ( (TLorentzVector)tau1Fit->getFit4Vector()  );
    m_map_fittedTau1[m_hypos[i]] = fittedTau1;
    TLorentzVector fittedTau2 =  ( (TLorentzVector)tau2Fit->getFit4Vector()  );
    m_map_fittedTau2[m_hypos[i]] = fittedTau2;
    TLorentzVector fittedB1 =  ( (TLorentzVector)b1Fit->getFit4Vector() ) ;
    m_map_fittedB1[m_hypos[i]] = fittedB1;
    TLorentzVector fittedB2 =  ( (TLorentzVector)b2Fit->getFit4Vector() ) ;
    m_map_fittedB2[m_hypos[i]]= fittedB2;
    
    /*
      if( entry_convergence_full.second == 0 ){
      	Chi2Map chi2Map = advancedfitter.CreateChi2Map(15, 100);
	
	TCanvas* c1 = new TCanvas("canvas1");
	TGraph2D* graph2d = new TGraph2D( chi2Map.size() );
	
	int pointN = 0;
	for(Chi2Map::iterator iter = chi2Map.begin(); iter != chi2Map.end(); iter++){
	  graph2d->SetPoint(pointN, iter->first.first, iter->first.second, iter->second);
	  pointN++;
	}
	
	std::stringstream fileNameStream;
	TString fileName;
	fileNameStream << "MinB1_" << m_bjet1_fitted.E() << "_MinTau1_" << m_tau1_fitted.E() << std::endl;
	fileNameStream >> fileName;
	
	graph2d->Draw("Cont3COLZ");
	c1->SaveAs(fileName + ".pdf");
	
	/*
	Chi2Map chi2MapAroundStartvalues = advancedfitter.CreateChi2MapAroundStartvalues(15, 10);
	
	TCanvas* c2 = new TCanvas("canvas2");
	TGraph2D* graph2dAroundStartvalues = new TGraph2D( chi2MapAroundStartvalues.size() );
	
	pointN = 0;
	for(Chi2Map::iterator iter = chi2MapAroundStartvalues.begin(); iter != chi2MapAroundStartvalues.end(); iter++){
	  graph2dAroundStartvalues->SetPoint(pointN, iter->first.first, iter->first.second, iter->second);
	  pointN++;
	}
	
	graph2dAroundStartvalues->Draw("Cont3COLZ");
	c2->SaveAs(fileName + "AroundStartvalues.pdf");
	
	Chi2Map chi2MapTauValues = advancedfitter.CreateChi2Map(0, 100);

	TCanvas* c3 = new TCanvas("canvas3");
	TGraph* graph1dTauValues = new TGraph( chi2MapTauValues.size() );
	
	pointN = 0;
	for(Chi2Map::iterator iter = chi2MapTauValues.begin(); iter != chi2MapTauValues.end(); iter++){
	  graph1dTauValues->SetPoint(pointN, iter->first.second, iter->second);
	  pointN++;
	}

	graph1dTauValues->Draw("AC*");
	c3->SaveAs(fileName + "TauValues.pdf");
	
	
	delete graph1dTauValues;
	//delete graph2dAroundStartvalues;
	//delete c2;
	
	delete c3;
	
	delete c1;
	delete graph2d;
	}
    */

    delete c_invmh1;
    delete c_invmh2;
    delete c_b1;
    delete c_b2;
    delete c_balance;
    
    delete fitObject;
   
    tau1Fit->reset();
    tau2Fit->reset();
    b1Fit->reset();
    b2Fit->reset();
    metFit->reset();
    higgs->reset();
    higgs1->reset();
    higgs2->reset();
  }
  
  delete tau1Fit;
  delete tau2Fit;
  delete b1Fit;
  delete b2Fit;
  delete metFit;
  delete higgs;
  delete higgs1;
  delete higgs2;
}

void HHKinFit2::HHKinFitMasterHeavyHiggs::addHypo(HHFitHypothesisHeavyHiggs hypo)
{
  m_hypos.push_back(hypo);
}

void HHKinFit2::HHKinFitMasterHeavyHiggs::addHypo(int mh1, int mh2)
{
  HHFitHypothesisHeavyHiggs hypo(mh1, mh2);
  m_hypos.push_back(hypo);
}

HHKinFit2::HHKinFitMasterHeavyHiggs::HHKinFitMasterHeavyHiggs(
                                                    TLorentzVector bjet1, 
						    double sigmaEbjet1,
						    TLorentzVector bjet2, 
						    double sigmaEbjet2,
						    TLorentzVector tauvis1,
						    TLorentzVector tauvis2,
						    TVector2 met, 
						    TMatrixD met_cov, 
						    bool istruth, 
						    TLorentzVector* heavyhiggsgen)
  :m_MET_COV(TMatrixD(4,4))
{
  
  m_bjet1 = HHLorentzVector(bjet1.Px(), bjet1.Py(), bjet1.Pz(), bjet1.E());
  m_bjet2 = HHLorentzVector(bjet2.Px(), bjet2.Py(), bjet2.Pz(), bjet2.E());
  m_tauvis1 = HHLorentzVector(tauvis1.Px(), tauvis1.Py(), tauvis1.Pz(), tauvis1.E());  
  m_tauvis2 = HHLorentzVector(tauvis2.Px(), tauvis2.Py(), tauvis2.Pz(), tauvis2.E());
 
  m_tauvis1.SetMkeepE(1.77682);
  m_tauvis2.SetMkeepE(1.77682);
   
  m_MET = met;
  m_MET_COV = met_cov;

  m_sigma_bjet1 = sigmaEbjet1;
  m_sigma_bjet2 = sigmaEbjet2;

  m_chi2_best = 99999.0;
  m_mH_best = -1.0; 
  m_bestHypo = HHFitHypothesisHeavyHiggs(0,0);

  if (istruth){
    TRandom3 r(0);

    Double_t bjet1_E  = r.Gaus(bjet1.E(), sigmaEbjet1);
    Double_t bjet1_P  = sqrt(pow(bjet1_E,2) - pow(bjet1.M(),2));
    Double_t bjet1_Pt = sin(bjet1.Theta())*bjet1_P;

    m_bjet1.SetPtEtaPhiE(bjet1_Pt, bjet1.Eta(), bjet1.Phi(), bjet1_E);

    TMatrixD bjet1Cov(2,2);
    // error propagation p=sqrt(e^2-m^2)
    Double_t bjet1_dpt = sin(bjet1.Theta())*bjet1.E()/bjet1.P()*sigmaEbjet1;
    bjet1Cov(0,0) = pow(cos(bjet1.Phi())*bjet1_dpt,2);                           
    bjet1Cov(0,1) = sin(bjet1.Phi())*cos(bjet1.Phi())*bjet1_dpt*bjet1_dpt;
    bjet1Cov(1,0) = sin(bjet1.Phi())*cos(bjet1.Phi())*bjet1_dpt*bjet1_dpt; 
    bjet1Cov(1,1) = pow(sin(bjet1.Phi())*bjet1_dpt,2);

    Double_t bjet2_E  = r.Gaus(bjet2.E(), sigmaEbjet2);
    Double_t bjet2_P  = sqrt(pow(bjet2_E,2) - pow(bjet2.M(),2));
    Double_t bjet2_Pt = sin(bjet2.Theta())*bjet2_P;

    m_bjet2.SetPtEtaPhiE(bjet2_Pt, bjet2.Eta(), bjet2.Phi(), bjet2_E);

    TMatrixD bjet2Cov(2,2);
    // error propagation p=sqrt(e^2-m^2)
    Double_t bjet2_dpt = sin(bjet2.Theta())*bjet2.E()/bjet2.P()*sigmaEbjet2;  
    bjet2Cov(0,0) = pow(cos(bjet2.Phi())*bjet2_dpt,2);  
    bjet2Cov(0,1) = sin(bjet2.Phi())*cos(bjet2.Phi())*bjet2_dpt*bjet2_dpt;
    bjet2Cov(1,0) = sin(bjet2.Phi())*cos(bjet2.Phi())*bjet2_dpt*bjet2_dpt;    
    bjet2Cov(1,1) = pow(sin(bjet2.Phi())*bjet2_dpt,2);


    HHLorentzVector recoil;
    if(heavyhiggsgen != NULL){
       Double_t pxRecoil = r.Gaus(-(heavyhiggsgen->Px() ), 10.0);
       Double_t pyRecoil = r.Gaus(-(heavyhiggsgen->Py() ), 10.0);

       recoil = HHLorentzVector(pxRecoil, pyRecoil, 0,
				sqrt(pxRecoil*pxRecoil+pyRecoil*pyRecoil));
    }
    else{
      recoil = HHLorentzVector(0,0,0,0);
      std::cout << "WARNING! Truthinput mode active but no Heavy Higgs gen-information given! Setting Recoil to Zero!" << std::endl;
    }

    TMatrixD recoilCov(2,2);
    recoilCov(0,0)=100;  recoilCov(0,1)=0;
    recoilCov(1,0)=0;    recoilCov(1,1)=100;

    HHLorentzVector recoHH = m_bjet1 + m_bjet2 + m_tauvis1 + m_tauvis2 + recoil;
    m_MET = TVector2(-recoHH.Px(), -recoHH.Py() );

    m_MET_COV = TMatrixD(2,2);
    m_MET_COV = recoilCov + bjet1Cov + bjet2Cov;

  }  
}


//Getters
HHKinFit2::HHFitHypothesisHeavyHiggs HHKinFit2::HHKinFitMasterHeavyHiggs::getBestHypothesis(){
  return(m_bestHypo);
}


//Getters for fit results

double HHKinFit2::HHKinFitMasterHeavyHiggs::getChi2(int mh1, int mh2){
HHFitHypothesisHeavyHiggs hypo(mh1, mh2);
return(m_map_chi2[hypo]);
}

  
double HHKinFit2::HHKinFitMasterHeavyHiggs::getChi2BJet1(int mh1, int mh2){
HHFitHypothesisHeavyHiggs hypo(mh1, mh2);
return(m_map_chi2BJet1[hypo]);
}



double HHKinFit2::HHKinFitMasterHeavyHiggs::getChi2BJet2(int mh1, int mh2){
HHFitHypothesisHeavyHiggs hypo(mh1, mh2);
return(m_map_chi2BJet2[hypo]);
}
    

double HHKinFit2::HHKinFitMasterHeavyHiggs::getChi2Balance(int mh1, int mh2){
HHFitHypothesisHeavyHiggs hypo(mh1, mh2);
return(m_map_chi2Balance[hypo]);
}
  
double HHKinFit2::HHKinFitMasterHeavyHiggs::getFitProb(int mh1, int mh2){
HHFitHypothesisHeavyHiggs hypo(mh1, mh2);
return( m_map_prob[hypo]);
}
  
double HHKinFit2::HHKinFitMasterHeavyHiggs::getMH(int mh1, int mh2){
HHFitHypothesisHeavyHiggs hypo(mh1, mh2);
return(m_map_mH[hypo]);
}


int HHKinFit2::HHKinFitMasterHeavyHiggs::getConvergence(int mh1, int mh2){
HHFitHypothesisHeavyHiggs hypo(mh1, mh2);
return(m_map_convergence[hypo]);
}

TLorentzVector HHKinFit2::HHKinFitMasterHeavyHiggs::getFittedTau1(int mh1, int mh2){
HHFitHypothesisHeavyHiggs hypo(mh1, mh2);
return(m_map_fittedTau1[hypo]);
}


TLorentzVector HHKinFit2::HHKinFitMasterHeavyHiggs::getFittedTau2(int mh1, int mh2){
HHFitHypothesisHeavyHiggs hypo(mh1, mh2);
return(m_map_fittedTau2[hypo]);
}

TLorentzVector HHKinFit2::HHKinFitMasterHeavyHiggs::getFittedBJet1(int mh1, int mh2){
HHFitHypothesisHeavyHiggs hypo(mh1, mh2);
return(m_map_fittedB1[hypo]);
}
  

TLorentzVector HHKinFit2::HHKinFitMasterHeavyHiggs::getFittedBJet2(int mh1, int mh2){
HHFitHypothesisHeavyHiggs hypo(mh1, mh2);
return(m_map_fittedB2[hypo]);
}
    
