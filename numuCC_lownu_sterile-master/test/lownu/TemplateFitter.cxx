#include "TemplateFitter.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TCanvas.h"
#include "THStack.h"
TemplateFitter::TemplateFitter(TH1D * templates_m[280], TH1D * templates_e[280], TH1D * target)
{
  for( int i = 0; i < 280; ++i ) {
  m_templates[i] = templates_m[i];
  e_templates[i] = templates_e[i];
  }
  
  m_target = target;
}

void TemplateFitter::setEnergyBins( double bins[281] )
{
  for( int i = 0; i < 281; ++ i ) {m_energy_bins[i] = bins[i];
  }
}

double TemplateFitter::getPmue( double energy, double Uee2, double Umm2, double dm2 )
{
  double L = 0.5;
  double del = 1.27*L*dm2/energy;
  double s2mue2 = 4 * Uee2 * Umm2;
  double prob = s2mue2  * pow(sin(del),2);
  return prob;
}
double TemplateFitter::getPee( double energy, double Uee2, double Umm2, double dm2 )
{
  double L = 0.5;
  double del = 1.27*L*dm2/energy;
  double s2ee2 = 4 * Uee2 * (1 - Uee2);
  double prob = 1.0 - s2ee2  * pow(sin(del),2);
  return prob;
}
double TemplateFitter::getPmm( double energy, double Uee2, double Umm2, double dm2 )
{
  double L = 0.5;
  double del = 1.27*L*dm2/energy;
  double s2mm2 = 4 * Umm2 * (1 - Umm2);
  double prob = 1.0 - s2mm2  * pow(sin(del),2);
  return prob;
}

// function whose return Minuit mimizes, must take const double* and return double
double TemplateFitter::getChi2( const double * par )
{ 
  // Create a histogram "temp" from the templates
  // Start with the intrinsic nu_e CC template, which doesn't change with oscillations
  TH1D * tp_m = (TH1D*) m_templates[0]->Clone();
  tp_m->Reset();
  TH1D *tp_em = (TH1D*) m_templates[0]->Clone();
  tp_em->Reset();
  TH1D *tp_mm = (TH1D*) m_templates[0]->Clone();
  tp_mm->Reset();

  // Add in oscillated neutrinos by taking the nu_mu CC templates and weighting by the oscillation probability
  for( int i = 10; i < 280; ++i ) {
    double mue = 0;
    double mm = 0;
    for(int j = 0; j<1000; j++){
      double e = m_energy_bins[i] + (m_energy_bins[i+1] - m_energy_bins[i])/1000.;
      //std::cout << m_energy_bins[i] << "\t";
      mue = mue + getPmue(e, par[0], par[1], par[2]); // par[0] = Uee2, par[1] = Umm2, par[2] = dm2
      mm  = mm  + getPmm(e, par[0], par[1], par[2]);
    } 
    double Pmue = mue/1000.0; 
    double Pmm  = mm/1000.0; 
    tp_em->Add(e_templates[i], Pmue);
    tp_mm->Add(m_templates[i], Pmm);
  }
  // now we have temp = intrinsic + oscillated nu_e CC
  //tp_m->Add(tp_em); 
  tp_m->Add(tp_mm);

  // calculate the chi2 with the "data" target
  double chi2 = 0.0;
  for( int bx = 1; bx <= m_target->GetNbinsX(); ++bx ) {
    double tgt = m_target->GetBinContent(bx);
    double diff = tp_m->GetBinContent(bx) - tgt;
    if( tgt > 0. ) chi2 += (diff*diff) / tgt;
  }

  TCanvas *c = new TCanvas("c","",1600,600);
  c->Divide(2,1);
  c->cd(1);
  m_target->Draw();
  tp_m->SetLineColor(2);
  tp_m->Draw("same");
  c->cd(2);
  m_target->SetLineColor(kBlack);
  m_target->Draw();
  THStack *h = new THStack("h","");
  tp_em->SetFillColor(kRed);
  tp_mm->SetFillColor(kBlue);
  h->Add(tp_em);
  h->Add(tp_mm);
  h->Draw("same");
  //tp_ee->Draw("same");
  c->SaveAs("m_fit_322.png");

  std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << chi2 << "\n";
 
  return chi2;
}



bool TemplateFitter::doFit( double &Uee2, double &Umm2, double &dm2 )
{
  // Make a Minuit fitter object
  ROOT::Math::Minimizer* fitter = ROOT::Math::Factory::CreateMinimizer("Minuit2"); 
  fitter->SetMaxFunctionCalls(1000000); // maximum number of times to try to find the minimum before failing
  fitter->SetMaxIterations(1000000);
  fitter->SetTolerance(0.1); // You might have to play with this -- how close to the correct value do you need to be?

  // The variables will be normalizations of the templates, we will start the with seed values of 1.0
  // fourth argument is step size, i.e. how much to change the normalization by at each step
  fitter->SetVariable( 0, "Uee2", 0.1, 0.0001 );
  fitter->SetVariable( 1, "Umm2", 0.1, 0.00001 );
  fitter->SetVariable( 2, "dm2", 10.1, 0.01 );
  fitter->SetVariableLowerLimit(0, 0.0);
  fitter->SetVariableLowerLimit(1, 0.0);
  fitter->SetVariableLowerLimit(2, 0.0);

  // 3 free parameters = theta, dm2
  ROOT::Math::Functor lf( this, &TemplateFitter::getChi2, 3 );
  ROOT::Math::Functor functor( lf, 3 );
  fitter->SetFunction( functor );

  // Go!
  fitter->Minimize();
/*
  if( fitter->Status() != 0 ) {
    std::cout << "Something bad happened" << std::endl;
    return false;
  }
*/
  const double *bestfit = fitter->X();
  Uee2 = bestfit[0];
  Umm2 = bestfit[1];
  dm2 = bestfit[2];
  double chi2 = fitter->MinValue();
  std::cout << "chi2= " << chi2 << "\n";
  return true;

}

