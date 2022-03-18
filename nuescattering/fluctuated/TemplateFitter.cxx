#include "TemplateFitter.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TCanvas.h"
#include "THStack.h"
TemplateFitter::TemplateFitter(TH1D * templates_m[480], TH1D * templates_m_w[480], TH1D * templates_e[480], TH1D * templates_e_w[480], TH1D * target)
{
  for( int i = 0; i < 480; ++i ) {
  m_templates[i] = templates_m[i];
  w_m_templates[i] = templates_m_w[i];
  e_templates[i] = templates_e[i];
  w_e_templates[i] = templates_e_w[i];
  }
  
  m_target = target;
}

void TemplateFitter::setEnergyBins( double bins[481] )
{
  for( int i = 0; i < 481; ++ i ) {m_energy_bins[i] = bins[i];
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
  TH1D * tp = (TH1D*) m_templates[0]->Clone();
  tp->Reset();
  TH1D *tp_os = (TH1D*) tp->Clone();
  TH1D *tp_unos = (TH1D*) tp->Clone();
  TH1D *tp_me = (TH1D*) tp->Clone();
  TH1D *tp_mm = (TH1D*) tp->Clone();
  TH1D *tp_em = (TH1D*) tp->Clone();
  TH1D *tp_ee = (TH1D*) tp->Clone();

  // Add in oscillated neutrinos by taking the nu_mu CC templates and weighting by the oscillation probability
  for( int i = 10; i < 480; ++i ) {
    double mue = 0;
    double ee = 0;
    double mm = 0;
    for(int j = 0; j<1001; j++){
      double e = m_energy_bins[i] + j*(m_energy_bins[i+1] - m_energy_bins[i])/1000.;
      //std::cout << m_energy_bins[i] << "\t";
      mue = mue + getPmue(e, par[0], par[1], par[2]); // par[0] = Uee2, par[1] = Umm2, par[2] = dm2
      ee  = ee  + getPee(e, par[0], par[1], par[2]);
      mm  = mm  + getPmm(e, par[0], par[1], par[2]);
    } 
    double Pmue = mue/1001.0; 
    double Pee  = ee/1001.0; 
    double Pmm  = mm/1001.0; 
    tp_me->Add(w_m_templates[i], Pmue);
    tp_mm->Add(m_templates[i], Pmm);
    tp_ee->Add(e_templates[i], Pee);
    tp_em->Add(w_e_templates[i], Pmue);
  }
  // now we have temp = intrinsic + oscillated nu_e CC
  tp_os->Add(tp_me); tp_os->Add(tp_em);
  tp_unos->Add(tp_mm); tp_unos->Add(tp_ee);
  tp->Add(tp_os); tp->Add(tp_unos);

  // calculate the chi2 with the "data" target
  double chi2 = 0.0;
  for( int bx = 1; bx <= m_target->GetNbinsX(); ++bx ) {
    double tgt = m_target->GetBinContent(bx);
    double diff = tp->GetBinContent(bx) - tgt;
    if( tgt > 0. ) chi2 += (diff*diff) / tgt;
  }

  TCanvas *c = new TCanvas("c","",1600,600);
  c->Divide(2,1);
  c->cd(1);
  //gPad->SetLogy();
  m_target->Draw();
  tp->SetLineColor(2);
  tp->Draw("same");
  c->cd(2);
  //gPad->SetLogy();
  m_target->SetLineColor(kBlack);
  m_target->Draw();
  THStack *h = new THStack("h","");
  tp_os->SetMarkerStyle(21);
  tp_os->SetMarkerSize(0.5);
  tp_os->SetMarkerColor(kRed);
  tp_unos->SetMarkerStyle(21);
  tp_unos->SetMarkerSize(0.5);
  tp_unos->SetMarkerColor(kBlue);
  h->Add(tp_os);
  h->Add(tp_unos);
  h->Draw("same");
  c->SaveAs("true_nue_fit_12_ft.png");

  TCanvas *c1 = new TCanvas("c1","",1600,600);
  c1->Divide(2,1);
  c1->cd(1);
  gPad->SetLogy();
  m_target->Draw();
  tp->SetLineColor(2);
  tp->Draw("same");
  c1->cd(2);
  gPad->SetLogy();
  m_target->SetLineColor(kBlack);
  m_target->Draw();
  THStack *h1 = new THStack("h1","");
  tp_os->SetMarkerStyle(21);
  tp_os->SetMarkerSize(0.5);
  tp_os->SetMarkerColor(kRed);
  tp_unos->SetMarkerStyle(21);
  tp_unos->SetMarkerSize(0.5);
  tp_unos->SetMarkerColor(kBlue);
  h1->Add(tp_os);
  h1->Add(tp_unos);
  h1->Draw("same");
  c1->SaveAs("true_nue_fit_12_ft_Log.png");

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
  fitter->SetVariable( 0, "Uee2", 0.1, 0.00001 );
  fitter->SetVariable( 1, "Umm2", 0.1, 0.00001 );
  fitter->SetVariable( 2, "dm2", 10.0, 0.01 );
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

void TemplateFitter::TrueDraw()
{ 
  double par[3];
  par[0]=0.01;
  par[1]=0.0016;
  par[2]=1.3;

  // Create a histogram "temp" from the templates
  // Start with the intrinsic nu_e CC template, which doesn't change with oscillations
  TH1D * tp = (TH1D*) m_templates[0]->Clone();
  tp->Reset();
  TH1D *tp_os = (TH1D*) tp->Clone();
  TH1D *tp_unos = (TH1D*) tp->Clone();
  TH1D *tp_me = (TH1D*) tp->Clone();
  TH1D *tp_mm = (TH1D*) tp->Clone();
  TH1D *tp_em = (TH1D*) tp->Clone();
  TH1D *tp_ee = (TH1D*) tp->Clone();

  // Add in oscillated neutrinos by taking the nu_mu CC templates and weighting by the oscillation probability
  for( int i = 10; i < 480; ++i ) {
    double mue = 0;
    double ee = 0;
    double mm = 0;
    for(int j = 0; j<1001; j++){
      double e = m_energy_bins[i] + j*(m_energy_bins[i+1] - m_energy_bins[i])/1000.;
      //std::cout << m_energy_bins[i] << "\t";
      mue = mue + getPmue(e, par[0], par[1], par[2]); // par[0] = Uee2, par[1] = Umm2, par[2] = dm2
      ee  = ee  + getPee(e, par[0], par[1], par[2]);
      mm  = mm  + getPmm(e, par[0], par[1], par[2]);
    } 
    double Pmue = mue/1001.0; 
    double Pee  = ee/1001.0; 
    double Pmm  = mm/1001.0; 
    tp_me->Add(w_m_templates[i], Pmue);
    tp_mm->Add(m_templates[i], Pmm);
    tp_ee->Add(e_templates[i], Pee);
    tp_em->Add(w_e_templates[i], Pmue);
  }
  // now we have temp = intrinsic + oscillated nu_e CC
  tp_os->Add(tp_me); tp_os->Add(tp_em);
  tp_unos->Add(tp_mm); tp_unos->Add(tp_ee);
  tp->Add(tp_os); tp->Add(tp_unos);

  TCanvas *c = new TCanvas("c","",1600,600);
  c->Divide(2,1);
  c->cd(1);
  //gPad->SetLogy();
  m_target->Draw();
  tp->SetLineColor(2);
  tp->Draw("same");
  c->cd(2);
  //gPad->SetLogy();
  m_target->SetLineColor(kBlack);
  m_target->Draw();
  THStack *h = new THStack("h","");
  tp_os->SetMarkerStyle(21);
  tp_os->SetMarkerSize(0.5);
  tp_os->SetMarkerColor(kRed);
  tp_unos->SetMarkerStyle(21);
  tp_unos->SetMarkerSize(0.5);
  tp_unos->SetMarkerColor(kBlue);
  h->Add(tp_os);
  h->Add(tp_unos);
  h->Draw("same");
  c->SaveAs("true_nue_TrueDraw_1_ft.png");

  TCanvas *c1 = new TCanvas("c1","",1600,600);
  c1->Divide(2,1);
  c1->cd(1);
  gPad->SetLogy();
  m_target->Draw();
  tp->SetLineColor(2);
  tp->Draw("same");
  c1->cd(2);
  gPad->SetLogy();
  m_target->SetLineColor(kBlack);
  m_target->Draw();
  THStack *h1 = new THStack("h1","");
  tp_os->SetMarkerStyle(21);
  tp_os->SetMarkerSize(0.5);
  tp_os->SetMarkerColor(kRed);
  tp_unos->SetMarkerStyle(21);
  tp_unos->SetMarkerSize(0.5);
  tp_unos->SetMarkerColor(kBlue);
  h1->Add(tp_os);
  h1->Add(tp_unos);
  h1->Draw("same");
  c1->SaveAs("true_nue_TrueDraw_1_ft_Log.png");
}
