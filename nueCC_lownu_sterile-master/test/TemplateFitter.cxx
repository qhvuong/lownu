#include "TemplateFitter.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TStyle.h"

TemplateFitter::TemplateFitter(TH1D * templates_m[400], TH1D * templates_e[400], TH1D * target)
{
  for( int i = 0; i < 400; ++i ) {
  m_templates[i] = templates_m[i];
  e_templates[i] = templates_e[i];
  }
  
  //m_intrinsic = intrinsic;
  m_target = target;
}

void TemplateFitter::setEnergyBins( double bins[401] )
{
  for( int i = 0; i < 401; ++ i ) {m_energy_bins[i] = bins[i];
  //std::cout << m_energy_bins[i] << "\t";
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
double TemplateFitter::getChi2( double * par )
{ 
  // Create a histogram "temp" from the templates
  // Start with the intrinsic nu_e CC template, which doesn't change with oscillations
  TH1D * tp_e = (TH1D*) e_templates[0]->Clone();
  tp_e->Reset();
  
  TH1D *tp_me = (TH1D*) e_templates[0]->Clone();
  tp_me->Reset();
  TH1D *tp_ee = (TH1D*) e_templates[0]->Clone();
  tp_ee->Reset();

  // Add in oscillated neutrinos by taking the nu_mu CC templates and weighting by the oscillation probability
  for( int i = 0; i < 4; ++i ) {
    double mue = 0;
    double ee  = 0;
    for(int j=0; j<10; j++){
    double e = m_energy_bins[i] + j*(m_energy_bins[i+1] - m_energy_bins[i])/10.;
    std::cout << m_energy_bins[i] << "\t" << e << "\n";
    mue = mue + getPmue(e, par[0], par[1], par[2]); // par[0] = Uee2, par[1] = Umm2, par[2] = dm2
    ee  = ee  + getPee(e, par[0], par[1], par[2]); 
    }
    double Pmue = mue/10.;
    double Pee  = ee/10.;
    //if(0.0014<par[1] && par[1]<0.0018 && 1.2<par[2] && par[2]<1.4) std::cout << e << "\t" << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << Pmue << "\t" << Pee << "\t" << Pmm << "\n"; 
    tp_me->Add(m_templates[i], Pmue);
    tp_ee->Add(e_templates[i], Pee);
  }
  // now we have temp = intrinsic + oscillated nu_e CC
  tp_e->Add(tp_me); tp_e->Add(tp_ee);

  // calculate the chi2 with the "data" target
  double chi2 = 0.0;
  for( int bx = 1; bx <= m_target->GetNbinsX(); ++bx ) {
    double tgt = m_target->GetBinContent(bx);
    double diff = tp_e->GetBinContent(bx) - tgt;
    if( tgt > 0. ) chi2 += (diff*diff) / tgt;
  }

  //std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << chi2 << "\n";
 
  return chi2;
}


void TemplateFitter::Draw()
{
  TH2D *h = new TH2D("h","",100,0,0.08,100,0,0.02);
  double par[3];
  par[2] = 6.0;
  par[0] = 0.0;
  par[1] = 0.0;
  double chi2 = 0;
  for(int i=1; i<=100; i++)
  {
    par[0] = h->GetXaxis()->GetBinCenter(i);
    for(int j=1; j<=100; j++)
    {
      par[1] = h->GetYaxis()->GetBinCenter(j);
      chi2 = getChi2(par);
      std::cout << par[0] << "\t" << par[1] << "\t" << chi2 << "\n";
      h->Fill(par[0], par[1], chi2);
    }
  }

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
  TCanvas *c = new TCanvas("c","",800,600);
  c->SetLogz();
  h->SetTitle("Chi2 Surface");
  h->GetXaxis()->SetTitle("Uee2");
  h->GetYaxis()->SetTitle("Umm2");
  h->Draw("colz");
  c->SaveAs("e_chi2Surface_220_0.png");
}


void TemplateFitter::TrueDraw()
{
  double par[3];
  par[0]=0.01;
  par[1]=0.0016;
  par[2]=1.3;
  TH1D * tp_e = (TH1D*) e_templates[0]->Clone();
  tp_e->Reset();
  TH1D *tp_me = (TH1D*) e_templates[0]->Clone();
  tp_me->Reset();
  TH1D *tp_ee = (TH1D*) e_templates[0]->Clone();
  tp_ee->Reset();

  // Add in oscillated neutrinos by taking the nu_mu CC templates and weighting by the oscillation probability
  for( int i = 0; i < 400; ++i ) {
    double mue = 0;
    double ee  = 0;
    for(int j=0; j<100; j++){
      double e = m_energy_bins[i] + (m_energy_bins[i+1] - m_energy_bins[i])/100.;
      mue = mue + getPmue(e, par[0], par[1], par[2]); // par[0] = Uee2, par[1] = Umm2, par[2] = dm2
      ee  = ee  + getPee(e, par[0], par[1], par[2]); 
    }
    double Pmue = mue/100.;
    double Pee  = ee/100.;
    tp_me->Add(m_templates[i], Pmue);
    tp_ee->Add(e_templates[i], Pee);
  }
  // now we have temp = intrinsic + oscillated nu_e CC
  tp_e->Add(tp_me); tp_e->Add(tp_ee);

  // calculate the chi2 with the "data" target
  double chi2 = 0.0;
  for( int bx = 1; bx <= m_target->GetNbinsX(); ++bx ) {
    double tgt = m_target->GetBinContent(bx);
    double diff = tp_e->GetBinContent(bx) - tgt;
    if( tgt > 0. ) chi2 += (diff*diff) / tgt;
  }

  TCanvas *c = new TCanvas("c","",1600,600);
  c->Divide(2,1);
  c->cd(1);
  m_target->Draw();
  tp_e->SetLineColor(2);
  tp_e->Draw("same");
  c->cd(2);
  m_target->SetLineColor(kBlack);
  m_target->Draw();
  THStack *h = new THStack("h","");
  tp_me->SetFillColor(kRed);
  tp_ee->SetFillColor(kBlue);
  h->Add(tp_me);
  h->Add(tp_ee);
  h->Draw("same");
  c->SaveAs("TrueDraw_10_0.png"); 
}
