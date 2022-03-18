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

TemplateFitter::TemplateFitter(TH1D * templates_m[480], TH1D * templates_m_w[480], TH1D * templates_e[480], TH1D * templates_e_w[480], TH1D * target)
{
  for( int i = 0; i < 480; ++i ) {
  m_templates[i] = templates_m[i];
  w_m_templates[i] = templates_m_w[i];
  e_templates[i] = templates_e[i];
  w_e_templates[i] = templates_e_w[i];
  }
  
  //m_intrinsic = intrinsic;
  m_target = target;
}

void TemplateFitter::setEnergyBins( double bins[481] )
{
  for( int i = 0; i < 481; ++ i ) {m_energy_bins[i] = bins[i];
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
    double ee  = 0;
    double mm  = 0;
    for(int j=0; j<1001; j++){
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
    tp_em->Add(w_e_templates[i], Pmue);
    tp_ee->Add(e_templates[i], Pee);
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

  return chi2;
}

void TemplateFitter::Draw()
{
  TH2D *h = new TH2D("h","",100,0,0.02,100,0,12.0);
  double par[3];
  par[0] = 0.04;
  par[1] = 0.0;
  par[2] = 0.0;
  double chi2 = 0;
  for(int i=1; i<=100; i++)
  {
    par[1] = h->GetXaxis()->GetBinCenter(i);
    for(int j=1; j<=100; j++)
    {
      par[2] = h->GetYaxis()->GetBinCenter(j);
      chi2 = getChi2(par);
      std::cout << par[1] << "\t" << par[2] << "\t" << chi2 << "\n";
      h->Fill(par[1], par[2], chi2);
    }
  }

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
  TCanvas *c = new TCanvas("c","",800,600);
  c->SetLogz();
  h->SetTitle("Chi2 Surface");
  h->GetXaxis()->SetTitle("Umm2");
  h->GetYaxis()->SetTitle("dm2");
  h->GetZaxis()->SetRangeUser(1, 500);
  h->SetStats(0);
  h->Draw("colz");
  c->SaveAs("nue_test.png");
  //c->SaveAs("true_nue_chi2Surface_20.png");

}

