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

TemplateFitter::TemplateFitter(TH1D * CC_templates_m[480], TH1D * CC_templates_m_nc[480], TH1D * CC_templates_e[480], TH1D * nue_templates_m[480], TH1D * nue_templates_m_w[480], TH1D * nue_templates_e[480], TH1D * nue_templates_e_w[480], TH1D * CC_target_e, TH1D * CC_target_m, TH1D * nue_target)
{
  for( int i = 0; i < 480; ++i ) {
  CC_m_templates[i] = CC_templates_m[i];
  CC_nc_m_templates[i] = CC_templates_m_nc[i];
  CC_e_templates[i] = CC_templates_e[i];
  nue_m_templates[i] = nue_templates_m[i];
  nue_w_m_templates[i] = nue_templates_m_w[i];
  nue_e_templates[i] = nue_templates_e[i];
  nue_w_e_templates[i] = nue_templates_e_w[i];
  }
  
  CC_m_target = CC_target_m;
  CC_e_target = CC_target_e;
  nue_em_target = nue_target;
/*
  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
  TCanvas *c = new TCanvas("c","",800,600);
  CC_e_target->Draw();
  c->SaveAs("CC_e_hElep_w0.png");
  CC_m_target->Draw();
  c->SaveAs("CC_m_hElep_w0.png");
  nue_em_target->Draw();
  c->SaveAs("nue_hElep_w0.png");
*/
  //std::cout << m_target->GetNbinsX() << "\n";
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
double TemplateFitter::getChi2( double * par )
{ 
  // Create a histogram "temp" from the templates
  // Start with the intrinsic nu_e CC template, which doesn't change with oscillations
  TH1D * CC_tp_e = (TH1D*) CC_m_templates[0]->Clone();
  CC_tp_e->Reset();
  TH1D * CC_tp_me = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_me_nc = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_ee = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_m = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_em = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_mm = (TH1D*) CC_tp_e->Clone();

  TH1D * nue_tp = (TH1D*) nue_m_templates[0]->Clone();
  nue_tp->Reset();
  TH1D * nue_tp_em = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_ee = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_me = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_mm = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_os = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_unos = (TH1D*) nue_tp->Clone();  

  // Add in oscillated neutrinos by taking the nu_mu CC templates and weighting by the oscillation probability
  for( int i = 10; i < 480; ++i ) {
    double mue = 0;
    double ee = 0;
    double mm = 0;
    for(int j = 0; j<1001; j++){
      double e = m_energy_bins[i] + j*(m_energy_bins[i+1] - m_energy_bins[i])/1000.;
      mue = mue + getPmue(e, par[0], par[1], par[2]); // par[0] = Uee2, par[1] = Umm2, par[2] = dm2
      ee  = ee  + getPee(e, par[0], par[1], par[2]);
      mm  = mm  + getPmm(e, par[0], par[1], par[2]);
    } 
    double Pmue = mue/1001.0; 
    double Pee  = ee/1001.0; 
    double Pmm  = mm/1001.0; 

    //CC_tp_me->Add(CC_m_templates[i], Pmue);
    CC_tp_me_nc->Add(CC_nc_m_templates[i], Pmue);
    CC_tp_ee->Add(CC_e_templates[i], Pee);

    CC_tp_em->Add(CC_e_templates[i], Pmue);
    //CC_tp_me_nc->Add(CC_nc_m_templates[i], Pmue);
    CC_tp_mm->Add(CC_m_templates[i], Pmm);

    nue_tp_me->Add(nue_w_m_templates[i], Pmue);
    nue_tp_mm->Add(nue_m_templates[i], Pmm);
    nue_tp_em->Add(nue_w_e_templates[i], Pmue);
    nue_tp_ee->Add(nue_e_templates[i], Pee);

  }
  // Now we have nue temp = mu-->e (no reco cut) + e-->e (no reco cut)
  CC_tp_e->Add(CC_tp_me_nc); CC_tp_e->Add(CC_tp_ee);
  CC_tp_m->Add(CC_tp_em);    CC_tp_m->Add(CC_tp_mm);

  nue_tp_os->Add(nue_tp_me);    nue_tp_os->Add(nue_tp_em);
  nue_tp_unos->Add(nue_tp_ee);  nue_tp_unos->Add(nue_tp_mm);
  nue_tp->Add(nue_tp_os);       nue_tp->Add(nue_tp_unos);
   
  // calculate the chi2 with the "data" target
  double chi2 = 0.0;
  double CC_chi2_e = 0.0;
  double CC_chi2_m = 0.0;
  double nue_chi2  = 0.0;

  for( int bx = 1; bx <= CC_e_target->GetNbinsX(); ++bx ) {
    double tgt = CC_e_target->GetBinContent(bx);
    double diff = CC_tp_e->GetBinContent(bx) - tgt;
    if( tgt > 0. ) CC_chi2_e += (diff*diff) / tgt;
  }
  for( int bx = 1; bx <= CC_m_target->GetNbinsX(); ++bx ) {
    double tgt = CC_m_target->GetBinContent(bx);
    double diff = CC_tp_m->GetBinContent(bx) - tgt;
    if( tgt > 0. ) CC_chi2_m += (diff*diff) / tgt;
  }
  for( int bx = 1; bx <= nue_em_target->GetNbinsX(); ++bx ) {
    double tgt = nue_em_target->GetBinContent(bx);
    double diff = nue_tp->GetBinContent(bx) - tgt;
    if( tgt > 0. ) nue_chi2 += (diff*diff) / tgt;
  }

  chi2 = CC_chi2_e + CC_chi2_m + nue_chi2;
  return CC_chi2_e;
}


void TemplateFitter::Draw()
{
  TH2D *h = new TH2D("h","",100,0,0.0032,100,0,2.6);
  double par[3];
  par[0] = 0.01;
  par[1] = 0.0;
  par[2] = 0.0;
  double chi2 = 0.0;
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
  h->GetZaxis()->SetRangeUser(80,200);
  h->SetStats(0);
  h->Draw("colz");
  c->SaveAs("CC_e_chi2Surface_100.png");
}
