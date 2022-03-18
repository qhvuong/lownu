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
   
  //nue_tp_me->Scale(6.); nue_tp_em->Scale(1/6.);
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

  TCanvas *c_e = new TCanvas("c_e","",1600,600);
  c_e->Divide(2,1);
  c_e->cd(1);
  CC_e_target->Draw();
  CC_tp_e->SetLineColor(2);
  CC_tp_e->Draw("same");
  c_e->cd(2);
  CC_e_target->SetLineColor(kBlack);
  CC_e_target->Draw();
  THStack *h_e = new THStack("h","");
  CC_tp_me_nc->SetFillColor(kRed);
  CC_tp_ee->SetFillColor(kBlue);
  h_e->Add(CC_tp_me_nc);
  h_e->Add(CC_tp_ee);
  h_e->Draw("same");
  //tp_ee->Draw("same");
  c_e->SaveAs("CC_e_fit_132_ft.png");

  TCanvas *c_m = new TCanvas("c","",1600,600);
  c_m->Divide(2,1);
  c_m->cd(1);
  CC_m_target->Draw();
  CC_tp_m->SetLineColor(2);
  CC_tp_m->Draw("same");
  c_m->cd(2);
  CC_m_target->SetLineColor(kBlack);
  CC_m_target->Draw();
  THStack *h_m = new THStack("h","");
  CC_tp_em->SetFillColor(kRed);
  CC_tp_mm->SetFillColor(kBlue);
  h_m->Add(CC_tp_em);
  h_m->Add(CC_tp_mm);
  h_m->Draw("same");
  //tp_ee->Draw("same");
  c_m->SaveAs("CC_m_fit_132_ft.png");

  TCanvas *c = new TCanvas("c","",1600,600);
  c->Divide(2,1);
  c->cd(1);
  gPad->SetLogy();
  nue_em_target->Draw();
  nue_tp->SetLineColor(2);
  nue_tp->Draw("same");
  c->cd(2);
  gPad->SetLogy();
  nue_em_target->SetLineColor(kBlack);
  nue_em_target->Draw();
  THStack *h = new THStack("h","");
  //nue_tp_os->SetFillColor(kRed);
  //nue_tp_unos->SetFillColor(kBlue);
  nue_tp_os->SetMarkerStyle(21);
  nue_tp_os->SetMarkerColor(kRed);
  nue_tp_os->SetMarkerSize(0.5);
  nue_tp_unos->SetMarkerStyle(21);
  nue_tp_unos->SetMarkerColor(kBlue);
  nue_tp_unos->SetMarkerSize(0.5);
  h->Add(nue_tp_os);
  h->Add(nue_tp_unos);
  h->Draw("same");
  //tp_ee->Draw("same");
  c->SaveAs("nue_fit_132_ft_Log.png");


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

  return true;

}


void TemplateFitter::TrueDraw()
{
  double par[3];
  par[0]=0.01;
  par[1]=0.0016;
  par[2]=1.3;

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
   
  //nue_tp_me->Scale(6.); nue_tp_em->Scale(1/6.);
  nue_tp_os->Add(nue_tp_me);    nue_tp_os->Add(nue_tp_em);
  nue_tp_unos->Add(nue_tp_ee);  nue_tp_unos->Add(nue_tp_mm);
  nue_tp->Add(nue_tp_os);       nue_tp->Add(nue_tp_unos);

  TCanvas *c_e = new TCanvas("c_e","",1600,600);
  c_e->Divide(2,1);
  c_e->cd(1);
  CC_e_target->Draw();
  CC_tp_e->SetLineColor(2);
  CC_tp_e->Draw("same");
  c_e->cd(2);
  CC_e_target->SetLineColor(kBlack);
  CC_e_target->Draw();
  THStack *h_e = new THStack("h","");
  CC_tp_me_nc->SetFillColor(kRed);
  CC_tp_ee->SetFillColor(kBlue);
  h_e->Add(CC_tp_me_nc);
  h_e->Add(CC_tp_ee);
  h_e->Draw("same");
  //tp_ee->Draw("same");
  c_e->SaveAs("CC_e_TrueDraw_13_ft.png");

  TCanvas *c_m = new TCanvas("c","",1600,600);
  c_m->Divide(2,1);
  c_m->cd(1);
  CC_m_target->Draw();
  CC_tp_m->SetLineColor(2);
  CC_tp_m->Draw("same");
  c_m->cd(2);
  CC_m_target->SetLineColor(kBlack);
  CC_m_target->Draw();
  THStack *h_m = new THStack("h","");
  CC_tp_em->SetFillColor(kRed);
  CC_tp_mm->SetFillColor(kBlue);
  h_m->Add(CC_tp_em);
  h_m->Add(CC_tp_mm);
  h_m->Draw("same");
  //tp_ee->Draw("same");
  c_m->SaveAs("CC_m_TrueDraw_13_ft.png");

  TCanvas *c = new TCanvas("c","",1600,600);
  c->Divide(2,1);
  c->cd(1);
  gPad->SetLogy();
  nue_em_target->Draw();
  nue_tp->SetLineColor(2);
  nue_tp->Draw("same");
  c->cd(2);
  gPad->SetLogy();
  nue_em_target->SetLineColor(kBlack);
  nue_em_target->Draw();
  THStack *h = new THStack("h","");
  //nue_tp_os->SetFillColor(kRed);
  //nue_tp_unos->SetFillColor(kBlue);
  nue_tp_os->SetMarkerStyle(21);
  nue_tp_os->SetMarkerColor(kRed);
  nue_tp_os->SetMarkerSize(0.5);
  nue_tp_unos->SetMarkerStyle(21);
  nue_tp_unos->SetMarkerColor(kBlue);
  nue_tp_unos->SetMarkerSize(0.5);
  h->Add(nue_tp_os);
  h->Add(nue_tp_unos);
  h->Draw("same");
  //tp_ee->Draw("same");
  c->SaveAs("nue_TrueDraw_13_ft_Log.png");

}

