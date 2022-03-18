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
#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TGraph.h"
#include "TLegend.h"

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
  target_nue = nue_target;
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

TMatrixD covmtr(300,300);

void TemplateFitter::setCovmtr( double bincontent[301][301] )
{
  for(int i=0; i<300; i++) {
    for(int j=0; j<300; j++) {
      covmtr[i][j] = bincontent[i][j];
    }
  }
}

const int N = 40;
TH2D *h0 = new TH2D("h0","",N,0,0.1, N,0,12.0);
TH2D *h1 = new TH2D("h1","",N,0,12.0,N,0,0.1);
TH2D *h2 = new TH2D("h2","",N,0,0.1, N,0,0.1);

int para = 2;
int cutNo = 0;
int cutEv = 0;
int s = 0;
// function whose return Minuit mimizes, must take const double* and return double
double TemplateFitter::getChi2( const double * par )
{ 
  // Create a histogram "temp" from the templates
  // Start with the intrinsic nu_e CC template, which doesn't change with oscillations
   TH1D * CC_tp_e = (TH1D*) CC_m_templates[0]->Clone();
  CC_tp_e->Reset();
  TH1D * CC_tp_me = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_ee = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_m  = (TH1D*) CC_tp_e->Clone();
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

    CC_tp_me->Add(CC_nc_m_templates[i], Pmue);
    CC_tp_ee->Add(CC_e_templates[i], Pee);

    CC_tp_em->Add(CC_e_templates[i], Pmue);
    CC_tp_mm->Add(CC_m_templates[i], Pmm);

    nue_tp_me->Add(nue_w_m_templates[i], Pmue);
    nue_tp_mm->Add(nue_m_templates[i], Pmm);
    nue_tp_em->Add(nue_w_e_templates[i], Pmue);
    nue_tp_ee->Add(nue_e_templates[i], Pee);
  }

  // Now we have nue temp = mu-->e (no reco cut) + e-->e (no reco cut)
  CC_tp_e->Add(CC_tp_me); CC_tp_e->Add(CC_tp_ee);
  CC_tp_m->Add(CC_tp_em);    CC_tp_m->Add(CC_tp_mm);

  nue_tp_os->Add(nue_tp_me);    nue_tp_os->Add(nue_tp_em);
  nue_tp_unos->Add(nue_tp_ee);  nue_tp_unos->Add(nue_tp_mm);
  nue_tp->Add(nue_tp_os);       nue_tp->Add(nue_tp_unos);

  // calculate the chi2 with the "data" target
  const int nbins_Elep = 100;
  double chi2 = 0.0;

  TMatrixD target(3*nbins_Elep, 1);
  TMatrixD temp(3*nbins_Elep, 1);
  TMatrixD diff(3*nbins_Elep, 1);
  TMatrixD diff_T(1, 3*nbins_Elep);
  TMatrixD unc(3*nbins_Elep, 3*nbins_Elep);
  TMatrixD cov(3*nbins_Elep, 3*nbins_Elep);
  TMatrixD covmtr_tot(3*nbins_Elep, 3*nbins_Elep);

  for( int bx = 0; bx < nbins_Elep; bx++ ) {
    temp[bx][0] = CC_tp_m->GetBinContent(bx+1); 
    temp[bx+nbins_Elep][0] = CC_tp_e->GetBinContent(bx+1); 
    temp[bx+2*nbins_Elep][0] = nue_tp->GetBinContent(bx+1); 

    target[bx][0] = CC_m_target->GetBinContent(bx+1); 
    target[bx+nbins_Elep][0] = CC_e_target->GetBinContent(bx+1); 
    target[bx+2*nbins_Elep][0] = target_nue->GetBinContent(bx+1); 

    diff[bx][0] = temp[bx][0] - target[bx][0];
    diff[bx+nbins_Elep][0] = temp[bx+nbins_Elep][0] - target[bx+nbins_Elep][0];
    diff[bx+2*nbins_Elep][0] = temp[bx+2*nbins_Elep][0] - target[bx+2*nbins_Elep][0];
 
    diff_T[0][bx] = diff[bx][0]; 
    diff_T[0][bx+nbins_Elep] = diff[bx+nbins_Elep][0]; 
    diff_T[0][bx+2*nbins_Elep] = diff[bx+2*nbins_Elep][0]; 

    for( int by = 0; by < nbins_Elep; by++ ) {
      unc[bx][by] 	       = 0;
      unc[bx][by+nbins_Elep]   = 0;
      unc[bx][by+2*nbins_Elep] = 0;

      unc[bx+nbins_Elep][by]              = 0;
      unc[bx+nbins_Elep][by+nbins_Elep]   = 0;
      unc[bx+nbins_Elep][by+2*nbins_Elep] = 0;

      unc[bx+2*nbins_Elep][by]              = 0;
      unc[bx+2*nbins_Elep][by+nbins_Elep]   = 0;
      unc[bx+2*nbins_Elep][by+2*nbins_Elep] = 0;

      if(bx == by) {
        unc[bx][by] 			      = CC_m_target->GetBinContent(bx+1);
        unc[bx+nbins_Elep][by+nbins_Elep]     = CC_e_target->GetBinContent(bx+1);
        unc[bx+2*nbins_Elep][by+2*nbins_Elep] = target_nue->GetBinContent(bx+1);
      }
/*
      cov[bx][by]              = covmtr[bx][by]              * CC_m_target->GetBinContent(bx+1) * CC_m_target->GetBinContent(by+1);
      cov[bx][by+nbins_Elep]   = covmtr[bx][by+nbins_Elep]   * CC_m_target->GetBinContent(bx+1) * CC_e_target->GetBinContent(by+1);
      cov[bx][by+2*nbins_Elep] = covmtr[bx][by+2*nbins_Elep] * CC_m_target->GetBinContent(bx+1) * target_nue->GetBinContent(by+1);

      cov[bx+nbins_Elep][by]              = covmtr[bx+nbins_Elep][by]              * CC_e_target->GetBinContent(bx+1) * CC_m_target->GetBinContent(by+1);
      cov[bx+nbins_Elep][by+nbins_Elep]   = covmtr[bx+nbins_Elep][by+nbins_Elep]   * CC_e_target->GetBinContent(bx+1) * CC_e_target->GetBinContent(by+1);
      cov[bx+nbins_Elep][by+2*nbins_Elep] = covmtr[bx+nbins_Elep][by+2*nbins_Elep] * CC_e_target->GetBinContent(bx+1) * target_nue->GetBinContent(by+1);

      cov[bx+2*nbins_Elep][by]              = covmtr[bx+2*nbins_Elep][by]              * target_nue->GetBinContent(bx+1) * CC_m_target->GetBinContent(by+1);
      cov[bx+2*nbins_Elep][by+nbins_Elep]   = covmtr[bx+2*nbins_Elep][by+nbins_Elep]   * target_nue->GetBinContent(bx+1) * CC_e_target->GetBinContent(by+1);
      cov[bx+2*nbins_Elep][by+2*nbins_Elep] = covmtr[bx+2*nbins_Elep][by+2*nbins_Elep] * target_nue->GetBinContent(bx+1) * target_nue->GetBinContent(by+1);
*/
    }
  }

  covmtr_tot = unc;
  TDecompSVD svd(covmtr_tot);
  TMatrixD inv = svd.Invert();

  TMatrixD diff_cov(1, 3*nbins_Elep);
  for( int bx = 0; bx < 3*nbins_Elep; bx++ ) {
    for( int by = 0; by < 3*nbins_Elep; by++ ) {
      diff_cov[0][bx] += diff_T[0][by] * inv[by][bx];
    }
    chi2 += diff_cov[0][bx] * diff[bx][0];
  }

  std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << chi2 << "\n";

  h0->Fill(par[1], par[2]);
  h1->Fill(par[2], par[0]);
  h2->Fill(par[0], par[1]);
 
  return chi2;
}

double bf0[1], bf1[1], bf2[1];

bool TemplateFitter::doFit( double &Uee2, double &Umm2, double &dm2 )
{
  double seed1[3][3], seed2[3][3];
  seed1[0][0] = 0.01;
  seed1[0][1] = 0.0016;
  seed1[0][2] = 1.3;
  seed1[1][0] = 0.001;
  seed1[1][1] = 0.001;
  seed1[1][2] = 0.1;
  seed1[2][0] = 0.1;
  seed1[2][1] = 0.1;
  seed1[2][2] = 10.0;

  seed2[0][0] = 0.04;
  seed2[0][1] = 0.01;
  seed2[0][2] = 6.0;
  seed2[1][0] = 0.001;
  seed2[1][1] = 0.001;
  seed2[1][2] = 0.1;
  seed2[2][0] = 0.1;
  seed2[2][1] = 0.1;
  seed2[2][2] = 10.0;

  double p0, p1, p2;
  if(para == 1) {
    p0 = seed1[s][0];
    p1 = seed1[s][1];
    p2 = seed1[s][2];
  }
  else if(para == 2) {
    p0 = seed2[s][0];
    p1 = seed2[s][1];
    p2 = seed2[s][2];
  }
  // Make a Minuit fitter object
  ROOT::Math::Minimizer* fitter = ROOT::Math::Factory::CreateMinimizer("Minuit2"); 
  fitter->SetMaxFunctionCalls(1000000); // maximum number of times to try to find the minimum before failing
  fitter->SetMaxIterations(1000000);
  fitter->SetTolerance(0.1); // You might have to play with this -- how close to the correct value do you need to be?

  // The variables will be normalizations of the templates, we will start the with seed values of 1.0
  // fourth argument is step size, i.e. how much to change the normalization by at each step
  fitter->SetVariable( 0, "Uee2", p0, 0.00001 );
  fitter->SetVariable( 1, "Umm2", p1, 0.00001 );
  fitter->SetVariable( 2, "dm2",  p2, 0.01 );
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
  bf0[0] = bestfit[0];
  bf1[0] = bestfit[1];
  bf2[0] = bestfit[2];
  
  double chi2 = fitter->MinValue();
  //h->Fill(Uee2, Umm2, dm2);
  return true;

}

void TemplateFitter::Draw()
{
  double par[3];
  par[0] = bf0[0];
  par[1] = bf1[0];
  par[2] = bf2[0];

  TH1D * CC_tp_e = (TH1D*) CC_m_templates[0]->Clone();
  CC_tp_e->Reset();
  TH1D * CC_tp_me = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_ee = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_m  = (TH1D*) CC_tp_e->Clone();
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

    CC_tp_me->Add(CC_nc_m_templates[i], Pmue);
    CC_tp_ee->Add(CC_e_templates[i], Pee);

    CC_tp_em->Add(CC_e_templates[i], Pmue);
    CC_tp_mm->Add(CC_m_templates[i], Pmm);

    nue_tp_me->Add(nue_w_m_templates[i], Pmue);
    nue_tp_mm->Add(nue_m_templates[i], Pmm);
    nue_tp_em->Add(nue_w_e_templates[i], Pmue);
    nue_tp_ee->Add(nue_e_templates[i], Pee);
  }

  // Now we have nue temp = mu-->e (no reco cut) + e-->e (no reco cut)
  CC_tp_e->Add(CC_tp_me); CC_tp_e->Add(CC_tp_ee);
  CC_tp_m->Add(CC_tp_em);    CC_tp_m->Add(CC_tp_mm);

  nue_tp_os->Add(nue_tp_me);    nue_tp_os->Add(nue_tp_em);
  nue_tp_unos->Add(nue_tp_ee);  nue_tp_unos->Add(nue_tp_mm);
  nue_tp->Add(nue_tp_os);       nue_tp->Add(nue_tp_unos);

  CC_e_target->SetName("e_target");
  CC_e_target->SetLineColor(kBlack);
  CC_m_target->SetName("m_target");
  CC_m_target->SetLineColor(kBlack);
  target_nue->SetName("nue_target");
  target_nue->SetLineColor(kBlack);

  CC_tp_e->SetName("tp_e");
  CC_tp_e->SetMarkerStyle(21);
  CC_tp_e->SetMarkerColor(kRed);
  CC_tp_e->SetMarkerSize(0.5);
  CC_tp_me->SetName("tp_me");
  CC_tp_me->SetMarkerStyle(21);
  CC_tp_me->SetMarkerColor(kGreen);
  CC_tp_me->SetMarkerSize(0.5);
  CC_tp_ee->SetName("tp_ee");
  CC_tp_ee->SetMarkerStyle(21);
  CC_tp_ee->SetMarkerColor(kBlue);
  CC_tp_ee->SetMarkerSize(0.5);

  CC_tp_m->SetName("tp_m");
  CC_tp_m->SetMarkerStyle(21);
  CC_tp_m->SetMarkerColor(kRed);
  CC_tp_m->SetMarkerSize(0.5);
  CC_tp_em->SetName("tp_em");
  CC_tp_em->SetMarkerStyle(21);
  CC_tp_em->SetMarkerColor(kGreen);
  CC_tp_em->SetMarkerSize(0.5);
  CC_tp_mm->SetName("tp_mm");
  CC_tp_mm->SetMarkerStyle(21);
  CC_tp_mm->SetMarkerColor(kBlue);
  CC_tp_mm->SetMarkerSize(0.5);

  nue_tp->SetName("tp_nue");
  nue_tp->SetMarkerStyle(21);
  nue_tp->SetMarkerColor(kRed);
  nue_tp->SetMarkerSize(0.5);
  nue_tp_os->SetName("tp_nue_os");
  nue_tp_os->SetMarkerStyle(21);
  nue_tp_os->SetMarkerColor(kGreen);
  nue_tp_os->SetMarkerSize(0.5);
  nue_tp_unos->SetName("tp_nue_unos");
  nue_tp_unos->SetMarkerStyle(21);
  nue_tp_unos->SetMarkerColor(kBlue);
  nue_tp_unos->SetMarkerSize(0.5);

  TCanvas *c_e = new TCanvas("c_e","",900,700);
  gPad->SetLogy();
  CC_e_target->GetXaxis()->SetRangeUser(0,60);
  CC_e_target->Draw();
  CC_tp_e->Draw("same");
  CC_tp_me->Draw("same");
  CC_tp_ee->Draw("same");
  TLegend *legend_e = new TLegend(0.75,0.78,0.9,0.9);
  legend_e->AddEntry("e_target","nueCC target");
  legend_e->AddEntry("tp_e","universal template fit");
  legend_e->AddEntry("tp_me","oscillated mu->e");
  legend_e->AddEntry("tp_ee","unoscillated e->e");
  legend_e->Draw();
  c_e->SaveAs(Form("Ev_CC_e_fit_%d%d_%d%d_ft.png",para,cutNo,cutEv,s));

  TCanvas *c_m = new TCanvas("c_m","",900,700);
  gPad->SetLogy();
  CC_m_target->GetXaxis()->SetRangeUser(0,100);
  CC_m_target->Draw();
  CC_tp_m->Draw("same");
  CC_tp_em->Draw("same");
  CC_tp_mm->Draw("same");
  TLegend *legend_m = new TLegend(0.75,0.78,0.9,0.9);
  legend_m->AddEntry("m_target","numuCC target");
  legend_m->AddEntry("tp_m","universal template fit");
  legend_m->AddEntry("tp_em","oscillated mu->e");
  legend_m->AddEntry("tp_mm","unoscillated e->e");
  legend_m->Draw();
  c_m->SaveAs(Form("Ev_CC_m_fit_%d%d_%d%d_ft.png",para,cutNo,cutEv,s));

  TCanvas *c = new TCanvas("c","",900,700);
  gPad->SetLogy();
  target_nue->GetXaxis()->SetRangeUser(0,80);
  target_nue->Draw();
  nue_tp->Draw("same");
  nue_tp_os->Draw("same");
  nue_tp_unos->Draw("same");
  TLegend *legend_nue = new TLegend(0.75,0.78,0.9,0.9);
  legend_nue->AddEntry("nue_target","nu+e target");
  legend_nue->AddEntry("tp_nue","universal template fit");
  legend_nue->AddEntry("tp_nue_os","oscillated mu->e + e->mu");
  legend_nue->AddEntry("tp_nue_unos","unoscillated e->e + mu->mu");
  legend_nue->Draw();
  c->SaveAs(Form("Ev_nue_fit_%d%d_%d%d_ft_Log.png",para,cutNo,cutEv,s));

  TFile *f = new TFile(Form("/dune/app/users/qvuong/lownu_analysis/Ev_combine/chi2/Ev_chi2_%d%d_%d.root",para,cutNo,cutEv));
  TH2D *hc0 = (TH2D*)f->Get("h0");
  TH2D *hc1 = (TH2D*)f->Get("h1");
  TH2D *hc2 = (TH2D*)f->Get("h2");

  double seed1[3][3], seed2[3][3];
  seed1[0][0] = 0.01;
  seed1[0][1] = 0.0016;
  seed1[0][2] = 1.3;
  seed1[1][0] = 0.001;
  seed1[1][1] = 0.001;
  seed1[1][2] = 0.1;
  seed1[2][0] = 0.1;
  seed1[2][1] = 0.1;
  seed1[2][2] = 10.0;

  seed2[0][0] = 0.04;
  seed2[0][1] = 0.01;
  seed2[0][2] = 6.0;
  seed2[1][0] = 0.001;
  seed2[1][1] = 0.001;
  seed2[1][2] = 0.1;
  seed2[2][0] = 0.1;
  seed2[2][1] = 0.1;
  seed2[2][2] = 10.0;

  double p0[1], p1[1], p2[1]; //true values
  double s0[1], s1[1], s2[1];
  if(para == 1) {
    p0[0] = seed1[0][0];
    p1[0] = seed1[0][1];
    p2[0] = seed1[0][2];

    s0[0] = seed1[s][0];
    s1[0] = seed1[s][1];
    s2[0] = seed1[s][2];
  }
  else if(para == 2) {
    p0[0] = seed2[0][0];
    p1[0] = seed2[0][1];
    p2[0] = seed2[0][2];

    s0[0] = seed2[s][0];
    s1[0] = seed2[s][1];
    s2[0] = seed2[s][2];
  }

  int nG = 1;
  TGraph *g0   = new TGraph(nG,p1,p2);
  TGraph *g1   = new TGraph(nG,p2,p0);
  TGraph *g2   = new TGraph(nG,p0,p1);
  TGraph *g0s  = new TGraph(nG,s1,s2);
  TGraph *g1s  = new TGraph(nG,s2,s0);
  TGraph *g2s  = new TGraph(nG,s0,s1);
  TGraph *g0BF = new TGraph(nG,bf1,bf2);
  TGraph *g1BF = new TGraph(nG,bf2,bf0);
  TGraph *g2BF = new TGraph(nG,bf0,bf1);

  g0->SetName("g0");
  g0->SetMarkerStyle(29);
  g0->SetMarkerColor(kRed);
  g0->SetMarkerSize(3);
  g1->SetName("g1");
  g1->SetMarkerStyle(29);
  g1->SetMarkerColor(kRed);
  g1->SetMarkerSize(3);
  g2->SetName("g2");
  g2->SetMarkerStyle(29);
  g2->SetMarkerColor(kRed);
  g2->SetMarkerSize(3);

  g0s->SetName("g0s");
  g0s->SetMarkerStyle(3);
  g0s->SetMarkerColor(kGreen);
  g0s->SetMarkerSize(3);
  g1s->SetName("g1s");
  g1s->SetMarkerStyle(3);
  g1s->SetMarkerColor(kGreen);
  g1s->SetMarkerSize(3);
  g2s->SetName("g2s");
  g2s->SetMarkerStyle(3);
  g2s->SetMarkerColor(kGreen);
  g2s->SetMarkerSize(3);

  g0BF->SetName("g0BF");
  g0BF->SetMarkerStyle(8);
  g0BF->SetMarkerColor(kBlue);
  g0BF->SetMarkerSize(1.1);
  g1BF->SetName("g1BF");
  g1BF->SetMarkerStyle(8);
  g1BF->SetMarkerColor(kBlue);
  g1BF->SetMarkerSize(1.1);
  g2BF->SetName("g2BF");
  g2BF->SetMarkerStyle(8);
  g2BF->SetMarkerColor(kBlue);
  g2BF->SetMarkerSize(1.1);

  h0->GetXaxis()->SetTitle("Umm2");
  h0->GetYaxis()->SetTitle("dm2");
  h0->SetStats(0);
  h1->GetXaxis()->SetTitle("dm2");
  h1->GetYaxis()->SetTitle("Uee2");
  h1->SetStats(0);
  h2->GetXaxis()->SetTitle("Uee2");
  h2->GetYaxis()->SetTitle("Umm2");
  h2->SetStats(0);

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();

  TCanvas *c0 = new TCanvas("c0","",800,600);
  c0->SetLogz();
  hc0->SetMaximum(1E6);
  hc0->Draw("colz");
  h0->Draw("same");
  g0->Draw("same P");
  g0s->Draw("same P");
  g0BF->Draw("same P");
  TLegend *legend0 = new TLegend(0.75,0.78,0.9,0.9);
  legend0->AddEntry("g0","  True values");
  legend0->AddEntry("g0s","  Initial Seeds");
  legend0->AddEntry("g0BF","  Bestfit values");
  legend0->Draw();
  c0->SaveAs(Form("Ev_parDraw_%d%d0_%d%d.png",para,cutNo,cutEv,s));
  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->SetLogz();
  hc1->SetMaximum(1E6);
  hc1->Draw("colz");
  h1->Draw("same");
  g1->Draw("same P");
  g1s->Draw("same P");
  g1BF->Draw("same P");
  TLegend *legend1 = new TLegend(0.75,0.78,0.9,0.9);
  legend1->AddEntry("g1","  True values");
  legend1->AddEntry("g1s","  Initial Seeds");
  legend1->AddEntry("g1BF","  Bestfit values");
  legend1->Draw();
  c1->SaveAs(Form("Ev_parDraw_%d%d1_%d%d.png",para,cutNo,cutEv,s));
  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->SetLogz();
  hc2->SetMaximum(1E6);
  hc2->Draw("colz");
  h2->Draw("same");
  g2->Draw("same P");
  g2s->Draw("same P");
  g2BF->Draw("same P");
  TLegend *legend2 = new TLegend(0.75,0.78,0.9,0.9);
  legend2->AddEntry("g2","  True values");
  legend2->AddEntry("g2s","  Initial Seeds");
  legend2->AddEntry("g2BF","  Bestfit values");
  legend2->Draw();
  c2->SaveAs(Form("Ev_parDraw_%d%d2_%d%d.png",para,cutNo,cutEv,s));
}

void TemplateFitter::TrueDraw()
{
  double seed1[3][3], seed2[3][3];
  seed1[0][0] = 0.01;
  seed1[0][1] = 0.0016;
  seed1[0][2] = 1.3;
  seed2[0][0] = 0.04;
  seed2[0][1] = 0.01;
  seed2[0][2] = 6.0;
  double par[3];
  if(para == 1) {
    par[0] = seed1[0][0];
    par[1] = seed1[0][1];
    par[2] = seed1[0][2];
  }
  else if(para == 2) {
    par[0] = seed2[0][0];
    par[1] = seed2[0][1];
    par[2] = seed2[0][2];
  }

  TH1D * CC_tp_e = (TH1D*) CC_m_templates[0]->Clone();
  CC_tp_e->Reset();
  TH1D * CC_tp_me = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_ee = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_m  = (TH1D*) CC_tp_e->Clone();
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

    CC_tp_me->Add(CC_nc_m_templates[i], Pmue);
    CC_tp_ee->Add(CC_e_templates[i], Pee);

    CC_tp_em->Add(CC_e_templates[i], Pmue);
    CC_tp_mm->Add(CC_m_templates[i], Pmm);

    nue_tp_me->Add(nue_w_m_templates[i], Pmue);
    nue_tp_mm->Add(nue_m_templates[i], Pmm);
    nue_tp_em->Add(nue_w_e_templates[i], Pmue);
    nue_tp_ee->Add(nue_e_templates[i], Pee);
  }

  // Now we have nue temp = mu-->e (no reco cut) + e-->e (no reco cut)
  CC_tp_e->Add(CC_tp_me); CC_tp_e->Add(CC_tp_ee);
  CC_tp_m->Add(CC_tp_em);    CC_tp_m->Add(CC_tp_mm);

  nue_tp_os->Add(nue_tp_me);    nue_tp_os->Add(nue_tp_em);
  nue_tp_unos->Add(nue_tp_ee);  nue_tp_unos->Add(nue_tp_mm);
  nue_tp->Add(nue_tp_os);       nue_tp->Add(nue_tp_unos);

  CC_e_target->SetName("e_target");
  CC_e_target->SetLineColor(kBlack);
  CC_m_target->SetName("m_target");
  CC_m_target->SetLineColor(kBlack);
  target_nue->SetName("nue_target");
  target_nue->SetLineColor(kBlack);

  CC_tp_e->SetName("tp_e");
  CC_tp_e->SetMarkerStyle(21);
  CC_tp_e->SetMarkerColor(kRed);
  CC_tp_e->SetMarkerSize(0.5);
  CC_tp_me->SetName("tp_me");
  CC_tp_me->SetMarkerStyle(21);
  CC_tp_me->SetMarkerColor(kGreen);
  CC_tp_me->SetMarkerSize(0.5);
  CC_tp_ee->SetName("tp_ee");
  CC_tp_ee->SetMarkerStyle(21);
  CC_tp_ee->SetMarkerColor(kBlue);
  CC_tp_ee->SetMarkerSize(0.5);

  CC_tp_m->SetName("tp_m");
  CC_tp_m->SetMarkerStyle(21);
  CC_tp_m->SetMarkerColor(kRed);
  CC_tp_m->SetMarkerSize(0.5);
  CC_tp_em->SetName("tp_em");
  CC_tp_em->SetMarkerStyle(21);
  CC_tp_em->SetMarkerColor(kGreen);
  CC_tp_em->SetMarkerSize(0.5);
  CC_tp_mm->SetName("tp_mm");
  CC_tp_mm->SetMarkerStyle(21);
  CC_tp_mm->SetMarkerColor(kBlue);
  CC_tp_mm->SetMarkerSize(0.5);

  nue_tp->SetName("tp_nue");
  nue_tp->SetMarkerStyle(21);
  nue_tp->SetMarkerColor(kRed);
  nue_tp->SetMarkerSize(0.5);
  nue_tp_os->SetName("tp_nue_os");
  nue_tp_os->SetMarkerStyle(21);
  nue_tp_os->SetMarkerColor(kGreen);
  nue_tp_os->SetMarkerSize(0.5);
  nue_tp_unos->SetName("tp_nue_unos");
  nue_tp_unos->SetMarkerStyle(21);
  nue_tp_unos->SetMarkerColor(kBlue);
  nue_tp_unos->SetMarkerSize(0.5);

  TCanvas *c_e = new TCanvas("c_e","",900,700);
  gPad->SetLogy();
  CC_e_target->GetXaxis()->SetRangeUser(0,60);
  CC_e_target->Draw();
  CC_tp_e->Draw("same");
  CC_tp_me->Draw("same");
  CC_tp_ee->Draw("same");
  TLegend *legend_e = new TLegend(0.75,0.78,0.9,0.9);
  legend_e->AddEntry("e_target","nueCC target");
  legend_e->AddEntry("tp_e","universal template fit");
  legend_e->AddEntry("tp_me","oscillated mu->e");
  legend_e->AddEntry("tp_ee","unoscillated e->e");
  legend_e->Draw();
  c_e->SaveAs(Form("Ev_CC_e_TrueDraw_%d%d_%d_ft.png",para,cutNo,cutEv));

  TCanvas *c_m = new TCanvas("c_m","",900,700);
  gPad->SetLogy();
  CC_m_target->GetXaxis()->SetRangeUser(0,100);
  CC_m_target->Draw();
  CC_tp_m->Draw("same");
  CC_tp_em->Draw("same");
  CC_tp_mm->Draw("same");
  TLegend *legend_m = new TLegend(0.75,0.78,0.9,0.9);
  legend_m->AddEntry("m_target","numuCC target");
  legend_m->AddEntry("tp_m","universal template fit");
  legend_m->AddEntry("tp_em","oscillated mu->e");
  legend_m->AddEntry("tp_mm","unoscillated e->e");
  legend_m->Draw();
  c_m->SaveAs(Form("Ev_CC_m_TrueDraw_%d%d_%d_ft.png",para,cutNo,cutEv));

  TCanvas *c = new TCanvas("c","",900,700);
  gPad->SetLogy();
  target_nue->GetXaxis()->SetRangeUser(0,80);
  target_nue->Draw();
  nue_tp->Draw("same");
  nue_tp_os->Draw("same");
  nue_tp_unos->Draw("same");
  TLegend *legend_nue = new TLegend(0.75,0.78,0.9,0.9);
  legend_nue->AddEntry("nue_target","nu+e target");
  legend_nue->AddEntry("tp_nue","universal template fit");
  legend_nue->AddEntry("tp_nue_os","oscillated mu->e + e->mu");
  legend_nue->AddEntry("tp_nue_unos","unoscillated e->e + mu->mu");
  legend_nue->Draw();
  c->SaveAs(Form("Ev_nue_TrueDraw_%d%d_%d_ft_Log.png",para,cutNo,cutEv));
}
/*
  double CC_chi2_e = 0.0;
  double CC_chi2_m = 0.0;
  double nue_chi2 = 0.0;

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
*/

