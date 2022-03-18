#include "TemplateFitter.cxx"
#include "TFile.h"
#include "TH2.h"

int main()
{
  TFile *f = new TFile("/nashome/q/qvuong/gen_data/output_2.root","READ");
  Int_t cutNo = 0;
  TH2D* hm = (TH2D*)f->Get(Form("m_hElepVsEv%d",cutNo));
  TH2D* he = (TH2D*)f->Get(Form("e_hElepVsEv%d",cutNo));
  TH1D * templates_m[400];
  TH1D * templates_e[400];
  for(int i=0; i<400; i++) { 
    templates_m[i] = (TH1D*)hm->ProjectionY(Form("bin%d",i+1),i+1,i+1);
    templates_e[i] = (TH1D*)he->ProjectionY(Form("bin%d",i+1),i+1,i+1);
  }
  //TH1D * intrinsic = (TH1D*)f->Get(Form("hLepE_sm%d",cutNo));
  TH1D * target_e = (TH1D*)f->Get(Form("e_hElep_w%d",cutNo));
  TH1D * target_m = (TH1D*)f->Get(Form("m_hElep_w%d",cutNo));

  TemplateFitter tf_e( templates_m, templates_e, target_e );
  //TemplateFitter tf_m( templates_m, templates_e, target_m );
  double energy_bins[401];
  for( int b = 0; b <= 400; ++b ) {
    energy_bins[b] = he->GetXaxis()->GetBinLowEdge(b+1);
  }

  tf_e.setEnergyBins( energy_bins );
  //tf_m.setEnergyBins( energy_bins );
  double bf_dm2, bf_Uee2, bf_Umm2;
  bool isOK_e = tf_e.doFit( bf_Uee2 );
  printf( "nue Best-fit Uee2 = %f\n", bf_Uee2 );
  //bool isOK_m = tf_m.doFit( bf_dm2, bf_Uee2, bf_Umm2 );
  //printf( "numu Best-fit Uee2 = %f, Umm2 = %f, dm2 = %f\n", bf_Uee2, bf_Umm2, bf_dm2 );
  //double *par;
  //par[0] = bf_theta; par[1] = bf_dm2;  
  //double chi2 = tf.getChi2( par );
/*
  // __________________ nue ______________________________

  TFile *f_nue = new TFile("/nashome/q/qvuong/data/output_12_nue.root","READ");
  Int_t cutNo_nue = 0;
  TH2D* h_nue = (TH2D*)f_nue->Get(Form("nm_hElepVsEv%d",cutNo_nue));
  TH1D * templates_nue[100];
  for(int i=0; i<100; i++) { 
    templates_nue[i] = (TH1D*)h_nue->ProjectionY(Form("bin%d",i+1),i+1,i+1);
  }
  TH1D * intrinsic_nue = (TH1D*)f_nue->Get(Form("hLepE_sm%d",cutNo_nue));
  TH1D * target_nue = (TH1D*)f_nue->Get(Form("hElep_w%d",cutNo_nue));

  
  TemplateFitter tf_nue( templates_nue, intrinsic_nue, target_nue );
  double energy_bins_nue[101];
  for( int b = 0; b <= 100; ++b ) {
    energy_bins_nue[b] = h_nue->GetXaxis()->GetBinLowEdge(b+1);
  }


  tf_nue.setEnergyBins( energy_bins_nue );
  double bf_dm2_nue, bf_theta_nue;
  bool isOK_nue = tf_nue.doFit( bf_theta_nue, bf_dm2_nue );
*/
  //printf( "Best-fit theta_nue = %f, dm2_nue = %f\n", bf_theta_nue, bf_dm2_nue );


  //tf.Draw( bf_theta, bf_dm2 );

}






