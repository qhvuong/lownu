#include "TemplateFitter.cxx"
#include "TFile.h"
#include "TH2.h"

int main()
{
  TFile *f = new TFile("/nashome/q/qvuong/gen_data/output_1.root","READ");
  TFile *f_nue = new TFile("/nashome/q/qvuong/gen_data/nue_output_1.root","READ");
  Int_t cutNo = 0;
  TH2D *hm = (TH2D*)f->Get(Form("m_hElepVsEv%d",cutNo));
  TH2D *hm_nc = (TH2D*)f->Get(Form("nc_m_hElepVsEv%d",cutNo));
  TH2D *he = (TH2D*)f->Get(Form("e_hElepVsEv%d",cutNo));
  TH2D *nue_hm = (TH2D*)f_nue->Get(Form("m_hElepVsEv%d",cutNo));
  TH2D *nue_he = (TH2D*)f_nue->Get(Form("e_hElepVsEv%d",cutNo));
  
  TH1D *CC_templates_m[400];
  TH1D *CC_templates_m_nc[400];
  TH1D *CC_templates_e[400];
  TH1D *nue_templates_m[400];
  TH1D *nue_templates_e[400];
  
  for(int i=0; i<400; i++) { 
    CC_templates_m[i]    = (TH1D*)hm->ProjectionY(Form("m_bin%d",i+1),i+1,i+1);
    CC_templates_m_nc[i] = (TH1D*)hm_nc->ProjectionY(Form("nc_m_bin%d",i+1),i+1,i+1);
    CC_templates_e[i]    = (TH1D*)he->ProjectionY(Form("e_bin%d",i+1),i+1,i+1);
    
    nue_templates_m[i]   = (TH1D*)nue_hm->ProjectionY(Form("nue_m_bin%d",i+1),i+1,i+1);
    nue_templates_e[i]   = (TH1D*)nue_he->ProjectionY(Form("nue_e_bin%d",i+1),i+1,i+1);
  }
  
  //TH1D * intrinsic = (TH1D*)f->Get(Form("hLepE_sm%d",cutNo));
  TH1D *CC_target_e = (TH1D*)f->Get(Form("CC_e_hElep_w%d",cutNo));
  TH1D *CC_target_m = (TH1D*)f->Get(Form("CC_m_hElep_w%d",cutNo));
  TH1D *nue_target_e = (TH1D*)f_nue->Get(Form("nue_e_hElep_w%d",cutNo));
  TH1D *nue_target_m = (TH1D*)f_nue->Get(Form("nue_m_hElep_w%d",cutNo));

  TH1D *target_e = (TH1D*)f->Get(Form("e_hElep_w%d",cutNo));
  TH1D *target_m = (TH1D*)f->Get(Form("m_hElep_w%d",cutNo));
  target_e->Add(nue_target_e);
  target_m->Add(nue_target_m);

  TemplateFitter tf( CC_templates_m, CC_templates_m_nc, CC_templates_e, nue_templates_m, nue_templates_e, target_m, target_e );
  //TemplateFitter tf_m( templates_m, templates_e, target_m );
  double energy_bins[401];
  for( int b = 0; b <= 400; ++b ) {
    energy_bins[b] = he->GetXaxis()->GetBinLowEdge(b+1);
  }

  tf.setEnergyBins( energy_bins );
  double bf_dm2, bf_Uee2, bf_Umm2;
  bool isOK = tf.doFit( bf_Uee2, bf_Umm2 , bf_dm2);
  printf( "nue Best-fit Uee2 = %f, Umm2 = %f, dm2 = %f\n", bf_Uee2, bf_Umm2, bf_dm2 );
  //tf.TrueDraw();

}



