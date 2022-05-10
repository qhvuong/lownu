#include "TemplateFitter.cxx"
#include "TFile.h"
#include "TH2.h"
#include "TStyle.h"

int main()
{
  char var[20] = "Elep";
  int oscpar, nuCut, EvCut, seed;
  oscpar = 2;
  TFile *CC_f = new TFile(Form("/dune/app/users/qvuong/lownu/gen_data/CC/output_%d.root",oscpar),"READ");
  TFile *nue_f = new TFile(Form("/dune/app/users/qvuong/lownu/gen_data/nuescattering/nue_output_%d.root",oscpar),"READ");
  for(nuCut = 0; nuCut < 4; nuCut ++) {
  if(nuCut != 0 && nuCut != 3) continue;
  //nuCut = 0;
  //for(EvCut=0; EvCut<3; EvCut++) {
  EvCut = 0;
  for(seed=0; seed<3; seed++) {

  TH2D* CC_hm = (TH2D*)CC_f->Get(Form("m_h%sVsEv%d",var,nuCut));
  TH2D* CC_hm_nc = (TH2D*)CC_f->Get(Form("nc_m_h%sVsEv%d",var,nuCut));
  TH2D* CC_he = (TH2D*)CC_f->Get(Form("e_h%sVsEv%d",var,nuCut));
  TH2D* nue_hm = (TH2D*)nue_f->Get(Form("m_h%sVsEv%d",var,EvCut));
  TH2D* nue_hm_w = (TH2D*)nue_f->Get(Form("m_h%sVsEv%d_w",var,EvCut));
  TH2D* nue_he = (TH2D*)nue_f->Get(Form("e_h%sVsEv%d",var,EvCut));
  TH2D* nue_he_w = (TH2D*)nue_f->Get(Form("e_h%sVsEv%d_w",var,EvCut));

  TH1D * CC_templates_m[480];
  TH1D * CC_templates_m_nc[480];
  TH1D * CC_templates_e[480];
  TH1D * nue_templates_m[480];
  TH1D * nue_templates_m_w[480];
  TH1D * nue_templates_e[480];
  TH1D * nue_templates_e_w[480];

  for(int i=0; i<480; i++) {
    CC_templates_m[i]    = (TH1D*)CC_hm->ProjectionY(Form("CC_m_bin%d",i+1),i+1,i+1);
    CC_templates_m_nc[i] = (TH1D*)CC_hm_nc->ProjectionY(Form("CC_nc_m_bin%d",i+1),i+1,i+1);
    CC_templates_e[i]    = (TH1D*)CC_he->ProjectionY(Form("CC_e_bin%d",i+1),i+1,i+1);
    nue_templates_m[i]    = (TH1D*)nue_hm->ProjectionY(Form("nue_m_bin%d",i+1),i+1,i+1);
    nue_templates_m_w[i]  = (TH1D*)nue_hm_w->ProjectionY(Form("nue_w_m_bin%d",i+1),i+1,i+1);
    nue_templates_e[i]    = (TH1D*)nue_he->ProjectionY(Form("nue_e_bin%d",i+1),i+1,i+1);
    nue_templates_e_w[i]  = (TH1D*)nue_he_w->ProjectionY(Form("nue_w_e_bin%d",i+1),i+1,i+1);
  }

  TH1D * CC_target_e_ft = (TH1D*)CC_f->Get(Form("e_h%s_w%d_ft",var,nuCut));
  TH1D * CC_target_m_ft = (TH1D*)CC_f->Get(Form("m_h%s_w%d_ft",var,nuCut));
  TH1D * nue_target_ft  = (TH1D*)nue_f->Get(Form("h%s%d_w_ft",var,EvCut));

  TFile *cov_f = new TFile(Form("/dune/app/users/qvuong/lownu/cov_matrix/%sReco_covmtr%d%d%d_10000.root",var,oscpar,nuCut,EvCut),"READ");
  TH2D *cov = (TH2D*)cov_f->Get("hcv");

  double cov_bins[301][301];
  for(int i=0; i<300; i++) {
    for(int j=0; j<300; j++) {
      cov_bins[i][j] = cov->GetBinContent(i+1, j+1);
    }
  }    

  TemplateFitter tf( CC_templates_m, CC_templates_m_nc, CC_templates_e, nue_templates_m, nue_templates_m_w, nue_templates_e, nue_templates_e_w, CC_target_m_ft, CC_target_e_ft, nue_target_ft );

  double energy_bins[481];
  for( int b = 0; b <= 480; ++b ) {
    energy_bins[b] = CC_he->GetXaxis()->GetBinLowEdge(b+1);
  }

  tf.setEnergyBins( energy_bins );
  tf.setCovmtr( cov_bins );
  tf.setPara( var, oscpar, nuCut, EvCut, seed );

  double bf_dm2, bf_Uee2, bf_Umm2;
  bool isOK = tf.doFit( bf_Uee2, bf_Umm2 , bf_dm2);
  printf( "nue Best-fit Uee2 = %f, Umm2 = %f, dm2 = %f\n", bf_Uee2, bf_Umm2, bf_dm2 );
  tf.Draw();
  tf.TrueDraw();
  }
  }
}



