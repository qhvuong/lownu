#include "TemplateFitter.cxx"
#include "TFile.h"
#include "TH2.h"
#include "TStyle.h"

int main()
{
  int para[3];
  para[0] = 2;
  //for(para[0]=1; para[0]<3; para[0]++) {
  TFile *CC_f = new TFile(Form("/dune/app/users/qvuong/lownu/gen_data/CC/output_%d.root",para[0]),"READ");
  TFile *nue_f = new TFile(Form("/dune/app/users/qvuong/lownu/gen_data/nuescattering/nue_output_%d.root",para[0]),"READ");
  //for(para[1]=0; para[1]<4; para[1]++) {
  //if(para[1] != 0 && para[1] != 3) continue;
  para[1] = 3;

  //for(para[2]=2; para[2]<3; para[2]++) {
  para[2] = 0;
  Int_t cutNu = para[1];
  Int_t cutEv = para[2];
  std::cout << para[0] << "\t" << para[1] << "\t" << para[2] << "\n";

  TH2D* CC_hm = (TH2D*)CC_f->Get(Form("m_hElepVsEv%d",cutNu));
  TH2D* CC_hm_nc = (TH2D*)CC_f->Get(Form("nc_m_hElepVsEv%d",cutNu));
  TH2D* CC_he = (TH2D*)CC_f->Get(Form("e_hElepVsEv%d",cutNu));
  TH2D* nue_hm = (TH2D*)nue_f->Get(Form("m_hElepVsEv%d",cutEv));
  TH2D* nue_hm_w = (TH2D*)nue_f->Get(Form("m_hElepVsEv%d_w",cutEv));
  TH2D* nue_he = (TH2D*)nue_f->Get(Form("e_hElepVsEv%d",cutEv));
  TH2D* nue_he_w = (TH2D*)nue_f->Get(Form("e_hElepVsEv%d_w",cutEv));

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

  TH1D * CC_target_e = (TH1D*)CC_f->Get(Form("e_hElep_w%d",cutNu));
  TH1D * CC_target_m = (TH1D*)CC_f->Get(Form("m_hElep_w%d",cutNu));
  TH1D * nue_target  = (TH1D*)nue_f->Get(Form("hElep%d_w",cutEv));
  TH1D * CC_target_e_ft = (TH1D*)CC_f->Get(Form("e_hElep_w%d_ft",cutNu));
  TH1D * CC_target_m_ft = (TH1D*)CC_f->Get(Form("m_hElep_w%d_ft",cutNu));
  TH1D * nue_target_ft  = (TH1D*)nue_f->Get(Form("hElep%d_w_ft",cutEv));

  TemplateFitter tf( CC_templates_m, CC_templates_m_nc, CC_templates_e, nue_templates_m, nue_templates_m_w, nue_templates_e, nue_templates_e_w, CC_target_m_ft, CC_target_e_ft, nue_target_ft );

  double energy_bins[481];
  for( int b = 0; b <= 480; ++b ) {
    energy_bins[b] = CC_he->GetXaxis()->GetBinLowEdge(b+1);
  }

  tf.setEnergyBins( energy_bins );
  tf.setPara( para );

  tf.Draw();
  //}
  //}
  //}
}



