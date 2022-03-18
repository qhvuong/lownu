#include "TemplateFitter.cxx"
#include "TFile.h"
#include "TH2.h"

int main()
{
  TFile *f = new TFile("/nashome/q/qvuong/gen_data/nuescattering/nue_output_2.root","READ");
  Int_t cutNo = 0;
  TH2D* hm = (TH2D*)f->Get(Form("m_hElepVsEv%d",cutNo));
  TH2D* hm_w = (TH2D*)f->Get(Form("m_hElepVsEv%d_w",cutNo));
  TH2D* he = (TH2D*)f->Get(Form("e_hElepVsEv%d",cutNo));
  TH2D* he_w = (TH2D*)f->Get(Form("e_hElepVsEv%d_w",cutNo));
  TH1D * templates_m[480];
  TH1D * templates_m_w[480];
  TH1D * templates_e[480];
  TH1D * templates_e_w[480];
  for(int i=0; i<480; i++) {
    templates_m[i] = (TH1D*)hm->ProjectionY(Form("m_bin%d",i+1),i+1,i+1);
    templates_m_w[i] = (TH1D*)hm_w->ProjectionY(Form("w_m_bin%d",i+1),i+1,i+1);
    templates_e[i] = (TH1D*)he->ProjectionY(Form("e_bin%d",i+1),i+1,i+1);
    templates_e_w[i] = (TH1D*)he_w->ProjectionY(Form("w_e_bin%d",i+1),i+1,i+1);
  }
  //TH1D * intrinsic = (TH1D*)f->Get(Form("hLepE_sm%d",cutNo));
  TH1D * target = (TH1D*)f->Get(Form("hElep_w%d",cutNo));
  TH1D * target_ft = (TH1D*)f->Get(Form("hElep_w%d_ft",cutNo));

  TemplateFitter tf( templates_m, templates_m_w, templates_e, templates_e_w, target );
  //TemplateFitter tf_m( templates_m, templates_e, target_m );
  double energy_bins[481];
  for( int b = 0; b <= 480; ++b ) {
    energy_bins[b] = he->GetXaxis()->GetBinLowEdge(b+1);
  }

  tf.setEnergyBins( energy_bins );
  double bf_dm2, bf_Uee2, bf_Umm2;
  tf.Draw();
  //tf.TrueDraw();
}




