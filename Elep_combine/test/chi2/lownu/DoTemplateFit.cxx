#include "TemplateFitter.cxx"
#include "TFile.h"
#include "TH2.h"

int main()
{
  TFile *f = new TFile("/nashome/q/qvuong/gen_data/output_2.root","READ");
  Int_t cutNo = 3;
  TH2D* hm = (TH2D*)f->Get(Form("m_hElepVsEv%d",cutNo));
  TH2D* hm_nc = (TH2D*)f->Get(Form("nc_m_hElepVsEv%d",cutNo));
  TH2D* he = (TH2D*)f->Get(Form("e_hElepVsEv%d",cutNo));
  TH1D * templates_m[280];
  TH1D * templates_m_nc[280];
  TH1D * templates_e[280];
  for(int i=0; i<280; i++) {
    templates_m[i] = (TH1D*)hm->ProjectionY(Form("m_bin%d",i+1),i+1,i+1);
    templates_m_nc[i] = (TH1D*)hm_nc->ProjectionY(Form("nc_m_bin%d",i+1),i+1,i+1);
    templates_e[i] = (TH1D*)he->ProjectionY(Form("e_bin%d",i+1),i+1,i+1);
  }
  //TH1D * intrinsic = (TH1D*)f->Get(Form("hLepE_sm%d",cutNo));
  TH1D * target_e = (TH1D*)f->Get(Form("e_hElep_w%d",cutNo));
  TH1D * target_m = (TH1D*)f->Get(Form("m_hElep_w%d",cutNo));

  TemplateFitter tf_e( templates_m, templates_m_nc, templates_e, target_e );
  //TemplateFitter tf_m( templates_m, templates_e, target_m );
  double energy_bins[281];
  for( int b = 0; b <= 280; ++b ) {
  energy_bins[b] = he->GetXaxis()->GetBinLowEdge(b+1);
  }
  
  tf_e.setEnergyBins( energy_bins );
  double bf_dm2, bf_Uee2, bf_Umm2;
  tf_e.Draw();
  //tf_e.TrueDraw();
}

