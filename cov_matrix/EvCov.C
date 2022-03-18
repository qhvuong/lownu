static const int N = 10000; // number of universes
static const int nbins = 52;
static const int nbins_Elep = 100;

const int n_mu = 19; // number of muon bins in covariance
const int n_e = 7; // number of electron bins in covariance

// bin edges -- fix these to be whatever they actually are
const double mubins[20] = {0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,7.,8.,12.,16.,20.,40.,100.};
const double ebins[8] = {0.,2.,4.,6.,8.,10.,20.,100.};

void EvCov()
{

  // covariance matrix
  TFile * covfile = new TFile( "/dune/app/users/qvuong/lownu_analysis/cov_matrix/total_covariance_DUNE_opt.root", "OLD" );
  TH2D * hcovmx = (TH2D*) covfile->Get( "total_covariance" );

  // only need the ND FHC part, which is the first 52 bins probably
  TMatrixD covmx( nbins, nbins );
  for( int x = 0; x < nbins; ++x ) {
    for( int y = 0; y < nbins; ++y ) {
      covmx[x][y] = hcovmx->GetBinContent( x+1, y+1 );

      //std::cout << x << "\t" << y << "\t" << covmx[x][y] << "\n";
    }
  }

    // cholesky decomposition
  TDecompChol decomp( covmx );
  if( !decomp.Decompose() ) {
    printf( "Main covariance matrix failed Cholesky decomposition\n" );
    return;
  }
  const TMatrixD chol = decomp.GetU();

  TRandom3 * rand = new TRandom3(12345);

  // make random number matrix
  TMatrixD scales( N, nbins );
  for( int i = 0; i < N; ++i ) {
    double mean = 0.;
    for( int j = 0; j < nbins; ++j ) {
      double val = rand->Gaus( 0., 1. );
      scales[i][j] = val;
      mean += val;
    }

    mean /= N;
    // force each column mean to be exactly 0, this just eliminates tiny statistical fluctuations in the mean weight
    for( int j = 0; j < nbins; ++j ) {
      double old = scales[i][j];
      scales[i][j] = old - mean;
    }

  }

  TMatrixD scaleCovars( nbins, nbins ); // pairwise covariance of columns of random numbers
  for( int i = 0; i < nbins; ++i ) { // columns
    for( int j = 0; j < nbins; ++j ) { // columns

      // compute column covariance
      double covar = 0.;
      for( int k = 0; k < N; ++k ) { // rows
        covar += (scales[k][i] * scales[k][j]); // means are 0 already by construction
      }
      scaleCovars[i][j] = covar/N;

      //std::cout << i << "\t" << j << "\t" << scaleCovars[i][j] << "\n";
    }
  }

  TDecompChol scaleDecomp( scaleCovars );
  if( !scaleDecomp.Decompose() ) printf( "Scale matrix didn't decompolse\n" );
  TMatrixD toInvert = scaleDecomp.GetU();
  TMatrixD inverse = toInvert.Invert();
  scales *= inverse;

  scales *= chol;

  int cutEv = 2;
  TFile *f = new TFile("/dune/app/users/qvuong/lownu_analysis/gen_data/CC/output_2.root");
  TFile *f_nue = new TFile("/dune/app/users/qvuong/lownu_analysis/gen_data/nuescattering/nue_output_2.root");
  TH2D *CC_m0 = (TH2D*)f->Get("m_hEvRecoVsEv0_cov");
  TH2D *CC_m3 = (TH2D*)f->Get("m_hEvRecoVsEv3_cov");
  TH2D *CC_e0 = (TH2D*)f->Get("e_hEvRecoVsEv0_cov");
  TH2D *CC_e3 = (TH2D*)f->Get("e_hEvRecoVsEv3_cov");
  TH2D *nue_m = (TH2D*)f_nue->Get(Form("m_hEvRecoVsEv%d_cov",cutEv));
  TH2D *nue_e = (TH2D*)f_nue->Get(Form("e_hEvRecoVsEv%d_cov",cutEv));
  TH1D *tp_m0[n_mu];
  TH1D *tp_m3[n_mu];
  TH1D *tp_e0[n_e];
  TH1D *tp_e3[n_e];
  TH1D *tp_m_nue[n_mu];
  TH1D *tp_e_nue[n_e];

  std::cout << CC_m0->GetNbinsY() << "\n";

  for(int mb=0; mb<n_mu; mb++){
    tp_m0[mb] = (TH1D*)CC_m0->ProjectionY(Form("m_bin0%d",mb+1),mb+1,mb+1);
    tp_m3[mb] = (TH1D*)CC_m3->ProjectionY(Form("m_bin3%d",mb+1),mb+1,mb+1);
    tp_m_nue[mb] = (TH1D*)nue_m->ProjectionY(Form("m_bin_nue%d",mb+1),mb+1,mb+1);
  }
  for(int eb=0; eb<n_e; eb++){
    tp_e0[eb] = (TH1D*)CC_e0->ProjectionY(Form("e_bin0%d",eb+1),eb+1,eb+1);
    tp_e3[eb] = (TH1D*)CC_e3->ProjectionY(Form("e_bin3%d",eb+1),eb+1,eb+1);
    tp_e_nue[eb] = (TH1D*)nue_e->ProjectionY(Form("e_bin_nue%d",eb+1),eb+1,eb+1);
  }

  TH1D *m0 = new TH1D("m0","",100,0,120);
  TH1D *m3 = new TH1D("m3","",100,0,120);
  TH1D *e0 = new TH1D("e0","",100,0,120);
  TH1D *e3 = new TH1D("e3","",100,0,120);
  TH1D *m_nue = new TH1D("m_nue","",100,0,120);
  TH1D *e_nue = new TH1D("e_nue","",100,0,120);
  TH1D *nue = new TH1D("nue","",100,0,120);
  TH1D *CC_m0_nom = new TH1D("CC_m0_nom","",100,0,120);
  TH1D *CC_m3_nom = new TH1D("CC_m3_nom3","",100,0,120);
  TH1D *CC_e0_nom = new TH1D("CC_e0_nom","",100,0,120);
  TH1D *CC_e3_nom = new TH1D("CC_e3_nom","",100,0,120);
  TH1D *nue_m_nom = new TH1D("nue_m_nom","",100,0,120);
  TH1D *nue_e_nom = new TH1D("nue_e_nom","",100,0,120);
  TH1D *nue_nom = new TH1D("nue_nom","",100,0,120);

  TMatrixD Elep_m0(N, nbins_Elep);
  TMatrixD Elep_m3(N, nbins_Elep);
  TMatrixD Elep_e0(N, nbins_Elep);
  TMatrixD Elep_e3(N, nbins_Elep);
  TMatrixD Elep_m_nue(N, nbins_Elep);
  TMatrixD Elep_e_nue(N, nbins_Elep);
  TMatrixD Elep_nue(N, nbins_Elep);
/*
  TCanvas *c = new TCanvas("c","",800,800);
  //TCanvas *c_m3 = new TCanvas("c_m3","",800,800);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.4, 1, 1.0);
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.4);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
*/
  for( int u = 0; u < N; ++u ) {
    m0->Reset();
    m3->Reset();
    e0->Reset();
    e3->Reset();
    m_nue->Reset();
    e_nue->Reset();
    nue->Reset();
    CC_m0_nom->Reset();
    CC_m3_nom->Reset();
    CC_e0_nom->Reset();
    CC_e3_nom->Reset();
    nue_m_nom->Reset();
    nue_e_nom->Reset();
    nue_nom->Reset();
    
    for(int mb=0; mb<n_mu; mb++) {
      int fluxbin = mb+1;
      double evtwgt = scales[u][fluxbin];

      m0->Add(tp_m0[mb],1.+evtwgt);
      m3->Add(tp_m3[mb],1.+evtwgt);
      m_nue->Add(tp_m_nue[mb],1.+evtwgt);

      CC_m0_nom->Add(tp_m0[mb], 1.);
      CC_m3_nom->Add(tp_m3[mb], 1.);
      nue_m_nom->Add(tp_m_nue[mb],1.);
    }
    for(int eb=0; eb<n_e; eb++) {
      int fluxbin = 38+eb+1;
      double evtwgt = scales[u][fluxbin];

      e0->Add(tp_e0[eb],1.+evtwgt);
      e3->Add(tp_e3[eb],1.+evtwgt);
      e_nue->Add(tp_e_nue[eb],1.+evtwgt);

      CC_e0_nom->Add(tp_e0[eb], 1.);
      CC_e3_nom->Add(tp_e3[eb], 1.);
      nue_e_nom->Add(tp_e_nue[eb],1.);
    }
    
    nue->Add(m_nue); nue->Add(e_nue);
    nue_nom->Add(nue_m_nom); nue_nom->Add(nue_e_nom);
    //std::cout << m0->GetBinLowEdge(100) << "\n";

    for(int i=0; i<100; i++){
      Elep_m0[u][i] = m0->GetBinContent(i+1);
      Elep_m3[u][i] = m3->GetBinContent(i+1);
      Elep_e0[u][i] = e0->GetBinContent(i+1);
      Elep_e3[u][i] = e3->GetBinContent(i+1);
      Elep_m_nue[u][i] = m_nue->GetBinContent(i+1);
      Elep_e_nue[u][i] = e_nue->GetBinContent(i+1);
      Elep_nue[u][i] = nue->GetBinContent(i+1);
    }
    //std::cout << CC_m0_nom->GetBinContent(48) << "\t" << Elep_CC_m0_nom[u][47] << "\n";
/*
    if(u<20) {

      c->cd();
      pad1->cd();
      TH1D *CC_m0 = (TH1D*)nue->Clone();
      CC_m0->SetTitle("Nue Elep in 20 universes (nu<10.0)");
      CC_m0->GetYaxis()->SetTitle("Elep histograms in 20 universes");
      //CC_m0->GetYaxis()->SetTitleSize(20);
      CC_m0->SetLineColor(u+1);
      CC_m0->SetStats(0);
      CC_m0->Draw("same");
      
      c->cd();
      pad2->cd();
      TH1D *ratio_m0 = (TH1D*)nue->Clone("ratio");
      ratio_m0->SetLineColor(u+1);
      ratio_m0->SetMaximum(1.7);
      ratio_m0->SetMinimum(0.3);
      ratio_m0->Sumw2(); 
      ratio_m0->SetStats(0); 
      ratio_m0->Divide(nue_nom);
      ratio_m0->SetTitle("");
      ratio_m0->GetYaxis()->SetTitle("ratio Elep/Elep_nom");
      ratio_m0->GetXaxis()->SetTitle("Elep (GeV)");
      ratio_m0->Draw("same");

      c_m3->cd();
      pad1->cd();
      TH1D *CC_m3 = (TH1D*)m3->Clone();
      CC_m3->SetTitle("CC_m Elep in 20 universes (nu<0.3)");
      CC_m3->GetYaxis()->SetTitle("Elep histograms in 20 universes");
      //CC_m0->GetYaxis()->SetTitleSize(20);
      CC_m3->SetLineColor(u+1);
      CC_m3->SetStats(0);
      CC_m3->Draw("same");
      
      c_m3->cd();
      pad2->cd();
      TH1D *ratio_m3 = (TH1D*)m3->Clone("ratio");
      ratio_m3->SetLineColor(u+1);
      ratio_m3->SetMaximum(1.7);
      ratio_m3->SetMinimum(0.3);
      ratio_m3->Sumw2(); 
      ratio_m3->SetStats(0); 
      ratio_m3->Divide(CC_m3_nom);
      ratio_m3->SetTitle("");
      ratio_m3->GetYaxis()->SetTitle("ratio Elep/Elep_nom");
      ratio_m3->GetXaxis()->SetTitle("Elep (GeV)");
      ratio_m3->Draw("same");

    }
    else continue;
*/
  }

  //c->SaveAs("ratio_nue.png");
  //c_m3->SaveAs("ratio_CC_m3.png");

  TMatrixD ElepCovars_m0( nbins_Elep, nbins_Elep );
  TMatrixD ElepCovars_m3( nbins_Elep, nbins_Elep );
  TMatrixD ElepCovars_e0( nbins_Elep, nbins_Elep );
  TMatrixD ElepCovars_e3( nbins_Elep, nbins_Elep );
  TMatrixD ElepCovars_m_nue( nbins_Elep, nbins_Elep );
  TMatrixD ElepCovars_e_nue( nbins_Elep, nbins_Elep );
  TMatrixD ElepCovars_nue( nbins_Elep, nbins_Elep );

  TMatrixD ElepCovars_m0e0( nbins_Elep, nbins_Elep );
  TMatrixD ElepCovars_m0nue( nbins_Elep, nbins_Elep );
  TMatrixD ElepCovars_e0m0( nbins_Elep, nbins_Elep );
  TMatrixD ElepCovars_e0nue( nbins_Elep, nbins_Elep );
  TMatrixD ElepCovars_nuem0( nbins_Elep, nbins_Elep );
  TMatrixD ElepCovars_nuee0( nbins_Elep, nbins_Elep );

  TMatrixD ElepCovars_m3e3( nbins_Elep, nbins_Elep );
  TMatrixD ElepCovars_m3nue( nbins_Elep, nbins_Elep );
  TMatrixD ElepCovars_e3m3( nbins_Elep, nbins_Elep );
  TMatrixD ElepCovars_e3nue( nbins_Elep, nbins_Elep );
  TMatrixD ElepCovars_nuem3( nbins_Elep, nbins_Elep );
  TMatrixD ElepCovars_nuee3( nbins_Elep, nbins_Elep );

  TMatrixD ElepCorrel_m0( nbins_Elep, nbins_Elep );
  TMatrixD ElepCorrel_m3( nbins_Elep, nbins_Elep );
  TMatrixD ElepCorrel_e0( nbins_Elep, nbins_Elep );
  TMatrixD ElepCorrel_e3( nbins_Elep, nbins_Elep );
  TMatrixD ElepCorrel_nue( nbins_Elep, nbins_Elep );

  TMatrixD ElepCorrel_m0e0( nbins_Elep, nbins_Elep );
  TMatrixD ElepCorrel_m0nue( nbins_Elep, nbins_Elep );
  TMatrixD ElepCorrel_e0m0( nbins_Elep, nbins_Elep );
  TMatrixD ElepCorrel_e0nue( nbins_Elep, nbins_Elep );
  TMatrixD ElepCorrel_nuem0( nbins_Elep, nbins_Elep );
  TMatrixD ElepCorrel_nuee0( nbins_Elep, nbins_Elep );

  TMatrixD ElepCorrel_m3e3( nbins_Elep, nbins_Elep );
  TMatrixD ElepCorrel_m3nue( nbins_Elep, nbins_Elep );
  TMatrixD ElepCorrel_e3m3( nbins_Elep, nbins_Elep );
  TMatrixD ElepCorrel_e3nue( nbins_Elep, nbins_Elep );
  TMatrixD ElepCorrel_nuem3( nbins_Elep, nbins_Elep );
  TMatrixD ElepCorrel_nuee3( nbins_Elep, nbins_Elep );

  for( int i = 0; i < nbins_Elep; ++i ) { // columns
    for( int j = 0; j < nbins_Elep; ++j ) { // columns
      // compute column covariance
      double covar_m0 = 0.;
      double covar_m3 = 0.;
      double covar_e0 = 0.;
      double covar_e3 = 0.;
      double covar_m_nue = 0.;
      double covar_e_nue = 0.;
      double covar_nue = 0.;
      double covar_m0e0 = 0.;
      double covar_m0nue = 0.;
      double covar_e0m0 = 0.;
      double covar_e0nue = 0.;
      double covar_nuem0 = 0.;
      double covar_nuee0 = 0.;
      double covar_m3e3 = 0.;
      double covar_m3nue = 0.;
      double covar_e3m3 = 0.;
      double covar_e3nue = 0.;
      double covar_nuem3 = 0.;
      double covar_nuee3 = 0.;

      double var_m0_i = 0.;
      double var_m3_i = 0.;
      double var_m0_j = 0.;
      double var_m3_j = 0.;
      double var_e0_i = 0.;
      double var_e3_i = 0.;
      double var_e0_j = 0.;
      double var_e3_j = 0.;
      double var_nue_i = 0.;
      double var_nue_j = 0.;

      for( int k = 0; k < N; ++k ) { // rows
        covar_m0 += (Elep_m0[k][i] - CC_m0_nom->GetBinContent(i+1)) * (Elep_m0[k][j] - CC_m0_nom->GetBinContent(j+1));
        covar_m3 += (Elep_m3[k][i] - CC_m3_nom->GetBinContent(i+1)) * (Elep_m3[k][j] - CC_m3_nom->GetBinContent(j+1));
        covar_e0 += (Elep_e0[k][i] - CC_e0_nom->GetBinContent(i+1)) * (Elep_e0[k][j] - CC_e0_nom->GetBinContent(j+1)); 
        covar_e3 += (Elep_e3[k][i] - CC_e3_nom->GetBinContent(i+1)) * (Elep_e3[k][j] - CC_e3_nom->GetBinContent(j+1)); 
        covar_m_nue += (Elep_m_nue[k][i] - nue_m_nom->GetBinContent(i+1)) * (Elep_m_nue[k][j] - nue_m_nom->GetBinContent(j+1));
        covar_e_nue += (Elep_e_nue[k][i] - nue_e_nom->GetBinContent(i+1)) * (Elep_e_nue[k][j] - nue_e_nom->GetBinContent(j+1)); 
        covar_nue += (Elep_nue[k][i] - nue_nom->GetBinContent(i+1)) * (Elep_nue[k][j] - nue_nom->GetBinContent(j+1)); 

        covar_m0e0 += (Elep_m0[k][i] - CC_m0_nom->GetBinContent(i+1)) * (Elep_e0[k][j] - CC_e0_nom->GetBinContent(j+1));
        covar_m0nue += (Elep_m0[k][i] - CC_m0_nom->GetBinContent(i+1)) * (Elep_nue[k][j] - nue_nom->GetBinContent(j+1));
        covar_e0m0 += (Elep_e0[k][i] - CC_e0_nom->GetBinContent(i+1)) * (Elep_m0[k][j] - CC_m0_nom->GetBinContent(j+1)); 
        covar_e0nue += (Elep_e0[k][i] - CC_e0_nom->GetBinContent(i+1)) * (Elep_nue[k][j] - nue_nom->GetBinContent(j+1)); 
        covar_nuem0 += (Elep_nue[k][i] - nue_nom->GetBinContent(i+1)) * (Elep_m0[k][j] - CC_m0_nom->GetBinContent(j+1)); 
        covar_nuee0 += (Elep_nue[k][i] - nue_nom->GetBinContent(i+1)) * (Elep_e0[k][j] - CC_e0_nom->GetBinContent(j+1)); 

        covar_m3e3 += (Elep_m3[k][i] - CC_m3_nom->GetBinContent(i+1)) * (Elep_e3[k][j] - CC_e3_nom->GetBinContent(j+1));
        covar_m3nue += (Elep_m3[k][i] - CC_m3_nom->GetBinContent(i+1)) * (Elep_nue[k][j] - nue_nom->GetBinContent(j+1));
        covar_e3m3 += (Elep_e3[k][i] - CC_e3_nom->GetBinContent(i+1)) * (Elep_m3[k][j] - CC_m3_nom->GetBinContent(j+1)); 
        covar_e3nue += (Elep_e3[k][i] - CC_e3_nom->GetBinContent(i+1)) * (Elep_nue[k][j] - nue_nom->GetBinContent(j+1)); 
        covar_nuem3 += (Elep_nue[k][i] - nue_nom->GetBinContent(i+1)) * (Elep_m3[k][j] - CC_m3_nom->GetBinContent(j+1)); 
        covar_nuee3 += (Elep_nue[k][i] - nue_nom->GetBinContent(i+1)) * (Elep_e3[k][j] - CC_e3_nom->GetBinContent(j+1)); 

        
        var_m0_i += (Elep_m0[k][i] - CC_m0_nom->GetBinContent(i+1)) * (Elep_m0[k][i] - CC_m0_nom->GetBinContent(i+1));
        var_m3_i += (Elep_m3[k][i] - CC_m3_nom->GetBinContent(i+1)) * (Elep_m3[k][i] - CC_m3_nom->GetBinContent(i+1));
        var_m0_j += (Elep_m0[k][j] - CC_m0_nom->GetBinContent(j+1)) * (Elep_m0[k][j] - CC_m0_nom->GetBinContent(j+1));
        var_m3_j += (Elep_m3[k][j] - CC_m3_nom->GetBinContent(j+1)) * (Elep_m3[k][j] - CC_m3_nom->GetBinContent(j+1));

        var_e0_i += (Elep_e0[k][i] - CC_e0_nom->GetBinContent(i+1)) * (Elep_e0[k][i] - CC_e0_nom->GetBinContent(i+1));
        var_e3_i += (Elep_e3[k][i] - CC_e3_nom->GetBinContent(i+1)) * (Elep_e3[k][i] - CC_e3_nom->GetBinContent(i+1));
        var_e0_j += (Elep_e0[k][j] - CC_e0_nom->GetBinContent(j+1)) * (Elep_e0[k][j] - CC_e0_nom->GetBinContent(j+1));
        var_e3_j += (Elep_e3[k][j] - CC_e3_nom->GetBinContent(j+1)) * (Elep_e3[k][j] - CC_e3_nom->GetBinContent(j+1));

        var_nue_i += (Elep_nue[k][i] - nue_nom->GetBinContent(i+1)) * (Elep_nue[k][i] - nue_nom->GetBinContent(i+1));
        var_nue_j += (Elep_nue[k][j] - nue_nom->GetBinContent(j+1)) * (Elep_nue[k][j] - nue_nom->GetBinContent(j+1));
      }

      ElepCovars_m0[i][j] = covar_m0/(N * CC_m0_nom->GetBinContent(i+1) * CC_m0_nom->GetBinContent(j+1));
      ElepCovars_m3[i][j] = covar_m3/(N * CC_m3_nom->GetBinContent(i+1) * CC_m3_nom->GetBinContent(j+1));
      ElepCovars_e0[i][j] = covar_e0/(N * CC_e0_nom->GetBinContent(i+1) * CC_e0_nom->GetBinContent(j+1));
      ElepCovars_e3[i][j] = covar_e3/(N * CC_e3_nom->GetBinContent(i+1) * CC_e3_nom->GetBinContent(j+1));
      ElepCovars_m_nue[i][j] = covar_m_nue/(N * nue_m_nom->GetBinContent(i+1) * nue_m_nom->GetBinContent(j+1));
      ElepCovars_e_nue[i][j] = covar_e_nue/(N * nue_e_nom->GetBinContent(i+1) * nue_e_nom->GetBinContent(j+1));
      ElepCovars_nue[i][j] = covar_nue/(N * nue_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j+1));

      ElepCovars_m0e0[i][j] = covar_m0e0/(N * CC_m0_nom->GetBinContent(i+1) * CC_e0_nom->GetBinContent(j+1));
      ElepCovars_m0nue[i][j] = covar_m0nue/(N * CC_m0_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j+1));
      ElepCovars_e0m0[i][j] = covar_e0m0/(N * CC_e0_nom->GetBinContent(i+1) * CC_m0_nom->GetBinContent(j+1));
      ElepCovars_e0nue[i][j] = covar_e0nue/(N * CC_e0_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j+1));
      ElepCovars_nuem0[i][j] = covar_nuem0/(N * nue_nom->GetBinContent(i+1) * CC_m0_nom->GetBinContent(j+1));
      ElepCovars_nuee0[i][j] = covar_nuee0/(N * nue_nom->GetBinContent(i+1) * CC_e0_nom->GetBinContent(j+1));

      ElepCovars_m3e3[i][j] = covar_m3e3/(N * CC_m3_nom->GetBinContent(i+1) * CC_e3_nom->GetBinContent(j+1));
      ElepCovars_m3nue[i][j] = covar_m3nue/(N * CC_m3_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j+1));
      ElepCovars_e3m3[i][j] = covar_e3m3/(N * CC_e3_nom->GetBinContent(i+1) * CC_m3_nom->GetBinContent(j+1));
      ElepCovars_e3nue[i][j] = covar_e3nue/(N * CC_e3_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j+1));
      ElepCovars_nuem3[i][j] = covar_nuem3/(N * nue_nom->GetBinContent(i+1) * CC_m3_nom->GetBinContent(j+1));
      ElepCovars_nuee3[i][j] = covar_nuee3/(N * nue_nom->GetBinContent(i+1) * CC_e3_nom->GetBinContent(j+1));


      ElepCorrel_m0[i][j] = covar_m0/sqrt(var_m0_i * var_m0_j);
      ElepCorrel_m3[i][j] = covar_m3/sqrt(var_m3_i * var_m3_j);
      ElepCorrel_e0[i][j] = covar_e0/sqrt(var_e0_i * var_e0_j);
      ElepCorrel_e3[i][j] = covar_e3/sqrt(var_e3_i * var_e3_j);
      ElepCorrel_nue[i][j] = covar_nue/sqrt(var_nue_i * var_nue_j);

      ElepCorrel_m0e0[i][j] = covar_m0e0/sqrt(var_m0_i * var_e0_j);
      ElepCorrel_m0nue[i][j] = covar_m0nue/sqrt(var_m0_i * var_nue_j);
      ElepCorrel_e0m0[i][j] = covar_e0m0/sqrt(var_e0_i * var_m0_j);
      ElepCorrel_e0nue[i][j] = covar_e0nue/sqrt(var_e0_i * var_nue_j);
      ElepCorrel_nuem0[i][j] = covar_nuem0/sqrt(var_nue_i * var_m0_j);
      ElepCorrel_nuee0[i][j] = covar_nuee0/sqrt(var_nue_i * var_e0_j);

      ElepCorrel_m3e3[i][j] = covar_m3e3/sqrt(var_m3_i * var_e3_j);
      ElepCorrel_m3nue[i][j] = covar_m3nue/sqrt(var_m3_i * var_nue_j);
      ElepCorrel_e3m3[i][j] = covar_e3m3/sqrt(var_e3_i * var_m3_j);
      ElepCorrel_e3nue[i][j] = covar_e3nue/sqrt(var_e3_i * var_nue_j);
      ElepCorrel_nuem3[i][j] = covar_nuem3/sqrt(var_nue_i * var_m3_j);
      ElepCorrel_nuee3[i][j] = covar_nuee3/sqrt(var_nue_i * var_e3_j);
    }
  }

  TMatrixD ElepCovars0( 3*nbins_Elep, 3*nbins_Elep );
  TMatrixD ElepCovars3( 3*nbins_Elep, 3*nbins_Elep );
  TMatrixD ElepCorrel0( 3*nbins_Elep, 3*nbins_Elep );
  TMatrixD ElepCorrel3( 3*nbins_Elep, 3*nbins_Elep );

  for(int i = 0; i < nbins_Elep; i++) {
    for(int j = 0; j < nbins_Elep; j++) {
      ElepCovars0[i][j] = ElepCovars_m0[i][j];
      ElepCovars0[i][j+nbins_Elep] = ElepCovars_m0e0[i][j];
      ElepCovars0[i][j+2*nbins_Elep] = ElepCovars_m0nue[i][j];
      ElepCovars3[i][j] = ElepCovars_m3[i][j];
      ElepCovars3[i][j+nbins_Elep] = ElepCovars_m3e3[i][j];
      ElepCovars3[i][j+2*nbins_Elep] = ElepCovars_m3nue[i][j];

      ElepCovars0[i+nbins_Elep][j] = ElepCovars_e0m0[i][j];
      ElepCovars0[i+nbins_Elep][j+nbins_Elep] = ElepCovars_e0[i][j];
      ElepCovars0[i+nbins_Elep][j+2*nbins_Elep] = ElepCovars_e0nue[i][j];
      ElepCovars3[i+nbins_Elep][j] = ElepCovars_e3m3[i][j];
      ElepCovars3[i+nbins_Elep][j+nbins_Elep] = ElepCovars_e3[i][j];
      ElepCovars3[i+nbins_Elep][j+2*nbins_Elep] = ElepCovars_e3nue[i][j];

      ElepCovars0[i+2*nbins_Elep][j] = ElepCovars_nuem0[i][j];
      ElepCovars0[i+2*nbins_Elep][j+nbins_Elep] = ElepCovars_nuee0[i][j];
      ElepCovars0[i+2*nbins_Elep][j+2*nbins_Elep] = ElepCovars_nue[i][j];
      ElepCovars3[i+2*nbins_Elep][j] = ElepCovars_nuem3[i][j];
      ElepCovars3[i+2*nbins_Elep][j+nbins_Elep] = ElepCovars_nuee3[i][j];
      ElepCovars3[i+2*nbins_Elep][j+2*nbins_Elep] = ElepCovars_nue[i][j];


      ElepCorrel0[i][j] = ElepCorrel_m0[i][j];
      ElepCorrel0[i][j+nbins_Elep] = ElepCorrel_m0e0[i][j];
      ElepCorrel0[i][j+2*nbins_Elep] = ElepCorrel_m0nue[i][j];
      ElepCorrel3[i][j] = ElepCorrel_m3[i][j];
      ElepCorrel3[i][j+nbins_Elep] = ElepCorrel_m3e3[i][j];
      ElepCorrel3[i][j+2*nbins_Elep] = ElepCorrel_m3nue[i][j];

      ElepCorrel0[i+nbins_Elep][j] = ElepCorrel_e0m0[i][j];
      ElepCorrel0[i+nbins_Elep][j+nbins_Elep] = ElepCorrel_e0[i][j];
      ElepCorrel0[i+nbins_Elep][j+2*nbins_Elep] = ElepCorrel_e0nue[i][j];
      ElepCorrel3[i+nbins_Elep][j] = ElepCorrel_e3m3[i][j];
      ElepCorrel3[i+nbins_Elep][j+nbins_Elep] = ElepCorrel_e3[i][j];
      ElepCorrel3[i+nbins_Elep][j+2*nbins_Elep] = ElepCorrel_e3nue[i][j];

      ElepCorrel0[i+2*nbins_Elep][j] = ElepCorrel_nuem0[i][j];
      ElepCorrel0[i+2*nbins_Elep][j+nbins_Elep] = ElepCorrel_nuee0[i][j];
      ElepCorrel0[i+2*nbins_Elep][j+2*nbins_Elep] = ElepCorrel_nue[i][j];
      ElepCorrel3[i+2*nbins_Elep][j] = ElepCorrel_nuem3[i][j];
      ElepCorrel3[i+2*nbins_Elep][j+nbins_Elep] = ElepCorrel_nuee3[i][j];
      ElepCorrel3[i+2*nbins_Elep][j+2*nbins_Elep] = ElepCorrel_nue[i][j];
    }
  }

      
  TH2D *h_m0 = new TH2D(ElepCovars_m0);
  TH2D *h_m3 = new TH2D(ElepCovars_m3);
  TH2D *h_e0 = new TH2D(ElepCovars_e0);
  TH2D *h_e3 = new TH2D(ElepCovars_e3);
  TH2D *h_m_nue = new TH2D(ElepCovars_m_nue);
  TH2D *h_e_nue = new TH2D(ElepCovars_e_nue);
  TH2D *h_nue = new TH2D(ElepCovars_nue);
  //TH2D *h0 = new TH2D(ElepCovars0);
  //TH2D *h3 = new TH2D(ElepCovars3);

  TH2D *hcr_m0 = new TH2D(ElepCorrel_m0);
  TH2D *hcr_m3 = new TH2D(ElepCorrel_m3);
  TH2D *hcr_e0 = new TH2D(ElepCorrel_e0);
  TH2D *hcr_e3 = new TH2D(ElepCorrel_e3);
  TH2D *hcr_nue = new TH2D(ElepCorrel_nue);
  //TH2D *hcr0 = new TH2D(ElepCorrel0);
  //TH2D *hcr3 = new TH2D(ElepCorrel3);

  TH2D *h0 = new TH2D("h0","",300,0,300,300,0,300);
  TH2D *h3 = new TH2D("h3","",300,0,300,300,0,300);
  TH2D *hcr0 = new TH2D("hcr0","",300,0,300,300,0,300);
  TH2D *hcr3 = new TH2D("hcr3","",300,0,300,300,0,300);
  for(int i=0; i<300; i++) {
    for(int j=0; j<300; j++) {
      h0->SetBinContent(i+1, j+1, ElepCovars0[i][j]);
      h3->SetBinContent(i+1, j+1, ElepCovars3[i][j]);

      hcr0->SetBinContent(i+1, j+1, ElepCorrel0[i][j]);
      hcr3->SetBinContent(i+1, j+1, ElepCorrel3[i][j]);
    }
  }
 
  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();

  h0->SetStats(0);
  h0->SetTitle("EvReco Covariance (nu<10.0)");
  h3->SetStats(0);
  h3->SetTitle("EvReco Covariance (nu<0.3)");

  hcr0->SetStats(0);
  hcr0->SetTitle("EvReco Correlation (nu<10.0)");
  hcr3->SetStats(0);
  hcr3->SetTitle("EvReco Correlation (nu<0.3)");

  TCanvas *c1 = new TCanvas("c1","",1800,1350);
  c1->Divide(2,1);
  c1->cd(1);
  h0->Draw("colz");
  c1->cd(2);
  h3->Draw("colz");
  c1->SaveAs(Form("Ev_total_Cov%d_%d.png",cutEv,N));

  const Int_t Number = 3;
  Double_t Red[Number]    = { 0.00, 1.00, 1.00};
  Double_t Green[Number]  = { 0.00, 1.00, 0.00};
  Double_t Blue[Number]   = { 1.00, 1.00, 0.00};
  Double_t Length[Number] = { 0.00, 0.50, 1.00 };
  Int_t nb=50;
  TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);

  TCanvas *ccr1 = new TCanvas("ccr1","",1800,1350);
  ccr1->Divide(2,1);
  ccr1->cd(1);
  hcr0->GetZaxis()->SetRangeUser(-1., 1.);
  hcr0->Draw("colz");
  ccr1->cd(2);
  hcr3->GetZaxis()->SetRangeUser(-1., 1.);
  hcr3->Draw("colz");
  ccr1->SaveAs(Form("Ev_total_Cor%d_%d.png",cutEv,N));

  TFile *out = new TFile(Form("/dune/app/users/qvuong/lownu_analysis/cov_matrix/Ev_cov%d_%d.root",cutEv,N),"RECREATE");
  h0->Write();
  h3->Write();
  hcr0->Write();
  hcr3->Write();
  out->Close();

}

/*
  TDecompSVD svd0(ElepCovars0);
  //svd0.Decompose();
  TMatrixD inv0 = svd0.Invert();
  TDecompSVD svd3(ElepCovars3);
  //svd3.Decompose();
  TMatrixD inv3 = svd3.Invert();

  TH2D *test0 = new TH2D(inv0);
  TH2D *test3 = new TH2D(inv3);
  TCanvas *c2 = new TCanvas("c2","",1800,1350);
  c2->Divide(2,1);
  c2->cd(1);
  test0->Draw("colz");
  c2->cd(2);
  test3->Draw("colz");
  c2->SaveAs(Form("test1_total_Cov_%d.png",N));
*/
