void test()
{
  TMatrixD a(2,1);
  TMatrixD b(2,2);

  a[0][0] = 0;
  a[0][1] = 1;
  //a[0][1] = 1.;

  b[0][0] = 1.;
  b[0][1] = 0.;
  b[1][0] = -1.;
  b[1][1] = 2.;
  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();

  TH2D *hb = new TH2D(b);
  hb->Draw("colz");
/*
  TMatrixD c = b.Invert(1);

  //TMatrixD c(a, TMatrixD(a, TMatrixD::kTransposeMult,b));

  std::cout << c[0][0] << "\t" << c[0][1] << "\t" << c[1][0] << "\t" << c[1][1] << "\n";

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();

  TH2D *h = new TH2D(c);
  h->Draw("colz");
*/
}

