#ifndef TEMPLATEFITTER_H 
#define TEMPLATEFITTER_H 

#include "Math/Minimizer.h"
#include "TH1.h"

class TemplateFitter {

  public:

    TemplateFitter(TH1D * templates_m[480], TH1D * templates_m_w[480], TH1D * templates_e[480], TH1D * templates_e_w[480], TH1D * target);
    ~TemplateFitter(){};
    void setEnergyBins(double bins[481]);
    //bool doFit(double &theta, double &dm2);
    void Draw();
    //void TrueDraw();
    //void Test(double theta);
    double getChi2(double * par);

  private:

    double getPmue(double energy, double Uee2, double Umm2, double dm2);
    double getPee(double energy, double Uee2, double Umm2, double dm2);
    double getPmm(double energy, double Uee2, double Umm2, double dm2);

    // The templates are reconstructed lepton energy, in a slice of true neutrino energy
    TH1D * m_templates[480];
    TH1D * w_m_templates[480];
    // The template for the intrinsic nue
    TH1D * e_templates[480]; 
    TH1D * w_e_templates[480];
    // define the energy bins used in the template
    double m_energy_bins[481];
    // this is the thing you are trying to fit to, i.e. the data distribution
    TH1D * m_target;
    
};

#endif
