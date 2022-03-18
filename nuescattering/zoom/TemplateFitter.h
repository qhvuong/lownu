#ifndef TEMPLATEFITTER_H 
#define TEMPLATEFITTER_H 

#include "Math/Minimizer.h"
#include "TH1.h"

class TemplateFitter {

  public:

    TemplateFitter(TH1D * templates_m[280], TH1D * templates_m_w[280], TH1D * templates_e[280], TH1D * templates_e_w[280], TH1D * target);
    ~TemplateFitter(){};
    void setEnergyBins(double bins[281]);
    bool doFit(double &Uee2, double &Umm2, double &dm2);
    void Draw(double Uee2, double Umm2, double dm2);
    void TrueDraw();

  private:

    double getPmue(double energy, double Uee2, double Umm2, double dm2);
    double getPee(double energy, double Uee2, double Umm2, double dm2);
    double getPmm(double energy, double Uee2, double Umm2, double dm2);
    double getChi2(const double * par);

    // The templates are reconstructed lepton energy, in a slice of true neutrino energy
    TH1D * m_templates[280];
    TH1D * w_m_templates[280];
    // The template for the intrinsic nue
    TH1D * e_templates[280]; 
    TH1D * w_e_templates[280]; 
    // define the energy bins used in the template
    double m_energy_bins[281];
    // this is the thing you are trying to fit to, i.e. the data distribution
    TH1D * m_target;
    
};

#endif
