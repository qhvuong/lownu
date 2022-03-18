#ifndef TEMPLATEFITTER_H 
#define TEMPLATEFITTER_H 

#include "Math/Minimizer.h"
#include "TH1.h"

class TemplateFitter {

  public:

    TemplateFitter(TH1D * CC_templates_m[280], TH1D * CC_templates_m_nc[280], TH1D * CC_templates_e[280], TH1D * nue_templates_m[280], TH1D * nue_templates_e[280], TH1D * CC_target_m, TH1D * CC_target_e, TH1D * nue_target);
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
    TH1D * CC_m_templates[280];
    // The templates are reconstructed lepton energy, in a slice of true neutrino energy no cuts
    TH1D * CC_nc_m_templates[280];
    // The template for the intrinsic nue
    TH1D * CC_e_templates[280]; 
    TH1D * nue_m_templates[280];
    TH1D * nue_e_templates[280];
    // define the energy bins used in the template
    double m_energy_bins[281];
    // this is the thing you are trying to fit to, i.e. the data distribution
    TH1D * CC_m_target;
    TH1D * CC_e_target;
    TH1D * nue_em_target;
    
};
/*
class TemplateFitter_nue {

  public:

    TemplateFitter_nue(TH1D * templates_nue[100], TH1D * intrinsic_nue, TH1D * target_nue);
    ~TemplateFitter_nue(){};
    void setEnergyBins(double bins[101]);
    bool doFit(double &theta, double &dm2);

  private:

    double getOscProb(double energy, double theta, double dm2);
    double getChi2(const double * par);

    // The templates are reconstructed lepton energy, in a slice of true neutrino energy
    TH1D * m_templates_nue[100];
    // The template for the intrinsic nue
    TH1D * m_intrinsic_nue; 
    // define the energy bins used in the template
    double m_energy_bins_nue[101];
    // this is the thing you are trying to fit to, i.e. the data distribution
    TH1D * m_target_nue;
};*/

#endif
