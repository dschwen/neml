#include "walker.h"

namespace neml {

WalkerKremplSwitchRule::WalkerKremplSwitchRule(
    std::shared_ptr<LinearElasticModel> elastic,
    std::shared_ptr<ViscoPlasticFlowRule> flow,
    std::shared_ptr<Interpolate> lambda,
    double eps0) :
      elastic_(elastic), flow_(flow), lambda_(lambda), eps0_(eps0)
{

}

std::string WalkerKremplSwitchRule::type()
{
  return "WalkerKremplSwitchRule";
}

ParameterSet WalkerKremplSwitchRule::parameters()
{
  ParameterSet pset(WalkerKremplSwitchRule::type());

  pset.add_parameter<NEMLObject>("elastic");
  pset.add_parameter<NEMLObject>("flow");
  pset.add_parameter<NEMLObject>("lambda");
  pset.add_parameter<double>("eps_ref");

  return pset;
}

std::unique_ptr<NEMLObject> WalkerKremplSwitchRule::initialize(ParameterSet & params)
{
  return neml::make_unique<WalkerKremplSwitchRule>(
      params.get_object_parameter<LinearElasticModel>("elastic"),
      params.get_object_parameter<ViscoPlasticFlowRule>("flow"),
      params.get_object_parameter<Interpolate>("lambda"),
      params.get_parameter<double>("eps_ref")
      ); 
}


size_t WalkerKremplSwitchRule::nhist() const
{
  return flow_->nhist();
}

int WalkerKremplSwitchRule::init_hist(double * const h)
{
  return flow_->init_hist(h);
}

int WalkerKremplSwitchRule::s(const double * const s, const double * const alpha,
              const double * const edot, double T,
              double Tdot,
              double * const sdot)
{
  double erate[6];
  std::copy(edot, edot+6, erate);

  double temp[6];
  double yv;
  int ier = flow_->g(s, alpha, T, temp);
  if (ier != SUCCESS) return ier;
  ier = flow_->y(s, alpha, T, yv);
  if (ier != SUCCESS) return ier;
  
  double kap;
  ier = kappa(edot, T, kap);
  if (ier != SUCCESS) return ier;

  for (int i=0; i<6; i++) {
    erate[i] -= yv * kap * temp[i];
  }

  double C[36];
  elastic_->C(T, C);

  mat_vec(C, 6, erate, 6, sdot);

  return 0;

}

int WalkerKremplSwitchRule::ds_ds(const double * const s, const double * const alpha,
              const double * const edot, double T,
              double Tdot,
              double * const d_sdot)
{
  double yv;
  int ier = flow_->y(s, alpha, T, yv);
  if (ier != SUCCESS) return ier;

  double kap;
  ier = kappa(edot, T, kap);
  if (ier != SUCCESS) return ier;

  double work[36];
  ier = flow_->dg_ds(s, alpha, T, work);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<36; i++) {
    work[i] *= -yv * kap;
  }

  double t1[6];
  ier = flow_->g(s, alpha, T, t1);
  if (ier != SUCCESS) return ier;
  double t2[6];
  ier = flow_->dy_ds(s, alpha, T, t2);
  if (ier != SUCCESS) return ier;
  for (size_t i = 0; i < 6;  i++) t2[i] *= kap;
  outer_update_minus(t1, 6, t2, 6, work);
  
  double t3[36];
  elastic_->C(T, t3);

  mat_mat(6,6,6, t3, work, d_sdot);

  return 0;
}

int WalkerKremplSwitchRule::ds_da(const double * const s, const double * const alpha,
              const double * const edot, double T,
              double Tdot,
              double * const d_sdot)
{
  double yv;
  int ier = flow_->y(s, alpha, T, yv);
  if (ier != SUCCESS) return ier;
 
  double kap;
  ier = kappa(edot, T, kap);
  if (ier != SUCCESS) return ier;

  int sz = 6 * nhist();
  
  std::vector<double> workv(sz);
  double * work = &workv[0];
  ier = flow_->dg_da(s, alpha, T, work);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<sz; i++) {
    work[i] *= -yv * kap;
  }

  double t1[6];
  ier = flow_->g(s, alpha, T, t1);
  if (ier != SUCCESS) return ier;
  std::vector<double> t2v(nhist());
  double * t2 = &t2v[0];
  ier = flow_->dy_da(s, alpha, T, t2);
  if (ier != SUCCESS) return ier;
  for (size_t i = 0; i < nhist(); i++) t2[i] *= kap;
  outer_update_minus(t1, 6, t2, nhist(), work);
  
  double C[36];
  elastic_->C(T, C);

  mat_mat(6, nhist(), 6, C, work, d_sdot);

  return 0;

}

int WalkerKremplSwitchRule::ds_de(const double * const s, const double * const alpha,
              const double * const edot, double T,
              double Tdot,
              double * const d_sdot)
{
  double work[36];
  
  double yv;
  int ier = flow_->y(s, alpha, T, yv);
  if (ier != SUCCESS) return ier;
 
  double dkap[6];
  ier = dkappa(edot, T, dkap);
  if (ier != SUCCESS) return ier;

  double g[6];
  ier = flow_->g(s, alpha, T, g);
  if (ier != SUCCESS) return ier;

  for (size_t i = 0; i < 6; i++) g[i] *= yv;

  std::fill(work, work+36, 0.0);
  for (size_t i = 0; i < 6; i++) work[CINDEX(i,i,6)] = 1.0;
  outer_update_minus(g, 6, dkap, 6, work);

  double C[36];
  ier = elastic_->C(T, C);
  if (ier != SUCCESS) return ier;

  mat_mat(6, 6, 6, C, work, d_sdot);

  return 0;
}

int WalkerKremplSwitchRule::a(const double * const s, const double * const alpha,
              const double * const edot, double T,
              double Tdot,
              double * const adot)
{
  double dg;
  int ier = flow_->y(s, alpha, T, dg);
  if (ier != SUCCESS) return 0;

  double kap;
  ier = kappa(edot, T, kap);
  if (ier != SUCCESS) return ier;

  ier = flow_->h(s, alpha, T, adot);
  if (ier != SUCCESS) return 0;
  for (size_t i=0; i<nhist(); i++) adot[i] *= (dg * kap);
  
  std::vector<double> tempv(nhist());
  double * temp = &tempv[0];
  ier = flow_->h_temp(s, alpha, T, temp);
  if (ier != SUCCESS) return ier;
  for (size_t i=0; i<nhist(); i++) adot[i] += temp[i] * Tdot;

  ier = flow_->h_time(s, alpha, T, temp);
  if (ier != SUCCESS) return ier;
  for (size_t i=0; i<nhist(); i++) adot[i] += (temp[i] * kap);

  return 0;

}

int WalkerKremplSwitchRule::da_ds(const double * const s, const double * const alpha,
              const double * const edot, double T,
              double Tdot,
              double * const d_adot)
{
  double dg;
  int ier = flow_->y(s, alpha, T, dg);
  if (ier != SUCCESS) return ier;

  double kap;
  ier = kappa(edot, T, kap);
  if (ier != SUCCESS) return ier;

  int sz = nhist() * 6;

  ier = flow_->dh_ds(s, alpha, T, d_adot);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<sz; i++) d_adot[i] *= (dg * kap);

  std::vector<double> t1v(nhist());
  double * t1 = &t1v[0];
  ier = flow_->h(s, alpha, T, t1);
  if (ier != SUCCESS) return ier;

  double t2[6];
  ier = flow_->dy_ds(s, alpha, T, t2);
  if (ier != SUCCESS) return ier;
  for (size_t i = 0; i < 6; i++) t2[i] *= kap;

  outer_update(t1, nhist(), t2, 6, d_adot);
  
  std::vector<double> t3v(sz);
  double * t3 = &t3v[0];
  ier = flow_->dh_ds_temp(s, alpha, T, t3);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<sz; i++) d_adot[i] += t3[i] * Tdot;

  ier = flow_->dh_ds_time(s, alpha, T, t3);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<sz; i++) d_adot[i] += t3[i] * kap;

  return 0;
  
}

int WalkerKremplSwitchRule::da_da(const double * const s, const double * const alpha,
              const double * const edot, double T,
              double Tdot,
              double * const d_adot)
{
  double dg;
  int ier = flow_->y(s, alpha, T, dg);
  if (ier != SUCCESS) return ier;

  double kap;
  ier = kappa(edot, T, kap);
  if (ier != SUCCESS) return ier;

  int nh = nhist();
  int sz = nh * nh;

  ier = flow_->dh_da(s, alpha, T, d_adot);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<sz; i++) d_adot[i] *= dg * kap;
  
  std::vector<double> t1v(nh);
  double * t1 = &t1v[0];
  ier = flow_->h(s, alpha, T, t1);
  if (ier != SUCCESS) return ier;
  
  std::vector<double> t2v(nh);
  double * t2 = &t2v[0];
  ier = flow_->dy_da(s, alpha, T, t2);
  if (ier != SUCCESS) return ier;

  for (int i = 0 ; i < nh; i++) t2[i] *= kap;

  outer_update(t1, nh, t2, nh, d_adot);
  
  std::vector<double> t3v(sz);
  double * t3 = &t3v[0];
  ier = flow_->dh_da_temp(s, alpha, T, t3);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<sz; i++) d_adot[i] += t3[i] * Tdot;

  ier = flow_->dh_da_time(s, alpha, T, t3);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<sz; i++) d_adot[i] += t3[i] * kap;

  return 0;
}

int WalkerKremplSwitchRule::da_de(const double * const s, const double * const alpha,
              const double * const edot, double T,
              double Tdot,
              double * const d_adot)
{
  double dg;
  int ier = flow_->y(s, alpha, T, dg);
  if (ier != SUCCESS) return ier;

  double dkap[6];
  ier = dkappa(edot, T, dkap);
  if (ier != SUCCESS) return ier;

  int nh = nhist();
  double * hr = new double [nh];
  ier = flow_->h(s, alpha, T, hr);
  if (ier != SUCCESS) return 0;
  for (int i = 0; i < nh; i++) hr[i] *= dg;

  outer_vec(hr, nh, dkap, 6, d_adot);

  ier = flow_->h_time(s, alpha, T, hr);
  if (ier != SUCCESS) return 0;
  outer_update(hr, nh, dkap, 6, d_adot);

  delete [] hr;

  return 0;
}

int WalkerKremplSwitchRule::work_rate(const double * const s,
                                    const double * const alpha,
                                    const double * const edot, double T,
                                    double Tdot, double & p_dot)
{
  double erate[6];
  std::fill(erate, erate+6, 0.0);

  double kap;
  int ier = kappa(edot, T, kap);
  if (ier != SUCCESS) return ier;

  double temp[6];
  double yv;
  ier = flow_->g(s, alpha, T, temp);
  if (ier != SUCCESS) return ier;
  ier = flow_->y(s, alpha, T, yv);
  if (ier != SUCCESS) return ier;

  for (int i=0; i<6; i++) {
    erate[i] += yv * kap * temp[i];
  }

  ier = flow_->g_temp(s, alpha, T, temp);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<6; i++) {
    erate[i] += Tdot * temp[i];
  }

  ier = flow_->g_time(s, alpha, T, temp);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<6; i++) {
    erate[i] += temp[i];
  }

  p_dot = dot_vec(s, erate, 6);
  
  return 0;
}

int WalkerKremplSwitchRule::elastic_strains(const double * const s_np1, double T_np1,
                                 double * const e_np1) const
{
  double S[36];
  elastic_->S(T_np1, S);
  mat_vec(S, 6, s_np1, 6, e_np1);

  return 0;
}

int WalkerKremplSwitchRule::set_elastic_model(std::shared_ptr<LinearElasticModel> emodel)
{
  elastic_ = emodel;
  return 0;
}

int WalkerKremplSwitchRule::kappa(const double * const edot, double T, double &
                                  kap)
{
  double edev[6];
  std::copy(edot, edot+6, edev);
  dev_vec(edev);

  double de = std::sqrt(2.0/3.0) * norm2_vec(edev, 6);

  kap = 1.0 - lambda_->value(T) + lambda_->value(T) * de / eps0_;

  return 0;
}

int WalkerKremplSwitchRule::dkappa(const double * const edot, double T,
                                   double * const dkap)
{
  std::copy(edot, edot+6, dkap);
  dev_vec(dkap);

  double fact = lambda_->value(T)  / eps0_ * std::sqrt(2.0/3.0) / norm2_vec(dkap, 6);

  for (size_t i = 0; i < 6; i++)
    dkap[i] *= fact;

  return 0;
}

SofteningModel::SofteningModel()
{

}

std::string SofteningModel::type()
{
  return "SofteningModel";
}

std::unique_ptr<NEMLObject> SofteningModel::initialize(ParameterSet & params)
{
  return neml::make_unique<SofteningModel>(); 
}

ParameterSet SofteningModel::parameters()
{
  ParameterSet pset(SofteningModel::type());

  return pset;
}

double SofteningModel::phi(double alpha, double T) const
{
  return 1.0;
}

double SofteningModel::dphi(double alpha, double T) const
{
  return 0.0;
}

WalkerSofteningModel::WalkerSofteningModel(std::shared_ptr<Interpolate> phi0,
                                           std::shared_ptr<Interpolate> phi1) :
    phi_0_(phi0), phi_1_(phi1)
{

}

std::string WalkerSofteningModel::type()
{
  return "WalkerSofteningModel";
}

std::unique_ptr<NEMLObject> WalkerSofteningModel::initialize(ParameterSet & params)
{
  return neml::make_unique<WalkerSofteningModel>(
      params.get_object_parameter<Interpolate>("phi_0"),
      params.get_object_parameter<Interpolate>("phi_1")
      ); 
}

ParameterSet WalkerSofteningModel::parameters()
{
  ParameterSet pset(WalkerSofteningModel::type());

  pset.add_parameter<NEMLObject>("phi_0");
  pset.add_parameter<NEMLObject>("phi_1");

  return pset;
}

double WalkerSofteningModel::phi(double alpha, double T) const
{
  return 1.0 + phi_0_->value(T) * std::pow(alpha, phi_1_->value(T));
}

double WalkerSofteningModel::dphi(double alpha, double T) const
{
  return phi_1_->value(T) * phi_0_->value(T) * std::pow(alpha, phi_1_->value(T)
                                                        - 1.0);
}

ThermalScaling::ThermalScaling()
{

}

std::string ThermalScaling::type()
{
  return "ThermalScaling";
}

std::unique_ptr<NEMLObject> ThermalScaling::initialize(ParameterSet & params)
{
  return neml::make_unique<ThermalScaling>(); 
}

ParameterSet ThermalScaling::parameters()
{
  ParameterSet pset(ThermalScaling::type());

  return pset;
}

double ThermalScaling::value(double T) const
{
  return 1.0;
}

ArrheniusThermalScaling::ArrheniusThermalScaling(std::shared_ptr<Interpolate> Q,
                                                 double R, double T_ref) :
    Q_(Q), R_(R), T_ref_(T_ref)
{

}

std::string ArrheniusThermalScaling::type()
{
  return "ArrheniusThermalScaling";
}

std::unique_ptr<NEMLObject> ArrheniusThermalScaling::initialize(ParameterSet & params)
{
  return neml::make_unique<ArrheniusThermalScaling>(
      params.get_object_parameter<Interpolate>("Q"),
      params.get_parameter<double>("R"),
      params.get_parameter<double>("T_ref")
      ); 
}

ParameterSet ArrheniusThermalScaling::parameters()
{
  ParameterSet pset(ArrheniusThermalScaling::type());

  pset.add_parameter<NEMLObject>("Q");
  pset.add_parameter<double>("R");
  pset.add_parameter<double>("T_ref");

  return pset;
}

double ArrheniusThermalScaling::value(double T) const
{
  return arr_(T) / arr_(T_ref_);
}

double ArrheniusThermalScaling::arr_(double T) const
{
  return std::exp(-Q_->value(T) / (R_ * T));
}

IsotropicHardening::IsotropicHardening(std::string name, 
                                       std::shared_ptr<ThermalScaling> scale) :
    ScalarInternalVariable(name), scale_(scale)
{}

/// Return zero for time rate by default 
double IsotropicHardening::ratet(VariableState & state)
{
  return 0;
}

/// Return zero for the time rate derivatives by default
double IsotropicHardening::d_ratet_d_h(VariableState & state)
{
  return 0;
}

/// Return zero for the time rate derivatives by default
double IsotropicHardening::d_ratet_d_a(VariableState & state) 
{
  return 0;
}

/// Return zero for the time rate derivatives by default
double IsotropicHardening::d_ratet_d_adot(VariableState & state)
{
  return 0;
}

/// Return zero for the time rate derivatives by default
double IsotropicHardening::d_ratet_d_D(VariableState & state)
{
  return 0;
}

/// Return zero for the time rate derivatives by default
Symmetric IsotropicHardening::d_ratet_d_s(VariableState & state)
{
  return Symmetric::zero();
}

/// Return zero for the time rate derivatives by default
Symmetric IsotropicHardening::d_ratet_d_g(VariableState & state)
{
  return Symmetric::zero();
}

/// Return zero for temperature rate by default 
double IsotropicHardening::rateT(VariableState & state)
{
  return 0;
}

/// Return zero for the temperature rate derivatives by default
double IsotropicHardening::d_rateT_d_h(VariableState & state)
{
  return 0;
}

/// Return zero for the temperature rate derivatives by default
double IsotropicHardening::d_rateT_d_a(VariableState & state) 
{
  return 0;
}

/// Return zero for the temperature rate derivatives by default
double IsotropicHardening::d_rateT_d_adot(VariableState & state)
{
  return 0;
}

/// Return zero for the temperature rate derivatives by default
double IsotropicHardening::d_rateT_d_D(VariableState & state)
{
  return 0;
}

/// Return zero for the temperature rate derivatives by default
Symmetric IsotropicHardening::d_rateT_d_s(VariableState & state)
{
  return Symmetric::zero();
}

/// Return zero for the temperature rate derivatives by default
Symmetric IsotropicHardening::d_rateT_d_g(VariableState & state)
{
  return Symmetric::zero();
}


ConstantIsotropicHardening::ConstantIsotropicHardening(
    std::shared_ptr<ThermalScaling> scale) :
      IsotropicHardening("R", scale)
{

}

std::string ConstantIsotropicHardening::type()
{
  return "ConstantIsotropicHardening";
}

ParameterSet ConstantIsotropicHardening::parameters()
{
  ParameterSet pset(ConstantIsotropicHardening::type());

  pset.add_optional_parameter<NEMLObject>("scaling", 
                                          std::make_shared<ThermalScaling>());

  return pset;
}

std::unique_ptr<NEMLObject> ConstantIsotropicHardening::initialize(
    ParameterSet & params)
{
  return neml::make_unique<ConstantIsotropicHardening>(
      params.get_object_parameter<ThermalScaling>("scaling")
      ); 
}

double ConstantIsotropicHardening::initial_value()
{
  return 0;
}

double ConstantIsotropicHardening::ratep(VariableState & state)
{
  return 0;
}

double ConstantIsotropicHardening::d_ratep_d_h(VariableState & state)
{
  return 0;
}

double ConstantIsotropicHardening::d_ratep_d_a(VariableState & state)
{
  return 0;
}

double ConstantIsotropicHardening::d_ratep_d_adot(VariableState & state)
{
  return 0;
}

double ConstantIsotropicHardening::d_ratep_d_D(VariableState & state)
{
  return 0;
}

Symmetric ConstantIsotropicHardening::d_ratep_d_s(VariableState & state)
{
  return Symmetric();
}

Symmetric ConstantIsotropicHardening::d_ratep_d_g(VariableState & state)
{
  return Symmetric();
}


WalkerIsotropicHardening::WalkerIsotropicHardening(
    std::shared_ptr<Interpolate> r0, std::shared_ptr<Interpolate> Rinf,
    std::shared_ptr<Interpolate> R0, std::shared_ptr<Interpolate> r1,
    std::shared_ptr<Interpolate> r2, std::shared_ptr<ThermalScaling> scale) :
      IsotropicHardening("R", scale), r0_(r0), Rinf_(Rinf), R0_(R0), r1_(r1), 
      r2_(r2)
{

}

std::string WalkerIsotropicHardening::type()
{
  return "WalkerIsotropicHardening";
}

ParameterSet WalkerIsotropicHardening::parameters()
{
  ParameterSet pset(WalkerIsotropicHardening::type());

  pset.add_parameter<NEMLObject>("r0");
  pset.add_parameter<NEMLObject>("Rinf");
  pset.add_parameter<NEMLObject>("R0");
  pset.add_parameter<NEMLObject>("r1");
  pset.add_parameter<NEMLObject>("r2");
  pset.add_optional_parameter<NEMLObject>("scaling", 
                                          std::make_shared<ThermalScaling>());

  return pset;
}

std::unique_ptr<NEMLObject> WalkerIsotropicHardening::initialize(ParameterSet & params)
{
  return neml::make_unique<WalkerIsotropicHardening>(
      params.get_object_parameter<Interpolate>("r0"),
      params.get_object_parameter<Interpolate>("Rinf"),
      params.get_object_parameter<Interpolate>("R0"),
      params.get_object_parameter<Interpolate>("r1"),
      params.get_object_parameter<Interpolate>("r2"),
      params.get_object_parameter<ThermalScaling>("scaling")
      ); 
}

double WalkerIsotropicHardening::initial_value()
{
  return 0;
}

double WalkerIsotropicHardening::ratep(VariableState & state)
{
  return r0_->value(state.T) * (Rinf_->value(state.T) - state.h);
}

double WalkerIsotropicHardening::d_ratep_d_h(VariableState & state)
{
  return -r0_->value(state.T);
}

double WalkerIsotropicHardening::d_ratep_d_a(VariableState & state)
{
  return 0;
}

double WalkerIsotropicHardening::d_ratep_d_adot(VariableState & state)
{
  return 0;
}

double WalkerIsotropicHardening::d_ratep_d_D(VariableState & state)
{
  return 0;
}

Symmetric WalkerIsotropicHardening::d_ratep_d_s(VariableState & state)
{
  return Symmetric();
}

Symmetric WalkerIsotropicHardening::d_ratep_d_g(VariableState & state)
{
  return Symmetric();
}

double WalkerIsotropicHardening::ratet(VariableState & state)
{
  double d = R0_->value(state.T) - state.h;
  return r1_->value(state.T) * d * std::pow(std::fabs(d), r2_->value(state.T) - 1.0);
}

double WalkerIsotropicHardening::d_ratet_d_h(VariableState & state)
{
  return -r1_->value(state.T) * r2_->value(state.T) * 
      std::pow(std::fabs(R0_->value(state.T) - state.h), r2_->value(state.T) - 1.0);
}

double WalkerIsotropicHardening::d_ratet_d_a(VariableState & state)
{
  return 0;
}

double WalkerIsotropicHardening::d_ratet_d_adot(VariableState & state)
{
  return 0;
}

double WalkerIsotropicHardening::d_ratet_d_D(VariableState & state)
{
  return 0;
}

Symmetric WalkerIsotropicHardening::d_ratet_d_s(VariableState & state)
{
  return Symmetric();
}

Symmetric WalkerIsotropicHardening::d_ratet_d_g(VariableState & state)
{
  return Symmetric();
}

DragStress::DragStress(std::string name, 
                                       std::shared_ptr<ThermalScaling> scale) :
    ScalarInternalVariable(name), scale_(scale)
{}

/// Makes no sense in this context
double DragStress::d_ratep_d_D(VariableState & state)
{
  return 0;
}

/// Return zero for time rate by default 
double DragStress::ratet(VariableState & state)
{
  return 0;
}

/// Return zero for the time rate derivatives by default
double DragStress::d_ratet_d_h(VariableState & state)
{
  return 0;
}

/// Return zero for the time rate derivatives by default
double DragStress::d_ratet_d_a(VariableState & state) 
{
  return 0;
}

/// Return zero for the time rate derivatives by default
double DragStress::d_ratet_d_adot(VariableState & state)
{
  return 0;
}

/// Return zero for the time rate derivatives by default
double DragStress::d_ratet_d_D(VariableState & state)
{
  return 0;
}

/// Return zero for the time rate derivatives by default
Symmetric DragStress::d_ratet_d_s(VariableState & state)
{
  return Symmetric::zero();
}

/// Return zero for the time rate derivatives by default
Symmetric DragStress::d_ratet_d_g(VariableState & state)
{
  return Symmetric::zero();
}

/// Return zero for temperature rate by default 
double DragStress::rateT(VariableState & state)
{
  return 0;
}

/// Return zero for the temperature rate derivatives by default
double DragStress::d_rateT_d_h(VariableState & state)
{
  return 0;
}

/// Return zero for the temperature rate derivatives by default
double DragStress::d_rateT_d_a(VariableState & state) 
{
  return 0;
}

/// Return zero for the temperature rate derivatives by default
double DragStress::d_rateT_d_adot(VariableState & state)
{
  return 0;
}

/// Return zero for the temperature rate derivatives by default
double DragStress::d_rateT_d_D(VariableState & state)
{
  return 0;
}

/// Return zero for the temperature rate derivatives by default
Symmetric DragStress::d_rateT_d_s(VariableState & state)
{
  return Symmetric::zero();
}

/// Return zero for the temperature rate derivatives by default
Symmetric DragStress::d_rateT_d_g(VariableState & state)
{
  return Symmetric::zero();
}

ConstantDragStress::ConstantDragStress(double value,
    std::shared_ptr<ThermalScaling> scale) :
      DragStress("D", scale), value_(value)
{

}

std::string ConstantDragStress::type()
{
  return "ConstantDragStress";
}

ParameterSet ConstantDragStress::parameters()
{
  ParameterSet pset(ConstantDragStress::type());
  
  pset.add_parameter<double>("value");
  pset.add_optional_parameter<NEMLObject>("scaling", 
                                          std::make_shared<ThermalScaling>());

  return pset;
}

std::unique_ptr<NEMLObject> ConstantDragStress::initialize(
    ParameterSet & params)
{
  return neml::make_unique<ConstantDragStress>(
      params.get_parameter<double>("value"),
      params.get_object_parameter<ThermalScaling>("scaling")
      ); 
}

double ConstantDragStress::initial_value()
{
  return value_;
}

double ConstantDragStress::D_xi(double T)
{
  return 1; // Can be an arbitrary value for this model
}

double ConstantDragStress::D_0(double T)
{
  return value_; // Straightforward!
}

double ConstantDragStress::ratep(VariableState & state)
{
  return 0;
}

double ConstantDragStress::d_ratep_d_h(VariableState & state)
{
  return 0;
}

double ConstantDragStress::d_ratep_d_a(VariableState & state)
{
  return 0;
}

double ConstantDragStress::d_ratep_d_adot(VariableState & state)
{
  return 0;
}

Symmetric ConstantDragStress::d_ratep_d_s(VariableState & state)
{
  return Symmetric();
}

Symmetric ConstantDragStress::d_ratep_d_g(VariableState & state)
{
  return Symmetric();
}

WalkerDragStress::WalkerDragStress(
    std::shared_ptr<Interpolate> d0, std::shared_ptr<Interpolate> d1,
    std::shared_ptr<Interpolate> d2, std::shared_ptr<Interpolate> D_xi,
    double D_0, std::shared_ptr<SofteningModel> softening,
    std::shared_ptr<ThermalScaling> scale) :
      DragStress("D", scale), d0_(d0), d1_(d1), d2_(d2),
      D_xi_(D_xi), D_0_(D_0), softening_(softening)
{

}

std::string WalkerDragStress::type()
{
  return "WalkerDragStress";
}

ParameterSet WalkerDragStress::parameters()
{
  ParameterSet pset(WalkerDragStress::type());
  
  pset.add_parameter<NEMLObject>("d0");
  pset.add_parameter<NEMLObject>("d1");
  pset.add_parameter<NEMLObject>("d2");
  pset.add_parameter<NEMLObject>("D_xi");
  pset.add_parameter<double>("D_0");
  pset.add_parameter<NEMLObject>("softening");
  pset.add_optional_parameter<NEMLObject>("scaling", 
                                          std::make_shared<ThermalScaling>());

  return pset;
}

std::unique_ptr<NEMLObject> WalkerDragStress::initialize(
    ParameterSet & params)
{
  return neml::make_unique<WalkerDragStress>(
      params.get_object_parameter<Interpolate>("d0"),
      params.get_object_parameter<Interpolate>("d1"),
      params.get_object_parameter<Interpolate>("d2"),
      params.get_object_parameter<Interpolate>("D_xi"),
      params.get_parameter<double>("D_0"),
      params.get_object_parameter<SofteningModel>("softening"),
      params.get_object_parameter<ThermalScaling>("scaling")
      ); 
}

double WalkerDragStress::initial_value()
{
  return D_0_;
}

double WalkerDragStress::D_xi(double T)
{
  return D_xi_->value(T); // Explicit parameter
}

double WalkerDragStress::D_0(double T)
{
  return D_0_; // Straightforward!
}

double WalkerDragStress::ratep(VariableState & state)
{
  return d0_->value(state.T) * (1.0 - state.h / D_xi_->value(state.T));
}

double WalkerDragStress::d_ratep_d_h(VariableState & state)
{
  return -d0_->value(state.T) / D_xi_->value(state.T);
}

double WalkerDragStress::d_ratep_d_a(VariableState & state)
{
  return 0;
}

double WalkerDragStress::d_ratep_d_adot(VariableState & state)
{
  return 0;
}

Symmetric WalkerDragStress::d_ratep_d_s(VariableState & state)
{
  return Symmetric();
}

Symmetric WalkerDragStress::d_ratep_d_g(VariableState & state)
{
  return Symmetric();
}

double WalkerDragStress::ratet(VariableState & state)
{
  return -scale_->value(state.T) * softening_->phi(state.a, state.T) * 
      d1_->value(state.T) * std::pow(state.h, d2_->value(state.T));
}

double WalkerDragStress::d_ratet_d_h(VariableState & state)
{
  return -d2_->value(state.T) * scale_->value(state.T) * softening_->phi(state.a, state.T) * 
      d1_->value(state.T) * std::pow(state.h, d2_->value(state.T) - 1.0);
}

double WalkerDragStress::d_ratet_d_a(VariableState & state)
{
  return -scale_->value(state.T) * softening_->dphi(state.a, state.T) * 
      d1_->value(state.T) * std::pow(state.h, d2_->value(state.T));
}

double WalkerDragStress::d_ratet_d_adot(VariableState & state)
{
  return 0;
}

Symmetric WalkerDragStress::d_ratet_d_s(VariableState & state)
{
  return Symmetric();
}

Symmetric WalkerDragStress::d_ratet_d_g(VariableState & state)
{
  return Symmetric();
}

KinematicHardening::KinematicHardening(std::string name, 
                                       std::shared_ptr<ThermalScaling> scale) :
    SymmetricInternalVariable(name), scale_(scale)
{}

/// Return zero for time rate by default 
Symmetric KinematicHardening::ratet(VariableState & state)
{
  return Symmetric::zero();
}

/// Return zero for the time rate derivatives by default
SymSymR4 KinematicHardening::d_ratet_d_h(VariableState & state)
{
  return SymSymR4::zero();
}

/// Return zero for the time rate derivatives by default
Symmetric KinematicHardening::d_ratet_d_a(VariableState & state) 
{
  return Symmetric::zero();
}

/// Return zero for the time rate derivatives by default
Symmetric KinematicHardening::d_ratet_d_adot(VariableState & state)
{
  return Symmetric::zero();
}

/// Return zero for the time rate derivatives by default
Symmetric KinematicHardening::d_ratet_d_D(VariableState & state)
{
  return Symmetric::zero();
}

/// Return zero for the time rate derivatives by default
SymSymR4 KinematicHardening::d_ratet_d_s(VariableState & state)
{
  return SymSymR4::zero();
}

/// Return zero for the time rate derivatives by default
SymSymR4 KinematicHardening::d_ratet_d_g(VariableState & state)
{
  return SymSymR4::zero();
}

/// Return zero for temperature rate by default 
Symmetric KinematicHardening::rateT(VariableState & state)
{
  return Symmetric::zero();
}

/// Return zero for the temperature rate derivatives by default
SymSymR4 KinematicHardening::d_rateT_d_h(VariableState & state)
{
  return SymSymR4::zero();
}

/// Return zero for the temperature rate derivatives by default
Symmetric KinematicHardening::d_rateT_d_a(VariableState & state) 
{
  return Symmetric::zero();
}

/// Return zero for the temperature rate derivatives by default
Symmetric KinematicHardening::d_rateT_d_adot(VariableState & state)
{
  return Symmetric::zero();
}

/// Return zero for the temperature rate derivatives by default
Symmetric KinematicHardening::d_rateT_d_D(VariableState & state)
{
  return Symmetric::zero();
}

/// Return zero for the temperature rate derivatives by default
SymSymR4 KinematicHardening::d_rateT_d_s(VariableState & state)
{
  return SymSymR4::zero();
}

/// Return zero for the temperature rate derivatives by default
SymSymR4 KinematicHardening::d_rateT_d_g(VariableState & state)
{
  return SymSymR4::zero();
}

FAKinematicHardening::FAKinematicHardening(std::shared_ptr<Interpolate> c,
                                           std::shared_ptr<Interpolate> g,
                                           std::shared_ptr<ThermalScaling>
                                           scale) : 
    KinematicHardening("X", scale), c_(c), g_(g)
{

}

std::string FAKinematicHardening::type()
{
  return "FAKinematicHardening";
}

ParameterSet FAKinematicHardening::parameters()
{
  ParameterSet pset(FAKinematicHardening::type());
  
  pset.add_parameter<NEMLObject>("c");
  pset.add_parameter<NEMLObject>("g");
  pset.add_optional_parameter<NEMLObject>("scaling", 
                                          std::make_shared<ThermalScaling>());

  return pset;
}

std::unique_ptr<NEMLObject> FAKinematicHardening::initialize(
    ParameterSet & params)
{
  return neml::make_unique<FAKinematicHardening>(
      params.get_object_parameter<Interpolate>("c"),
      params.get_object_parameter<Interpolate>("g"),
      params.get_object_parameter<ThermalScaling>("scaling")
      ); 
}

Symmetric FAKinematicHardening::initial_value()
{
  return Symmetric::zero();
}

Symmetric FAKinematicHardening::ratep(VariableState & state)
{
  return 2.0/3.0 * c_->value(state.T) * state.g - g_->value(state.T) * state.h;
}

SymSymR4 FAKinematicHardening::d_ratep_d_h(VariableState & state)
{
  return -g_->value(state.T) * SymSymR4::id();
}

Symmetric FAKinematicHardening::d_ratep_d_a(VariableState & state) 
{
  return Symmetric::zero();
}

Symmetric FAKinematicHardening::d_ratep_d_adot(VariableState & state)
{
  return Symmetric::zero();
}

Symmetric FAKinematicHardening::d_ratep_d_D(VariableState & state)
{
  return Symmetric::zero();
}

SymSymR4 FAKinematicHardening::d_ratep_d_s(VariableState & state)
{
  return SymSymR4::zero();
}

SymSymR4 FAKinematicHardening::d_ratep_d_g(VariableState & state)
{
  return 2.0/3.0 * c_->value(state.T) * SymSymR4::id();
}


WalkerKinematicHardening::WalkerKinematicHardening(
    std::shared_ptr<Interpolate> c0, std::shared_ptr<Interpolate> c1, 
    std::shared_ptr<Interpolate> c2, std::shared_ptr<Interpolate> l0,
    std::shared_ptr<Interpolate> l1, std::shared_ptr<Interpolate> l,
    std::shared_ptr<Interpolate> b0, std::shared_ptr<Interpolate> x0,
    std::shared_ptr<Interpolate> x1, std::shared_ptr<SofteningModel> softening,
    std::shared_ptr<ThermalScaling> scale) : 
      KinematicHardening("X", scale), c0_(c0), c1_(c1), c2_(c2), l0_(l0), l1_(l1),
      l_(l), b0_(b0), x0_(x0), x1_(x1), softening_(softening)
{

}

std::string WalkerKinematicHardening::type()
{
  return "WalkerKinematicHardening";
}

ParameterSet WalkerKinematicHardening::parameters()
{
  ParameterSet pset(WalkerKinematicHardening::type());
  
  pset.add_parameter<NEMLObject>("c0");
  pset.add_parameter<NEMLObject>("c1");
  pset.add_parameter<NEMLObject>("c2");
  pset.add_parameter<NEMLObject>("l0");
  pset.add_parameter<NEMLObject>("l1");
  pset.add_parameter<NEMLObject>("l");
  pset.add_parameter<NEMLObject>("b0");
  pset.add_parameter<NEMLObject>("x0");
  pset.add_parameter<NEMLObject>("x1");
  pset.add_parameter<NEMLObject>("softening");

  pset.add_optional_parameter<NEMLObject>("scaling", 
                                          std::make_shared<ThermalScaling>());

  return pset;
}

std::unique_ptr<NEMLObject> WalkerKinematicHardening::initialize(
    ParameterSet & params)
{
  return neml::make_unique<WalkerKinematicHardening>(
      params.get_object_parameter<Interpolate>("c0"),
      params.get_object_parameter<Interpolate>("c1"),
      params.get_object_parameter<Interpolate>("c2"),
      params.get_object_parameter<Interpolate>("l0"),
      params.get_object_parameter<Interpolate>("l1"),
      params.get_object_parameter<Interpolate>("l"),
      params.get_object_parameter<Interpolate>("b0"),
      params.get_object_parameter<Interpolate>("x0"),
      params.get_object_parameter<Interpolate>("x1"),
      params.get_object_parameter<SofteningModel>("softening"),
      params.get_object_parameter<ThermalScaling>("scaling")
      ); 
}

Symmetric WalkerKinematicHardening::initial_value()
{
  return Symmetric::zero();
}

Symmetric WalkerKinematicHardening::ratep(VariableState & state)
{
  return c_(state) * (2.0/3.0 * state.g - b_(state) / L_(state));
}

SymSymR4 WalkerKinematicHardening::d_ratep_d_h(VariableState & state)
{
  return -c_(state) * db_dx_(state) / L_(state);
}

Symmetric WalkerKinematicHardening::d_ratep_d_a(VariableState & state) 
{
  return c_(state) * b_(state) / std::pow(L_(state), 2.0) * dL_(state);
}

Symmetric WalkerKinematicHardening::d_ratep_d_adot(VariableState & state)
{
  return dc_(state) * (2.0/3.0 * state.g - b_(state) / L_(state));;
}

Symmetric WalkerKinematicHardening::d_ratep_d_D(VariableState & state)
{
  return Symmetric::zero();
}

SymSymR4 WalkerKinematicHardening::d_ratep_d_s(VariableState & state)
{
  return -c_(state) * db_ds_(state) / L_(state);
}

SymSymR4 WalkerKinematicHardening::d_ratep_d_g(VariableState & state)
{
  return c_(state) * 2.0/3.0 * SymSymR4::id();
}

Symmetric WalkerKinematicHardening::ratet(VariableState & state)
{
  return -scale_->value(state.T) * x0_->value(state.T) * 
      softening_->phi(state.a, state.T) * 
      std::pow(std::sqrt(3.0/2) * state.h.norm() / state.D, x1_->value(state.T))
      * state.h / (state.h.norm() * std::sqrt(3.0/2));
}

SymSymR4 WalkerKinematicHardening::d_ratet_d_h(VariableState & state)
{
  Symmetric d = state.h / state.h.norm();
  SymSymR4 dd = 1.0 / (std::sqrt(3.0/2.0) * state.h.norm()) * (SymSymR4::id() - 
                                                               douter(d, d));
  return -scale_->value(state.T) * x0_->value(state.T) * 
      softening_->phi(state.a, state.T) * 
      (x1_->value(state.T)/state.D * std::pow(std::sqrt(3.0/2) * state.h.norm() / state.D,
                               x1_->value(state.T) - 1.0)  * douter(d,d)
       + std::pow(std::sqrt(3.0/2) * state.h.norm() / state.D,
                  x1_->value(state.T)) * dd);
}

Symmetric WalkerKinematicHardening::d_ratet_d_a(VariableState & state) 
{
  return -scale_->value(state.T) * x0_->value(state.T) * 
      softening_->dphi(state.a, state.T) * 
      std::pow(std::sqrt(3.0/2) * state.h.norm() / state.D, x1_->value(state.T))
      * state.h / (state.h.norm() * std::sqrt(3.0/2));
}

Symmetric WalkerKinematicHardening::d_ratet_d_adot(VariableState & state)
{
  return Symmetric::zero();
}

Symmetric WalkerKinematicHardening::d_ratet_d_D(VariableState & state)
{
  return scale_->value(state.T) * x0_->value(state.T) * 
      softening_->phi(state.a, state.T) * x1_->value(state.T) *
      std::pow(std::sqrt(3.0/2) * state.h.norm() / state.D,
               x1_->value(state.T)-1.0)
      * state.h / (state.h.norm() * std::sqrt(3.0/2))
      * std::sqrt(3.0/2) * state.h.norm() / (state.D*state.D);
}

SymSymR4 WalkerKinematicHardening::d_ratet_d_s(VariableState & state)
{
  return SymSymR4::zero();
}

SymSymR4 WalkerKinematicHardening::d_ratet_d_g(VariableState & state)
{
  return SymSymR4::zero();
}

double WalkerKinematicHardening::c_(VariableState & state)
{
  return c0_->value(state.T) + c1_->value(state.T) * 
      std::pow(state.adot, 1.0/c2_->value(state.T));
}

double WalkerKinematicHardening::dc_(VariableState & state)
{
  return c1_->value(state.T)/c2_->value(state.T) * 
      std::pow(state.adot, 1.0/c2_->value(state.T)-1.0);
}

double WalkerKinematicHardening::L_(VariableState & state)
{
  return l_->value(state.T) * (l1_->value(state.T) + (1.0 - l1_->value(state.T))
                               * std::exp(-l0_->value(state.T) * state.a));
}

double WalkerKinematicHardening::dL_(VariableState & state)
{
  return -l0_->value(state.T) * l_->value(state.T) * (1.0 - l1_->value(state.T))
                               * std::exp(-l0_->value(state.T) * state.a);
}

Symmetric WalkerKinematicHardening::n_(VariableState & state)
{
  return (3.0/2.0) * (state.s.dev() - state.h) / 
      (std::sqrt(3.0/2.0) * (state.s.dev() - state.h).norm());
}

Symmetric WalkerKinematicHardening::b_(VariableState & state)
{
  Symmetric n = n_(state);
  return (1.0 - b0_->value(state.T)) * state.h + 2.0/3.0 * b0_->value(state.T) *
      douter(n,n).dot(state.h);
}

SymSymR4 WalkerKinematicHardening::db_ds_(VariableState & state)
{
  SymSymR4 Nbar = dN_(state).dot(SymSymR4::id_dev());
  Symmetric n = n_(state);

  return 2.0/3.0 * b0_->value(state.T) * (Nbar * n.contract(state.h) + 
                                          douter(n,
                                                 Nbar.dot(state.h).transpose()));
}

SymSymR4 WalkerKinematicHardening::db_dx_(VariableState & state)
{
  Symmetric n = n_(state);
  SymSymR4 N = dN_(state);
  return (1.0 - b0_->value(state.T)) * SymSymR4::id() 
      + 2.0/3.0 * b0_->value(state.T) * douter(n,n)
      - 2.0/3.0 * b0_->value(state.T) * (N * n.contract(state.h) + 
                                         douter(n, 
                                                N.dot(state.h).transpose()));
}

SymSymR4 WalkerKinematicHardening::dN_(VariableState & state)
{
  Symmetric d = state.s.dev() - state.h;
  double dn = d.norm();
  return std::sqrt(3.0/2.0) / dn * (SymSymR4::id() - douter(d/dn,d/dn));
}

WrappedViscoPlasticFlowRule::WrappedViscoPlasticFlowRule() :
    stored_hist_(false)
{

}

History WrappedViscoPlasticFlowRule::blank_hist_() const
{
  return stored_hist_;
}

History WrappedViscoPlasticFlowRule::create_blank_hist_() const
{
  History h;
  populate_hist(h);
  return h;
}

History WrappedViscoPlasticFlowRule::gather_hist_(double * const h) const
{
  History hv = blank_hist_();
  hv.set_data(h);
  return hv;
}

History WrappedViscoPlasticFlowRule::gather_hist_(const double * const h) const
{
  History hv = blank_hist_();
  hv.set_data(const_cast<double*>(h));
  return hv;
}

State WrappedViscoPlasticFlowRule::make_state_(const double * const s, const double *
                                               const alpha, double T) const
{
  return State(Symmetric(s), gather_hist_(alpha), T);
}

size_t WrappedViscoPlasticFlowRule::nhist() const
{
  return blank_hist_().size();
}

int WrappedViscoPlasticFlowRule::init_hist(double * const h) const
{
  // Pointless memory error
  std::fill(h, h+nhist(), 0.0);
  // Actual stuff
  History hv = gather_hist_(h);
  initialize_hist(hv);
  return 0;
}

// Rate rule
int WrappedViscoPlasticFlowRule::y(const double* const s, const double* const alpha, double T,
              double & yv) const
{
  y(make_state_(s, alpha, T), yv);
  return 0;
}

int WrappedViscoPlasticFlowRule::dy_ds(const double* const s, const double* const alpha, double T,
              double * const dyv) const
{
  Symmetric res(dyv);
  dy_ds(make_state_(s, alpha, T), res);
  return 0;
}

int WrappedViscoPlasticFlowRule::dy_da(const double* const s, const double* const alpha, double T,
              double * const dyv) const
{
  // This is transposed, but that doesn't matter as it's flat
  History res = gather_derivative_<double>(dyv);
  dy_da(make_state_(s, alpha, T), res);
  return 0;
}

// Flow rule
int WrappedViscoPlasticFlowRule::g(const double * const s, const double * const alpha, double T,
              double * const gv) const
{
  Symmetric res(gv);
  g(make_state_(s, alpha, T), res);
  return 0;
}

int WrappedViscoPlasticFlowRule::dg_ds(const double * const s, const double * const alpha, double T,
              double * const dgv) const
{
  SymSymR4 res(dgv);
  dg_ds(make_state_(s, alpha, T), res);
  return 0;
}

int WrappedViscoPlasticFlowRule::dg_da(const double * const s, const double * const alpha, double T,
             double * const dgv) const
{
  // This is transposed and it does matter
  double * temp = new double [nhist() * 6];
  History res = gather_derivative_<Symmetric>(temp);
  dg_da(make_state_(s, alpha, T), res);
  
  for (size_t i = 0; i < nhist(); i++)
    for (size_t j = 0; j < 6; j++)
      dgv[CINDEX(j,i,nhist())] = temp[CINDEX(i,j,6)];
  delete [] temp;

  return 0;
}

// Hardening rule
int WrappedViscoPlasticFlowRule::h(const double * const s, const double * const alpha, double T,
              double * const hv) const
{
  History res = gather_hist_(hv);
  h(make_state_(s, alpha, T), res);
  return 0;
}

int WrappedViscoPlasticFlowRule::dh_ds(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  History res = gather_derivative_<Symmetric>(dhv);
  dh_ds(make_state_(s, alpha, T), res);
  return 0;
}

int WrappedViscoPlasticFlowRule::dh_da(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  History res = gather_derivative_<History>(dhv);
  dh_da(make_state_(s, alpha, T), res);
  return 0;
}

// Hardening rule wrt time
int WrappedViscoPlasticFlowRule::h_time(const double * const s, const double * const alpha, double T,
              double * const hv) const
{
  History res = gather_hist_(hv);
  h_time(make_state_(s, alpha, T), res);
  return 0;
}

void WrappedViscoPlasticFlowRule::h_time(const State & state, History & res) const
{
  res.zero();
}

int WrappedViscoPlasticFlowRule::dh_ds_time(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  History res = gather_derivative_<Symmetric>(dhv);
  dh_ds_time(make_state_(s, alpha, T), res);
  return 0;
}

void WrappedViscoPlasticFlowRule::dh_ds_time(const State & state, History & res) const
{
  res.zero();
}

int WrappedViscoPlasticFlowRule::dh_da_time(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  History res = gather_derivative_<History>(dhv);
  dh_da_time(make_state_(s, alpha, T), res);
  return 0;
}

void WrappedViscoPlasticFlowRule::dh_da_time(const State & state, History & res) const
{
  res.zero();
}

// Hardening rule wrt temperature
int WrappedViscoPlasticFlowRule::h_temp(const double * const s, const double * const alpha, double T,
              double * const hv) const
{
  History res = gather_hist_(hv);
  h_temp(make_state_(s, alpha, T), res);
  return 0;
}

void WrappedViscoPlasticFlowRule::h_temp(const State & state, History & res) const
{
  res.zero();
}

int WrappedViscoPlasticFlowRule::dh_ds_temp(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  History res = gather_derivative_<Symmetric>(dhv);
  dh_ds_temp(make_state_(s, alpha, T), res);
  return 0;
}

void WrappedViscoPlasticFlowRule::dh_ds_temp(const State & state, History & res) const
{
  res.zero();
}

int WrappedViscoPlasticFlowRule::dh_da_temp(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  History res = gather_derivative_<History>(dhv);
  dh_da_temp(make_state_(s, alpha, T), res);
  return 0;
}

void WrappedViscoPlasticFlowRule::dh_da_temp(const State & state, History & res) const
{
  res.zero();
}

TestFlowRule::TestFlowRule(double eps0, double D, double n, double s0, double K)
  : WrappedViscoPlasticFlowRule(), eps0_(eps0), D_(D), n_(n), s0_(s0), K_(K)
{
  populate_hist(stored_hist_);
}

std::string TestFlowRule::type()
{
  return "TestFlowRule";
}

std::unique_ptr<NEMLObject> TestFlowRule::initialize(ParameterSet & params)
{
  return neml::make_unique<TestFlowRule>(
      params.get_parameter<double>("eps0"),
      params.get_parameter<double>("D"),
      params.get_parameter<double>("n"),
      params.get_parameter<double>("s0"),
      params.get_parameter<double>("K")
      );
}

ParameterSet TestFlowRule::parameters()
{
  ParameterSet pset(TestFlowRule::type());

  pset.add_parameter<double>("eps0");
  pset.add_parameter<double>("D");
  pset.add_parameter<double>("n");
  pset.add_parameter<double>("s0");
  pset.add_parameter<double>("K");

  return pset;
}

void TestFlowRule::populate_hist(History & h) const
{
  h.add<double>("alpha");
  h.add<double>("iso");
}

void TestFlowRule::initialize_hist(History & h) const
{
  h.get<double>("alpha") = 0.0;
  h.get<double>("iso") = s0_;
}

void TestFlowRule::y(const State & state, double & res) const
{
  double h = (std::sqrt(3.0/2.0) * state.S.dev().norm() - state.h.get<double>("iso"))
      / D_;
  if (h > 0.0) {
    res = eps0_ * std::pow(h, n_);
  }
  else {
    res = 0;
  }
}

void TestFlowRule::dy_ds(const State & state, Symmetric & res) const
{
  double h = (std::sqrt(3.0/2.0) * state.S.dev().norm() - state.h.get<double>("iso"))
      / D_;
  if (h > 0.0) {
    res = eps0_ * n_ * std::pow(h, n_-1.0) * std::sqrt(3.0/2.0) * state.S.dev()
        / state.S.dev().norm() / D_;
  }
  else {
    res = Symmetric::zero();
  }
}

void TestFlowRule::dy_da(const State & state, History & res) const
{
  double h = (std::sqrt(3.0/2.0) * state.S.dev().norm() - state.h.get<double>("iso"))
      / D_;
  res.zero();
  if (h > 0.0) {
    res.get<double>("iso") = -eps0_ * n_ * std::pow(h, n_-1.0) / D_;
  }
}

void TestFlowRule::g(const State & state, Symmetric & res) const
{
  double ns = state.S.dev().norm();
  if (ns > 0)
    res = std::sqrt(3.0 / 2.0) * state.S.dev() / ns;
  else
    res = Symmetric::zero();
}

void TestFlowRule::dg_ds(const State & state, SymSymR4 & res) const
{
  Symmetric s = state.S.dev();
  double ns = s.norm();
  Symmetric sn = s / ns;
  res = std::sqrt(3.0 / 2.0) / ns * (SymSymR4::id_dev() - douter(sn, sn));
}

void TestFlowRule::dg_da(const State & state, History & res) const
{
  res.zero();
}

void TestFlowRule::h(const State & state, History & res) const
{
  res.get<double>("alpha") = 1.0;
  res.get<double>("iso") = K_;
}

void TestFlowRule::dh_ds(const State & state, History & res) const
{
  res.zero();
}

void TestFlowRule::dh_da(const State & state, History & res) const
{
  res.zero();
}

}
