#include "H2MuAnalyzer/MakeHistos/interface/KinematicAngles.h"

KinAngles::KinAngles() {
  mu1_vec.SetPxPyPzE(0,0,0,0);
  mu2_vec.SetPxPyPzE(0,0,0,0);
  lep1_vec.SetPxPyPzE(0,0,0,0);
  lep2_vec.SetPxPyPzE(0,0,0,0);
  H_vec.SetPxPyPzE(0,0,0,0);
  Z1_vec.SetPxPyPzE(0,0,0,0);
  Z2_vec.SetPxPyPzE(0,0,0,0);

  Z1_boost = Z1_vec.BoostVector();
}

KinAngles::KinAngles(TLorentzVector mu1, TLorentzVector mu2, TLorentzVector lep1, TLorentzVector lep2) {
  mu1_vec = mu1;
  mu2_vec = mu2;
  lep1_vec = lep1;
  lep2_vec = lep2;

  H_vec  = mu1_vec + mu2_vec;
  Z2_vec = lep1_vec + lep2_vec;
  Z1_vec = H_vec + Z2_vec;

  Z1_boost = Z1_vec.BoostVector();

  mu1_vec.Boost( -Z1_boost );
  mu2_vec.Boost( -Z1_boost );
  lep1_vec.Boost( -Z1_boost );
  lep2_vec.Boost( -Z1_boost );
  H_vec.Boost( -Z1_boost );
  Z2_vec.Boost( -Z1_boost );
}


float KinAngles::MELA_Cos_Theta1() {
  float cos_theta_1;

  cos_theta_1 = H_vec.CosTheta();
  return cos_theta_1;
}

float KinAngles::MELA_Cos_PhiH() {
  TVector3  z_axis, norm_z, norm_H;
  float cos_phi_H;

  z_axis.SetXYZ(0,0,1);
  norm_z = z_axis.Cross( H_vec.BoostVector() );
  norm_H = mu1_vec.BoostVector().Cross( mu2_vec.BoostVector() );
  cos_phi_H = norm_z * norm_H / (norm_z.Mag() * norm_H.Mag());
  return cos_phi_H;
}

float KinAngles::MELA_Cos_Phi1() {
  TVector3  norm_1, norm_H;
  float cos_phi_1;

  norm_1 = mu1_vec.BoostVector().Cross( mu2_vec.BoostVector() );
  norm_H = lep1_vec.BoostVector().Cross( lep2_vec.BoostVector() );
  cos_phi_1 = norm_1 * norm_H / (norm_1.Mag() * norm_H.Mag());
  return cos_phi_1;
}
