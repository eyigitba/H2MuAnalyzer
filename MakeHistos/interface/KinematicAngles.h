# ifndef KINEMATIC_ANGLES
# define KINEMATIC_ANGLES

#include "TLorentzVector.h"
#include "TVector3.h"

class KinAngles{
  TLorentzVector mu1_vec, mu2_vec, lep1_vec, lep2_vec, H_vec, Z2_vec, Z1_vec;
  TVector3 Z1_boost;
  public:
    KinAngles();
    KinAngles(TLorentzVector mu1, TLorentzVector mu2, TLorentzVector lep1, TLorentzVector lep2);

//    void SetMu1Vec(TLorentzVector mu1);
//    void SetMu2Vec(TLorentzVector mu2);
//    void SetLep1Vec(TLorentzVector lep1);
//    void SetLep2Vec(TLorentzVector lep2);

    float MELA_Cos_Theta1();
    float MELA_Cos_PhiH();
    float MELA_Cos_Phi1();
};


# endif
