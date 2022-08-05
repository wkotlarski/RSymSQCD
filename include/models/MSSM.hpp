#ifndef MSSM_H_
#define MSSM_H_

#include <array>

struct MSSMParameters {
   double MassTop;
   double MassGlu;
   double MassSq;
   double eta_sign;
   double delta;
   double WidthGlu;
};

class MSSM {

public:
   MSSM() = delete;
   MSSM(MSSMParameters const& params);

   double matrixMSSMTree_uubar_suLsuRdagger(double, double, double) const;
   double matrixMSSMVirt_uu_suLsuR(double, double, double, double, double, int, double) const;

private:
      // gauge vector for DR ME
      std::array<double, 4> eta;

      // particle masses
      const double MassGlu;
      const double MassTop;
      const double MassSq;
      const double WidthGlu;
};
#endif // MSSM_H_
