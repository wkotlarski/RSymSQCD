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

   double matrixTree_uu_suLsuR(double, double, double) const;
   double sigmaTree_uu_suLsuR(double, double, double = 0.) const;
   double matrixVirt_uu_suLsuR(double, double, double, double, double, int, double) const;
   double matrixSoft_uu_suLsuRg(double, double, double, double, double) const;
   double matrixHard_uu_suLsuRg(double, std::array<std::array<double, 4>, 5> const &) const;

   double matrixTree_uubar_suLsuLdagger(double, double, double) const;
   double sigmaTree_uubar_suLsuLdagger(double, double, double = 0.) const;
   double matrixVirt_uubar_suLsuLdagger(double, double, double, double, double, int, double) const;
   double matrixSoft_uubar_suLsuLdaggerg(double, double, double, double, double) const;
   double matrixHard_uubar_suLsuLdaggerg(double, std::array<std::array<double, 4>, 5> const &) const;

   double matrixTree_ddbar_suLsuLdagger(double, double, double) const;
   double sigmaTree_ddbar_suLsuLdagger(double, double, double = 0.) const;
   double matrixVirt_ddbar_suLsuLdagger(double, double, double, double, double, int, double) const;
   double matrixSoft_ddbar_suLsuLdaggerg(double, double, double, double, double) const;
   double matrixHard_ddbar_suLsuLdaggerg(double, std::array<std::array<double, 4>, 5> const &) const;

   double matrixTree_gg_suLsuLdagger(double, double, double) const;
   double sigmaTree_gg_suLsuLdagger(double, double, double = 0.) const;
   double matrixVirt_gg_suLsuLdagger(double, double, double, double, double, int, double) const;
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
