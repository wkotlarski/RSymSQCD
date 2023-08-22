#ifndef MRSSM_H_
#define MRSSM_H_

#include <array>

struct MRSSMParameters {
   double MassTop;
   double MassGlu;
   double MasssigmaO;
   double MassSq;
   double eta_sign;
   double delta;
   double WidthGlu;
};

class MRSSM {

public:
   MRSSM() = delete;
   MRSSM(MRSSMParameters const& params);

   // pp -> suL suR
   // ------------------------------------------------------------------------------------

   // Born
   double sigmaTree_uu_suLsuR(double, double, double = 0.) const;
   double matrixTree_uu_suLsuR(double, double, double, int) const;

   // virtual
   double matrixVirt_uu_suLsuR(double, double, double, double, int, int, double) const;
   double matrixVirt_ud_suLsdR(double, double, double, double, int, int, double) const;

   // soft
   double matrixSoft_uu_suLsuRg_dp(double, double, double) const;
   double matrixSoft_uu_suLsuRg_sp(double, double, double, double, double) const;
   double matrixSoft_uu_suLsuRg_finite(double, double, double, double, double) const;

   // hard
   double matrixHard_uu_suLsuRg(
      double, std::array<std::array<double, 4>, 5> const &) const;
   double matrixHard_gu_suLsuRubar(
      double, std::array<std::array<double, 4>, 5> const &) const;
   double matrixHard_gu_suLsuRubar_DR(
      double, std::array<std::array<double, 4>, 5> const &) const;
   double matrixHard_gu_suLsuRubar_DS(double, std::array<std::array<double, 4>, 5> const&) const;
   double matrixHard_gu_suLsuRubar_DS_unsimp(double, std::array<std::array<double, 4>, 5> const&) const;
   // double matrixHard_gu_suLsuRubar_DR_wEtaDep(double, std::array<std::array<double, 4>, 5> const&) const;

   // pp > suL suL*
   // ------------------------------------------------------------------------------------
   double matrixVirt_uubar_suLsuLdagger(double, double, double, double, int, int, double) const;
   double matrixVirt_ddbar_suLsuLdagger(double, double, double, double, int, int, double) const;
   double matrixVirt_gg_suLsuLdagger(double, double, double, double, int, int, double) const;
   double matrixTree_ddbar_suLsuLdagger(double, double, double, int) const;
   double matrixTree_uubar_suLsuLdagger(double, double, double, int) const;
   double matrixTree_udbar_suLsdLdagger(double, double, double, int) const;
   double matrixTree_gg_suLsuLdagger(double, double, double, int) const;

   double sigmaTree_uubar_suLsuLdagger(double, double, double = 0.) const;
   double sigmaTree_ddbar_suLsuLdagger(double, double, double = 0.) const;
   double sigmaTree_gg_suLsuLdagger(double, double, double = 0.) const;
   double matrixSoft_gg_suLsuLdaggerg_dp(double, double, double) const;
   double matrixSoft_gg_suLsuLdaggerg_sp(double, double, double, double, double) const;
   double matrixSoft_gg_suLsuLdaggerg_finite(double, double, double, double, double) const;
   double matrixSoft_ddbar_suLsuLdaggerg_dp(double, double, double) const;
   double matrixSoft_ddbar_suLsuLdaggerg_sp(double, double, double, double, double) const;
   double matrixSoft_ddbar_suLsuLdaggerg_finite(double, double, double, double, double) const;
   double matrixSoft_uubar_suLsuLdaggerg_dp(double, double, double) const;
   double matrixSoft_uubar_suLsuLdaggerg_sp(double, double, double, double, double) const;
   double matrixSoft_uubar_suLsuLdaggerg_finite(double, double, double, double, double) const;
   double matrixHard_gg_suLsuLdaggerg(double, std::array<std::array<double, 4>, 5> const&) const;
   double matrixHard_uubar_suLsuLdaggerg(double, std::array<std::array<double, 4>, 5> const&) const;
   double matrixHard_ddbar_suLsuLdaggerg(double, std::array<std::array<double, 4>, 5> const&) const;
   double matrixHard_gd_suLsuLdaggerd(double, std::array<std::array<double, 4>, 5> const&) const;
   double matrixHard_gu_suLsuLdaggeru_DR(double, std::array<std::array<double, 4>, 5> const&) const;
   // double matrixHard_gu_suLsuLdaggeru_DR_wEta(double, std::array<std::array<double, 4>, 5> const&) const;
   double matrixHard_gu_suLsuLdaggeru_DS(double, std::array<std::array<double, 4>, 5> const&) const;
   double matrixHard_gu_suLsuLdaggeru(double, std::array<std::array<double, 4>, 5> const&) const;
   double matrixHard_gubar_suLsuLdaggerubar(double, std::array<std::array<double, 4>, 5> const&) const;
   double matrixHard_gubar_suLsuLdaggerubar_DR(double, std::array<std::array<double, 4>, 5> const&) const;
   double matrixHard_gdbar_suLsuLdaggerdbar(double, std::array<std::array<double, 4>, 5> const&) const;

   // pp -> glglbar
   // ------------------------------------------------------------------------------------
   double matrixTree_uubar_glglbar(double, double, double, int) const;
   double matrixTree_gg_glglbar(double, double, double, int) const;

private:
      // gauge vector for DR ME
      const std::array<double, 4> eta;

      // particle masses
      const double MasssigmaO;
      const double MassphiO;
      const double MassGlu;
      const double MassTop;
      const double MassSq;
      const double WidthGlu;
};

#endif // MRSSM_H_
