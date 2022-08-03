#ifndef PROCESS_H_
#define PROCESS_H_

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
   double sigmaMRSSMTree_uu_suLsuR(double, double );
   double matrixMRSSMTree_uu_suLsuR(double, double, double) const;

   // virtual
   double matrixMRSSMVirt_uu_suLsuR(double, double, double, double, double, int, double);
   double matrixMRSSMVirt_ud_suLsdR(double, double, double, double, double, int, double);

   // soft
   double matrixMRSSMSoft_uu_suLsuRg(double, double, double, double, double);

   // hard
   double matrixMRSSMHard_uu_suLsuRg(double, std::array<std::array<double, 4>, 5> const&) const noexcept;
   double matrixMRSSMHard_gu_suLsuRubar(double, std::array<std::array<double, 4>, 5> const&) const noexcept;
   double matrixMRSSMHard_gu_suLsuRubar_DR(double, std::array<std::array<double, 4>, 5> const&) const noexcept;
   double matrixMRSSMHard_gu_suLsuRubar_DS(double, std::array<std::array<double, 4>, 5> const&) const;
   double matrixMRSSMHard_gu_suLsuRubar_DS_unsimp(double, std::array<std::array<double, 4>, 5> const&) const;
   // double matrixMRSSMHard_gu_suLsuRubar_DR_wEtaDep(double, std::array<std::array<double, 4>, 5> const&) const;

   // pp > suL suL*
   // ------------------------------------------------------------------------------------
   double matrixMRSSMVirt_uubar_suLsuLdagger(double, double, double, double, double, int, double);
   double matrixMRSSMVirt_ddbar_suLsuLdagger(double, double, double, double, double, int, double);
   double matrixMRSSMVirt_GG_suLsuLdagger(double, double, double, double, double, int, double);
   double matrixMRSSMTree_ddbar_suLsuLdagger(double, double, double) const;
   double matrixMRSSMTree_uubar_suLsuLdagger(double, double, double) const;
   double matrixMRSSMTree_udbar_suLsdLdagger(double, double, double) const;
   double matrixMRSSMTree_GG_suLsuLdagger(double, double, double) const;

   double sigmaMRSSMTree_uubar_suLsuLdagger(double, double);
   double sigmaMRSSMTree_ddbar_suLsuLdagger(double, double);
   double sigmaMRSSMTree_gg_suLsuLdagger(double, double );
   double matrixMRSSMSoft_gg_suLsuLdaggerg(double, double, double, double, double);
   double matrixMRSSMSoft_ddbar_suLsuLdaggerg(double, double, double, double, double);
   double matrixMRSSMSoft_uubar_suLsuLdaggerg(double, double, double, double, double);
   double matrixMRSSMHard_gg_suLsuLdaggerg(double, std::array<std::array<double, 4>, 5> const&) const;
   double matrixMRSSMHard_uubar_suLsuLdaggerg(double, std::array<std::array<double, 4>, 5> const&) const;
   double matrixMRSSMHard_ddbar_suLsuLdaggerg(double, std::array<std::array<double, 4>, 5> const&) const;
   double matrixMRSSMHard_gd_suLsuLdaggerd(double, std::array<std::array<double, 4>, 5> const&) const;
   double matrixMRSSMHard_gu_suLsuLdaggeru_DR(double, std::array<std::array<double, 4>, 5> const&) const;
   // double matrixMRSSMHard_gu_suLsuLdaggeru_DR_wEta(double, std::array<std::array<double, 4>, 5> const&) const;
   double matrixMRSSMHard_gu_suLsuLdaggeru_DS(double, std::array<std::array<double, 4>, 5> const&) const;
   double matrixMRSSMHard_gu_suLsuLdaggeru(double, std::array<std::array<double, 4>, 5> const&) const;

   // pp -> glglbar
   // ------------------------------------------------------------------------------------
   double matrixMRSSMTree_uubar_glglbar(double, double, double) const;
   double matrixMRSSMTree_gg_glglbar(double, double, double) const;

private:
      // gauge vector for DR ME
      std::array<double, 4> eta;

      // particle masses
      const double MasssigmaO;
      const double MassphiO;
      const double MassGlu;
      const double MassTop;
      const double MassSq;
      const double WidthGlu;
};
#endif /* PROCESS_H_ */
