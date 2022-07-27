#ifndef SGLUONS_H_
#define SGLUONS_H_

#include <boost/property_tree/ptree.hpp>

class Sgluons {
public:
   Sgluons(boost::property_tree::ptree const&);
      // double matrixSimplifiedSoft_uubar_OOg( double, double );
      // double matrixSimplifiedHard_uubar_OOg(std::array<std::array<double, 4>, 5> const&) const;
      // double matrixSimplifiedSoft_gg_OOg( double, double );
      // double matrixSimplifiedHard_gg_OOg(std::array<std::array<double, 4>, 5> const&) const;
      // pp > OO
      //double matrixMRSSMSoft_qqbar_OOg( double, double );
      //double matrixMRSSMSoft_gg_OOg( double, double );
      // double sigmaMRSSMTree_uubar_OO(double, double );
      double matrixSgluonsTree_qqbar_OO(double, double, double) {return 0.;}
      double matrixSgluonsTree_gg_OO(double, double, double);

private:
   double mO;
};

#endif
