/*
 *
 */

#include "XSection_Tree_MSSM.h"
#include "MSSM_Tree_uu_su1su4_EqualFlavor.h"

XSection_Tree_MSSM::XSection_Tree_MSSM() {

  // TODO Auto-generated constructor stub

}

XSection_Tree_MSSM::~XSection_Tree_MSSM() {
  // TODO Auto-generated destructor stub
}



int XSection_Tree_MSSM::integrand(const int *ndim, const cubareal xx[],
                             const int *ncomp, cubareal ff[], void *userdata) {


    static double MassSq = squark_mass.at(0).at(0);
    static double mu = MassSq;
    double x1 = xx[1];
    double x2 = xx[2];
    double s = S * x1 * x2;     //partonic 
    double Tmin = pow(MassSq,2) - s/2 - sqrt(pow(s,2)/4 -
                 pow(MassSq,2)*s);
    double Tmax = pow(MassSq,2) - s/2 + sqrt(pow(s,2)/4 -
                  pow(MassSq,2)*s);
    double T = xx[0]*(Tmax-Tmin) + Tmin;
    double jacobian = (Tmax-Tmin);
    double U = 2*pow(MassSq,2) - s - T;
 
    if(s < 4*pow(MassSq,2) || T*U < pow(MassSq,4))
    {
        ff[0] = 0;
        return 0;
    }
    // no prefactors
    double SquaredMReal = MsquaredMSSMTree_uu_suLsuR_EqualFlavor(pdf_nlo->alphasQ(MassSq), 
                          MassSq, gluino_mass, s, T, U);
    double dSigmaPart = SquaredMReal*4.*M_PI/(pow(4.*M_PI,2))/
                         (4.*9)/(pow(s,2));

    double dSigmaHad = dSigmaPart
                     * pdf_nlo->xfxQ(2,x1,mu)/x1
                     * pdf_nlo->xfxQ(2,x2,mu)/x2;

    ff[0] = dSigmaHad*jacobian*3.89379*pow(10,11);   // in femto barn
    return 1;
}


 std::array<double, 3> XSection_Tree_MSSM::integrate() {
    constexpr int ndim = 3;
    constexpr int ncomp = 1;
    constexpr int nvec = 1;
    constexpr double accuracy_rel = 1e-4;
    constexpr double accuracy_abs = 1e-12;
    constexpr int eval_min = 1000;
    constexpr int eval_max = 1000000;
    constexpr int verbose = 0;        // adjust shown output 0 ... 3
    constexpr int last = 4;
    constexpr int key = 0;
    int nregions, neval, fail;
    cubareal integral[ncomp], error[ncomp], prob[ncomp];

      Cuhre(ndim, ncomp, integrand, NULL, nvec,
      accuracy_rel, accuracy_abs, verbose | last,
      eval_min, eval_max, key, NULL, NULL,
      &nregions, &neval, &fail, integral, error, prob);

    std::array <double, 3> result{integral[0], error[0], prob[0]}; 

  //  int CuhreNeval = neval;
  //  int CuhreFail = fail;
  //  double CuhreResult = (double)integral[0];
  //  double CuhreError = (double)error[0];
  //  double CuhreProb = (double)prob[0];
  //  //double CuhreTime = duration(timeNow()-t0);
  //  printf("Cuhre Result:\tneval %d\tfail %d\n", CuhreNeval, 
  //         CuhreFail);
  //  printf("Cuhre Result:\t%.8e +- %.8e\tp = %.3e\n",
  //         CuhreResult, CuhreError, CuhreProb);
  //  //printf("Suave Time:\t %f s \n",SuaveTime);
  //  printf("\n");


    return result;
}
