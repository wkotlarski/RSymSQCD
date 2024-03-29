#include "constants.hpp"
#include "mathematica_wrapper.hpp"

#include <Eigen/Dense>
#include "clooptools.h"

#include <cmath>
#include <iostream>

Process::Process(boost::property_tree::ptree const& pt) {

   MasssigmaO = pt.get<double>("masses.pseudoscalar_sgluon");
   MassphiO  = sqrt( pow(MasssigmaO,2) + 4.0 * pow(MassGlu, 2) );
   double eta_sign = pt.get<double>("technical parameters.eta_sign");
   double delta = pt.get<double>("technical parameters.delta");
   WidthGlu = pt.get<double>("technical parameters.WidthOverMass") * MassGlu;

   // choose a gage vector \eta for DR matrix elements
   eta = {std::sqrt(1.+Sqr(delta)), 0., delta, eta_sign};

/* -------------------- Squark-squark production---------------------------------*/

   // if(processID == "MRSSM,uu_suLsuR") {
      //matrixelementTree = &Process::matrixMSSMTree_uu_suLsuR; // same as in MSSM
      // // // sigmaPartTree1 = &Process::sigmaMRSSMTree_uu_suLsuR;
      // // splitting_kernel1 = &Pqq;
      // sigmaPartTree2 = &Process::sigmaMRSSMTree_uu_suLsuR;
      // splitting_kernel2 = &Pqq;
      //k = 2.*2*3*3;
      //h = 2.*2;
   // }
   // else if(processID == "MRSSM,ud_suLsdR") { // same as in MSSM
      // matrixelmatrixMRSSMTree_uubar_suLsuLdaggerementTree = &Process::matrixMSSMTree_ud_suLsdR;
      // matrixelementVirt = &Process::matrixMRSSMVirt_ud_suLsdR;
      // m1 = MassSquarks;
      // m2 = MassSquarks;
      // result doubled up, as there is ud and du initial state
      // flav.push_back({2, 1, 2});
      // k = 2.*2*3*3;
      // h = 2.*2;
   // }
   // else if (processID == "MRSSM,qq_sqLsqR") {
      // matrixelementTree = &Process::matrixMSSMTree_uu_suLsuR; // same as in MSSM
      // matrixelementVirt = &Process::matrixMRSSMVirt_uu_suLsuR;
      // matrixelementReal_SC = &Process::matrixMRSSMSoft_uu_suLsuRg;
      // sigmaPartTree1 = &Process::sigmaMRSSMTree_uu_suLsuR;
      // splitting_kernel1 = &Pqq;
      // sigmaPartTree2 = &Process::sigmaMRSSMTree_uu_suLsuR;
      // splitting_kernel2 = &Pqq;
      // matrixelementReal_HnonC = &Process::matrixMRSSMHard_uu_suLsuRg;
      // m1 = MassSquarks;
      // m2 = MassSquarks;
      // for (int i : {1, 2, 3, 4, 5}) {
         // flav.push_back({ i,  i, 1});
         // flav.push_back({-i, -i, 1});
      // }
      // k = 2.*2*3*3;
      // h = 2.*2;
   // }
   // else if(processID == "MRSSM,qq'_sqLsqR'") {
      // matrixelementTree = &Process::matrixMSSMTree_uu_suLsuR;
      // matrixelementVirt = &Process::matrixMRSSMVirt_ud_suLsdR;
      // matrixelementReal_SC = &Process::matrixMRSSMSoft_uu_suLsuRg;
      // sigmaPartTree1 = &Process::sigmaMRSSMTree_uu_suLsuR;
      // splitting_kernel1 = &Pqq;
      // sigmaPartTree2 = &Process::sigmaMRSSMTree_uu_suLsuR;
      // splitting_kernel2 = &Pqq;
      // matrixelementReal_HnonC = &Process::matrixMRSSMHard_uu_suLsuRg;
      // m1 = MassSquarks;
      // m2 = MassSquarks;
      // result doubled up, as there is ud and du initial state
      // for (int i : {1, 2, 3, 4, 5}) {
         // for (int j : {1, 2, 3, 4, 5}) {
            // if (j <= i) continue;
            // flav.push_back({ i,  j, 2});
            // flav.push_back({-i, -j, 2});
         // }
      // }
      // k = 2.*2*3*3;
      // h = 2.*2;
   // }
    /* squark production, MSSM */
   // else if(processID == "MSSM,uu_suLsuR") {
      // matrixelementTree = &Process::matrixMSSMTree_uu_suLsuR;
      // matrixelementVirt = &Process::matrixMSSMVirt_uu_suLsuR;
      //matrixelementReal_SC = &Process::matrixMRSSMSoft_uu_suLsuRg;
      // matrixelementReal_SC = &Process::matrix_soft_stub;
      // flav.push_back({2, 2, 1});
      // m1 = MassSquarks;
      // m2 = MassSquarks;
      // k = 2.*2*3*3;
      // h = 2.*2;
   // }
   // else if(processID == "MSSM,uu_suLsuL") {
      // matrixelementTree = &Process::matrixMSSMTree_uu_suLsuL;
      //matrixelementVirt = &Process::matrixMSSMVirt_uu_suLsuL;
      //matrixelementReal_SC = &Process::matrixMRSSMSoft_uu_suLsuLg;
      // matrixelementReal_SC = &Process::matrix_soft_stub;
      // m1 = MassSquarks;
      // m2 = MassSquarks;
      // k = 2.*2*3*3;
      // h = 2.*2;
   // }
   // else if(processID == "MSSM,ud_suLsdR") {
      // matrixelementTree = &Process::matrixMSSMTree_ud_suLsdR;
      // matrixelementVirt = &Process::matrixMSSMVirt_ud_suLsdR;
      // m1 = MassSquarks;
      // m2 = MassSquarks;
      // k = 2.*2*3*3;
      // h = 2.*2;
   // }
   // else if(processID == "MSSM,ud_suLsdL") {
      // matrixelementTree = &Process::matrixMSSMTree_ud_suLsdL;
      //matrixelementVirt = &matrixMSSMVirt_ud_suLsdL;
      // m1 = MassSquarks;
      // m2 = MassSquarks;
      // k = 2.*2*3*3;
      // h = 2.*2;
   // }

/* -------------------- Squark-antisquark production ---------------------*/

   // else if(processID == "MRSSM,uubar_suLsuLdagger") {
      // sigmaPartTree1 = &Process::sigmaMRSSMTree_uubar_suLsuLdagger;
      // sigmaPartTree2 = &Process::sigmaMRSSMTree_uubar_suLsuLdagger;
      // splitting_kernel1 = &Pqq;
      // splitting_kernel2 = &Pqq;
   // }
   // else if(processID == "MRSSM,ddbar_suLsuLdagger") {
      // sigmaPartTree1 = &Process::sigmaMRSSMTree_ddbar_suLsuLdagger;
      // sigmaPartTree2 = &Process::sigmaMRSSMTree_ddbar_suLsuLdagger;
      // splitting_kernel1 = &Pqq;
      // splitting_kernel2 = &Pqq;

   // }
   // else if(processID == "MRSSM,udbar_suLsdLdagger") {
      // matrixelementTree = &Process::matrixMSSMTree_ud_suLsdR;  // matrix elements are identical
      // m1 = MassSquarks;
      // m2 = MassSquarks;
      // k = 2.*2*3*3;
      // h = 2.*2;
   // }
   // else if(processID == "MRSSM,GG_suLsuLdagger") {
      // sigmaPartTree1 = &Process::sigmaMRSSMTree_gg_suLsuLdagger;
      // sigmaPartTree2 = &Process::sigmaMRSSMTree_gg_suLsuLdagger;
      // splitting_kernel1 = &Pgg;
      // splitting_kernel2 = &Pgg;
   // }
   // else if(processID == "MRSSM,gq_suLsuLdagger") {
      // sigmaPartTree1 = &Process::sigmaMRSSMTree_ddbar_suLsuLdagger;
      // splitting_kernel1 = &Pqg;
      // sigmaPartTree2 = &Process::sigmaMRSSMTree_gg_suLsuLdagger;
      // splitting_kernel2 = &Pgq;
   // }
   // else if(processID == "MRSSM,gu_suLsuLdagger") {
      // sigmaPartTree1 = &Process::sigmaMRSSMTree_uubar_suLsuLdagger;
      // splitting_kernel1 = &Pqg;

      // sigmaPartTree2 = &Process::sigmaMRSSMTree_gg_suLsuLdagger;
      // splitting_kernel2 = &Pgq;

      // matrixelementReal_SC = &Process::matrix_soft_stub;
      // m1 = MassSquarks;
      // m2 = MassSquarks;
      // if gluino cannot be on-shell use full real ME
      // if( pt.get<double>("collider setup.sqrt_S") < m1 +  MassGlu ||
         //  MassGlu < m2 ) {
         // matrixelementReal_HnonC = &Process::matrixMRSSMHard_gu_suLsuLdaggeru;
         // std::cout << "\nINFO: Not using any subtraction for the gu_suLsuLdagger channel.\n";
      // }
      // else {
         // otherwise, if WidthGlu < 0 use DR
         // if (WidthGlu < 0) {
            // matrixelementReal_HnonC = &Process::matrixMRSSMHard_gu_suLsuLdaggeru_DR;//_wEta;
            // std::cout << "\nINFO: Using diagram removal for the gu_suLsuLdagger channel.\n";
         // }
         // else use DS
         // else {
            // std::cout << "\nINFO: Using diagram subtraction for the gu_suLsuLdagger channel with Γ/M = "
                     //  << pt.get<double>("technical parameters.WidthOverMass") << ".\n";
            // matrixelementReal_HnonC = &Process::matrixMRSSMHard_gu_suLsuLdaggeru_DS;
         // }
      // }
   // }
   // else if(processID == "MRSSM,gu_suLsuR") {
      // sigmaPartTree1 = &Process::sigmaMRSSMTree_uu_suLsuR;
      // splitting_kernel1 = &Pqg;

      // sigmaPartTree2 = &Process::matrix_xsec_stub;
      // splitting_kernel2 = &Pgq;

      // matrixelementReal_SC = &Process::matrix_soft_stub;

      // m1 = MassSquarks;
      // m2 = MassSquarks;
      // if gluino cannot be on-shell use full real ME
      // if( pt.get<double>("collider setup.sqrt_S") < std::min(m1,m2) +  MassGlu ||
         //  MassGlu < m2 ) {
         // matrixelementReal_HnonC = &Process::matrixMRSSMHard_gu_suLsuRubar;
         // std::cout << "\nINFO: Not using any subtraction for the gu_suLsuRubar channel.\n";
      // }
      // else {
         // otherwise, if WidthGlu < 0 use DR
         // if( WidthGlu < 0) {
            // matrixelementReal_HnonC = &Process::matrixMRSSMHard_gu_suLsuRubar_DR_wEtaDep;
            // std::cout << "\nINFO: Using diagram removal for the gu_suLsuRubar channel.\n";
         // }
         // else use DS
         // else {
            // std::cout << "\nINFO: Using diagram subtraction for the gu_suLsuRubar channel with Γ/M = "
                     //  << pt.get<double>("technical parameters.WidthOverMass") << ".\n";
            // matrixelementReal_HnonC = &Process::matrixMRSSMHard_gu_suLsuRubar_DS;
         // }
      // }
      //flav.push_back({21, 2, 2});
   // }
   // else if(processID == "MRSSM,gd_sdLsdR") {
      // sigmaPartTree1 = &Process::sigmaMRSSMTree_uu_suLsuR;
      // splitting_kernel1 = &Pqg;

      // sigmaPartTree2 = &Process::matrix_xsec_stub;
      // splitting_kernel2 = &Pgq;

      // matrixelementReal_SC = &Process::matrix_soft_stub;

      // m1 = MassSquarks;
      // m2 = MassSquarks;
      // if gluino cannot be on-shell use full real ME
      // if( pt.get<double>("collider setup.sqrt_S") < std::min(m1,m2) +  MassGlu ||
         //  MassGlu < m2 ) {
         // matrixelementReal_HnonC = &Process::matrixMRSSMHard_gu_suLsuRubar;
         // std::cout << "\nINFO: Not using any subtraction for the gu_suLsuRubar channel.\n";
      // }
      // else {
         // otherwise, if WidthGlu < 0 use DR
         // if( WidthGlu < 0) {
            // matrixelementReal_HnonC = &Process::matrixMRSSMHard_gu_suLsuRubar_DR_wEtaDep;
            // std::cout << "\nINFO: Using diagram removal for the gu_suLsuRubar channel.\n";
         // }
         // else use DS
         // else {
            // std::cout << "\nINFO: Using diagram subtraction for the gu_suLsuRubar channel with Γ/M = "
                     //  << pt.get<double>("technical parameters.WidthOverMass") << ".\n";
            // matrixelementReal_HnonC = &Process::matrixMRSSMHard_gu_suLsuRubar_DS;
         // }
      // }
      // for (int el : {1, 2, 3, 4, 5}) {
         // flav.push_back({21, el, 20});
      // }
   // }
   // else if(processID == "MSSM,uubar_suLsuLdagger") {
      // matrixelementTree = &Process::matrixMSSMTree_uu_suLsuR;  // matrix elements are identical
      //std::vector<int> row {2, -2};
      // flav.push_back({2, -2, 2});
      // m1 = MassSquarks;
      // m2 = MassSquarks;
      // k = 2.*2*3*3;
      // h = 2.*2;
   // }
   // else if(processID == "MSSM,uubar_suLsuRdagger") {
      // matrixelementTree = &Process::matrixMSSMTree_uubar_suLsuRdagger;
      // flav.push_back({2, -2, 2} );
      // m1 = MassSquarks;
      // m2 = MassSquarks;
      // k = 2.*2*3*3;
      // h = 2.*2;
   // }
   // else if(processID == "MSSM,udbar_suLsdLdagger") {
      // matrixelementTree = &Process::matrixMSSMTree_ud_suLsdR;  // matrix elements are identical
      // m1 = MassSquarks;
      // m2 = MassSquarks;
      // k = 2.*2*3*3;
      // h = 2.*2;
   // }
   // else if(processID == "MSSM,udbar_suLsdRdagger") {
      // matrixelementTree = &Process::matrixMSSMTree_uubar_suLsuRdagger; // matrix elements are identical
      // m1 = MassSquarks;
      // m2 = MassSquarks;
      // k = 2.*2*3*3;
      // h = 2.*2;
   // }
   // else if(processID == "MSSM,ddbar_suLsuLdagger") {
      // matrixelementTree = &Process::matrixMRSSMTree_ddbar_suLsuLdagger; // matrix elements are identical
      // m1 = MassSquarks;
      // m2 = MassSquarks;
      // k = 2.*2*3*3;
      // h = 2.*2;
   // }
/*      else if(processID == "MSSM,GG_suLsuLdagger") {
      matrixelementTree = &Process::matrixMRSSMTree_GG_suLsuLdagger; // matrix elements are identical
      m1 = MassSquarks;
      m2 = MassSquarks;
      f1 = 0;
      f2 = 0;
      k = 2.*2*8*8;
      h = 1.;
   }
*/
}

#include "matrix_elements_and_xsections/MRSSM/mrssm_born_mes.cpp"
#include "matrix_elements_and_xsections/MRSSM/mrssm_born_xsecs.cpp"

/* -------------------- Squark-squark production: q+q > sq+sq ---------------------------------*/
// #include "matrix_elements_and_xsections/MSSM/mssm_uu_suLsuR_tree_matrix.cpp"
// #include "matrix_elements_and_xsections/MSSM/mssm_uu_suLsuL_tree_matrix.cpp"
// #include "matrix_elements_and_xsections/MSSM/mssm_ud_suLsdR_tree_matrix.cpp"
// #include "matrix_elements_and_xsections/MSSM/mssm_ud_suLsdL_tree_matrix.cpp"

/* -------------------- Squark-antisquark production: q+q^bar > sq+sq^dagger ---------------------*/
// #include "matrix_elements_and_xsections/MSSM/mssm_uubar_suLsuRdagger_tree_matrix.cpp"



/* -----------------------------------------------------------------------------------------------*/
/* ------------------------------------------ Virtual --------------------------------------------*/
/* -----------------------------------------------------------------------------------------------*/

/* ----------------------- Squark-squark production: q+q > sq+sq ---------------------------------*/
// mssm
// #include "matrix_elements_and_xsections/MSSM/mssm_uu_suLsuR_virt_matrix.cpp"
// #include "matrix_elements_and_xsections/MSSM/mssm_ud_suLsdR_virt_matrix.cpp"

//mrssm
#include "matrix_elements_and_xsections/MRSSM/mrssm_ud_suLsdR_virt_matrix.cpp"

// wEta_noSimplify is faster than wEta
#include "matrix_elements_and_xsections/MRSSM/mrssm_gu_suLsuLdaggeru_hard-DR_wEta_noSimplify.cpp"


#include "matrix_elements_and_xsections/MRSSM/mrssm_gu_suLsuRubar_hard-DR_wEtaDep.cpp"
