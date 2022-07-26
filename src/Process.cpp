#include "Process.hpp"
#include "splitting_kernels.hpp"

#include <Eigen/Dense>
#include "clooptools.h"

#include <cmath>
#include <iostream>

Process::Process(std::string const& processID, boost::property_tree::ptree const& pt) {

   MassTop = pt.get<double>("masses.top");
   MassGlu = pt.get<double>("masses.gluino");
   MasssigmaO = pt.get<double>("masses.pseudoscalar_sgluon");
   MassphiO  = sqrt( pow(MasssigmaO,2) + 4.0 * pow(MassGlu, 2) );
   MassSquarks = pt.get<double>("masses.squarks");
   double eta_sign = pt.get<double>("technical parameters.eta_sign");
   double delta = pt.get<double>("technical parameters.delta");
   WidthGlu = pt.get<double>("technical parameters.WidthOverMass") * MassGlu;

   // choose a gage vector \eta for DR matrix elements
   eta = {std::sqrt(1.+Sqr(delta)), 0., delta, eta_sign};

   // @todo remove
   MassSq = MassSquarks;
   partonic = false; // calculate sigma_had with |M|^2 and not sigma_part

/* -------------------- Squark-squark production---------------------------------*/

   if(processID == "MRSSM,uu_suLsuR") {
      //matrixelementTree = &Process::matrixMSSMTree_uu_suLsuR; // same as in MSSM
      //matrixelementVirt = &Process::matrixMRSSMVirt_uu_suLsuR;
      //matrixelementReal_SC = &Process::matrixMRSSMSoft_uu_suLsuRg;
      sigmaPartTree1 = &Process::sigmaMRSSMTree_uu_suLsuR;
      splitting_kernel1 = &Pqq;
      sigmaPartTree2 = &Process::sigmaMRSSMTree_uu_suLsuR;
      splitting_kernel2 = &Pqq;
      //m1 = MassSquarks;
      //m2 = MassSquarks;
      //flav.push_back({2, 2, 1});
      //k = 2.*2*3*3;
      //h = 2.*2;
   }
   else if(processID == "MRSSM,ud_suLsdR") { // same as in MSSM
      matrixelementTree = &Process::matrixMSSMTree_ud_suLsdR;
      matrixelementVirt = &Process::matrixMRSSMVirt_ud_suLsdR;
      m1 = MassSquarks;
      m2 = MassSquarks;
      // result doubled up, as there is ud and du initial state
      // flav.push_back({2, 1, 2});
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if (processID == "MRSSM,qq_sqLsqR") {
      matrixelementTree = &Process::matrixMSSMTree_uu_suLsuR; // same as in MSSM
      matrixelementVirt = &Process::matrixMRSSMVirt_uu_suLsuR;
      matrixelementReal_SC = &Process::matrixMRSSMSoft_uu_suLsuRg;
      sigmaPartTree1 = &Process::sigmaMRSSMTree_uu_suLsuR;
      splitting_kernel1 = &Pqq;
      sigmaPartTree2 = &Process::sigmaMRSSMTree_uu_suLsuR;
      splitting_kernel2 = &Pqq;
      // matrixelementReal_HnonC = &Process::matrixMRSSMHard_uu_suLsuRg;
      m1 = MassSquarks;
      m2 = MassSquarks;
      for (int i : {1, 2, 3, 4, 5}) {
         // flav.push_back({ i,  i, 1});
         // flav.push_back({-i, -i, 1});
      }
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MRSSM,qq'_sqLsqR'") {
      matrixelementTree = &Process::matrixMSSMTree_uu_suLsuR;
      matrixelementVirt = &Process::matrixMRSSMVirt_ud_suLsdR;
      matrixelementReal_SC = &Process::matrixMRSSMSoft_uu_suLsuRg;
      sigmaPartTree1 = &Process::sigmaMRSSMTree_uu_suLsuR;
      splitting_kernel1 = &Pqq;
      sigmaPartTree2 = &Process::sigmaMRSSMTree_uu_suLsuR;
      splitting_kernel2 = &Pqq;
      // matrixelementReal_HnonC = &Process::matrixMRSSMHard_uu_suLsuRg;
      m1 = MassSquarks;
      m2 = MassSquarks;
      // result doubled up, as there is ud and du initial state
      for (int i : {1, 2, 3, 4, 5}) {
         for (int j : {1, 2, 3, 4, 5}) {
            if (j <= i) continue;
            // flav.push_back({ i,  j, 2});
            // flav.push_back({-i, -j, 2});
         }
      }
      k = 2.*2*3*3;
      h = 2.*2;
   }
    /* squark production, MSSM */
   else if(processID == "MSSM,uu_suLsuR") {
      matrixelementTree = &Process::matrixMSSMTree_uu_suLsuR;
      matrixelementVirt = &Process::matrixMSSMVirt_uu_suLsuR;
      //matrixelementReal_SC = &Process::matrixMRSSMSoft_uu_suLsuRg;
      matrixelementReal_SC = &Process::matrix_soft_stub;
      // flav.push_back({2, 2, 1});
      m1 = MassSquarks;
      m2 = MassSquarks;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MSSM,uu_suLsuL") {
      matrixelementTree = &Process::matrixMSSMTree_uu_suLsuL;
      //matrixelementVirt = &Process::matrixMSSMVirt_uu_suLsuL;
      //matrixelementReal_SC = &Process::matrixMRSSMSoft_uu_suLsuLg;
      matrixelementReal_SC = &Process::matrix_soft_stub;
      m1 = MassSquarks;
      m2 = MassSquarks;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MSSM,ud_suLsdR") {
      matrixelementTree = &Process::matrixMSSMTree_ud_suLsdR;
      matrixelementVirt = &Process::matrixMSSMVirt_ud_suLsdR;
      m1 = MassSquarks;
      m2 = MassSquarks;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MSSM,ud_suLsdL") {
      matrixelementTree = &Process::matrixMSSMTree_ud_suLsdL;
      //matrixelementVirt = &matrixMSSMVirt_ud_suLsdL;
      m1 = MassSquarks;
      m2 = MassSquarks;
      k = 2.*2*3*3;
      h = 2.*2;
   }

/* -------------------- Squark-antisquark production ---------------------*/

   else if(processID == "MRSSM,uubar_suLsuLdagger") {
      sigmaPartTree1 = &Process::sigmaMRSSMTree_uubar_suLsuLdagger;
      sigmaPartTree2 = &Process::sigmaMRSSMTree_uubar_suLsuLdagger;
      matrixelementTree = &Process::matrixMRSSMTree_uubar_suLsuLdagger;
      matrixelementVirt = &Process::matrixMRSSMVirt_uubar_suLsuLdagger;
      matrixelementReal_SC = &Process::matrixMRSSMSoft_uubar_suLsuLdaggerg;
      splitting_kernel1 = &Pqq;
      splitting_kernel2 = &Pqq;
      // matrixelementReal_HnonC = &Process::matrixMRSSMHard_uubar_suLsuLdaggerg;
      m1 = MassSquarks;
      m2 = MassSquarks;
      // flav.push_back({2, -2, 2});
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MRSSM,ddbar_suLsuLdagger") {
      sigmaPartTree1 = &Process::sigmaMRSSMTree_ddbar_suLsuLdagger;
      sigmaPartTree2 = &Process::sigmaMRSSMTree_ddbar_suLsuLdagger;
      matrixelementTree = &Process::matrixMRSSMTree_ddbar_suLsuLdagger;
      matrixelementVirt = &Process::matrixMRSSMVirt_ddbar_suLsuLdagger;
      matrixelementReal_SC = &Process::matrixMRSSMSoft_ddbar_suLsuLdaggerg;
      splitting_kernel1 = &Pqq;
      splitting_kernel2 = &Pqq;
      // matrixelementReal_HnonC = &Process::matrixMRSSMHard_ddbar_suLsuLdaggerg;
      m1 = MassSquarks;
      m2 = MassSquarks;
      // flav.push_back({1, -1, 2});
      // flav.push_back({3, -3, 2});
      // flav.push_back({4, -4, 2});
      // flav.push_back({5, -5, 2});
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MRSSM,udbar_suLsdLdagger") {
      matrixelementTree = &Process::matrixMSSMTree_ud_suLsdR;  // matrix elements are identical
      m1 = MassSquarks;
      m2 = MassSquarks;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MRSSM,GG_suLsuLdagger") {
      sigmaPartTree1 = &Process::sigmaMRSSMTree_gg_suLsuLdagger;
      sigmaPartTree2 = &Process::sigmaMRSSMTree_gg_suLsuLdagger;
      matrixelementTree = &Process::matrixMRSSMTree_GG_suLsuLdagger; 
      matrixelementVirt = &Process::matrixMRSSMVirt_GG_suLsuLdagger;
      matrixelementReal_SC = &Process::matrixMRSSMSoft_gg_suLsuLdaggerg;
      splitting_kernel1 = &Pgg;
      splitting_kernel2 = &Pgg;
      // matrixelementReal_HnonC = &Process::matrixMRSSMHard_gg_suLsuLdaggerg;
      m1 = MassSquarks;
      m2 = MassSquarks;
      // flav.push_back({21, 21, 1});
      k = 2.*2*8*8;
      h = 1.;
   }
   else if(processID == "MRSSM,gq_suLsuLdagger") {
      sigmaPartTree1 = &Process::sigmaMRSSMTree_ddbar_suLsuLdagger;
      splitting_kernel1 = &Pqg;

      sigmaPartTree2 = &Process::sigmaMRSSMTree_gg_suLsuLdagger;
      splitting_kernel2 = &Pgq;

      matrixelementReal_SC = &Process::matrix_soft_stub;
      // matrixelementReal_HnonC = &Process::matrixMRSSMHard_gd_suLsuLdaggerd;

      m1 = MassSquarks;
      m2 = MassSquarks;
      // for( int el : { 1, -1, 3, -3, 4, -4, 5, -5}) flav.push_back({21, el, 2});
   }
   else if(processID == "MRSSM,gu_suLsuLdagger") {
      sigmaPartTree1 = &Process::sigmaMRSSMTree_uubar_suLsuLdagger;
      splitting_kernel1 = &Pqg;

      sigmaPartTree2 = &Process::sigmaMRSSMTree_gg_suLsuLdagger;
      splitting_kernel2 = &Pgq;

      matrixelementReal_SC = &Process::matrix_soft_stub;
      m1 = MassSquarks;
      m2 = MassSquarks;
      // if gluino cannot be on-shell use full real ME
      if( pt.get<double>("collider setup.sqrt_S") < m1 +  MassGlu ||
          MassGlu < m2 ) {
         // matrixelementReal_HnonC = &Process::matrixMRSSMHard_gu_suLsuLdaggeru;
         std::cout << "\nINFO: Not using any subtraction for the gu_suLsuLdagger channel.\n";
      }
      else {
         // otherwise, if WidthGlu < 0 use DR
         if (WidthGlu < 0) {
            // matrixelementReal_HnonC = &Process::matrixMRSSMHard_gu_suLsuLdaggeru_DR;//_wEta;
            std::cout << "\nINFO: Using diagram removal for the gu_suLsuLdagger channel.\n";
         }
         // else use DS
         else {
            std::cout << "\nINFO: Using diagram subtraction for the gu_suLsuLdagger channel with Γ/M = "
                      << pt.get<double>("technical parameters.WidthOverMass") << ".\n";
            // matrixelementReal_HnonC = &Process::matrixMRSSMHard_gu_suLsuLdaggeru_DS;
         }
      }
      // for( int el : { 2, -2 }) flav.push_back({21, el, 2});
   }
   else if(processID == "MRSSM,gu_suLsuR") {
      sigmaPartTree1 = &Process::sigmaMRSSMTree_uu_suLsuR;
      splitting_kernel1 = &Pqg;

      sigmaPartTree2 = &Process::matrix_xsec_stub;
      splitting_kernel2 = &Pgq;

      // matrixelementReal_SC = &Process::matrix_soft_stub;

      m1 = MassSquarks;
      m2 = MassSquarks;
      // if gluino cannot be on-shell use full real ME
      if( pt.get<double>("collider setup.sqrt_S") < std::min(m1,m2) +  MassGlu ||
          MassGlu < m2 ) {
         // matrixelementReal_HnonC = &Process::matrixMRSSMHard_gu_suLsuRubar;
         std::cout << "\nINFO: Not using any subtraction for the gu_suLsuRubar channel.\n";
      }
      else {
         // otherwise, if WidthGlu < 0 use DR
         if( WidthGlu < 0) {
            // matrixelementReal_HnonC = &Process::matrixMRSSMHard_gu_suLsuRubar_DR_wEtaDep;
            std::cout << "\nINFO: Using diagram removal for the gu_suLsuRubar channel.\n";
         }
         // else use DS
         else {
            std::cout << "\nINFO: Using diagram subtraction for the gu_suLsuRubar channel with Γ/M = "
                      << pt.get<double>("technical parameters.WidthOverMass") << ".\n";
            // matrixelementReal_HnonC = &Process::matrixMRSSMHard_gu_suLsuRubar_DS;
         }
      }
      //flav.push_back({21, 2, 2});
   }
   else if(processID == "MRSSM,gd_sdLsdR") {
      sigmaPartTree1 = &Process::sigmaMRSSMTree_uu_suLsuR;
      splitting_kernel1 = &Pqg;

      sigmaPartTree2 = &Process::matrix_xsec_stub;
      splitting_kernel2 = &Pgq;

      matrixelementReal_SC = &Process::matrix_soft_stub;

      m1 = MassSquarks;
      m2 = MassSquarks;
      // if gluino cannot be on-shell use full real ME
      if( pt.get<double>("collider setup.sqrt_S") < std::min(m1,m2) +  MassGlu ||
          MassGlu < m2 ) {
         // matrixelementReal_HnonC = &Process::matrixMRSSMHard_gu_suLsuRubar;
         std::cout << "\nINFO: Not using any subtraction for the gu_suLsuRubar channel.\n";
      }
      else {
         // otherwise, if WidthGlu < 0 use DR
         if( WidthGlu < 0) {
            // matrixelementReal_HnonC = &Process::matrixMRSSMHard_gu_suLsuRubar_DR_wEtaDep;
            std::cout << "\nINFO: Using diagram removal for the gu_suLsuRubar channel.\n";
         }
         // else use DS
         else {
            std::cout << "\nINFO: Using diagram subtraction for the gu_suLsuRubar channel with Γ/M = "
                      << pt.get<double>("technical parameters.WidthOverMass") << ".\n";
            // matrixelementReal_HnonC = &Process::matrixMRSSMHard_gu_suLsuRubar_DS;
         }
      }
      for (int el : {1, 2, 3, 4, 5}) {
         // flav.push_back({21, el, 20});
      }
   }
   else if(processID == "MSSM,uubar_suLsuLdagger") {
      matrixelementTree = &Process::matrixMSSMTree_uu_suLsuR;  // matrix elements are identical
      //std::vector<int> row {2, -2};
      // flav.push_back({2, -2, 2});
      m1 = MassSquarks;
      m2 = MassSquarks;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MSSM,uubar_suLsuRdagger") {
      matrixelementTree = &Process::matrixMSSMTree_uubar_suLsuRdagger;
      // flav.push_back({2, -2, 2} );
      m1 = MassSquarks;
      m2 = MassSquarks;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MSSM,udbar_suLsdLdagger") {
      matrixelementTree = &Process::matrixMSSMTree_ud_suLsdR;  // matrix elements are identical
      m1 = MassSquarks;
      m2 = MassSquarks;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MSSM,udbar_suLsdRdagger") {
      matrixelementTree = &Process::matrixMSSMTree_uubar_suLsuRdagger; // matrix elements are identical
      m1 = MassSquarks;
      m2 = MassSquarks;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MSSM,ddbar_suLsuLdagger") {
      matrixelementTree = &Process::matrixMRSSMTree_ddbar_suLsuLdagger; // matrix elements are identical
      m1 = MassSquarks;
      m2 = MassSquarks;
      k = 2.*2*3*3;
      h = 2.*2;
   }
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
/* -------------------- simpliefied models ---------------------*/

   else if(processID == "Simplified,uubar_OO") {
      matrixelementTree = &Process::matrix_tree_stub;
      matrixelementVirt = &Process::matrix_virt_stub;
      // matrixelementReal_SC = &Process::matrixSimplifiedSoft_uubar_OOg;
      sigmaPartTree1 = &Process::sigmaMRSSMTree_uubar_OO;
      splitting_kernel1 = &Pgg;
      sigmaPartTree2 = &Process::sigmaMRSSMTree_uubar_OO;
      splitting_kernel2 = &Pgg;
      // matrixelementReal_HnonC = &Process::matrixSimplifiedHard_uubar_OOg;
      m1 = MasssigmaO;
      m2 = MasssigmaO;
      for (int i : {1, 2, 3, 4, 5}) {
         // flav.push_back({i, -i, 2});
      }
   }
   else if(processID == "Simplified,gg_OO") {
      sigmaPartTree1 = &Process::matrix_xsec_stub;
      matrixelementTree = &Process::matrix_tree_stub;
      matrixelementVirt = &Process::matrix_virt_stub;
      // matrixelementReal_SC = &Process::matrixSimplifiedSoft_gg_OOg;
      // matrixelementReal_HnonC = &Process::matrixSimplifiedHard_gg_OOg;
      m1 = MasssigmaO;
      m2 = MasssigmaO;
      // flav.push_back({21, 21, 1});
   }
   else {
      std::cout << "Error! Subprocess " << processID << " not implemented.\n";
   }
}

/* --------------------------------------------------------------------------------------------*/
/* ------------------------------------------ Tree --------------------------------------------*/
/* --------------------------------------------------------------------------------------------*/


/* -------------------- Squark-squark production: q+q > sq+sq ---------------------------------*/
#include "matrix_elements_and_xsections/MSSM/mssm_uu_suLsuR_tree_matrix.cpp"
#include "matrix_elements_and_xsections/MSSM/mssm_uu_suLsuL_tree_matrix.cpp"
#include "matrix_elements_and_xsections/MSSM/mssm_ud_suLsdR_tree_matrix.cpp"
#include "matrix_elements_and_xsections/MSSM/mssm_ud_suLsdL_tree_matrix.cpp"

/* -------------------- Squark-antisquark production: q+q^bar > sq+sq^dagger ---------------------*/
#include "matrix_elements_and_xsections/MRSSM/mrssm_uubar_suLsuLdagger_tree_matrix.cpp"
#include "matrix_elements_and_xsections/MRSSM/mrssm_ddbar_suLsuLdagger_tree_matrix.cpp"
#include "matrix_elements_and_xsections/MRSSM/mrssm_udbar_suLsdLdagger_tree_matrix.cpp"
#include "matrix_elements_and_xsections/MSSM/mssm_uubar_suLsuRdagger_tree_matrix.cpp"

/* -------------------- Squark-antisquark production: G+G > sq+sq^dagger ---------------------*/
#include "matrix_elements_and_xsections/MRSSM/mrssm_gg_suLsuLdagger_tree_matrix.cpp"



/* /////////////////////////////////// partonic xsections /////////////////////////////////////// */
//#include "matrix_elements_and_xsections/mrssm_uu_suLsuR_tree_xsec.cpp"
#include "matrix_elements_and_xsections/MRSSM/mrssm_partonic_xsections.cpp"

/* -----------------------------------------------------------------------------------------------*/
/* ------------------------------------------ Virtual --------------------------------------------*/
/* -----------------------------------------------------------------------------------------------*/

/* ----------------------- Squark-squark production: q+q > sq+sq ---------------------------------*/
// mssm
#include "matrix_elements_and_xsections/MSSM/mssm_uu_suLsuR_virt_matrix.cpp"
#include "matrix_elements_and_xsections/MSSM/mssm_ud_suLsdR_virt_matrix.cpp"

//mrssm
#include "matrix_elements_and_xsections/MRSSM/mrssm_uu_suLsuR_virt_matrix.cpp"
#include "matrix_elements_and_xsections/MRSSM/mrssm_ud_suLsdR_virt_matrix.cpp"

/* --------------------- Squark-antisquark production: G+G > sq+sq^dagger ------------------------*/
//mrssm
#include "matrix_elements_and_xsections/MRSSM/mrssm_gg_suLsuLdagger_virt_matrix.cpp"

/* --------------------- Squark-antisquark production: q+qbar > sq+sq^dagger ------------------------*/
//mrssm
#include "matrix_elements_and_xsections/MRSSM/mrssm_ddbar_suLsuLdagger_virt_matrix.cpp"
#include "matrix_elements_and_xsections/MRSSM/mrssm_uubar_suLsuLdagger_virt_matrix.cpp"

/*

   Real emissions

 */
#include "matrix_elements_and_xsections/MRSSM/mrssm_uubar_suLsuLdaggerg_soft.cpp"
#include "matrix_elements_and_xsections/MRSSM/mrssm_uubar_suLsuLdaggerg_hard.cpp"

#include "matrix_elements_and_xsections/MRSSM/mrssm_ddbar_suLsuLdaggerg_soft.cpp"
#include "matrix_elements_and_xsections/MRSSM/mrssm_ddbar_suLsuLdaggerg_hard.cpp"

#include "matrix_elements_and_xsections/MRSSM/mrssm_gd_suLsuLdaggerd_hard.cpp"
// wEta_noSimplify is faster than wEta
#include "matrix_elements_and_xsections/MRSSM/mrssm_gu_suLsuLdaggeru_hard-DS.cpp"
#include "matrix_elements_and_xsections/MRSSM/mrssm_gu_suLsuLdaggeru_hard-DR.cpp"
#include "matrix_elements_and_xsections/MRSSM/mrssm_gu_suLsuLdaggeru_hard-DR_wEta_noSimplify.cpp"
#include "matrix_elements_and_xsections/MRSSM/mrssm_gu_suLsuLdaggeru_hard.cpp"

#include "matrix_elements_and_xsections/MRSSM/mrssm_gg_suLsuLdaggerg_hard.cpp"
#include "matrix_elements_and_xsections/MRSSM/mrssm_gg_suLsuLdaggerg_soft.cpp"

#include "matrix_elements_and_xsections/MRSSM/mrssm_uu_suLsuRg_soft.cpp"
#include "matrix_elements_and_xsections/MRSSM/mrssm_uu_suLsuRg_hard.cpp"
#include "matrix_elements_and_xsections/MRSSM/mrssm_gu_suLsuRubar_hard.cpp"
#include "matrix_elements_and_xsections/MRSSM/mrssm_gu_suLsuRubar_hard-DR.cpp"
#include "matrix_elements_and_xsections/MRSSM/mrssm_gu_suLsuRubar_hard-DR_wEtaDep.cpp"
#include "matrix_elements_and_xsections/MRSSM/mrssm_gu_suLsuRubar_hard-DS.cpp"

// #include "matrix_elements_and_xsections/simplified_uubar_OOg_soft.cpp"
// #include "matrix_elements_and_xsections/simplified_uubar_OOg_hard.cpp"

// #include "matrix_elements_and_xsections/simplified_gg_OOg_soft.cpp"
// #include "matrix_elements_and_xsections/simplified_gg_OOg_hard.cpp"
