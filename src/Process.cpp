#include <Eigen/Dense>
#include "Process.hpp"

Process::Process(std::string processID, boost::property_tree::ptree pt) {
   
   MassTop = pt.get<double>("masses.top");
   MassGlu = pt.get<double>("masses.gluino");
   MasssigmaO = pt.get<double>("masses.pseudoscalar_sgluon");
   MassphiO  = sqrt( pow(MasssigmaO,2) + 4.0 * pow(MassGlu, 2) );
   MassSuL = pt.get<double>("masses.suL");
   MassSuR = pt.get<double>("masses.suR");
   MassSdL = pt.get<double>("masses.sdL");
   MassSdR = pt.get<double>("masses.sdR");
   MassSsL = pt.get<double>("masses.ssL");
   MassSsR = pt.get<double>("masses.ssR");
   MassScL = pt.get<double>("masses.scL");
   MassScR = pt.get<double>("masses.scR");
   MassSbL = pt.get<double>("masses.sbL");
   MassSbR = pt.get<double>("masses.sbR");
   MassStL = pt.get<double>("masses.stL");
   MassStR = pt.get<double>("masses.stR");
   mu_r = pt.get<double>("collider setup.mu_r");
   mu_f = pt.get<double>("collider setup.mu_f");
   dS = pt.get<double>("technical parameters.dS");
   double eta_sign = pt.get<double>("technical parameters.eta_sign");
   double delta = pt.get<double>("technical parameters.delta");
   WidthGlu = pt.get<double>("technical parameters.WidthOverMass") * MassGlu;
   pdf = LHAPDF::mkPDF( pt.get<std::string>("collider setup.pdf") , 0);

   // choose a gage vector \eta for DR matrix elements
   eta = {sqrt(1.+delta*delta), 0., delta, eta_sign};
   
   // @todo remove
   MassSq = MassSuL;
   partonic = false; // calculate sigma_had with |M|^2 and not sigma_part
   
/* -------------------- Squark-squark production---------------------------------*/

   if(processID == "MRSSM,uu_suLsuR") {
      matrixelementTree = &Process::matrixMSSMTree_uu_suLsuR; // same as in MSSM 
      matrixelementVirt = &Process::matrixMRSSMVirt_uu_suLsuR;
      matrixelementReal_SC = &Process::matrixMRSSMSoft_uu_suLsuRg;
      sigmaPartTree1 = &Process::sigmaMRSSMTree_uu_suLsuR;
      splitting_kernel1 = &Process::Pqq;
      sigmaPartTree2 = &Process::sigmaMRSSMTree_uu_suLsuR;
      splitting_kernel2 = &Process::Pqq;
      // matrixelementReal_HnonC = &Process::matrixMRSSMHard_uu_suLsuRg;
      m1 = MassSuL;
      m2 = MassSuR;
      flav.push_back( std::vector<int> {2, 2, 1} );
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MRSSM,ud_suLsdR") { // same as in MSSM
      matrixelementTree = &Process::matrixMSSMTree_ud_suLsdR;
      matrixelementVirt = &Process::matrixMRSSMVirt_ud_suLsdR;
      m1 = MassSuL;
      m2 = MassSuR;
      // result doubled up, as there is ud and du initial state
      flav.push_back( std::vector<int> {2, 1, 2} );
      k = 2.*2*3*3;
      h = 2.*2;
   }
    /* squark production, MSSM */
   else if(processID == "MSSM,uu_suLsuR") {
      matrixelementTree = &Process::matrixMSSMTree_uu_suLsuR;
      matrixelementVirt = &Process::matrixMSSMVirt_uu_suLsuR;
      //matrixelementReal_SC = &Process::matrixMRSSMSoft_uu_suLsuRg;
      matrixelementReal_SC = &Process::matrix_soft_stub;
      flav.push_back( std::vector<int> {2, 2, 1} );
      m1 = MassSuL;
      m2 = MassSuR;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MSSM,uu_suLsuL") {
      matrixelementTree = &Process::matrixMSSMTree_uu_suLsuL;
      //matrixelementVirt = &Process::matrixMSSMVirt_uu_suLsuL;
      //matrixelementReal_SC = &Process::matrixMRSSMSoft_uu_suLsuLg;
      matrixelementReal_SC = &Process::matrix_soft_stub;
      m1 = MassSuL;
      m2 = MassSuR;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MSSM,ud_suLsdR") {
      matrixelementTree = &Process::matrixMSSMTree_ud_suLsdR;
      matrixelementVirt = &Process::matrixMSSMVirt_ud_suLsdR;
      m1 = MassSuL;
      m2 = MassSuR;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MSSM,ud_suLsdL") {
      matrixelementTree = &Process::matrixMSSMTree_ud_suLsdL;
      //matrixelementVirt = &matrixMSSMVirt_ud_suLsdL;
      m1 = MassSuL;
      m2 = MassSuR;
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
      splitting_kernel1 = &Process::Pqq;
      splitting_kernel2 = &Process::Pqq;
      // matrixelementReal_HnonC = &Process::matrixMRSSMHard_uubar_suLsuLdaggerg;
      m1 = MassSuL;
      m2 = MassSuL;
      flav.push_back( std::vector<int> {2, -2, 2} );
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MRSSM,ddbar_suLsuLdagger") {
      sigmaPartTree1 = &Process::sigmaMRSSMTree_ddbar_suLsuLdagger;
      sigmaPartTree2 = &Process::sigmaMRSSMTree_ddbar_suLsuLdagger;
      matrixelementTree = &Process::matrixMRSSMTree_ddbar_suLsuLdagger;
      matrixelementVirt = &Process::matrixMRSSMVirt_ddbar_suLsuLdagger;
      matrixelementReal_SC = &Process::matrixMRSSMSoft_ddbar_suLsuLdaggerg;
      splitting_kernel1 = &Process::Pqq;
      splitting_kernel2 = &Process::Pqq;
      // matrixelementReal_HnonC = &Process::matrixMRSSMHard_ddbar_suLsuLdaggerg;
      m1 = MassSuL;
      m2 = MassSuL;
      flav.push_back( std::vector<int> {1, -1, 2} );
      flav.push_back( std::vector<int> {3, -3, 2} );
      flav.push_back( std::vector<int> {4, -4, 2} );
      flav.push_back( std::vector<int> {5, -5, 2} );
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MRSSM,udbar_suLsdLdagger") {
      matrixelementTree = &Process::matrixMSSMTree_ud_suLsdR;  // matrix elements are identical
      m1 = MassSuL;
      m2 = MassSuR;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MRSSM,GG_suLsuLdagger") {
      sigmaPartTree1 = &Process::sigmaMRSSMTree_gg_suLsuLdagger;
      sigmaPartTree2 = &Process::sigmaMRSSMTree_gg_suLsuLdagger;
      matrixelementTree = &Process::matrixMRSSMTree_GG_suLsuLdagger; 
      matrixelementVirt = &Process::matrixMRSSMVirt_GG_suLsuLdagger;
      matrixelementReal_SC = &Process::matrixMRSSMSoft_gg_suLsuLdaggerg;
      splitting_kernel1 = &Process::Pgg;
      splitting_kernel2 = &Process::Pgg;
      // matrixelementReal_HnonC = &Process::matrixMRSSMHard_gg_suLsuLdaggerg;
      m1 = MassSuL;
      m2 = MassSuL;
      flav.push_back( std::vector<int> {21, 21, 1} );
      k = 2.*2*8*8;
      h = 1.;
   }
   else if(processID == "MRSSM,gq_suLsuLdagger") {
      sigmaPartTree1 = &Process::sigmaMRSSMTree_ddbar_suLsuLdagger;
      splitting_kernel1 = &Process::Pqg;
      
      sigmaPartTree2 = &Process::sigmaMRSSMTree_gg_suLsuLdagger;
      splitting_kernel2 = &Process::Pgq;
      
      matrixelementReal_SC = &Process::matrix_soft_stub;
      // matrixelementReal_HnonC = &Process::matrixMRSSMHard_gd_suLsuLdaggerd;
      
      m1 = MassSuL;
      m2 = MassSuL;
      std::array<int, 8> partons = { 1, -1, 3, -3, 4, -4, 5, -5 };
      for( int el : partons) flav.push_back( std::vector<int> {21, el, 2} );
   }   
   else if(processID == "MRSSM,gu_suLsuLdagger") {
      sigmaPartTree1 = &Process::sigmaMRSSMTree_uubar_suLsuLdagger;
      splitting_kernel1 = &Process::Pqg;
      
      sigmaPartTree2 = &Process::sigmaMRSSMTree_gg_suLsuLdagger;
      splitting_kernel2 = &Process::Pgq;
      
      matrixelementReal_SC = &Process::matrix_soft_stub;
      m1 = MassSuL;
      m2 = MassSuL;
      // if gluino cannot be on-shell use full real ME
      if( pt.get<double>("collider setup.sqrt_S") < m1 +  MassGlu ||
          MassGlu < m2 ) {
         std::cout << "\nINFO: Not using any subtraction for the gu_suLsuLdagger channel.\n";
         matrixelementReal_HnonC = &Process::matrixMRSSMHard_gu_suLsuLdaggeru;
         matrixelementReal_HnonC_CSub1 = nullptr;
         matrixelementReal_HnonC_CSub2 = nullptr;
      }
      else {
         // otherwise, if WidthGlu < 0 use DR
         if( WidthGlu < 0) {
            matrixelementReal_HnonC = &Process::matrixMRSSMHard_gu_suLsuLdaggeru_DR2;//_wEta;
            matrixelementReal_HnonC_CSub1 = nullptr;
            matrixelementReal_HnonC_CSub2 = nullptr;
            std::cout << "\nINFO: Using diagram removal for the gu_suLsuLdagger channel.\n";
         }
         // else use DS
         else {
            std::cout << "\nINFO: Using diagram subtraction for the gu_suLsuLdagger channel with width/mass = "
               << pt.get<double>("technical parameters.WidthOverMass") << ".\n";
            /*
            matrixelementReal_HnonC = &Process::matrixMRSSMHard_gu_suLsuLdaggeru_DS;
            matrixelementReal_HnonC_CSub1 = nullptr;
            matrixelementReal_HnonC_CSub2 = &Process::matrixMRSSMHard_gu_suLsuLdaggeru_DS_CSub2;
           */

            matrixelementReal_HnonC = &Process::matrixMRSSMHard_gubar_suLsuLdaggerubar_DS;
            matrixelementReal_HnonC_CSub1 = &Process::matrixMRSSMHard_gubar_suLsuLdaggerubar_DS_CSub1;
            matrixelementReal_HnonC_CSub2 = nullptr;
         }
      }
      std::array<int, 2> partons { 2, -2 };
      for( int el : partons) flav.push_back( std::vector<int> {21, el, 2} );
   }
   else if(processID == "MRSSM,gu_suLsuR") {
      sigmaPartTree1 = &Process::sigmaMRSSMTree_uu_suLsuR;
      splitting_kernel1 = &Process::Pqg;
      
      sigmaPartTree2 = &Process::matrix_xsec_stub;
      splitting_kernel2 = &Process::Pgq;
      
      matrixelementReal_SC = &Process::matrix_soft_stub;
      
      m1 = MassSuL;
      m2 = MassSuL;
      // if gluino cannot be on-shell use full real ME
      if( pt.get<double>("collider setup.sqrt_S") < std::min(m1,m2) +  MassGlu ||
          MassGlu < m2 ) {
         std::cout << "\nINFO: Not using any subtraction for the gu_suLsuR channel.\n";
         matrixelementReal_HnonC = &Process::matrixMRSSMHard_gu_suLsuRubar;
         matrixelementReal_HnonC_CSub1 = nullptr;
         matrixelementReal_HnonC_CSub2 = nullptr;
      }
      else {
         // otherwise, if WidthGlu < 0 use DR
         if( WidthGlu < 0) {
            std::cout << "\nINFO: Using diagram removal for the gu_suLsuR channel.\n";
            matrixelementReal_HnonC = &Process::matrixMRSSMHard_gu_suLsuRubar_DR_wEtaDep;
            matrixelementReal_HnonC_CSub1 = nullptr;
            matrixelementReal_HnonC_CSub2 = nullptr;
         }
         // else use DS
         else {
            std::cout << "\nINFO: Using diagram subtraction for the gu_suLsuR channel with WoM = "
               << pt.get<double>("technical parameters.WidthOverMass") << ".\n";
            matrixelementReal_HnonC = &Process::matrixMRSSMHard_gu_suLsuRubar_DS;
            matrixelementReal_HnonC_CSub1 = &Process::matrixMRSSMHard_gu_suLsuRubar_DS_CSub1;
            matrixelementReal_HnonC_CSub2 = &Process::matrixMRSSMHard_gu_suLsuRubar_DS_CSub2;
         }
      }
      std::array<int, 1> partons = { 2 };
      for( int el : partons) flav.push_back( std::vector<int> {21, el, 2} );
   } 
   else if(processID == "MSSM,uubar_suLsuLdagger") {
      matrixelementTree = &Process::matrixMSSMTree_uu_suLsuR;  // matrix elements are identical
      //std::vector<int> row {2, -2};
      flav.push_back( std::vector<int> {2, -2, 2} );
      m1 = MassSuL;
      m2 = MassSuL;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MSSM,uubar_suLsuRdagger") {
      matrixelementTree = &Process::matrixMSSMTree_uubar_suLsuRdagger;
      flav.push_back( std::vector<int> {2, -2, 2} );
      m1 = MassSuL;
      m2 = MassSuR;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MSSM,udbar_suLsdLdagger") {
      matrixelementTree = &Process::matrixMSSMTree_ud_suLsdR;  // matrix elements are identical
      m1 = MassSuL;
      m2 = MassSuR;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MSSM,udbar_suLsdRdagger") {
      matrixelementTree = &Process::matrixMSSMTree_uubar_suLsuRdagger; // matrix elements are identical
      m1 = MassSuL;
      m2 = MassSuR;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MSSM,ddbar_suLsuLdagger") {
      matrixelementTree = &Process::matrixMRSSMTree_ddbar_suLsuLdagger; // matrix elements are identical
      m1 = MassSuL;
      m2 = MassSuR;
      k = 2.*2*3*3;
      h = 2.*2;
   }
/*      else if(processID == "MSSM,GG_suLsuLdagger") {
      matrixelementTree = &Process::matrixMRSSMTree_GG_suLsuLdagger; // matrix elements are identical
      m1 = MassSuL;
      m2 = MassSuL;
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
      matrixelementReal_SC = &Process::matrixSimplifiedSoft_uubar_OOg;
      sigmaPartTree1 = &Process::sigmaMRSSMTree_uubar_OO;
      splitting_kernel1 = &Process::Pgg;
      sigmaPartTree2 = &Process::sigmaMRSSMTree_uubar_OO;
      splitting_kernel2 = &Process::Pgg;      
      // matrixelementReal_HnonC = &Process::matrixSimplifiedHard_uubar_OOg;
      m1 = MasssigmaO;
      m2 = MasssigmaO;
      flav.push_back( std::vector<int> {1, -1, 2} );
      flav.push_back( std::vector<int> {2, -2, 2} );      
      flav.push_back( std::vector<int> {3, -3, 2} );
      flav.push_back( std::vector<int> {4, -4, 2} );
      flav.push_back( std::vector<int> {5, -5, 2} );
   }
   else if(processID == "Simplified,gg_OO") {
      sigmaPartTree1 = &Process::matrix_xsec_stub;
      matrixelementTree = &Process::matrix_tree_stub;
      matrixelementVirt = &Process::matrix_virt_stub;
      matrixelementReal_SC = &Process::matrixSimplifiedSoft_gg_OOg;
      // matrixelementReal_HnonC = &Process::matrixSimplifiedHard_gg_OOg;
      m1 = MasssigmaO;
      m2 = MasssigmaO;
      flav.push_back( std::vector<int> {21, 21, 0} );
   }
   else {
      std::cout << "Error! Subprocess " << processID << " not implemented.\n";
   }
}

/* --------------------------------------------------------------------------------------------*/
/* ------------------------------------------ Tree --------------------------------------------*/
/* --------------------------------------------------------------------------------------------*/


/* -------------------- Squark-squark production: q+q > sq+sq ---------------------------------*/
#include "matrix_elements_and_xsections/mssm_uu_suLsuR_tree_matrix.cpp"
#include "matrix_elements_and_xsections/mssm_uu_suLsuL_tree_matrix.cpp"
#include "matrix_elements_and_xsections/mssm_ud_suLsdR_tree_matrix.cpp"
#include "matrix_elements_and_xsections/mssm_ud_suLsdL_tree_matrix.cpp"                

/* -------------------- Squark-antisquark production: q+q^bar > sq+sq^dagger ---------------------*/
#include "matrix_elements_and_xsections/mrssm_uubar_suLsuLdagger_tree_matrix.cpp"
#include "matrix_elements_and_xsections/mrssm_ddbar_suLsuLdagger_tree_matrix.cpp"
#include "matrix_elements_and_xsections/mrssm_udbar_suLsdLdagger_tree_matrix.cpp"
#include "matrix_elements_and_xsections/mssm_uubar_suLsuRdagger_tree_matrix.cpp"

/* -------------------- Squark-antisquark production: G+G > sq+sq^dagger ---------------------*/
#include "matrix_elements_and_xsections/mrssm_gg_suLsuLdagger_tree_matrix.cpp"



/* /////////////////////////////////// partonic xsections /////////////////////////////////////// */
//#include "matrix_elements_and_xsections/mrssm_uu_suLsuR_tree_xsec.cpp"
#include "matrix_elements_and_xsections/mrssm_partonic_xsections.cpp"

/* -----------------------------------------------------------------------------------------------*/
/* ------------------------------------------ Virtual --------------------------------------------*/
/* -----------------------------------------------------------------------------------------------*/

/* ----------------------- Squark-squark production: q+q > sq+sq ---------------------------------*/
// mssm
#include "matrix_elements_and_xsections/mssm_uu_suLsuR_virt_matrix.cpp"
#include "matrix_elements_and_xsections/mssm_ud_suLsdR_virt_matrix.cpp"

//mrssm
#include "matrix_elements_and_xsections/mrssm_uu_suLsuR_virt_matrix.cpp"
#include "matrix_elements_and_xsections/mrssm_ud_suLsdR_virt_matrix.cpp"

/* --------------------- Squark-antisquark production: G+G > sq+sq^dagger ------------------------*/
//mrssm
#include "matrix_elements_and_xsections/mrssm_gg_suLsuLdagger_virt_matrix.cpp"

/* --------------------- Squark-antisquark production: q+qbar > sq+sq^dagger ------------------------*/
//mrssm
#include "matrix_elements_and_xsections/mrssm_ddbar_suLsuLdagger_virt_matrix.cpp"
#include "matrix_elements_and_xsections/mrssm_uubar_suLsuLdagger_virt_matrix.cpp"

/*
 
   Real emissions
 
 */
#include "matrix_elements_and_xsections/mrssm_uubar_suLsuLdaggerg_soft.cpp"
#include "matrix_elements_and_xsections/mrssm_uubar_suLsuLdaggerg_hard.cpp"

#include "matrix_elements_and_xsections/mrssm_ddbar_suLsuLdaggerg_soft.cpp"
#include "matrix_elements_and_xsections/mrssm_ddbar_suLsuLdaggerg_hard.cpp"

#include "matrix_elements_and_xsections/mrssm_gd_suLsuLdaggerd_hard.cpp"
#include "matrix_elements_and_xsections/mrssm_gu_suLsuLdaggeru_hard-DS.cpp"
#include "matrix_elements_and_xsections/mrssm_gu_suLsuLdaggeru_hard-DS_CSub2.cpp"
#include "matrix_elements_and_xsections/mrssm_gu_suLsuLdaggeru_hard-DR2.cpp"
#include "matrix_elements_and_xsections/mrssm_gu_suLsuLdaggeru_hard-DR_wEta_noSimplify.cpp"
#include "matrix_elements_and_xsections/mrssm_gu_suLsuLdaggeru_hard.cpp"

#include "matrix_elements_and_xsections/mrssm_gubar_suLsuLdaggerubar_hard-DS.cpp"
#include "matrix_elements_and_xsections/mrssm_gubar_suLsuLdaggerubar_hard-DS_CSub1.cpp"

#include "matrix_elements_and_xsections/mrssm_gg_suLsuLdaggerg_hard.cpp"
#include "matrix_elements_and_xsections/mrssm_gg_suLsuLdaggerg_soft.cpp"

#include "matrix_elements_and_xsections/mrssm_uu_suLsuRg_soft.cpp"
#include "matrix_elements_and_xsections/mrssm_uu_suLsuRg_hard.cpp"
#include "matrix_elements_and_xsections/mrssm_gu_suLsuRubar_hard.cpp"
#include "matrix_elements_and_xsections/mrssm_gu_suLsuRubar_hard-DR.cpp"
#include "matrix_elements_and_xsections/mrssm_gu_suLsuRubar_hard-DR_wEtaDep.cpp"
#include "matrix_elements_and_xsections/mrssm_gu_suLsuRubar_hard-DS.cpp"
#include "matrix_elements_and_xsections/mrssm_gu_suLsuRubar_hard-DS_CSub1.cpp"
#include "matrix_elements_and_xsections/mrssm_gu_suLsuRubar_hard-DS_CSub2.cpp"

#include "matrix_elements_and_xsections/simplified_uubar_OOg_soft.cpp"
#include "matrix_elements_and_xsections/simplified_uubar_OOg_hard.cpp"

#include "matrix_elements_and_xsections/simplified_gg_OOg_soft.cpp"
#include "matrix_elements_and_xsections/simplified_gg_OOg_hard.cpp"
