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
   pdf = LHAPDF::mkPDF( pt.get<std::string>("collider setup.pdf") , 0);
   
   // @todo remove 
   MassSq = MassSuL;
   partonic = false; // calculate sigma_had with |M|^2 and not sigma_part
   
/* -------------------- Squark-squark production---------------------------------*/

	if(processID == "MRSSM,uu_suLsuR") {
      sigmaPartTree = &Process::sigmaMSSMTree_uu_suLsuR; // same as in MSSM
		matrixelementTree = &Process::matrixMSSMTree_uu_suLsuR; // same as in MSSM 
	    
      matrixelementVirt = &Process::matrixMRSSMVirt_uu_suLsuR;
      matrixelementReal_SC = &Process::matrixMRSSMSoft_uu_suLsuRg;
      m1 = MassSuL;
      m2 = MassSuR;
      flav.push_back( std::vector<int> {2, 2, 1} );
      f1 = 2;
      f2 = 2;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MRSSM,ud_suLsdR") { // same as in MSSM
      matrixelementTree = &Process::matrixMSSMTree_ud_suLsdR;
      matrixelementVirt = &Process::matrixMRSSMVirt_ud_suLsdR;
      m1 = MassSuL;
      m2 = MassSuR;
      f1 = 2.;
      f2 = 1.;
      k = 2.*2*3*3;
      h = 2.*2;
   }
    /* squark production, MSSM */
   else if(processID == "MSSM,uu_suLsuR") {
      matrixelementTree = &Process::matrixMSSMTree_uu_suLsuR;
      matrixelementVirt = &Process::matrixMSSMVirt_uu_suLsuR;
      //matrixelementReal_SC = &Process::matrixMRSSMSoft_uu_suLsuRg;
      matrixelementReal_SC = &Process::g;
      flav.push_back( std::vector<int> {2, 2, 1} );
      m1 = MassSuL;
      m2 = MassSuR;
      f1 = 2.;
      f2 = 2.;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MSSM,uu_suLsuL") {
      matrixelementTree = &Process::matrixMSSMTree_uu_suLsuL;
      //matrixelementVirt = &Process::matrixMSSMVirt_uu_suLsuL;
      //matrixelementReal_SC = &Process::matrixMRSSMSoft_uu_suLsuLg;
      matrixelementReal_SC = &Process::g;
      m1 = MassSuL;
      m2 = MassSuR;
      f1 = 2.;
      f2 = 2.;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MSSM,ud_suLsdR") {
      matrixelementTree = &Process::matrixMSSMTree_ud_suLsdR;
      matrixelementVirt = &Process::matrixMSSMVirt_ud_suLsdR;
      m1 = MassSuL;
      m2 = MassSuR;
      f1 = 2.;
      f2 = 1.;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MSSM,ud_suLsdL") {
      matrixelementTree = &Process::matrixMSSMTree_ud_suLsdL;
      //matrixelementVirt = &matrixMSSMVirt_ud_suLsdL;
      m1 = MassSuL;
      m2 = MassSuR;
      f1 = 2.;
      f2 = 1.;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   
/* -------------------- Squark-antisquark production ---------------------*/

   else if(processID == "MRSSM,uubar_suLsuLdagger") {
      matrixelementTree = &Process::matrixMRSSMTree_uubar_suLsuLdagger;
      //matrixelementVirt = &matrixMRSSMVirt_uubar_suLsuLdagger;
      matrixelementReal_SC = &Process::matrixMRSSMSoft_uubar_suLsuLdaggerg;
      m1 = MassSuL;
      m2 = MassSuL;
      flav.push_back( std::vector<int> {2, -2, 2} );
      f1 = 2;
      f2 = -2;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MRSSM,ddbar_suLsuLdagger") {
      matrixelementTree = &Process::matrixMRSSMTree_ddbar_suLsuLdagger;
      //matrixelementVirt = &matrixMRSSMVirt_ddbar_suLsuLdagger;
      matrixelementReal_SC = &Process::matrixMRSSMSoft_ddbar_suLsuLdaggerg;
      m1 = MassSuL;
      m2 = MassSuL;
      flav.push_back( std::vector<int> {1, -1, 2} );
      f1 = 1;
      f2 = -1;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MRSSM,udbar_suLsdLdagger") {
      matrixelementTree = &Process::matrixMSSMTree_ud_suLsdR;  // matrix elements are identical
      m1 = MassSuL;
      m2 = MassSuR;
      f1 = 2;
      f2 = -1;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MRSSM,GG_suLsuLdagger") {
      matrixelementTree = &Process::matrixMRSSMTree_GG_suLsuLdagger; 
      //matrixelementVirt = &Process::matrixMRSSMVirt_uu_suLsuR;
      matrixelementReal_SC = &Process::matrixMRSSMSoft_gg_suLsuLdaggerg;
      m1 = MassSuL;
      m2 = MassSuL;
      flav.push_back( std::vector<int> {0, 0, 1} );
      f1 = 0;
      f2 = 0;
      k = 2.*2*8*8;
      h = 1.;
   }
   else if(processID == "MSSM,uubar_suLsuLdagger") {
      matrixelementTree = &Process::matrixMSSMTree_uu_suLsuR;  // matrix elements are identical
      //std::vector<int> row {2, -2};
      flav.push_back( std::vector<int> {2, -2} );
      m1 = MassSuL;
      m2 = MassSuL;
      f1 = 2;
      f2 = -2;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MSSM,uubar_suLsuRdagger") {
      matrixelementTree = &Process::matrixMSSMTree_uubar_suLsuRdagger;
      m1 = MassSuL;
      m2 = MassSuR;
      f1 = 2;
      f2 = -2;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MSSM,udbar_suLsdLdagger") {
      matrixelementTree = &Process::matrixMSSMTree_ud_suLsdR;  // matrix elements are identical
      m1 = MassSuL;
      m2 = MassSuR;
      f1 = 2;
      f2 = -1;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MSSM,udbar_suLsdRdagger") {
      matrixelementTree = &Process::matrixMSSMTree_uubar_suLsuRdagger; // matrix elements are identical
      m1 = MassSuL;
      m2 = MassSuR;
      f1 = 2;
      f2 = -1;
      k = 2.*2*3*3;
      h = 2.*2;
   }
      else if(processID == "MSSM,ddbar_suLsuLdagger") {
      matrixelementTree = &Process::matrixMRSSMTree_ddbar_suLsuLdagger; // matrix elements are identical
      m1 = MassSuL;
      m2 = MassSuR;
      f1 = 1;
      f2 = -1;
      k = 2.*2*3*3;
      h = 2.*2;
   }
      else if(processID == "MSSM,GG_suLsuLdagger") {
      matrixelementTree = &Process::matrixMRSSMTree_GG_suLsuLdagger; // matrix elements are identical
      m1 = MassSuL;
      m2 = MassSuL;
      f1 = 0;
      f2 = 0;
      k = 2.*2*8*8;
      h = 1.;
   }

/* -------------------- Sgluon production ---------------------*/

   else if( processID == "sgluons-qqbar_OO" ) {	  
	  partonic = true;
      sigmaPartTree = &Process::matrixSgluonTree_qqbar_OO;
      matrixelementVirt = &Process::f;
      m1 =  MasssigmaO;
      m2 =  MasssigmaO;
      k = 1;
      f1 = 69;
      f2 = 69;
      h=1;
   }
   else if( processID == "sgluons-gg_OO" ) {
	  partonic = true;
      sigmaPartTree = &Process::matrixSgluonTree_gg_OO;
      matrixelementVirt = &Process::f;
      m1 = m2 = MasssigmaO;
      k = 1;
      f1 = 0;
      f2 = 0;
      h=1;
   }
   else {
      std::cout << "Error! Subprocess " << processID << " not implemented.\n";
   }
}

inline double Process::f(double S, double T, double x, double y, int z) {
   return 0.;
}

inline double Process::g(double S, double T) {
   return 0.;
}
/* ///////////////////////////////////// matrix elements ///////////////////////////////////// */

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

/* -------------------- Squark-antisquark production: q+q^bar > sq+sq^dagger ---------------------*/
#include "matrix_elements_and_xsections/mrssm_gg_suLsuLdagger_tree_matrix.cpp"



/* /////////////////////////////////// partonic xsections /////////////////////////////////////// */
#include "matrix_elements_and_xsections/mrssm_uu_suLsuR_tree_xsec.cpp"

/* --------------------------------- Sgluon production -------------------------------------------*/
#include "matrix_elements_and_xsections/mrssm_uubar_OO_tree_xsec.cpp" // is identical for all initial state quarks
#include "matrix_elements_and_xsections/mrssm_gg_OO_tree_xsec.cpp"

/* ///////////////////////////////////// matrix elements ///////////////////////////////////// */


 


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


// real
#include "matrix_elements_and_xsections/mrssm_gg_suLsuLdaggerg.cpp"
#include "matrix_elements_and_xsections/mrssm_ddbar_suLsuLdaggerg.cpp"
#include "matrix_elements_and_xsections/mrssm_uubar_suLsuLdaggerg.cpp"
/* -----------------------------------------------------------------------------------------------*/
/* --------------------------------------------- Soft --------------------------------------------*/ 
/* -----------------------------------------------------------------------------------------------*/

double Process::matrixMRSSMSoft_uu_suLsuRg(double s12, double th) {
   
   double Alfas = pdf->alphasQ( mu_r );
   double Alfas2 = pow(Alfas, 2);
   double beta = sqrt(1. - 4. * pow(m1, 2)/s12);
   double MGl2 = pow(MassGlu, 2);
   
   std::complex<double> temp = ((-42.666666666666667*Alfas*Alfas2*beta*(-1. + beta*beta)*MGl2*s12*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) + 
   (42.666666666666667*Alfas*Alfas2*beta*MGl2*s12*(-1. + beta*Cos(th))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) + 
   (21.333333333333333*Alfas*Alfas2*beta*MGl2*s12*Power(-1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) + 
   (21.333333333333333*Alfas*Alfas2*beta*(s12*s12)*Power(-1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) + 
   (10.666666666666667*Alfas*Alfas2*beta*(s12*s12)*Power(-1. + beta*Cos(th),3)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) + 
   (21.333333333333333*Alfas*Alfas2*beta*MGl2*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) - 
   (10.666666666666667*Alfas*Alfas2*beta*(s12*s12)*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) - 
   (2.3703703703703704*Alfas*Alfas2*beta*MGl2*s12*(-1. + beta*Cos(th))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),3) - 
   (0.59259259259259259*Alfas*Alfas2*beta*(-1. + beta*beta)*(s12*s12)*(-1. + beta*Cos(th))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),3) + 
   (1.1851851851851852*Alfas*Alfas2*beta*(s12*s12)*Power(-1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),3) + 
   (0.59259259259259259*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.59259259259259259*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (1.1851851851851852*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.14814814814814815*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.2962962962962963*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),3)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.2962962962962963*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.14814814814814815*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.2962962962962963*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (42.666666666666667*Alfas*Alfas2*beta*(-1. + beta*beta)*MGl2*s12*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),4) + 
   (42.666666666666667*Alfas*Alfas2*beta*MGl2*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),4) + 
   (21.333333333333333*Alfas*Alfas2*beta*(s12*s12)*Power(1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),4) + 
   (10.666666666666667*Alfas*Alfas2*beta*(s12*s12)*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),4) - 
   (10.666666666666667*Alfas*Alfas2*beta*(s12*s12)*Power(1. + beta*Cos(th),3)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),4) + 
   (0.2962962962962963*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.074074074074074074*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.59259259259259259*Alfas*Alfas2*beta*s12*(1. + beta*Cos(th))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.2962962962962963*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.14814814814814815*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.59259259259259259*Alfas*Alfas2*beta*s12*Power(1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.074074074074074074*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.14814814814814815*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.14814814814814815*Alfas*Alfas2*beta*s12*Power(1. + beta*Cos(th),3)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.14814814814814815*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.037037037037037037*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*s12*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.2962962962962963*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.2962962962962963*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.037037037037037037*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.074074074074074074*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),3)*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.074074074074074074*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.037037037037037037*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.074074074074074074*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.14814814814814815*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.037037037037037037*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*s12*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.14814814814814815*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.037037037037037037*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.14814814814814815*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.074074074074074074*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.037037037037037037*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.074074074074074074*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (1.1851851851851852*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*(MGl2*MGl2)*s12*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (1.1851851851851852*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*MGl2*(s12*s12)*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.59259259259259259*Alfas*Alfas2*beta*Power(-1. + beta*beta,3)*MGl2*(s12*s12)*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.074074074074074074*Alfas*Alfas2*beta*Power(-1. + beta*beta,4)*Power(s12,3)*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.14814814814814815*Alfas*Alfas2*beta*Power(-1. + beta*beta,3)*Power(s12,3)*(-1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.14814814814814815*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*Power(s12,3)*Power(-1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.14814814814814815*Alfas*Alfas2*beta*Power(-1. + beta*beta,3)*Power(s12,3)*(1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.14814814814814815*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*Power(s12,3)*Power(1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.2962962962962963*Alfas*Alfas2*(-1. + beta*beta)*s12*(-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.2962962962962963*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)
     - (0.2962962962962963*Alfas*Alfas2*(-1. + beta*beta)*s12*(-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.2962962962962963*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)
     - (9.4814814814814815*Alfas*Alfas2*(-1. + beta*beta)*(MGl2*MGl2)*s12*(-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*
      Sin(th))/(Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (9.4814814814814815*Alfas*Alfas2*(-1. + beta*beta)*MGl2*(s12*s12)*(-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*
      Sin(th))/(Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (4.7407407407407407*Alfas*Alfas2*Power(-1. + beta*beta,2)*MGl2*(s12*s12)*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.59259259259259259*Alfas*Alfas2*Power(-1. + beta*beta,3)*Power(s12,3)*(-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*
      Sin(th))/(Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (1.1851851851851852*Alfas*Alfas2*Power(-1. + beta*beta,2)*Power(s12,3)*(-1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (1.1851851851851852*Alfas*Alfas2*(-1. + beta*beta)*Power(s12,3)*Power(-1. + beta*Cos(th),2)*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (1.1851851851851852*Alfas*Alfas2*Power(-1. + beta*beta,2)*Power(s12,3)*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (9.4814814814814815*Alfas*Alfas2*(MGl2*MGl2)*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (9.4814814814814815*Alfas*Alfas2*MGl2*(s12*s12)*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (4.7407407407407407*Alfas*Alfas2*(-1. + beta*beta)*MGl2*(s12*s12)*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.59259259259259259*Alfas*Alfas2*Power(-1. + beta*beta,2)*Power(s12,3)*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (1.1851851851851852*Alfas*Alfas2*(-1. + beta*beta)*Power(s12,3)*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (1.1851851851851852*Alfas*Alfas2*Power(s12,3)*Power(-1. + beta*Cos(th),3)*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (1.1851851851851852*Alfas*Alfas2*(-1. + beta*beta)*Power(s12,3)*Power(1. + beta*Cos(th),2)*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (1.1851851851851852*Alfas*Alfas2*(-1. + beta*beta)*Power(s12,3)*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (1.1851851851851852*Alfas*Alfas2*Power(s12,3)*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),3)*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.074074074074074074*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.037037037037037037*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*s12*(4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*
      Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.074074074074074074*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),2)*(4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*
      Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.037037037037037037*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),3)*(4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*
      Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.074074074074074074*Alfas*Alfas2*beta*s12*(1. + beta*Cos(th))*(4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.074074074074074074*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)\
    + (0.018518518518518519*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)\
    + (0.018518518518518519*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)\
    - (0.018518518518518519*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)\
    - (0.037037037037037037*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)\
    - (0.018518518518518519*Alfas*Alfas2*beta*s12*Power(1. + beta*Cos(th),3)*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)\
    + (0.037037037037037037*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*s12*(4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*
      Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.037037037037037037*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)\
    + (0.037037037037037037*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)\
    - (0.037037037037037037*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)\
    - (1.1851851851851852*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*(MGl2*MGl2)*s12*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (1.1851851851851852*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*MGl2*(s12*s12)*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.59259259259259259*Alfas*Alfas2*beta*Power(-1. + beta*beta,3)*MGl2*(s12*s12)*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.074074074074074074*Alfas*Alfas2*beta*Power(-1. + beta*beta,4)*Power(s12,3)*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.14814814814814815*Alfas*Alfas2*beta*Power(-1. + beta*beta,3)*Power(s12,3)*(-1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.14814814814814815*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*Power(s12,3)*Power(-1. + beta*Cos(th),2)*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (1.1851851851851852*Alfas*Alfas2*beta*(-1. + beta*beta)*(MGl2*MGl2)*s12*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (1.1851851851851852*Alfas*Alfas2*beta*(-1. + beta*beta)*MGl2*(s12*s12)*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.59259259259259259*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*MGl2*(s12*s12)*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.074074074074074074*Alfas*Alfas2*beta*Power(-1. + beta*beta,3)*Power(s12,3)*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.14814814814814815*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*Power(s12,3)*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.14814814814814815*Alfas*Alfas2*beta*(-1. + beta*beta)*Power(s12,3)*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.14814814814814815*Alfas*Alfas2*beta*(-1. + beta*beta)*Power(s12,3)*Power(1. + beta*Cos(th),3)*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.074074074074074074*Alfas*Alfas2*(-1. + beta*beta)*s12*(-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 
        2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.037037037037037037*Alfas*Alfas2*Power(-1. + beta*beta,2)*s12*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.074074074074074074*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*(-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 
        2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.037037037037037037*Alfas*Alfas2*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.018518518518518519*Alfas*Alfas2*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.018518518518518519*Alfas*Alfas2*s12*Power(-1. + beta*Cos(th),3)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.074074074074074074*Alfas*Alfas2*s12*(1. + beta*Cos(th))*(-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 
        2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.037037037037037037*Alfas*Alfas2*(-1. + beta*beta)*s12*(1. + beta*Cos(th))*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.018518518518518519*Alfas*Alfas2*s12*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.018518518518518519*Alfas*Alfas2*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.018518518518518519*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.018518518518518519*Alfas*Alfas2*s12*Power(1. + beta*Cos(th),3)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.074074074074074074*Alfas*Alfas2*(-1. + beta*beta)*s12*(-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 
        2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.037037037037037037*Alfas*Alfas2*Power(-1. + beta*beta,2)*s12*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.074074074074074074*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*(-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 
        2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.037037037037037037*Alfas*Alfas2*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.018518518518518519*Alfas*Alfas2*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.018518518518518519*Alfas*Alfas2*s12*Power(-1. + beta*Cos(th),3)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.074074074074074074*Alfas*Alfas2*s12*(1. + beta*Cos(th))*(-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 
        2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.037037037037037037*Alfas*Alfas2*(-1. + beta*beta)*s12*(1. + beta*Cos(th))*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.018518518518518519*Alfas*Alfas2*s12*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.018518518518518519*Alfas*Alfas2*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.018518518518518519*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.018518518518518519*Alfas*Alfas2*s12*Power(1. + beta*Cos(th),3)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.0023148148148148148*Alfas*Alfas2*beta*Power(1. - 1.*beta*Cos(th),2)*
      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) - 
        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.0011574074074074074*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) - 
        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. - 1.*beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) + 
   (0.0011574074074074074*Alfas*Alfas2*beta*Power(1. + beta*Cos(th),2)*
      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) - 
        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) - 
   (0.005787037037037037*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.0046296296296296296*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.0005787037037037037*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.0023148148148148148*Alfas*Alfas2*beta*(1. + beta*Cos(th))*(-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 
        0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 
        1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.0011574074074074074*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.0011574074074074074*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.0005787037037037037*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*Power(1. + beta*Cos(th),2)*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. - 1.*beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) + 
   (0.040509259259259259*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) - 
   (0.0040509259259259259*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) - 
   (0.0081018518518518519*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.0081018518518518519*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. - 1.*beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2)) - 
   (0.032407407407407407*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) - 
   (0.0040509259259259259*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*Power(1. + beta*Cos(th),2)*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. - 1.*beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2)) + 
   (0.0034722222222222222*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*
      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) - 
        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.0011574074074074074*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) - 
        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.0005787037037037037*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*Power(1. - 1.*beta*Cos(th),2)*
      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) - 
        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2)) + 
   (0.0005787037037037037*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) - 
        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.024305555555555556*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.016203703703703704*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) - 
   (0.016203703703703704*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) - 
   (0.032407407407407407*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
   (0.0081018518518518519*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
   (0.0040509259259259259*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*Power(1. - 1.*beta*Cos(th),2)*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
   (0.0081018518518518519*Alfas*Alfas2*beta*Power(1. - 1.*beta*Cos(th),3)*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
   (0.0040509259259259259*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) - 
   (0.0081018518518518519*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.0081018518518518519*Alfas*Alfas2*beta*Power(1. + beta*Cos(th),2)*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) - 
   (0.064814814814814815*Alfas*Alfas2*beta*(-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 
        0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 
        1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.016203703703703704*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.064814814814814815*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
   (0.016203703703703704*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2)) + 
   (0.032407407407407407*Alfas*Alfas2*beta*Power(1. - 1.*beta*Cos(th),2)*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
   (0.0081018518518518519*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.0081018518518518519*Alfas*Alfas2*beta*Power(1. + beta*Cos(th),2)*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.0011574074074074074*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2)) + 
   (0.0046296296296296296*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2));
   
   return temp.real();
}
