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
		matrixelementTree = &Process::matrixMSSMTree_uu_suLsuR; // same as in MSSM 
      matrixelementVirt = &Process::matrixMRSSMVirt_uu_suLsuR;
      sigmaPartTree = &Process::sigmaMRSSMTree_uu_suLsuR; // same as in MSSM
      matrixelementReal_SC = &Process::matrixMRSSMSoft_uu_suLsuRg;
      matrixelementReal_HnonC = &Process::matrixMRSSMHard_uu_suLsuRg;
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
      matrixelementReal_SC = &Process::matrix_soft_stub;
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
      matrixelementReal_SC = &Process::matrix_soft_stub;
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
      sigmaPartTree = &Process::sigmaMRSSMTree_uubar_suLsuLdagger;
      matrixelementTree = &Process::matrixMRSSMTree_uubar_suLsuLdagger;
      matrixelementVirt = &Process::matrixMRSSMVirt_uubar_suLsuLdagger;
      matrixelementReal_SC = &Process::matrixMRSSMSoft_uubar_suLsuLdaggerg;
      matrixelementReal_HnonC = &Process::matrixMRSSMHard_uubar_suLsuLdaggerg;
      m1 = MassSuL;
      m2 = MassSuL;
      flav.push_back( std::vector<int> {2, -2, 2} );
      f1 = 2;
      f2 = -2;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MRSSM,ddbar_suLsuLdagger") {
      sigmaPartTree = &Process::sigmaMRSSMTree_ddbar_suLsuLdagger;
      matrixelementTree = &Process::matrixMRSSMTree_ddbar_suLsuLdagger;
      matrixelementVirt = &Process::matrixMRSSMVirt_ddbar_suLsuLdagger;
      matrixelementReal_SC = &Process::matrixMRSSMSoft_ddbar_suLsuLdaggerg;
      matrixelementReal_HnonC = &Process::matrixMRSSMHard_ddbar_suLsuLdaggerg;
      m1 = MassSuL;
      m2 = MassSuL;
      flav.push_back( std::vector<int> {1, -1, 2} );
      flav.push_back( std::vector<int> {3, -3, 2} );
      flav.push_back( std::vector<int> {4, -4, 2} );
      flav.push_back( std::vector<int> {5, -5, 2} );
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
      sigmaPartTree = &Process::sigmaMRSSMTree_gg_suLsuLdagger;
      matrixelementTree = &Process::matrixMRSSMTree_GG_suLsuLdagger; 
      matrixelementVirt = &Process::matrixMRSSMVirt_GG_suLsuLdagger;
      matrixelementReal_SC = &Process::matrixMRSSMSoft_gg_suLsuLdaggerg;
      matrixelementReal_HnonC = &Process::matrixMRSSMHard_gg_suLsuLdaggerg;
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
      matrixelementVirt = &Process::matrix_virt_stub;
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
      matrixelementVirt = &Process::matrix_virt_stub;
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

inline double Process::det( std::vector< double* >& p, int i1, int i2, int i3, int i4 ) {
   boost::numeric::ublas::matrix<double> m (4, 4);
   for (unsigned i = 0; i < m.size1 (); ++ i) {
      m(0, i) = p[i1-1][i];
      m(1, i) = p[i2-1][i];
      m(2, i) = p[i3-1][i];
      m(3, i) = p[i4-1][i];
   }
        
   return 0*determinant(m);
}

int Process::determinant_sign(const boost::numeric::ublas::permutation_matrix<std::size_t>& pm)
{
    int pm_sign=1;
    std::size_t size = pm.size();
    for (std::size_t i = 0; i < size; ++i)
        if (i != pm(i))
            pm_sign *= -1.0; // swap_rows would swap a pair of rows here, so we change sign
    return pm_sign;
}

double Process::determinant( boost::numeric::ublas::matrix<double>& m ) {
    boost::numeric::ublas::permutation_matrix<std::size_t> pm(m.size1());
    double det = 1.0;
    if( boost::numeric::ublas::lu_factorize(m,pm) ) {
        det = 0.0;
    } else {
        for(int i = 0; i < m.size1(); i++) 
            det *= m(i,i); // multiply by elements on diagonal
        det = det * determinant_sign( pm );
    }
    return det;
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
//#include "matrix_elements_and_xsections/mrssm_uu_suLsuR_tree_xsec.cpp"
#include "matrix_elements_and_xsections/mrssm_partonic_xsections.cpp"
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
#include "matrix_elements_and_xsections/mrssm_uubar_suLsuLdaggerg_soft.cpp"
#include "matrix_elements_and_xsections/mrssm_uubar_suLsuLdaggerg_hard.cpp"
#include "matrix_elements_and_xsections/mrssm_gg_suLsuLdaggerg_hard.cpp"
#include "matrix_elements_and_xsections/mrssm_ddbar_suLsuLdaggerg_hard.cpp"

#include "matrix_elements_and_xsections/mrssm_uu_suLsuRg_soft.cpp"
#include "matrix_elements_and_xsections/mrssm_uu_suLsuRg_hard.cpp"