#include "Process.h"
#include "MRSSM_Tree_uu_su1su4.h"
#include "MSSM_Tree_uu_su1su4.h"
#include "MRSSM_Tree_ud_su1sd4.h"
#include "MSSM_Tree_ud_su1sd4.h"
#include "MRSSM_1L_uu_su1su4.h"
#include "MSSM_1L_uu_su1su4.h"
#include "MRSSM_1L_ud_su1sd4.h"
#include "MSSM_1L_ud_su1sd4.h"
#include <string>

Process::Process(std::string processID) {	        
	if(processID == "MRSSM,uu_suLsuR") {
        matrixelementTree = &matrixMRSSMTree_uu_suLsuR;
        matrixelementVirt = &matrixMRSSMVirt_uu_suLsuR;
//        f1 = 2.;
//        f2 = 2.;
    }
    else if(processID == "MSSM,uu_suLsuR")
    {
        matrixelementTree = &matrixMSSMTree_uu_suLsuR;
        matrixelementVirt = &matrixMSSMVirt_uu_suLsuR;
//        f1 = 2.;
//        f2 = 2.;
    }
    else if(processID == "MRSSM,ud_suLsdR")
    {
        matrixelementTree = &matrixMRSSMTree_ud_suLsdR;
        matrixelementVirt = &matrixMRSSMVirt_ud_suLsdR;
//        f1 = 2.;
//        f2 = 1.;
    }
    else if(processID == "MSSM,ud_suLsdR")
    {
        matrixelementTree = &matrixMSSMTree_ud_suLsdR;
        matrixelementVirt = &matrixMSSMVirt_ud_suLsdR;
//        f1 = 2.;
//        f2 = 1.;
    }
}
