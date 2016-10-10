#include "Process.h"
// squark production
#include "MRSSM_Tree_uu_su1su4.h"
#include "MSSM_Tree_uu_su1su4.h"
#include "MRSSM_Tree_ud_su1sd4.h"
#include "MSSM_Tree_ud_su1sd4.h"
#include "MRSSM_1L_uu_su1su4.h"
#include "MSSM_1L_uu_su1su4.h"
#include "MRSSM_1L_ud_su1sd4.h"
#include "MSSM_1L_ud_su1sd4.h"
// squark anti-squark production
#include "MRSSM_Tree_uubar_su1su1dagger.h"
#include "MRSSM_Tree_GG_su1su1dagger.h"
#include "MRSSM_1L_uubar_su1su1dagger.h"
#include "MRSSM_1L_GG_su1su1dagger.h"

#include <string>

Process::Process(std::string processID) {	        
    /* squark production, MRSSM */
	if(processID == "MRSSM,uu_suLsuR") {
        matrixelementTree = &matrixMRSSMTree_uu_suLsuR;
        matrixelementVirt = &matrixMRSSMVirt_uu_suLsuR;
        f1 = 2.;
        f2 = 2.;
        k = 2.*2*3*3;
        h = 2.*2;
    }
    else if(processID == "MRSSM,ud_suLsdR")
    {
        matrixelementTree = &matrixMRSSMTree_ud_suLsdR;
        matrixelementVirt = &matrixMRSSMVirt_ud_suLsdR;
        f1 = 2.;
        f2 = 1.;
        k = 2.*2*3*3;
        h = 2.*2;
    }
    /* squark production, MSSM */
    else if(processID == "MSSM,uu_suLsuR")
    {
        matrixelementTree = &matrixMSSMTree_uu_suLsuR;
        matrixelementVirt = &matrixMSSMVirt_uu_suLsuR;
        f1 = 2.;
        f2 = 2.;
        k = 2.*2*3*3;
        h = 2.*2;
    }
    else if(processID == "MSSM,ud_suLsdR")
    {
        matrixelementTree = &matrixMSSMTree_ud_suLsdR;
        matrixelementVirt = &matrixMSSMVirt_ud_suLsdR;
        f1 = 2.;
        f2 = 1.;
        k = 2.*2*3*3;
        h = 2.*2;
    }
    /* squark anti-squark production, MRSSM */
    else if(processID == "MRSSM,uubar_suLsuLdagger")
    {
        matrixelementTree = &matrixMRSSMTree_uubar_suLsuLdagger;
        matrixelementVirt = &matrixMRSSMVirt_uubar_suLsuLdagger;
        f1 = 2.;
        f2 = -2.;
        k = 2.*2*3*3;
        h = 2.*2;
    }
    else if(processID == "MRSSM,GG_suLsuLdagger")
    {
        matrixelementTree = &matrixMRSSMTree_GG_suLsuLdagger;
        matrixelementVirt = &matrixMRSSMVirt_GG_suLsuLdagger;
        f1 = 0.;
        f2 = 0.;
        k = 2.*2*8*8;
        h = 1.;
    }
}
