#include "SM.h"
#include "MSSM.h"
#include "MRSSM.h"
#include "../../include/IMatrixElements.h"

IMatrixElements* IMatrixElements::create_process(boost::property_tree::ptree const& ptree) {
   if (ptree.get<std::string>("process.model") == "SM") {
      return new SM(ptree);
   }
   else if (ptree.get<std::string>("process.model") == "MSSM") {
      return new MSSM(ptree);
   }
   else if (ptree.get<std::string>("process.model") == "MRSSM") {
      return new MRSSM(ptree);
   }
}
