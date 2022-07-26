#include "XSection.hpp"

#include "LHAPDF/Info.h"
#include "LHAPDF/Config.h"

#include <string>

void XSection::init(boost::property_tree::ptree const& pt_in, boost::program_options::variables_map const& vm_in) {

   pt = pt_in;
   vm = vm_in;

   sqrtS_= pt.get<double>("collider setup.sqrt_S");
}
