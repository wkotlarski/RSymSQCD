#include "XSection.hpp"

#include "LHAPDF/Info.h"
#include "LHAPDF/Config.h"

#include <string>

void XSection::init(Process const& process_, boost::property_tree::ptree const& pt_in, boost::program_options::variables_map const& vm_in) {

   pt = pt_in;
   vm = vm_in;

   processID = process_;

   S_sqrt = pt.get<double>("collider setup.sqrt_S");
   S = Sqr(S_sqrt);
   mu_r = pt.get<double>("collider setup.mu_r");
   mu_f = pt.get<double>("collider setup.mu_f");
   // disable printint LHAPDF header since it gets printed twice
   LHAPDF::Info& cfg = LHAPDF::getConfig();
   cfg.set_entry("Verbosity", 0);
   pdf = LHAPDF::mkPDF(pt.get<std::string>("collider setup.pdf"), 0);
   LHAPDF::setVerbosity(0);

   //technical
   dS = pt.get<double>("technical parameters.dS");
   dC = pt.get<double>("technical parameters.dC");
}
