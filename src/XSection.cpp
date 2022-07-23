#include "XSection.hpp"

#include <string>

void XSection::init(Process&& process_, boost::property_tree::ptree const& pt_in, boost::program_options::variables_map const& vm_in) {

   pt = pt_in;
   vm = vm_in;

   processID = process_;
   m1 = processID.m1;
   m2 = processID.m2;

   S_sqrt = pt.get<double>("collider setup.sqrt_S");
   S = Sqr(S_sqrt);
   mu_r = pt.get<double>("collider setup.mu_r");
   mu_f = pt.get<double>("collider setup.mu_f");
   pdf = LHAPDF::mkPDF(pt.get<std::string>("collider setup.pdf"), 0);
   LHAPDF::setVerbosity(0);

   //technical
   dS = pt.get<double>("technical parameters.dS");
   dC = pt.get<double>("technical parameters.dC");
}
