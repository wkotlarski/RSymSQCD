#include <string>

#include "XSection.hpp"

void XSection::init(Process *processID_init, boost::property_tree::ptree pt_in, double x, double y, double z) {
    
   pt = pt_in;
   prec_virt = x;
   prec_sc = y;
   prec_hnc = z;

   processID = processID_init;
   m1 = processID->m1;
   m2 = processID->m2;
    
   S_sqrt = pt.get<double>("collider setup.sqrt_S");
   S = pow(S_sqrt, 2);
   mu_r = pt.get<double>("collider setup.mu_r");
   mu_f = pt.get<double>("collider setup.mu_f");
   pdf = LHAPDF::mkPDF( pt.get<std::string>("collider setup.pdf") , 0);
   LHAPDF::setVerbosity(0);
    
   //technical
   dS = pt.get<double>("technical parameters.dS");
   dC = pt.get<double>("technical parameters.dC");
}
