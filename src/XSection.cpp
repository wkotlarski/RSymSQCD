
#include "XSection.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
 
//XSection::XSection(Process process_init) {
//    process = process_init;
//}


void XSection::init(Process *processID_init, boost::property_tree::ptree pt, double x, double y, double z) {
    
   ptr = pt;
   prec_virt = x;
   prec_sc = y;
   prec_hnc = z;
   double m_sq = pt.get<double>("masses.suL");
   squark_mass = {{
          {m_sq, m_sq},
          {m_sq, m_sq},
          {m_sq, m_sq},
          {m_sq, m_sq},
          {m_sq, m_sq},
          {m_sq, m_sq}
    }};
    gluino_mass = pt.get<double>("masses.gluino");
    processID = processID_init;
    
    S_sqrt = ptr.get<double>("collider setup.sqrt_S");
    S = pow(S_sqrt, 2);
    muR = pt.get<double>("collider setup.mu_r");
    muF = pt.get<double>("collider setup.mu_f");
    pdf_nlo = LHAPDF::mkPDF( pt.get<std::string>("collider setup.pdf_nlo") , 0);
    LHAPDF::setVerbosity(0);
    
    //technical
    dS = ptr.get<double>("technical parameters.dS");
    dC = ptr.get<double>("technical parameters.dC");
}
