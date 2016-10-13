
#include "XSection.h"
#include <iostream>
#include <vector>
#include <string>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
 
//XSection::XSection(Process process_init) {
//    process = process_init;
//}


void XSection::init(Process *processID_init, boost::property_tree::ptree pt, double x, double y, double z) {
    
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
    muR = pt.get<double>("collider setup.mu_r");
    muF = pt.get<double>("collider setup.mu_f");
    pdf_nlo = LHAPDF::mkPDF( pt.get<std::string>("collider setup.pdf_nlo") , 0);
    LHAPDF::setVerbosity(0);
}
