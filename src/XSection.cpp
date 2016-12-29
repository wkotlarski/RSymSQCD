
#include "XSection.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

void XSection::init(Process *processID_init, boost::property_tree::ptree pt_in, double x, double y, double z) {
    
   pt = pt_in;
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

    processID = processID_init;
    m1 = processID->m1;
    m2 = processID->m2;
    
    S_sqrt = pt.get<double>("collider setup.sqrt_S");
    S = pow(S_sqrt, 2);
    mu_r = pt.get<double>("collider setup.mu_r");
    mu_f = pt.get<double>("collider setup.mu_f");
    //mu_r = muR;
    //mu_f = muF;
    pdf = LHAPDF::mkPDF( pt.get<std::string>("collider setup.pdf") , 0);
    LHAPDF::setVerbosity(0);
    
    //technical
    dS = pt.get<double>("technical parameters.dS");
    std::cout << dS << '\n';
    dC = pt.get<double>("technical parameters.dC");
}
