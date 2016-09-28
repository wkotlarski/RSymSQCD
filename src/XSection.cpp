
#include "XSection.h"
#include <iostream>
#include <vector>
#include <string>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
 
//XSection::XSection(Process process_init) {
//    process = process_init;
//}


void XSection::init(double m_sq, double m_gluino, double m_sgluon, Process *processID_init) {
    
    squark_mass = {{
          {m_sq, m_sq},
          {m_sq, m_sq},
          {m_sq, m_sq},
          {m_sq, m_sq},
          {m_sq, m_sq},
          {m_sq, m_sq}
    }};
    gluino_mass = m_gluino;
    sgluon_mass = m_sgluon;
    processID = processID_init;
    muR = squark_mass.at(0).at(0);
    muF = squark_mass.at(0).at(0);
    //XSection::squark_mass = 
    std::vector <std::vector<std::string>> read_card;
    read_card.clear();
    std::string temp;
    std::ifstream file("run_parameters.dat");
    if ( file.is_open() )  {
        while ( getline(file, temp) ) {
            
            std::vector<std::string> tokens;
            split(tokens, temp, boost::algorithm::is_any_of(" ,\t"));
            
            tokens.erase(
                std::remove_if(
                    tokens.begin(), 
                    tokens.end(), 
                    [](std::string x){
                        return std::all_of(x.begin(),x.end(),[](char c){ 
                  return std::isspace(static_cast<unsigned char>(c));
               });
                    }
                ), 
                tokens.end()
            );
                if ( tokens.at(0) == "pdf_nlo") {
                    pdf_nlo = LHAPDF::mkPDF(tokens.at(1), 0);
                    LHAPDF::setVerbosity(0);
                }
            read_card.push_back(tokens);
        }
    }
    std::cout << std::scientific << std::setprecision(18) << 
            pdf_nlo->alphasQ(1000.) << " " << pdf_nlo->xfxQ(0, 0.1, 1000.) << '\n';
    //for (auto const &element: read_card)
        //std::cout << element[1] << '\n';
        //std::cout << element.at(0) <<'\n';
}
