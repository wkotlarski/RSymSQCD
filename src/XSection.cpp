
#include "XSection.h"
#include <iostream>
#include <vector>
#include <string>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

void XSection::init(void) {
    
    squark_mass = {{
          {2000, 2000},
          {2000, 2000},
          {1500, 1500},
          {1500, 1500},
          {1500, 1500},
          {1500, 1500}
    }};
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
                if ( tokens.at(0) == "PDF_NLO") {
                    pdf_nlo = LHAPDF::mkPDF(tokens.at(1), 0);
                }
            read_card.push_back(tokens);
        }
    }
    std::cout << pdf_nlo->alphasQ(2000.) << '\n';
    for (auto const &element: read_card)
        //std::cout << element[1] << '\n';
        std::cout << element.at(0) <<'\n';
}
