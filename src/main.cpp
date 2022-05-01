#include <chrono>
using namespace std::chrono_literals;

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "XSection_Tree.hpp"
#include "XSection_Virt.hpp"
#include "XSection_SC.hpp"
#include "XSection_HnonC.hpp"

#include "LHAPDF/Info.h"
#include "LHAPDF/Config.h"

#include <nlohmann/json.hpp>

using namespace std;
using json = nlohmann::json;

/*
 *    XSection_* return results in the form xsection, its error, erros's pvalue
 *    overload + operator to add xsection and errors in quadrature
 *    pvalues are just added
 */
inline array<double,3> operator+(array<double,3> x, array<double,3> y) {
   return array<double,3> {
      x.at(0) + y.at(0), sqrt(pow(x.at(1),2)+pow(y.at(1),2)), x.at(2) + y.at(2)
   };
}


void xsec_to_json(json& j, string str, array<double,3> tree, array<double,3> virt, array<double,3> soft, array<double,3> hard) {
   j["cross sections"][str] = {
      {"tree", {{"res", tree.at(0)}, {"err", tree.at(1)}, {"p-val", tree.at(2)}}},
      {"virtual", {{"res", virt.at(0)}, {"err", virt.at(1)}, {"p-val", virt.at(2)}}},
      {"SC", {{"res", soft.at(0)}, {"err", soft.at(1)}, {"p-val", soft.at(2)}}},
      {"HnonC", {{"res", hard.at(0)}, {"err", hard.at(1)}, {"p-val", hard.at(2)}}}
   };
}

// why do I have to write this?
// why isn't init() enough?
Process *XSection::processID;
double XSection::mu_r;
double XSection::mu_f;
double XSection::m1;
double XSection::m2;
double XSection::S;
double XSection::S_sqrt;
double XSection::dC;
double XSection::dS;
boost::property_tree::ptree XSection::pt;
boost::program_options::variables_map XSection::vm;
const LHAPDF::PDF* XSection::pdf;
void print( string str, array<double,3> tree, array<double,3> virt, array<double,3> soft, array<double,3> hard) {
   cout << "\nResults for subprocess " << str << '\n';
   cout << scientific;
   // print out LO run statistics
   cout << "---------------------------------------------------------------" << endl;
   cout << setprecision(5);
   cout << setw(12) << "tree:" << setw(13) << tree.at(0)
         << " +/- " << setprecision(1) << tree.at(1)
         << " fb ( p-value = " << setw(8) << tree.at(2) << " )\n";
   cout << setprecision(5);
   cout << setw(12) << "virtual:" << setw(13) << virt.at(0) << " +/- "
           << setprecision(1) << virt.at(1) << " fb ( p-value = "
           << setw(8) << virt.at(2) << " )\n";

      cout << setprecision(5);
      cout << setw(12) << "real (soft):" << setw(13) << soft.at(0) << " +/- " << setprecision(1) << soft.at(1)
           << " fb ( p-value = " << setw(8) << soft.at(2) << " )\n";
      cout << setprecision(5);
      cout << setw(12) << "real (hard):" << setw(13) << hard.at(0) << " +/- " << setprecision(1) << hard.at(1)
           << " fb ( p-value = " << setw(8) << hard.at(2) << " )\n";
      cout << "---------------------------------------------------------------" << endl;
      cout << setprecision(5);
      cout << setw(12) << "sum:" << setw(13)
           << tree.at(0) + virt.at(0) + hard.at(0) + soft.at(0)
           << " +/- " << setprecision(1) << sqrt(pow(tree.at(1),2) +
           pow(virt.at(1),2) + pow(hard.at(1),2) +
           pow(soft.at(1),2)) << " fb" << endl;
}

void print( string str, array<double,3> tree) {
   cout << "\nResults for subprocess " << str << '\n';
      cout << scientific;
   //print out LO run statistics
   cout << "---------------------------------------------------------------" << endl;
   cout << setprecision(5);
   cout << setw(12) << "tree:" << setw(13) << tree.at(0)
         << " +/- " << setprecision(1) << tree.at(1)
         << " fb ( p-value = " << setw(8) << tree.at(2) << " )\n";
}

int main(int argc, char* argv[]) {

   // disable printint LHAPDF header since it gets printed twice
   LHAPDF::Info& cfg = LHAPDF::getConfig();
   cfg.set_entry("Verbosity", 0);

   //
   boost::program_options::positional_options_description p;
   p.add("card", -1);

   // program options
   boost::program_options::options_description desc("Allowed options");
   desc.add_options()
      ("help", "produce help message")
      ("precision-virt", po::value<int>() -> default_value(3), "")
      // gu_suLsuLdaggeru with SC precision 5 for BMP2 gives p-value 1
      ("precision-sc",   po::value<int>() -> default_value(6), "")
      // hard precision 3 gives p-valus of ~0.3 but on a home PC precision of 4
      // takes to long for a 'normal' run
      ("precision-hard", po::value<int>() -> default_value(3), "")
      ("enable-born", po::value<bool>() -> default_value(true), "")
      ("enable-virt", po::value<bool>() -> default_value(true), "")
      ("enable-sc",   po::value<bool>() -> default_value(true), "")
      ("enable-hard", po::value<bool>() -> default_value(true), "")
      ("card", po::value<string>(), "path to a run card")
      ("subprocess", po::value<string>() -> default_value(""), "")
   ;

   // verbosity of the integration routines
   boost::program_options::options_description verbosity("Integration verbosity");
   verbosity.add_options()
      ("verbosity-born", po::value<int>() -> default_value(0), "")
      ("verbosity-virt", po::value<int>() -> default_value(0), "")
      ("verbosity-sc",   po::value<int>() -> default_value(0), "")
      ("verbosity-hard", po::value<int>() -> default_value(0), "verbosity of hard collinear integration")
   ;
   desc.add(verbosity);

   boost::program_options::variables_map vm;
   boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
   boost::program_options::notify(vm);

   if (vm.count("help")) {
      cout << desc << '\n';
      return 1;
   }

   bool enable_born = vm["enable-born"].as<bool>();
   bool enable_virt = vm["enable-virt"].as<bool>();
   bool enable_sc = vm["enable-sc"].as<bool>();
   bool enable_hard = vm["enable-hard"].as<bool>();

   string card = vm["card"].as<string>();
   string subprocess = vm["subprocess"].as<string>();

   boost::property_tree::ptree pt;
   boost::property_tree::ini_parser::read_ini( card, pt );

   cout << "INFO: Using phase space slicing parameters dS = " << pt.get<double>("technical parameters.dS")
        << ", dC = " <<  pt.get<double>("technical parameters.dC") << '\n';

   // local arrays are not aumatically initialized to 0
   // need to use {}
   array<double,3> temp {}, xsection_tree {}, xsection_virt {}, xsection_SC {}, xsection_HnonC {},
           xsection_tree1 {}, xsection_virt1 {}, xsection_SC1 {}, xsection_HnonC1 {},
           xsection_tree2 {}, xsection_virt2 {}, xsection_SC2 {}, xsection_HnonC2 {},
           xsection_tree3 {}, xsection_virt3 {}, xsection_SC3 {}, xsection_HnonC3 {},
           xsection_tree4 {}, xsection_virt4 {}, xsection_SC4 {}, xsection_HnonC4 {},
           xsection_tree5 {}, xsection_virt5 {}, xsection_SC5 {}, xsection_HnonC5 {},
           xsection_tree_total {}, xsection_virt_total {}, xsection_SC_total {}, xsection_HnonC_total {};

   enum Model {
       MRSSM,
       MSSM,
       Simplified,
       no_model
   };
   enum Channel {
       pp_OO,
       pp_OsOs,
       pp_suLsuR,
       pp_suLsuL,
       pp_suLsdR,
       pp_suLsdL,
       pp_suLsuLdagger,
       pp_suLsuRdagger,
       pp_suLsdLdagger,
       pp_suLsdRdagger,
       no_channel
   };

   Model model = no_model;
   Channel channel = no_channel;

   if ( pt.get<string>("process.model") == "MRSSM" ) {
      model = MRSSM;
   } else if ( string(argv[1]) == "MSSM" ) {
      model = MSSM;
   } else if ( string(argv[1]) == "Simplified" ) {
      model = Simplified;
   } else {
	   cout << "\n Model not implemented! \n\n";
   }

   if ( pt.get<string>("process.process") == "pp_OsOs" ) {
      channel = pp_OsOs;
   } else if ( pt.get<string>("process.process") == "pp_suLsuR" ) {
      channel = pp_suLsuR;
   } else if ( pt.get<string>("process.process") == "pp_suLsuL" ) {
      channel = pp_suLsuL;
   } else if ( pt.get<string>("process.process") == "pp_suLsdR" ) {
      channel = pp_suLsdR;
   } else if ( pt.get<string>("process.process") == "pp_suLsdL" ) {
      channel = pp_suLsdL;
   } else if ( pt.get<string>("process.process") == "pp_suLsuLdagger" ) {
      channel = pp_suLsuLdagger;
   } else if ( pt.get<string>("process.process") == "pp_suLsuRdagger" ) {
      channel = pp_suLsuRdagger;
   } else if ( pt.get<string>("process.process") == "pp_suLsdLdagger" ) {
      channel = pp_suLsdLdagger;
   } else if ( pt.get<string>("process.process") == "pp_suLsdRdagger" ) {
      channel = pp_suLsdRdagger;
         } else if ( pt.get<string>("process.process") == "pp_OO" ) {
      channel = pp_OO;
   } else {
	   cout << "\n Process not implemented! \n\n";
   }

   auto start = chrono::steady_clock::now();

   json j;
   j["process"] = pt.get<string>("process.process");
   j["sqrt(S)"] = pt.get<double>("collider setup.sqrt_S");
   j["mu_r"] = pt.get<double>("collider setup.mu_r");
   j["mu_f"] = pt.get<double>("collider setup.mu_f");
   j["pdf"] = pt.get<string>("collider setup.pdf");
   j["masses"] = {
      {"gluino", pt.get<double>("masses.gluino")},
      {"pseudoscalar sgluon", pt.get<double>("masses.pseudoscalar_sgluon")},
      {"top", pt.get<double>("masses.top")},
      {"squarks", pt.get<double>("masses.squarks")}
   };

   if( pt.get<string>("process.order") == "LO" ) {
	  switch(model) {
         case MRSSM:
            switch(channel) {
               case pp_OsOs:
                  {
                  Process process1("sgluons-gg_OO", pt);
                  XSection::init( &process1, pt, vm );
                  XSection_Tree tree;
                  temp = tree.integrate();
                  Process process2("sgluons-qqbar_OO", pt);
                  XSection::init( &process2, pt, vm );
                  xsection_tree = tree.integrate() + temp;
                  print("pp > OO", xsection_tree);
                  break;
			      }
               case pp_suLsuR:
                  {                                                     // checked with MadGraph and Philip
                  Process process1("MRSSM,uu_suLsuR", pt);
                  XSection::init( &process1, pt, vm );
                  XSection_Tree tree;
                  xsection_tree = tree.integrate();
                  print("uu > suLsuR", xsection_tree);
                  break;
			      }
			   case pp_suLsdR:
                  {
				                                                        // checked with MadGraph and Philip
                  Process process1("MRSSM,ud_suLsdR", pt);
                  XSection::init( &process1, pt, vm );
                  XSection_Tree tree;
                  xsection_tree = tree.integrate();
                  print("uu > suLsdR", xsection_tree);
                  break;
			      }

               case pp_suLsuLdagger:
                  {
                  Process process1("MRSSM,GG_suLsuLdagger", pt);
                  XSection::init( &process1, pt, vm );
                  XSection_Tree tree;
                  xsection_tree1 = tree.integrate();
                  print("GG > suLsdLdagger", xsection_tree1);

                  Process process2("MRSSM,uubar_suLsuLdagger", pt);
                  XSection::init( &process2, pt, vm );
                  xsection_tree2 = tree.integrate();
                  print("uubar > suLsdLdagger", xsection_tree2);

                  Process process3("MRSSM,ddbar_suLsuLdagger", pt);
                  XSection::init( &process3, pt, vm );
                  xsection_tree3 = tree.integrate();
                  print("qqbar > suLsuLdagger", xsection_tree3);

                  xsection_tree_total = xsection_tree1 + xsection_tree2 + xsection_tree3;
                  print("pp > suLsdLdagger", xsection_tree_total);
                  break;
                  }
			   default:
			      {
			      xsection_tree = {0,0,0};
			      break;
			      }
            }
            break;
         case MSSM:
            switch(channel) {
			   //case pp_suLsuR:
			   //case pp_suLsuLdagger:
			   //default:
            }
         break;
      }
   }

   else if( pt.get<string>("process.order") == "NLO" ) {
	  switch(model) {
         case MRSSM:
            switch(channel) {
               case pp_OsOs:
                  {
                  // todo
                  break;
			      }
               case pp_suLsuR:
                  {
                  XSection_Tree tree;
                  XSection_Virt virt;
                  XSection_SC sc;
                  XSection_HnonC hc;

                  // uu > suL suR (+g) process
		            if( subprocess == "" ) {
                     Process process1("MRSSM,uu_suLsuR", pt);
	                  XSection::init( &process1, pt, vm );
                     if(enable_born) xsection_tree1 = tree.integrate();
                     if(enable_virt) xsection_virt1 = virt.integrate();
                     if(enable_sc) xsection_SC1 = sc.integrate();
                     if(enable_hard) xsection_HnonC1 = hc.integrate();
                     print( "uu > suLsuR(+X)", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
                     xsec_to_json(j, "uu->suLsuR(+X)", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
		            }

                  // the matrix element is regular in the limit dS -> 0 but the phase space parametrization
                  // fails if we are exactly on the threshold
                  pt.put( "technical parameters.dS", 1e-10 );

                  // gu > suL suR ubar process
		            if( subprocess == "gu_suLsuRubar" || subprocess == "" ) {
                     Process process2( "MRSSM,gu_suLsuR", pt);
                     XSection::init( &process2, pt, vm );
                     if(enable_sc) xsection_SC2 = sc.integrate();
                     if(enable_hard) xsection_HnonC2 = hc.integrate();
                     print( "gu > suLsuR(+X)", xsection_tree2, xsection_virt2, xsection_SC2, xsection_HnonC2 );
                     xsec_to_json(j, "gu->suLsuR(+X)", xsection_tree2, xsection_virt2, xsection_SC2, xsection_HnonC2);
		            }

                  xsection_tree_total = xsection_tree1;
                  xsection_virt_total = xsection_virt1;
                  xsection_SC_total = xsection_SC1 + xsection_SC2;
                  xsection_HnonC_total = xsection_HnonC1 + xsection_HnonC2;
                  print( "sum", xsection_tree1, xsection_virt1, xsection_SC_total, xsection_HnonC_total );
                  break;
			      }
               case pp_suLsuLdagger:
                  {
                  XSection_Tree tree;
                  XSection_Virt virt;
                  XSection_SC sc;
                  XSection_HnonC hc;

		            if( subprocess == "" || subprocess == "uubar_suLsuLdagger" ) {
                     Process process1("MRSSM,uubar_suLsuLdagger", pt);
                     XSection::init( &process1, pt, vm );
                     if(enable_born) xsection_tree1 = tree.integrate();
                     if(enable_virt) xsection_virt1 = virt.integrate();
                     if(enable_sc) xsection_SC1 = sc.integrate();
                     if(enable_hard) xsection_HnonC1 = hc.integrate();
                     print( "uubar > suLsuL*", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
                     xsec_to_json(j, "uubar->suLsuL*", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
	    	         }

                  if( subprocess == "") {
                     Process process2("MRSSM,ddbar_suLsuLdagger", pt);
                     XSection::init( &process2, pt, vm );
                     if(enable_born) xsection_tree2 = tree.integrate();
                     if(enable_virt) xsection_virt2 = virt.integrate();
                     if(enable_sc) xsection_SC2 = sc.integrate();
                     if(enable_hard) xsection_HnonC2 = hc.integrate();
                     print( "ddbar->suLsuL*", xsection_tree2, xsection_virt2, xsection_SC2, xsection_HnonC2);
                     xsec_to_json(j, "ddbar->suLsuL*", xsection_tree2, xsection_virt2, xsection_SC2, xsection_HnonC2);
                  }

                  if( subprocess == "") {
                     Process process3("MRSSM,GG_suLsuLdagger", pt);
                     XSection::init( &process3, pt, vm );
                     if(enable_born) xsection_tree3 = tree.integrate();
                     if(enable_virt) xsection_virt3 = virt.integrate();
                     if(enable_sc) xsection_SC3 = sc.integrate();
                     if(enable_hard) xsection_HnonC3 = hc.integrate();
                     print( "gg > suLsuL*", xsection_tree3, xsection_virt3, xsection_SC3, xsection_HnonC3);
                     xsec_to_json(j, "gg->suLsuL*", xsection_tree3, xsection_virt3, xsection_SC3, xsection_HnonC3);
                  }

                  // the matrix element is regular in the limit dS -> 0 but the phase space parametrization
                  // fails if we are exactly on the threshold
                  pt.put( "technical parameters.dS", 1e-10 );

                  if( subprocess == "") {
                     Process process4("MRSSM,gq_suLsuLdagger", pt);
                     XSection::init( &process4, pt, vm );
                     if(enable_sc) xsection_SC4 = sc.integrate();
                     if(enable_hard) xsection_HnonC4 = hc.integrate();
                     print( "gq > suLsuL*(+X)", xsection_tree4, xsection_virt4, xsection_SC4, xsection_HnonC4);
                     xsec_to_json(j, "gq->suLsuL*(+X)", xsection_tree4, xsection_virt4, xsection_SC4, xsection_HnonC4);
                  }

                  // g u > suL suLdagger 
                  pt.put( "technical parameters.dS", 1e-9 );
                  if( subprocess == "gu_suLsuLdaggeru" || subprocess == "" ) {
                     Process process5("MRSSM,gu_suLsuLdagger", pt);
                     XSection::init( &process5, pt, vm );
                     if(enable_sc) xsection_SC5 = sc.integrate();
                     if(enable_hard) xsection_HnonC5 = hc.integrate();
                     print( "gu > suLsuL*(+X)", xsection_tree5, xsection_virt5, xsection_SC5, xsection_HnonC5 );
                     xsec_to_json(j, "gu->suLsuL*(+X)", xsection_tree5, xsection_virt5, xsection_SC5, xsection_HnonC5);
                  }

                  xsection_tree_total = xsection_tree1 + xsection_tree2 + xsection_tree3;
                  xsection_virt_total = xsection_virt1 + xsection_virt2 + xsection_virt3;
                  xsection_SC_total = xsection_SC1 + xsection_SC2 + xsection_SC3 + xsection_SC4 + xsection_SC5;
                  xsection_HnonC_total = xsection_HnonC1 + xsection_HnonC2 + xsection_HnonC3
                          + xsection_HnonC4 + xsection_HnonC5;
                  print( "total", xsection_tree_total, xsection_virt_total, xsection_SC_total, xsection_HnonC_total);
                  break;
                  }
               default:
                  cout << "NLO process not implemented\n";
            }
     case Simplified:
        switch(channel) {
           case pp_OO:
           {
              {
           }
        }
        }
     }
   }

   auto end = chrono::steady_clock::now();
   cout << '\n';
   cout << "INFO: Calculation ended after ";
   if (end - start > 1h)
      cout << chrono::duration_cast<chrono::hours>(end-start).count() << " hour(s), ";
   if (end - start > 1min)
      cout << chrono::duration_cast<chrono::minutes>(end-start).count() %  60 << " minute(s) and ";
   cout << chrono::duration_cast<chrono::seconds>(end-start).count() % 60 << " second(s)\n";

   std::ofstream o(
      pt.get<string>("process.process") + "_" +
      to_string(pt.get<double>("masses.squarks")) + "_" +
      to_string(pt.get<double>("masses.gluino")) + "_" +
      to_string(pt.get<double>("masses.pseudoscalar_sgluon")) + "_" +
      to_string(pt.get<double>("collider setup.sqrt_S")) + "_" +
      to_string(pt.get<double>("collider setup.mu_r")) + "_" +
      to_string(pt.get<double>("collider setup.mu_f")) + "_" +
      pt.get<string>("collider setup.pdf") +
      ".json");
   o << std::setw(3) << j << std::endl;

   return 0;
}
