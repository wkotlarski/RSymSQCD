#include <chrono>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "XSection_Tree.hpp"
#include "XSection_Virt.hpp"
#include "XSection_SC.hpp"
#include "XSection_HnonC.hpp"

using namespace std;

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

// why do I have to write this?
// why isn't init() enough?
Process *XSection::processID; 
double XSection::mu_r;
double XSection::mu_f;
double XSection::m1;
double XSection::m2;
double XSection::prec_virt;
double XSection::prec_sc;
double XSection::prec_hnc;
double XSection::S;
double XSection::S_sqrt;
double XSection::dC;
double XSection::dS;
boost::property_tree::ptree XSection::pt;
const LHAPDF::PDF* XSection::pdf;

void print( string str, array<double,3> tree, array<double,3> virt, array<double,3> soft, array<double,3> hard) {
   cout << "\nResults for subprocess " << str << '\n';
      cout << scientific;
   //print out LO run statistics
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
   
   // program options
   boost::program_options::options_description desc("Allowed options");
   desc.add_options()
      ("help", "produce help message")
      ("compression", po::value<double>(), "set compression level")
      ("precision-virt", po::value<int>() -> default_value(4), "")
      ("precision-sc",   po::value<int>() -> default_value(4), "")
      ("precision-hard", po::value<int>() -> default_value(4), "")
      ("enable-born", po::value<bool>() -> default_value(true), "")
      ("enable-virt", po::value<bool>() -> default_value(true), "")
      ("enable-sc",   po::value<bool>() -> default_value(true), "")
      ("enable-hard", po::value<bool>() -> default_value(true), "")
      ("card", po::value<string>(), "path to a run card")
      ("subprocess", po::value<string>() -> default_value(""), "")
   ;

   boost::program_options::variables_map vm;        
   boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
   boost::program_options::notify(vm);

   if (vm.count("help")) {
      cout << desc << "\n";
      return 0;
   }

   bool enable_born = vm["enable-born"].as<bool>();
   bool enable_virt = vm["enable-virt"].as<bool>();
   bool enable_sc = vm["enable-sc"].as<bool>();
   bool enable_hard = vm["enable-hard"].as<bool>();

   double prec_virt = pow( 10., -vm["precision-virt"].as<int>() );
   double prec_sc   = pow( 10., -vm["precision-sc"].as<int>() );
   double prec_hard = pow( 10., -vm["precision-hard"].as<int>() );
   cout << prec_virt << ' ' << prec_sc << ' ' << prec_hard << '\n';

   string card = vm["card"].as<string>();
   string subprocess = vm["subprocess"].as<string>();

/* invoke programm like 
   "./RSymSQCD MSSM pp_suLsuR NLO 1 1 1 run.ini" for  NLO calculation  
       with last three numbers giving the desired accuracy of the virtual,
       soft-collinear and hard-noncollinear part, respectively 
   "./RSymSQCD MSSM pp_suLsuR LO" for  LO calculation  */
   
   boost::property_tree::ptree pt;
   boost::property_tree::ini_parser::read_ini( card, pt );   

   cout << "INFO: using dS = " << pt.get<double>("technical parameters.dS")
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

   if( pt.get<string>("process.order") == "LO" ) {
	  switch(model) {
         case MRSSM:
            switch(channel) {
               case pp_OsOs:
                  {
                  Process process1("sgluons-gg_OO", pt);  
                  XSection::init( &process1, pt, 1, 1, 1 );
                  XSection_Tree tree;
                  temp = tree.integrate();     
                  Process process2("sgluons-qqbar_OO", pt);
                  XSection::init( &process2, pt, 1, 1, 1 );
                  xsection_tree = tree.integrate() + temp;   
                  print("pp > OO", xsection_tree);
                  break;
			      }                  
               case pp_suLsuR: 
                  {                                                     // checked with MadGraph and Philip
                  Process process1("MRSSM,uu_suLsuR", pt);
                  XSection::init( &process1, pt, 1, 1, 1 );
                  XSection_Tree tree;
                  xsection_tree = tree.integrate();
                  print("uu > suLsuR", xsection_tree);
                  break;      	
			      }	
			   case pp_suLsdR: 
                  {  
				                                                        // checked with MadGraph and Philip 
                  Process process1("MRSSM,ud_suLsdR", pt);
                  XSection::init( &process1, pt, 1, 1, 1 );
                  XSection_Tree tree;
                  xsection_tree = tree.integrate();
                  print("uu > suLsdR", xsection_tree);
                  break;      	
			      }	

               case pp_suLsuLdagger:
                  {
                  Process process1("MRSSM,GG_suLsuLdagger", pt);
                  XSection::init( &process1, pt, 1, 1, 1 );
                  XSection_Tree tree;
                  xsection_tree1 = tree.integrate();
                  print("GG > suLsdLdagger", xsection_tree1);
                  
                  Process process2("MRSSM,uubar_suLsuLdagger", pt);
                  XSection::init( &process2, pt, 1, 1, 1 );
                  xsection_tree2 = tree.integrate();
                  print("uubar > suLsdLdagger", xsection_tree2);
                  
                  Process process3("MRSSM,ddbar_suLsuLdagger", pt);
                  XSection::init( &process3, pt, 1, 1, 1 );
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
		            if( atoi(argv[5]) == 1 || subprocess == "" ) {
                     Process process1("MRSSM,uu_suLsuR", pt);
	                  XSection::init( &process1, pt, prec_virt, prec_sc, prec_hard );
                     if(enable_born) xsection_tree1 = tree.integrate();      
                     if(enable_virt) xsection_virt1 = virt.integrate();
                     if(enable_sc) xsection_SC1 = sc.integrate();
                     if(enable_hard) xsection_HnonC1 = hc.integrate();
                     print( "uu > suLsuR(+X)", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
		            }
                  
                  // the matrix element is regular in the limit dS -> 0 but the phase space parametrization
                  // fails if we are exactly on the threshold
                  pt.put( "technical parameters.dS", 1e-10 );                  

                  // gu > suL suR ubar process
		            if( atoi(argv[5]) == 2 || subprocess == "" ) {
                     Process process2( "MRSSM,gu_suLsuR", pt);
                     XSection::init( &process2, pt, prec_virt, prec_sc, prec_hard );      
                     if(enable_sc) xsection_SC2 = sc.integrate();
                     if(enable_hard) xsection_HnonC2 = hc.integrate();
                     print( "gu > suLsuR(+X)", xsection_tree2, xsection_virt2, xsection_SC2, xsection_HnonC2 );
		            }
                  
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

		            if( atoi(argv[5]) == 1 || subprocess == "" ) {
                     Process process1("MRSSM,uubar_suLsuLdagger", pt);
                     XSection::init( &process1, pt, prec_virt, prec_sc, prec_hard );                 
                     if(enable_born) xsection_tree1 = tree.integrate();      
                     if(enable_virt) xsection_virt1 = virt.integrate();
                     if(enable_sc) xsection_SC1 = sc.integrate();
                     if(enable_hard) xsection_HnonC1 = hc.integrate();
                     print( "uubar > suLsuL*", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
	    	         }
         
		            if( atoi(argv[5]) == 2 || subprocess == "") {
                     Process process2("MRSSM,ddbar_suLsuLdagger", pt);
                     XSection::init( &process2, pt, prec_virt, prec_sc, prec_hard );                    
                     if(enable_born) xsection_tree2 = tree.integrate();
                     if(enable_virt) xsection_virt2 = virt.integrate();
                     if(enable_sc) xsection_SC2 = sc.integrate();
                     if(enable_hard) xsection_HnonC2 = hc.integrate();
                     print( "ddbar > suLsuL*", xsection_tree2, xsection_virt2, xsection_SC2, xsection_HnonC2);
		            }

		            if( atoi(argv[5]) == 3 || subprocess == "") {
                     Process process3("MRSSM,GG_suLsuLdagger", pt);
                     XSection::init( &process3, pt, prec_virt, prec_sc, prec_hard );
                     if(enable_born) xsection_tree3 = tree.integrate();
                     if(enable_virt) xsection_virt3 = virt.integrate();
                     if(enable_sc) xsection_SC3 = sc.integrate();
                     if(enable_hard) xsection_HnonC3 = hc.integrate();
                     print( "gg > suLsuL*", xsection_tree3, xsection_virt3, xsection_SC3, xsection_HnonC3);
		            }
                  
                  // the matrix element is regular in the limit dS -> 0 but the phase space parametrization
                  // fails if we are exactly on the threshold
                  pt.put( "technical parameters.dS", 1e-10 );         

  		            if( atoi(argv[5]) == 4 || subprocess == "") {
                     Process process4("MRSSM,gq_suLsuLdagger", pt);
                     XSection::init( &process4, pt, prec_virt, prec_sc, prec_hard );
                     if(enable_sc) xsection_SC4 = sc.integrate();
                     if(enable_hard) xsection_HnonC4 = hc.integrate();
                     print( "gq > suLsuL*(+X)", xsection_tree4, xsection_virt4, xsection_SC4, xsection_HnonC4);
                  }

                  // g u > suL suLdagger 
		            if( atoi(argv[5]) == 5 || subprocess == "" ) {
                     Process process5("MRSSM,gu_suLsuLdagger", pt);
                     XSection::init( &process5, pt, prec_virt, prec_sc, prec_hard );
                     if(enable_sc) xsection_SC5 = sc.integrate();
                     if(enable_hard) xsection_HnonC5 = hc.integrate();
                     print( "gu > suLsuL*(+X)", xsection_tree5, xsection_virt5, xsection_SC5, xsection_HnonC5 );
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
   } else {
      cout << "Third command line argument must be 'LO' or 'NLO'." << endl;
      return 0;
   }

   auto end = chrono::steady_clock::now();
   cout << "Calculation ended after " 
      << chrono::duration_cast<chrono::minutes>(end-start).count()
      << " minutes\n";   
   
   return 0;
}
