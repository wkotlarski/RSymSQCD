#include "version.hpp"
#include "XSection_Tree.hpp"
#include "XSection_Virt.hpp"
#include "XSection_SC.hpp"
#include "XSection_HnonC.hpp"
#include "splitting_kernels.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/program_options.hpp>
#include "LHAPDF/LHAPDF.h"
#include <nlohmann/json.hpp>

#include <chrono>
#include <functional>

using namespace std;
using namespace std::chrono_literals;
namespace po = boost::program_options;
using json = nlohmann::json;
using namespace std::placeholders;

/*
 *    XSection_* return results in the form xsection, its error, erros's pvalue
 *    overload + operator to add xsection and errors in quadrature
 *    pvalues are just added
 */
inline array<double, 3> operator+(array<double, 3> const& x, array<double, 3> const& y) {
   return {x.at(0) + y.at(0), std::hypot(x.at(1), y.at(1)), std::max(x.at(2), y.at(2))};
}

void xsec_to_json(json& j, std::string const& str, array<double, 3> const& tree, array<double, 3> const& virt, array<double, 3> const& soft, array<double, 3> const& hard) {
   j["cross sections"][str] = {
      {"tree", {{"res", tree.at(0)}, {"err", tree.at(1)}, {"p-val", tree.at(2)}}},
      {"virtual", {{"res", virt.at(0)}, {"err", virt.at(1)}, {"p-val", virt.at(2)}}},
      {"SC", {{"res", soft.at(0)}, {"err", soft.at(1)}, {"p-val", soft.at(2)}}},
      {"HnonC", {{"res", hard.at(0)}, {"err", hard.at(1)}, {"p-val", hard.at(2)}}}
   };
}

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

   boost::program_options::positional_options_description p;
   p.add("card", -1);

   // program options
   boost::program_options::options_description desc("Allowed options");
   desc.add_options()
      ("help,h",    "produce help message")
      ("version,v", "display the version number")
      ("json-outputfile-name", po::value<string>(), "name of output file in JSON format")
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
      return 0;
   }
   if (vm.count("version")) {
      std::cout << RSymSQCD_VERSION << std::endl;
      return 0;
   }

   bool enable_born = vm["enable-born"].as<bool>();
   bool enable_virt = vm["enable-virt"].as<bool>();
   bool enable_sc = vm["enable-sc"].as<bool>();
   bool enable_hard = vm["enable-hard"].as<bool>();

   string card = vm["card"].as<string>();
   string subprocess = vm["subprocess"].as<string>();

   boost::property_tree::ptree pt;
   boost::property_tree::ini_parser::read_ini(card, pt);

   cout << "\nINFO: Using phase space slicing parameters dS = " << pt.get<double>("technical parameters.dS")
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

   enum class Model {
       MRSSM,
       MSSM,
       Simplified,
       no_model
   };

   Model model = Model::no_model;
   if (pt.get<string>("process.model") == "MRSSM") {
      model = Model::MRSSM;
   }
   else if (pt.get<string>("process.model") == "MSSM") {
      model = Model::MSSM;
   }
   else if (pt.get<string>("process.model") == "Simplified") {
       model = Model::Simplified;
   }
   else {
      std::cout << "\nModel not implemented!\n\n";
      return 1;
   }

   enum class Channel {
       pp_OO,
       pp_OsOs,
       pp_suLsuR,
       pp_suLsuL,
       pp_suLsdR,
       pp_suLsdL,
       pp_sqLsqR,
       pp_suLsuLdagger,
       pp_suLsuRdagger,
       pp_suLsdLdagger,
       pp_suLsdRdagger,
       no_channel
   };

   Channel channel = Channel::no_channel;
   if ( pt.get<string>("process.process") == "pp_OsOs" ) {
      channel = Channel::pp_OsOs;
   } else if ( pt.get<string>("process.process") == "pp_suLsuR" ) {
      channel = Channel::pp_suLsuR;
   } else if ( pt.get<string>("process.process") == "pp_sqLsqR" ) {
      channel = Channel::pp_sqLsqR;
   } else if ( pt.get<string>("process.process") == "pp_suLsuL" ) {
      channel = Channel::pp_suLsuL;
   } else if ( pt.get<string>("process.process") == "pp_suLsdR" ) {
      channel = Channel::pp_suLsdR;
   } else if ( pt.get<string>("process.process") == "pp_suLsdL" ) {
      channel = Channel::pp_suLsdL;
   } else if ( pt.get<string>("process.process") == "pp_suLsuLdagger" ) {
      channel = Channel::pp_suLsuLdagger;
   } else if ( pt.get<string>("process.process") == "pp_suLsuRdagger" ) {
      channel = Channel::pp_suLsuRdagger;
   } else if ( pt.get<string>("process.process") == "pp_suLsdLdagger" ) {
      channel = Channel::pp_suLsdLdagger;
   } else if ( pt.get<string>("process.process") == "pp_suLsdRdagger" ) {
      channel = Channel::pp_suLsdRdagger;
         } else if ( pt.get<string>("process.process") == "pp_OO" ) {
      channel = Channel::pp_OO;
   } else {
	   cout << "\n Process not implemented! \n\n";
   }

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
   j["technical parameters"] = {
      {"dS", pt.get<double>("technical parameters.dS")},
      {"dC", pt.get<double>("technical parameters.dC")},
      {"WidthOverMass", pt.get<double>("technical parameters.WidthOverMass")},
      {"eta_sign", pt.get<double>("technical parameters.eta_sign")},
      {"delta", pt.get<double>("technical parameters.delta")}
   };

   auto start = chrono::steady_clock::now();

   const double dS = pt.get<double>("technical parameters.dS");
   const double dC = pt.get<double>("technical parameters.dC");
   const double muR = pt.get<double>("collider setup.mu_r");
   const double muF = pt.get<double>("collider setup.mu_f");

   // set PDFs
   std::cout << '\n';
   std::unique_ptr<LHAPDF::PDF> pdf(LHAPDF::mkPDF(pt.get<std::string>("collider setup.pdf"), 0));
   LHAPDF::Info& cfg = LHAPDF::getConfig();
   cfg.set_entry("Verbosity", 0);
   LHAPDF::setVerbosity(0);

   if( pt.get<string>("process.order") == "LO" ) {
      /*
      switch(model) {
         case Model::MRSSM:
            switch(channel) {
               case Channel::pp_OsOs:
               {
                  Process process1("sgluons-gg_OO", pt);
                  XSection::init(std::move(process1), pt, vm );
                  XSection_Tree tree;
                  temp = tree.integrate();
                  Process process2("sgluons-qqbar_OO", pt);
                  XSection::init(std::move(process2), pt, vm );
                  xsection_tree = tree.integrate() + temp;
                  print("pp > OO", xsection_tree);
                  break;
			      }
               case Channel::pp_suLsuR:
               {                                                     // checked with MadGraph and Philip
                  Process process("MRSSM,uu_suLsuR", pt);
                  XSection::init(std::move(process), pt, vm );
                  XSection_Tree tree;
                  xsection_tree = tree.integrate();
                  print("uu > suLsuR", xsection_tree);
                  break;
			      }
               case Channel::pp_suLsdR:
               {                                                     // checked with MadGraph and Philip
                  Process process("MRSSM,ud_suLsdR", pt);
                  XSection::init(std::move(process), pt, vm );
                  XSection_Tree tree;
                  xsection_tree = tree.integrate();
                  print("uu > suLsdR", xsection_tree);
                  break;
			      }
               case Channel::pp_suLsuLdagger:
               {
                  Process process1("MRSSM,GG_suLsuLdagger", pt);
                  XSection::init(std::move(process1), pt, vm );
                  XSection_Tree tree;
                  xsection_tree1 = tree.integrate();
                  print("GG > suLsdLdagger", xsection_tree1);

                  Process process2("MRSSM,uubar_suLsuLdagger", pt);
                  XSection::init(std::move(process2), pt, vm );
                  xsection_tree2 = tree.integrate();
                  print("uubar > suLsdLdagger", xsection_tree2);

                  Process process3("MRSSM,ddbar_suLsuLdagger", pt);
                  XSection::init(std::move(process3), pt, vm );
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
         case Model::MSSM:
            switch(channel) {
               //case pp_suLsuR:
               //case pp_suLsuLdagger:
               //default:
            }
         break;
      }
      */
   }
   else if (pt.get<string>("process.order") == "NLO") {
	   switch(model) {
         case Model::MRSSM:
         {
            Process mrssm(pt);
            switch(channel) {
               case Channel::pp_suLsuR:
               {
                  const double m1 = pt.get<double>("masses.squarks");
                  const double m2 = pt.get<double>("masses.squarks");
                  // uu > suL suR (+g) process
		            if( subprocess == "" ) {
                     std::vector<std::array<int, 3>> flav {{2,2,1}};
                     XSection_Tree tree(
                        m1, m2,
                        // same ME as in the MSSM
                        std::bind(&Process::matrixMRSSMTree_uu_suLsuR, mrssm, _1, _2, _3),
                        flav, muR, muF, pdf.get()
                     );
                     XSection_Virt virt(
                        m1, m2,
                        std::bind(&Process::matrixMRSSMVirt_uu_suLsuR, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav, muR, muF, pdf.get()
                     );
                     XSection_SC sc(
                        m1, m2,
                        std::bind(&Process::matrixMRSSMSoft_uu_suLsuRg, mrssm, _1, _2, _3, _4, _5),
                        dS, dC, muR, muF,
                        flav, pdf.get(),
                        {
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqq, std::bind(&Process::sigmaMRSSMTree_uu_suLsuR, mrssm, _1, _2)},
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqq, std::bind(&Process::sigmaMRSSMTree_uu_suLsuR, mrssm, _1, _2)}
                        }
                     );
                     XSection_HnonC hc(
                        m1, m2,
                        std::bind(&Process::matrixMRSSMHard_uu_suLsuRg, mrssm, _1, _2),
                        dS, dC, muR, muF,
                        flav, pdf.get()
                     );
	                  tree.init(pt, vm );
	                  virt.init(pt, vm );
	                  sc.init(pt, vm );
	                  hc.init(pt, vm );
                     if(enable_born) xsection_tree1 = tree.integrate();
                     if(enable_virt) xsection_virt1 = virt.integrate();
                     if(enable_sc) xsection_SC1 = sc.integrate();
                     if(enable_hard) xsection_HnonC1 = hc.integrate();
                     print( "uu > suLsuR(+X)", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
                     xsec_to_json(j, "uu->suLsuR(+X)", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
		            }

                  // the matrix element is regular in the limit dS -> 0 but the phase space parametrization
                  // fails if we are exactly on the threshold

                  // gu > suL suR ubar process
		            if (subprocess == "gu_suLsuRubar" || subprocess == "" ) {
                     std::vector<std::array<int, 3>> flav {};
                     for(int el : {2, -2}) flav.push_back({21, el, 2});
                     XSection_SC sc(
                        pt.get<double>("masses.squarks"), pt.get<double>("masses.squarks"),
                        std::bind(&Process::matrix_soft_stub, mrssm, _1, _2, _3, _4, _5),
                        0. /*dS*/, dC, muR, muF,
                        flav, pdf.get(),
                        {
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqg, std::bind(&Process::sigmaMRSSMTree_uu_suLsuR, mrssm, _1, _2)},
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pgq, std::bind(&Process::matrix_xsec_stub, mrssm, _1, _2)}
                        }

                     );
                     XSection_HnonC hc(
                        pt.get<double>("masses.squarks"), pt.get<double>("masses.squarks"),
                        std::bind(&Process::matrixMRSSMHard_gu_suLsuRubar, mrssm, _1, _2),
                        0. /*dS*/, dC, muR, muF,
                        flav, pdf.get()
                     );
	                  sc.init(pt, vm );
	                  hc.init(pt, vm );
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
               case Channel::pp_sqLsqR:
               {
                  // qq > sqL sqR (+g) process
		            if( subprocess == "" ) {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {1, 2, 3, 4, 5}) {
                        flav.push_back({ i,  i, 1});
                        flav.push_back({-i, -i, 1});
                     }
                     XSection_Tree tree(pt.get<double>("masses.squarks"), pt.get<double>("masses.squarks"), std::bind(&Process::matrixMRSSMTree_uu_suLsuR, mrssm, _1, _2, _3), flav,
                        pt.get<double>("collider setup.mu_r"),
                        pt.get<double>("collider setup.mu_f")
                        , pdf.get()
);
                     XSection_Virt virt(
                        pt.get<double>("masses.squarks"), pt.get<double>("masses.squarks"),
                        std::bind(&Process::matrixMRSSMVirt_uu_suLsuR, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav,
                        pt.get<double>("collider setup.mu_r"),
                        pt.get<double>("collider setup.mu_f")
                        , pdf.get()

                     );
                     XSection_SC sc(
                        pt.get<double>("masses.squarks"), pt.get<double>("masses.squarks"),
                        std::bind(&Process::matrixMRSSMSoft_uu_suLsuRg, mrssm, _1, _2, _3, _4, _5),
                        dS, dC, muR, muF,
                        flav
                        , pdf.get(),
                        {
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqq, std::bind(&Process::sigmaMRSSMTree_uu_suLsuR, mrssm, _1, _2)},
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqq, std::bind(&Process::sigmaMRSSMTree_uu_suLsuR, mrssm, _1, _2)}
                        }

                     );

                     XSection_HnonC hc(pt.get<double>("masses.squarks"), pt.get<double>("masses.squarks"), std::bind(&Process::matrixMRSSMHard_uu_suLsuRg, mrssm, _1, _2),
                        dS, dC, muR, muF,
                           flav
                        , pdf.get()
                           
                           );
	                  tree.init(pt, vm );
	                  virt.init(pt, vm );
	                  sc.init(pt, vm );
	                  hc.init(pt, vm );
                     if(enable_born) xsection_tree1 = tree.integrate();
                     if(enable_virt) xsection_virt1 = virt.integrate();
                     if(enable_sc) xsection_SC1 = sc.integrate();
                     if(enable_hard) xsection_HnonC1 = hc.integrate();
                     print( "qq > sqLsqR(+X)", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
                     xsec_to_json(j, "qq->sqLsqR(+X)", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
		            }
                  // qq' > sqL sq'R (+g) process
		            if( subprocess == "" ) {
                     XSection_Tree tree(pt.get<double>("masses.squarks"), pt.get<double>("masses.squarks"), std::bind(&Process::matrixMRSSMTree_uu_suLsuR, mrssm, _1, _2, _3), {},
                        pt.get<double>("collider setup.mu_r"),
                        pt.get<double>("collider setup.mu_f"),
                        pdf.get()
);
                     XSection_Virt virt(
                        pt.get<double>("masses.squarks"), pt.get<double>("masses.squarks"),
                        std::bind(&Process::matrixMRSSMVirt_uu_suLsuR, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        {},
                        pt.get<double>("collider setup.mu_r"),
                        pt.get<double>("collider setup.mu_f"),
                        pdf.get()

                     );
                     //XSection_HnonC hc;
	                  tree.init(pt, vm );
	                  virt.init(pt, vm );
	                  //sc.init(pt, vm );
	                  //hc.init(pt, vm );
                     //if(enable_born) xsection_tree2 = tree.integrate();
                     if(enable_virt) xsection_virt2 = virt.integrate();
                     //if(enable_sc) xsection_SC2 = sc.integrate();
                     //if(enable_hard) xsection_HnonC2 = hc.integrate();
                     print( "qq' > sqLsqR'(+X)", xsection_tree2, xsection_virt2, xsection_SC2, xsection_HnonC2);
                     xsec_to_json(j, "qq->sqLsqR'(+X)", xsection_tree2, xsection_virt2, xsection_SC2, xsection_HnonC2);
		            }

                  // the matrix element is regular in the limit dS -> 0 but the phase space parametrization
                  // fails if we are exactly on the threshold
                  pt.put( "technical parameters.dS", 1e-10 );

                  // gu > suL suR ubar process
		            if( subprocess == "gq_sqLsqRqbar" || subprocess == "" ) {
                     /*
                     Process process( "MRSSM,gd_sdLsdR", pt);
                     XSection::init(std::move(process), pt, vm );
                     if(enable_sc) xsection_SC3 = sc.integrate();
                     if(enable_hard) xsection_HnonC3 = hc.integrate();
                     print( "gq > sqLsqR(+X)", xsection_tree3, xsection_virt3, xsection_SC3, xsection_HnonC3 );
                     xsec_to_json(j, "gq->sqLsqR(+X)", xsection_tree3, xsection_virt3, xsection_SC3, xsection_HnonC3);
                     */
		            }

                  xsection_tree_total = xsection_tree1 + xsection_tree2;
                  xsection_virt_total = xsection_virt1 + xsection_virt2;
                  xsection_SC_total = xsection_SC1 + xsection_SC2 + xsection_SC3;
                  xsection_HnonC_total = xsection_HnonC1 + xsection_HnonC2 + xsection_HnonC3;
                  print( "sum", xsection_tree_total, xsection_virt_total, xsection_SC_total, xsection_HnonC_total );
                  break;
			      }
               case Channel::pp_suLsuLdagger:
               {
                  const double m1 = pt.get<double>("masses.squarks");
                  const double m2 = pt.get<double>("masses.squarks");
		            if (subprocess == "" || subprocess == "uubar_suLsuLdagger") {
                     std::vector<std::array<int, 3>> flav {{2, -2, 2}};
                     XSection_Tree tree(
                        m1, m2, 
                        std::bind(&Process::matrixMRSSMTree_uubar_suLsuLdagger, mrssm, _1, _2, _3), flav,
                        muR, muF, pdf.get()
                     );
                     XSection_Virt virt(
                        m1, m2, 
                        std::bind(&Process::matrixMRSSMVirt_uubar_suLsuLdagger, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav,
                        muR, muF, pdf.get()
                     );
                     XSection_SC sc(
                        m1, m2, 
                        std::bind(&Process::matrixMRSSMSoft_uubar_suLsuLdaggerg, mrssm, _1, _2, _3, _4, _5),
                        dS, dC, muR, muF,
                        flav, pdf.get(),
                        {
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqq, std::bind(&Process::sigmaMRSSMTree_uubar_suLsuLdagger, mrssm, _1, _2)},
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqq, std::bind(&Process::sigmaMRSSMTree_uubar_suLsuLdagger, mrssm, _1, _2)}
                        }

                     );
                     XSection_HnonC hc(
                        m1, m2, 
                        std::bind(&Process::matrixMRSSMHard_uubar_suLsuLdaggerg, mrssm, _1, _2), 
                        dS, dC, muR, muF,
                        flav, pdf.get()
                     );
	                  tree.init(pt, vm );
	                  virt.init(pt, vm );
	                  sc.init(pt, vm );
	                  hc.init(pt, vm );
                     if(enable_born) xsection_tree1 = tree.integrate();
                     if(enable_virt) xsection_virt1 = virt.integrate();
                     if(enable_sc) xsection_SC1 = sc.integrate();
                     if(enable_hard) xsection_HnonC1 = hc.integrate();
                     print( "uubar > suLsuL*", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
                     xsec_to_json(j, "uubar->suLsuL*", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
	    	         }

                  if( subprocess == "") {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {1, 3, 4, 5}) {
                        flav.push_back({i,-i, 2});
                     }
                     XSection_Tree tree(
                        m1, m2,
                        std::bind(&Process::matrixMRSSMTree_ddbar_suLsuLdagger, mrssm, _1, _2, _3), flav,
                        muR, muF, pdf.get()
                     );
                     XSection_Virt virt(
                        m1, m2,
                        std::bind(&Process::matrixMRSSMVirt_ddbar_suLsuLdagger, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav,
                        muR, muF, pdf.get()

                     );
                     XSection_SC sc(
                        m1, m2,
                        std::bind(&Process::matrixMRSSMSoft_ddbar_suLsuLdaggerg, mrssm, _1, _2, _3, _4, _5),
                        dS, dC, muR, muF,
                        flav
                        , pdf.get(),
                        {
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqq, std::bind(&Process::sigmaMRSSMTree_ddbar_suLsuLdagger, mrssm, _1, _2)},
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqq, std::bind(&Process::sigmaMRSSMTree_ddbar_suLsuLdagger, mrssm, _1, _2)}
                        }

                     );
                     XSection_HnonC hc(
                        m1, m2,
                        std::bind(&Process::matrixMRSSMHard_ddbar_suLsuLdaggerg, mrssm, _1, _2), 
                        dS, dC, muR, muF,
                           flav
                        , pdf.get()
                           );
	                  tree.init(pt, vm );
	                  virt.init(pt, vm );
	                  sc.init(pt, vm );
	                  hc.init(pt, vm );
                     if(enable_born) xsection_tree2 = tree.integrate();
                     if(enable_virt) xsection_virt2 = virt.integrate();
                     if(enable_sc) xsection_SC2 = sc.integrate();
                     if(enable_hard) xsection_HnonC2 = hc.integrate();
                     print( "ddbar->suLsuL*", xsection_tree2, xsection_virt2, xsection_SC2, xsection_HnonC2);
                     xsec_to_json(j, "ddbar->suLsuL*", xsection_tree2, xsection_virt2, xsection_SC2, xsection_HnonC2);
                  }

                  if( subprocess == "") {
                     std::vector<std::array<int, 3>> flav {{21, 21, 1}};
                     XSection_Tree tree(
                        m1, m2,
                        std::bind(&Process::matrixMRSSMTree_GG_suLsuLdagger, mrssm, _1, _2, _3), flav,
                        muR, muF,
                        pdf.get()
                     );
                     XSection_Virt virt(
                        m1, m2,
                        std::bind(&Process::matrixMRSSMVirt_GG_suLsuLdagger, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav,
                        muR, muF,
                        pdf.get()
                     );
                     XSection_SC sc(
                        m1, m2,
                        std::bind(&Process::matrixMRSSMSoft_gg_suLsuLdaggerg, mrssm, _1, _2, _3, _4, _5),
                        dS, dC, muR, muF,
                        flav, pdf.get(),
                        {
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pgg, std::bind(&Process::sigmaMRSSMTree_gg_suLsuLdagger, mrssm, _1, _2)},
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pgg, std::bind(&Process::sigmaMRSSMTree_gg_suLsuLdagger, mrssm, _1, _2)}
                        }

                     );
                     XSection_HnonC hc(
                        m1, m2,
                        std::bind(&Process::matrixMRSSMHard_gg_suLsuLdaggerg, mrssm, _1, _2), 
                        dS, dC, muR, muF,
                        flav, pdf.get()
                     );
	                  tree.init(pt, vm );
	                  virt.init(pt, vm );
	                  sc.init(pt, vm );
	                  hc.init(pt, vm );
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
                     std::vector<std::array<int, 3>> flav {};
                     for( int el : { 1, -1, 3, -3, 4, -4, 5, -5}) flav.push_back({21, el, 2});
                     XSection_SC sc(
                        m1, m2,
                        std::bind(&Process::matrix_soft_stub, mrssm, _1, _2, _3, _4, _5),
                        dS, dC, muR, muF,
                        flav
                        , pdf.get(),
                        {
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqg, std::bind(&Process::sigmaMRSSMTree_ddbar_suLsuLdagger, mrssm, _1, _2)},
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pgq, std::bind(&Process::sigmaMRSSMTree_gg_suLsuLdagger, mrssm, _1, _2)}
                        }

                     );
                     XSection_HnonC hc(
                        m1, m2,
                        std::bind(&Process::matrixMRSSMHard_gd_suLsuLdaggerd, mrssm, _1, _2), 
                        dS, dC, muR, muF,
                           flav
                        , pdf.get()
                           );
	                  sc.init(pt, vm );
	                  hc.init(pt, vm );
                     if(enable_sc) xsection_SC4 = sc.integrate();
                     if(enable_hard) xsection_HnonC4 = hc.integrate();
                     print( "gq > suLsuL*(+X)", xsection_tree4, xsection_virt4, xsection_SC4, xsection_HnonC4);
                     xsec_to_json(j, "gq->suLsuL*(+X)", xsection_tree4, xsection_virt4, xsection_SC4, xsection_HnonC4);
                  }

                  // g u > suL suLdagger
                  if( subprocess == "gu_suLsuLdaggeru" || subprocess == "" ) {
                     std::vector<std::array<int, 3>> flav {};
                     for( int el : { 2, -2 }) flav.push_back({21, el, 2});
                     XSection_SC sc(
                        m1, m2,
                        std::bind(&Process::matrix_soft_stub, mrssm, _1, _2, _3, _4, _5),
                        0. /*dS*/, dC, muR, muF,
                        flav
                        , pdf.get(),
                        {
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqg, std::bind(&Process::sigmaMRSSMTree_ddbar_suLsuLdagger, mrssm, _1, _2)},
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pgq, std::bind(&Process::sigmaMRSSMTree_gg_suLsuLdagger, mrssm, _1, _2)}
                        }

                     );
                     XSection_HnonC hc(
                        m1, m2,
                        std::bind(&Process::matrixMRSSMHard_gu_suLsuLdaggeru, mrssm, _1, _2), 
                        dS, dC, muR, muF,
                           flav
                        , pdf.get()
                           );
	                  sc.init(pt, vm );
	                  hc.init(pt, vm );
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
                            }
         case Model::Simplified:
            switch(channel) {
               case Channel::pp_OO:
               {
                  {
                  }
               }
            }
      }
   }

   // print out time statistics
   auto end = chrono::steady_clock::now();
   cout << "\nINFO: Calculation ended after ";
   if (end - start > 1h)
      cout << chrono::duration_cast<chrono::hours>(end-start).count() << " hour(s), ";
   if (end - start > 1min)
      cout << chrono::duration_cast<chrono::minutes>(end-start).count() %  60 << " minute(s) and ";
   cout << chrono::duration_cast<chrono::seconds>(end-start).count() % 60 << " second(s)\n";

   // write results to JSON file
   const string json_outputfile_name =
      vm.count("json-outputfile-name")
         ? vm["json-outputfile-name"].as<string>()
         : pt.get<string>("process.process") + "_"
           + to_string(pt.get<double>("masses.squarks")) + "_"
           + to_string(pt.get<double>("masses.gluino")) + "_"
           + to_string(pt.get<double>("masses.pseudoscalar_sgluon")) + "_"
           + to_string(pt.get<double>("collider setup.sqrt_S")) + "_"
           + to_string(pt.get<double>("collider setup.mu_r")) + "_"
           + to_string(pt.get<double>("collider setup.mu_f")) + "_"
           + pt.get<string>("collider setup.pdf")
           + ".json";
   std::ofstream o(json_outputfile_name);
   o << std::setw(3) << j << std::endl;

   return 0;
}
