#include "version.hpp"
#include "XSection_Tree.hpp"
#include "XSection_Virt.hpp"
#include "XSection_SC.hpp"
#include "XSection_HnonC.hpp"
#include "splitting_kernels.hpp"
#include "utils.hpp"

// models
#include "models/MSSM.hpp"
#include "models/MRSSM.hpp"
#include "models/Sgluons.hpp"

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
inline array<double, 3>& operator+=(array<double, 3>& x, array<double, 3> const& y) {
   x.at(0) += y.at(0);
   x.at(1) = std::hypot(x.at(1), y.at(1));
   x.at(2) = std::max(x.at(2), y.at(2));
   return x;
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
      ("precision-born", po::value<int>() -> default_value(6), "")
      ("precision-virt", po::value<int>() -> default_value(3), "")
      // gu_suLsuLdaggeru with SC precision 5 for BMP2 gives p-value 1
      ("precision-sc",   po::value<int>() -> default_value(6), "")
      ("precision-hard", po::value<int>() -> default_value(4), "")
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
   try {
      boost::property_tree::ini_parser::read_ini(card, pt);
   }
   catch (const std::exception& e) {
      std::cout << "Error while parsing trying to parse " << card << " (" << e.what() << ")\n";
      return 1;
   }

   // local arrays are not aumatically initialized to 0
   // need to use {}
   array<double,3> temp {}, xsection_tree {}, xsection_virt {}, xsection_SC {}, xsection_HnonC {},
           xsection_tree1 {}, xsection_virt1 {}, xsection_SC1 {}, xsection_HnonC1 {},
           xsection_tree2 {}, xsection_virt2 {}, xsection_SC2 {}, xsection_HnonC2 {},
           xsection_tree3 {}, xsection_virt3 {}, xsection_SC3 {}, xsection_HnonC3 {},
           xsection_tree4 {}, xsection_virt4 {}, xsection_SC4 {}, xsection_HnonC4 {},
           xsection_tree5 {}, xsection_virt5 {}, xsection_SC5 {}, xsection_HnonC5 {},
           xsection_tree6 {}, xsection_virt6 {}, xsection_SC6 {}, xsection_HnonC6 {},
           xsection_tree_total {}, xsection_virt_total {}, xsection_SC_total {}, xsection_HnonC_total {};

   enum class Order {LO, NLO/*, NLO+NLL*/};
   Order order;
   if (pt.get<string>("process.order") == "LO") {
      order = Order::LO;
   }
   else if (pt.get<string>("process.order") == "NLO") {
      order = Order::NLO;
   }
   else {
      std::cerr << "Error: Unknown order \"" << pt.get<string>("process.order") << "\". Allowed values are LO and NLO. Please check the input card.\n";
      return 1;
   }

   enum class Model {
       MRSSM,
       MSSM,
       Sgluons
   };

   MRSSMParameters mrssm_params;
   SgluonParameters sgluon_params;
   Model model;
   if (pt.get<string>("process.model") == "MRSSM") {
      model = Model::MRSSM;
      mrssm_params.MassTop = pt.get<double>("masses.top");
      mrssm_params.MassGlu = pt.get<double>("masses.gluino");
      mrssm_params.MasssigmaO = pt.get<double>("masses.pseudoscalar_sgluon");
      mrssm_params.MassSq = pt.get<double>("masses.squarks");
      mrssm_params.eta_sign = pt.get<int>("technical parameters.eta_sign", -1);
      mrssm_params.delta = pt.get<double>("technical parameters.delta", 0.);
      mrssm_params.WidthGlu = pt.get<double>("technical parameters.WidthOverMass", -1.) * mrssm_params.MassGlu;
   }
   else if (pt.get<string>("process.model") == "MSSM") {
      model = Model::MSSM;
   }
   else if (pt.get<string>("process.model") == "Sgluons") {
      sgluon_params.mO = pt.get<double>("masses.sgluons");
      model = Model::Sgluons;
   }
   else {
      std::cout << "\nModel not implemented!\n\n";
      return 1;
   }

   enum class Channel {
       pp_OO,
       pp_suLsuR,
       pp_suLsuL,
       pp_suLsdR,
       pp_sqLsqR,
       pp_sqLsqR_w_cc,
       pp_suLsuLdagger,
       pp_sqsqdagger,
       pp_glglbar
   };

   Channel channel;
   if (pt.get<string>("process.process") == "pp_suLsuR") {
      channel = Channel::pp_suLsuR;
   }
   else if (pt.get<string>("process.process") == "pp_sqLsqR") {
      channel = Channel::pp_sqLsqR;
   }
   else if (pt.get<string>("process.process") == "pp_sqLsqR+cc") {
      channel = Channel::pp_sqLsqR_w_cc;
   }
   else if (pt.get<string>("process.process") == "pp_suLsuLdagger") {
      channel = Channel::pp_suLsuLdagger;
   }
   else if (pt.get<string>("process.process") == "pp_sqsqdagger") {
      channel = Channel::pp_sqsqdagger;
   }
   else if (pt.get<string>("process.process") == "pp_OO") {
      channel = Channel::pp_OO;
   }
   else if (pt.get<string>("process.process") == "pp_glglbar") {
      channel = Channel::pp_glglbar;
   }
   else {
      cout << "\n Process not implemented! \n\n";
      return 1.;
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

   // set PDFs
   std::cout << '\n';
   std::unique_ptr<LHAPDF::PDF> pdf(LHAPDF::mkPDF(pt.get<std::string>("collider setup.pdf")));
   LHAPDF::setVerbosity(0);

   // parameters used by both LO and NLO calculations
   const double muR = pt.get<double>("collider setup.mu_r");
   const double muF = pt.get<double>("collider setup.mu_f");

   XSectionParameters parameters;
   parameters.sqrtS = pt.get<double>("collider setup.sqrt_S");
   parameters.muR =  muR;
   parameters.muF =  muF;
   parameters.pdf = pdf.get();

   const int born_verbosity = vm["verbosity-born"].as<int>();
   const int born_precision = vm["precision-born"].as<int>();
   const int virt_verbosity = vm["verbosity-virt"].as<int>();
   const int virt_precision = vm["precision-virt"].as<int>();
   const int sc_verbosity = vm["verbosity-sc"].as<int>();
   const int sc_precision = vm["precision-sc"].as<int>();
   const int hard_verbosity = vm["verbosity-hard"].as<int>();
   const int hard_precision = vm["precision-hard"].as<int>();

   auto start = chrono::steady_clock::now();

   if (pt.get<string>("process.order") == "LO") {
      switch(model) {
         case Model::MRSSM:
            MRSSM mrssm(mrssm_params);
            switch(channel) {
               case Channel::pp_suLsuR:
               {
                  const double m1 = pt.get<double>("masses.squarks");
                  const double m2 = pt.get<double>("masses.squarks");
                  {
                     std::vector<std::array<int, 3>> flav {{2,2,1}};
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMTree_uu_suLsuR, mrssm, _1, _2, _3),
                        flav,
                        born_precision, born_verbosity
                     );
                     xsection_tree = tree.integrate();
                     xsec_to_json(j, "uu->suLsuR", xsection_tree);
                     print_to_terminal("uu > suLsuR", xsection_tree);
                  }
                  break;
               }
               case Channel::pp_suLsuLdagger:
               {
                  break;
               }
               case Channel::pp_sqLsqR_w_cc:
               {
                  const double m1 = pt.get<double>("masses.squarks");
                  const double m2 = pt.get<double>("masses.squarks");
                  {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {1, 2, 3, 4, 5}) {
                        flav.push_back({ i,  i, 1});
                        flav.push_back({-i, -i, 1});
                        for (int j : {1, 2, 3, 4, 5}) {
                           if (j>=i) continue;
                           flav.push_back({ i,  j, 2});
                           flav.push_back({-i, -j, 2});
                        }
                     }
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMTree_uu_suLsuR, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     auto result = tree.integrate();
                     print_to_terminal("qq > sqLsqR + cc", result);
                     xsec_to_json(j, "qq->sqLsqR+cc", result);
                  }
                  break;
               }
               case Channel::pp_sqLsqR:
               {
                  const double m1 = pt.get<double>("masses.squarks");
                  const double m2 = pt.get<double>("masses.squarks");
                  {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {1, 2, 3, 4, 5}) {
                        flav.push_back({ i,  i, 1});
                        for (int j : {1, 2, 3, 4, 5}) {
                           if (j>=i) continue;
                           flav.push_back({ i,  j, 2});
                        }
                     }
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMTree_uu_suLsuR, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     auto result = tree.integrate();
                     print_to_terminal("qq > sqLsqR", result);
                     xsec_to_json(j, "qq->sqLsqR", result);
                  }
                  break;
               }
               case Channel::pp_sqsqdagger:
               {
                  const double m1 = pt.get<double>("masses.squarks");
                  const double m2 = pt.get<double>("masses.squarks");
                  std::array<double, 3> result {0., 0., 0.};
                  {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {1, 2, 3, 4, 5}) {
                        flav.push_back({i, -i, 4});
                     }
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMTree_uubar_suLsuLdagger, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     auto chan_res = tree.integrate();
                     print_to_terminal("qqbar > sqsq*", chan_res);
                     xsec_to_json(j, "qqbar->sqsq*", chan_res);
                     result = result + chan_res;
                  }
                  {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {1, 2, 3, 4, 5}) {
                        for (int j : {1, 2, 3, 4, 5}) {
                           if (j==i) continue;
                           flav.push_back({i, -j, 2});
                        }
                     }
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMTree_uu_suLsuR, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     auto chan_res = tree.integrate();
                     print_to_terminal("qqbar > sqsq*", chan_res);
                     xsec_to_json(j, "qqbar->sqsq*", chan_res);
                     result = result + chan_res;
                  }
                  {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {1, 2, 3, 4, 5}) {
                        // 4 (squark flavours) * 2 (pp symmetry) * 2 (L and R squarks)
                        flav.push_back({i,-i, 2*4*2});
                     }
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMTree_ddbar_suLsuLdagger, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     auto chan_res = tree.integrate();
                     print_to_terminal( "ddbar->suLsuL*", chan_res);
                     xsec_to_json(j, "ddbar->suLsuL*", chan_res);
                     result = result + chan_res;
                  }
                  {
                     // 5 squark flavours * L and R
                     std::vector<std::array<int, 3>> flav {{21, 21, 2*5}};
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMTree_GG_suLsuLdagger, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     auto chan_res = tree.integrate();
                     print_to_terminal( "gg > suLsuL*", chan_res);
                     xsec_to_json(j, "gg->suLsuL*", chan_res);
                     result = result + chan_res;
                  }
                  print_to_terminal("total", result);
                  break;
               }
               case Channel::pp_glglbar:
               {
                  const double m = pt.get<double>("masses.gluino");
                  std::array<double, 3> result {0., 0., 0.};
                  {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {1, 2, 3, 4, 5}) {
                        flav.push_back({i, -i, 2});
                     }
                     XSection_Tree tree(
                        parameters, m, m,
                        std::bind(&MRSSM::matrixMRSSMTree_uubar_glglbar, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     auto chan_res = tree.integrate();
                     print_to_terminal("qqbar > gluglubar", chan_res);
                     xsec_to_json(j, "qqbar->gluglubar", chan_res);
                     result = result + chan_res;
                  }
                  {
                     std::vector<std::array<int, 3>> flav {{21, 21, 1}};
                     XSection_Tree tree(
                        parameters, m, m,
                        std::bind(&MRSSM::matrixMRSSMTree_gg_glglbar, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     auto chan_res = tree.integrate();
                     print_to_terminal("gg > gluglubar", chan_res);
                     xsec_to_json(j, "gg->gluglubar", chan_res);
                     result = result + chan_res;
                  }
                  print_to_terminal("total", result);
               }
               default:
                  break;
            } // end of process block
         break;
      } // end of model block
   } // end of LO block
   else if (pt.get<string>("process.order") == "NLO") {
      const double dS = pt.get<double>("technical parameters.dS", 1e-5);
      // for the matrix elements that are regular in the limit dS -> 0 because the phase space parametrization
      // fails if we are exactly on the threshold
      constexpr double dS0 = 1e-10;
      const double dC = pt.get<double>("technical parameters.dC", 1e-6);
      cout << "\nINFO: Using phase space slicing parameters δS=" << dS
           << " and δC=" << dC << '\n';
      if (dC > dS) {
         cout << "Warning: δC should be always << than δS\n";
      }
      const double delta = pt.get<double>("technical parameters.delta", 0.);
      const int eta_sign = pt.get<double>("technical parameters.eta_sign", -1);
      switch(model) {
         case Model::MRSSM:
         {
            MRSSM mrssm(mrssm_params);
            switch(channel) {
               case Channel::pp_suLsuR:
               {
                  const double m1 = pt.get<double>("masses.squarks");
                  const double m2 = pt.get<double>("masses.squarks");
                  // uu > suL suR (+g) process
                  if( subprocess == "" ) {
                     const std::vector<std::array<int, 3>> flav {{2, 2, 1}};
                     XSection_Tree tree(
                        parameters, m1, m2,
                        // same ME as in the MSSM
                        std::bind(&MRSSM::matrixMRSSMTree_uu_suLsuR, mrssm, _1, _2, _3),
                        flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMVirt_uu_suLsuR, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMSoft_uu_suLsuRg, mrssm, _1, _2, _3, _4, _5),
                        dS, dC,
                        flav,
                        {
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_uu_suLsuR, mrssm, _1, _2)},
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_uu_suLsuR, mrssm, _1, _2)}
                        },
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_uu_suLsuRg, mrssm, _1, _2),
                        dS, dC,
                        flav,
                        hard_precision, hard_verbosity
                     );
                     if(enable_born) xsection_tree1 = tree.integrate();
                     if(enable_virt) xsection_virt1 = virt.integrate();
                     if(enable_sc) xsection_SC1 = sc.integrate();
                     if(enable_hard) xsection_HnonC1 = hc.integrate();
                     print_to_terminal( "uu > suLsuR(+X)", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
                     xsec_to_json(j, "uu->suLsuR(+X)", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
                  }

                  // gu > suL suR ubar process
                  if (subprocess == "gu_suLsuRubar" || subprocess == "" ) {
                     const std::vector<std::array<int, 3>> flav {{21, 2, 2}};
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrix_soft_stub, mrssm, _1, _2, _3, _4, _5),
                        dS0, dC,
                        flav,
                        {
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqg, std::bind(&MRSSM::sigmaMRSSMTree_uu_suLsuR, mrssm, _1, _2)},
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pgq, std::bind(&MRSSM::matrix_xsec_stub, mrssm, _1, _2)}
                        },
                        sc_precision, sc_verbosity
                     );
                     auto f =
                        mrssm_params.MassGlu < mrssm_params.MassSq
                           ? &MRSSM::matrixMRSSMHard_gu_suLsuRubar
                              : (mrssm_params.WidthGlu < 0 ? &MRSSM::matrixMRSSMHard_gu_suLsuRubar_DR : &MRSSM::matrixMRSSMHard_gu_suLsuRubar_DS);
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(f, mrssm, _1, _2),
                        dS0, dC,
                        flav,
                        hard_precision, hard_verbosity
                     );
                     if(enable_sc) xsection_SC2 = sc.integrate();
                     if(enable_hard) xsection_HnonC2 = hc.integrate();
                     print_to_terminal( "gu > suLsuR(+X)", xsection_tree2, xsection_virt2, xsection_SC2, xsection_HnonC2 );
                     xsec_to_json(j, "gu->suLsuR(+X)", xsection_tree2, xsection_virt2, xsection_SC2, xsection_HnonC2);
                  }

                  xsection_tree_total = xsection_tree1;
                  xsection_virt_total = xsection_virt1;
                  xsection_SC_total = xsection_SC1 + xsection_SC2;
                  xsection_HnonC_total = xsection_HnonC1 + xsection_HnonC2;
                  print_to_terminal( "sum", xsection_tree1, xsection_virt1, xsection_SC_total, xsection_HnonC_total );
                  break;
               }
               case Channel::pp_sqLsqR:
               {
                  // qq > sqL sqR (+g) process
                  const double m1 = pt.get<double>("masses.squarks");
                  const double m2 = pt.get<double>("masses.squarks");
                  if( subprocess == "" ) {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {1, 2, 3, 4, 5}) {
                        flav.push_back({i,  i, 1});
                        for (int j : {1, 2, 3, 4, 5}) {
                           if (j>=i) continue;
                           flav.push_back({i,  j, 2});
                        }
                     }
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMTree_uu_suLsuR, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMVirt_uu_suLsuR, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMSoft_uu_suLsuRg, mrssm, _1, _2, _3, _4, _5),
                        dS, dC,
                        flav,
                        {
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_uu_suLsuR, mrssm, _1, _2)},
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_uu_suLsuR, mrssm, _1, _2)}
                        },
                        sc_precision, sc_verbosity
                     );

                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_uu_suLsuRg, mrssm, _1, _2),
                        dS, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     if(enable_born) xsection_tree1 = tree.integrate();
                     if(enable_virt) xsection_virt1 = virt.integrate();
                     if(enable_sc) xsection_SC1 = sc.integrate();
                     if(enable_hard) xsection_HnonC1 = hc.integrate();
                     print_to_terminal( "qq > sqLsqR(+X)", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
                     xsec_to_json(j, "qq->sqLsqR(+X)", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
                  }
                  // gu > suL suR ubar process
                  if (subprocess == "gu_suLsuRubar" || subprocess == "" ) {
                     std::vector<std::array<int, 3>> flav {};
                     for(int el : {1, 2, 3, 4, 5}) flav.push_back({21, el, 5});
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrix_soft_stub, mrssm, _1, _2, _3, _4, _5),
                        dS0, dC,
                        flav,
                        {
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqg, std::bind(&MRSSM::sigmaMRSSMTree_uu_suLsuR, mrssm, _1, _2)},
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pgq, std::bind(&MRSSM::matrix_xsec_stub, mrssm, _1, _2)}
                        },
                        sc_precision, sc_verbosity
                     );
                     auto f =
                        mrssm_params.MassGlu < mrssm_params.MassSq
                           ? &MRSSM::matrixMRSSMHard_gu_suLsuRubar
                              : (mrssm_params.WidthGlu < 0 ? &MRSSM::matrixMRSSMHard_gu_suLsuRubar_DR : &MRSSM::matrixMRSSMHard_gu_suLsuRubar_DS);
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(f, mrssm, _1, _2),
                        dS0, dC,
                        flav,
                        hard_precision, hard_verbosity
                     );
                     if(enable_sc) xsection_SC2 = sc.integrate();
                     if(enable_hard) xsection_HnonC2 = hc.integrate();
                     print_to_terminal( "gu > suLsuR(+X)", xsection_tree2, xsection_virt2, xsection_SC2, xsection_HnonC2 );
                     xsec_to_json(j, "gu->suLsuR(+X)", xsection_tree2, xsection_virt2, xsection_SC2, xsection_HnonC2);
                  }

                  xsection_tree_total = xsection_tree1;
                  xsection_virt_total = xsection_virt1;
                  xsection_SC_total = xsection_SC1 + xsection_SC2;
                  xsection_HnonC_total = xsection_HnonC1 + xsection_HnonC2;
                  print_to_terminal( "sum", xsection_tree1, xsection_virt1, xsection_SC_total, xsection_HnonC_total );
                  break;
               }
               case Channel::pp_sqLsqR_w_cc:
               {
                  // qq > sqL sqR (+g) process
                  const double m1 = pt.get<double>("masses.squarks");
                  const double m2 = pt.get<double>("masses.squarks");
                  if( subprocess == "" ) {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {1, 2, 3, 4, 5}) {
                        flav.push_back({ i,  i, 1});
                        flav.push_back({-i, -i, 1});
                        for (int j : {1, 2, 3, 4, 5}) {
                           if (j>=i) continue;
                           flav.push_back({ i,  j, 2});
                           flav.push_back({-i, -j, 2});
                        }
                     }
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMTree_uu_suLsuR, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMVirt_uu_suLsuR, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMSoft_uu_suLsuRg, mrssm, _1, _2, _3, _4, _5),
                        dS, dC,
                        flav,
                        {
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_uu_suLsuR, mrssm, _1, _2)},
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_uu_suLsuR, mrssm, _1, _2)}
                        },
                        sc_precision, sc_verbosity
                     );

                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_uu_suLsuRg, mrssm, _1, _2),
                        dS, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     if(enable_born) xsection_tree1 = tree.integrate();
                     if(enable_virt) xsection_virt1 = virt.integrate();
                     if(enable_sc) xsection_SC1 = sc.integrate();
                     if(enable_hard) xsection_HnonC1 = hc.integrate();
                     print_to_terminal( "qq > sqLsqR(+X)", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
                     xsec_to_json(j, "qq->sqLsqR(+X)", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
                  }
                  // gu > suL suR ubar process
                  if (subprocess == "gu_suLsuRubar" || subprocess == "" ) {
                     std::vector<std::array<int, 3>> flav {};
                     for(int el : {1, -1, 2, -2, 3, -3, 4, -4, 5, -5}) flav.push_back({21, el, 2*5});
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrix_soft_stub, mrssm, _1, _2, _3, _4, _5),
                        dS0, dC,
                        flav,
                        {
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqg, std::bind(&MRSSM::sigmaMRSSMTree_uu_suLsuR, mrssm, _1, _2)},
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pgq, std::bind(&MRSSM::matrix_xsec_stub, mrssm, _1, _2)}
                        },
                        sc_precision, sc_verbosity
                     );
                     auto f =
                        mrssm_params.MassGlu < mrssm_params.MassSq
                           ? &MRSSM::matrixMRSSMHard_gu_suLsuRubar
                              : (mrssm_params.WidthGlu < 0 ? &MRSSM::matrixMRSSMHard_gu_suLsuRubar_DR : &MRSSM::matrixMRSSMHard_gu_suLsuRubar_DS);
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(f, mrssm, _1, _2),
                        dS0, dC,
                        flav,
                        hard_precision, hard_verbosity
                     );
                     if(enable_sc) xsection_SC2 = sc.integrate();
                     if(enable_hard) xsection_HnonC2 = hc.integrate();
                     print_to_terminal( "gu > suLsuR(+X)", xsection_tree2, xsection_virt2, xsection_SC2, xsection_HnonC2 );
                     xsec_to_json(j, "gu->suLsuR(+X)", xsection_tree2, xsection_virt2, xsection_SC2, xsection_HnonC2);
                  }

                  xsection_tree_total = xsection_tree1;
                  xsection_virt_total = xsection_virt1;
                  xsection_SC_total = xsection_SC1 + xsection_SC2;
                  xsection_HnonC_total = xsection_HnonC1 + xsection_HnonC2;
                  print_to_terminal( "sum", xsection_tree1, xsection_virt1, xsection_SC_total, xsection_HnonC_total );
                  break;
               }
               case Channel::pp_suLsuLdagger:
               {
                  const double m1 = pt.get<double>("masses.squarks");
                  const double m2 = pt.get<double>("masses.squarks");
                  if (subprocess == "" || subprocess == "uubar_suLsuLdagger") {
                     std::vector<std::array<int, 3>> flav {{2, -2, 2}};
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMTree_uubar_suLsuLdagger, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMVirt_uubar_suLsuLdagger, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMSoft_uubar_suLsuLdaggerg, mrssm, _1, _2, _3, _4, _5),
                        dS, dC,
                        flav,
                        {
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_uubar_suLsuLdagger, mrssm, _1, _2)},
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_uubar_suLsuLdagger, mrssm, _1, _2)}
                        },
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_uubar_suLsuLdaggerg, mrssm, _1, _2),
                        dS, dC,
                        flav,
                        hard_precision, hard_verbosity
                     );
                     if(enable_born) xsection_tree1 = tree.integrate();
                     if(enable_virt) xsection_virt1 = virt.integrate();
                     if(enable_sc) xsection_SC1 = sc.integrate();
                     if(enable_hard) xsection_HnonC1 = hc.integrate();
                     print_to_terminal( "uubar > suLsuL*", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
                     xsec_to_json(j, "uubar->suLsuL*", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
                  }

                  if( subprocess == "") {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {1, 3, 4, 5}) {
                        flav.push_back({i,-i, 2});
                     }
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMTree_ddbar_suLsuLdagger, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMVirt_ddbar_suLsuLdagger, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMSoft_ddbar_suLsuLdaggerg, mrssm, _1, _2, _3, _4, _5),
                        dS, dC,
                        flav,
                        {
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_ddbar_suLsuLdagger, mrssm, _1, _2)},
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_ddbar_suLsuLdagger, mrssm, _1, _2)}
                        },
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_ddbar_suLsuLdaggerg, mrssm, _1, _2),
                        dS, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     if(enable_born) xsection_tree2 = tree.integrate();
                     if(enable_virt) xsection_virt2 = virt.integrate();
                     if(enable_sc) xsection_SC2 = sc.integrate();
                     if(enable_hard) xsection_HnonC2 = hc.integrate();
                     print_to_terminal( "ddbar->suLsuL*", xsection_tree2, xsection_virt2, xsection_SC2, xsection_HnonC2);
                     xsec_to_json(j, "ddbar->suLsuL*", xsection_tree2, xsection_virt2, xsection_SC2, xsection_HnonC2);
                  }

                  if( subprocess == "") {
                     std::vector<std::array<int, 3>> flav {{21, 21, 1}};
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMTree_GG_suLsuLdagger, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMVirt_GG_suLsuLdagger, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMSoft_gg_suLsuLdaggerg, mrssm, _1, _2, _3, _4, _5),
                        dS, dC,
                        flav,
                        {
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pgg, std::bind(&MRSSM::sigmaMRSSMTree_gg_suLsuLdagger, mrssm, _1, _2)},
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pgg, std::bind(&MRSSM::sigmaMRSSMTree_gg_suLsuLdagger, mrssm, _1, _2)}
                        },
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_gg_suLsuLdaggerg, mrssm, _1, _2),
                        dS, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     if(enable_born) xsection_tree3 = tree.integrate();
                     if(enable_virt) xsection_virt3 = virt.integrate();
                     if(enable_sc) xsection_SC3 = sc.integrate();
                     if(enable_hard) xsection_HnonC3 = hc.integrate();
                     print_to_terminal( "gg > suLsuL*", xsection_tree3, xsection_virt3, xsection_SC3, xsection_HnonC3);
                     xsec_to_json(j, "gg->suLsuL*", xsection_tree3, xsection_virt3, xsection_SC3, xsection_HnonC3);
                  }

                  if( subprocess == "") {
                     std::vector<std::array<int, 3>> flav {};
                     for( int el : { 1, -1, 3, -3, 4, -4, 5, -5}) flav.push_back({21, el, 2});
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrix_soft_stub, mrssm, _1, _2, _3, _4, _5),
                        dS0, dC, flav,
                        {
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqg, std::bind(&MRSSM::sigmaMRSSMTree_ddbar_suLsuLdagger, mrssm, _1, _2)},
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pgq, std::bind(&MRSSM::sigmaMRSSMTree_gg_suLsuLdagger, mrssm, _1, _2)}
                        },
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_gd_suLsuLdaggerd, mrssm, _1, _2),
                        dS0, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     if(enable_sc) xsection_SC4 = sc.integrate();
                     if(enable_hard) xsection_HnonC4 = hc.integrate();
                     print_to_terminal( "gq > suLsuL*(+X)", xsection_tree4, xsection_virt4, xsection_SC4, xsection_HnonC4);
                     xsec_to_json(j, "gq->suLsuL*(+X)", xsection_tree4, xsection_virt4, xsection_SC4, xsection_HnonC4);
                  }

                  // g u > suL suLdagger
                  if( subprocess == "gu_suLsuLdaggeru" || subprocess == "" ) {
                     std::vector<std::array<int, 3>> flav {};
                     for( int el : { 2, -2 }) flav.push_back({21, el, 2});
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrix_soft_stub, mrssm, _1, _2, _3, _4, _5),
                        dS0, dC, flav,
                        {
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqg, std::bind(&MRSSM::sigmaMRSSMTree_uubar_suLsuLdagger, mrssm, _1, _2)},
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pgq, std::bind(&MRSSM::sigmaMRSSMTree_gg_suLsuLdagger, mrssm, _1, _2)}
                        },
                        sc_precision, sc_verbosity
                     );
                     auto f =
                        mrssm_params.MassGlu < mrssm_params.MassSq
                           ? &MRSSM::matrixMRSSMHard_gu_suLsuLdaggeru
                              : (mrssm_params.WidthGlu < 0 ? &MRSSM::matrixMRSSMHard_gu_suLsuLdaggeru_DR : &MRSSM::matrixMRSSMHard_gu_suLsuLdaggeru_DS);
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(f, mrssm, _1, _2),
                        dS0, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     if(enable_sc) xsection_SC5 = sc.integrate();
                     if(enable_hard) xsection_HnonC5 = hc.integrate();
                     print_to_terminal( "gu > suLsuL*(+X)", xsection_tree5, xsection_virt5, xsection_SC5, xsection_HnonC5 );
                     xsec_to_json(j, "gu->suLsuL*(+X)", xsection_tree5, xsection_virt5, xsection_SC5, xsection_HnonC5);
                  }

                  xsection_tree_total = xsection_tree1 + xsection_tree2 + xsection_tree3;
                  xsection_virt_total = xsection_virt1 + xsection_virt2 + xsection_virt3;
                  xsection_SC_total = xsection_SC1 + xsection_SC2 + xsection_SC3 + xsection_SC4 + xsection_SC5;
                  xsection_HnonC_total = xsection_HnonC1 + xsection_HnonC2 + xsection_HnonC3
                          + xsection_HnonC4 + xsection_HnonC5;
                  print_to_terminal( "total", xsection_tree_total, xsection_virt_total, xsection_SC_total, xsection_HnonC_total);
                  break;
               }
               case Channel::pp_sqsqdagger:
               {
                  const double m1 = pt.get<double>("masses.squarks");
                  const double m2 = pt.get<double>("masses.squarks");
                  if (subprocess == "" || subprocess == "uubar_suLsuLdagger") {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {1, 2, 3, 4, 5}) {
                        flav.push_back({i, -i, 4});
                     }
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMTree_uubar_suLsuLdagger, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMVirt_uubar_suLsuLdagger, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMSoft_uubar_suLsuLdaggerg, mrssm, _1, _2, _3, _4, _5),
                        dS, dC,
                        flav,
                        {
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_uubar_suLsuLdagger, mrssm, _1, _2)},
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_uubar_suLsuLdagger, mrssm, _1, _2)}
                        },
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_uubar_suLsuLdaggerg, mrssm, _1, _2),
                        dS, dC,
                        flav,
                        hard_precision, hard_verbosity
                     );
                     if(enable_born) xsection_tree1 = tree.integrate();
                     if(enable_virt) xsection_virt1 = virt.integrate();
                     if(enable_sc) xsection_SC1 = sc.integrate();
                     if(enable_hard) xsection_HnonC1 = hc.integrate();
                     print_to_terminal( "qqbar > sqsq*", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
                     xsec_to_json(j, "qqbar->sqsq*", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
                  }
                  if (subprocess == "" || subprocess == "uubar_suLsuLdagger") {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {1, 2, 3, 4, 5}) {
                        for (int j : {1, 2, 3, 4, 5}) {
                           if (j==i) continue;
                           flav.push_back({i,-j,2});
                        }
                     }
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMTree_uu_suLsuR, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMVirt_uubar_suLsuLdagger, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMSoft_uubar_suLsuLdaggerg, mrssm, _1, _2, _3, _4, _5),
                        dS, dC,
                        flav,
                        {
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_uubar_suLsuLdagger, mrssm, _1, _2)},
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_uubar_suLsuLdagger, mrssm, _1, _2)}
                        },
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_uubar_suLsuLdaggerg, mrssm, _1, _2),
                        dS, dC,
                        flav,
                        hard_precision, hard_verbosity
                     );
                     if(enable_born) xsection_tree2 = tree.integrate();
                     if(enable_virt) xsection_virt2 = virt.integrate();
                     if(enable_sc) xsection_SC2 = sc.integrate();
                     if(enable_hard) xsection_HnonC2 = hc.integrate();
                     print_to_terminal( "qqbar > sqsq*", xsection_tree2, xsection_virt2, xsection_SC2, xsection_HnonC2);
                     xsec_to_json(j, "qqbar->sqsq*", xsection_tree2, xsection_virt2, xsection_SC2, xsection_HnonC2);
                  }

                  if( subprocess == "") {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {1, 2, 3, 4, 5}) {
                        // 4 (squark flavours) * 2 (pp symmetry) * 2 (L and R squarks)
                        flav.push_back({i,-i, 2*4*2});
                     }
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMTree_ddbar_suLsuLdagger, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMVirt_ddbar_suLsuLdagger, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMSoft_ddbar_suLsuLdaggerg, mrssm, _1, _2, _3, _4, _5),
                        dS, dC,
                        flav,
                        {
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_ddbar_suLsuLdagger, mrssm, _1, _2)},
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_ddbar_suLsuLdagger, mrssm, _1, _2)}
                        },
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_ddbar_suLsuLdaggerg, mrssm, _1, _2),
                        dS, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     if(enable_born) xsection_tree3 = tree.integrate();
                     if(enable_virt) xsection_virt3 = virt.integrate();
                     if(enable_sc) xsection_SC3 = sc.integrate();
                     if(enable_hard) xsection_HnonC3 = hc.integrate();
                     print_to_terminal( "ddbar->suLsuL*", xsection_tree3, xsection_virt3, xsection_SC3, xsection_HnonC3);
                     xsec_to_json(j, "ddbar->suLsuL*", xsection_tree3, xsection_virt3, xsection_SC3, xsection_HnonC3);
                  }

                  if( subprocess == "") {
                     // 5 squark flavours * L and R
                     std::vector<std::array<int, 3>> flav {{21, 21, 2*5}};
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMTree_GG_suLsuLdagger, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMVirt_GG_suLsuLdagger, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMSoft_gg_suLsuLdaggerg, mrssm, _1, _2, _3, _4, _5),
                        dS, dC,
                        flav,
                        {
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pgg, std::bind(&MRSSM::sigmaMRSSMTree_gg_suLsuLdagger, mrssm, _1, _2)},
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pgg, std::bind(&MRSSM::sigmaMRSSMTree_gg_suLsuLdagger, mrssm, _1, _2)}
                        },
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_gg_suLsuLdaggerg, mrssm, _1, _2),
                        dS, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     if(enable_born) xsection_tree4 = tree.integrate();
                     if(enable_virt) xsection_virt4 = virt.integrate();
                     if(enable_sc) xsection_SC4 = sc.integrate();
                     if(enable_hard) xsection_HnonC4 = hc.integrate();
                     print_to_terminal( "gg > suLsuL*", xsection_tree4, xsection_virt4, xsection_SC4, xsection_HnonC4);
                     xsec_to_json(j, "gg->suLsuL*", xsection_tree4, xsection_virt4, xsection_SC4, xsection_HnonC4);
                  }

                  if( subprocess == "") {
                     std::vector<std::array<int, 3>> flav {};
                     for( int el : { 1, -1, 3, -3, 4, -4, 5, -5}) flav.push_back({21, el, 2});
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrix_soft_stub, mrssm, _1, _2, _3, _4, _5),
                        dS0, dC, flav,
                        {
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqg, std::bind(&MRSSM::sigmaMRSSMTree_ddbar_suLsuLdagger, mrssm, _1, _2)},
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pgq, std::bind(&MRSSM::sigmaMRSSMTree_gg_suLsuLdagger, mrssm, _1, _2)}
                        },
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_gd_suLsuLdaggerd, mrssm, _1, _2),
                        dS0, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     if(enable_sc) xsection_SC5 = sc.integrate();
                     if(enable_hard) xsection_HnonC5 = hc.integrate();
                     print_to_terminal( "gq > suLsuL*(+X)", xsection_tree5, xsection_virt5, xsection_SC5, xsection_HnonC5);
                     xsec_to_json(j, "gq->suLsuL*(+X)", xsection_tree5, xsection_virt5, xsection_SC5, xsection_HnonC5);
                  }

                  // g u > suL suLdagger
                  if( subprocess == "gu_suLsuLdaggeru" || subprocess == "" ) {
                     std::vector<std::array<int, 3>> flav {};
                     for( int el : { 2, -2 }) flav.push_back({21, el, 2});
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrix_soft_stub, mrssm, _1, _2, _3, _4, _5),
                        dS0, dC, flav,
                        {
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqg, std::bind(&MRSSM::sigmaMRSSMTree_ddbar_suLsuLdagger, mrssm, _1, _2)},
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pgq, std::bind(&MRSSM::sigmaMRSSMTree_gg_suLsuLdagger, mrssm, _1, _2)}
                        },
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_gu_suLsuLdaggeru, mrssm, _1, _2),
                        dS0, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     if(enable_sc) xsection_SC6 = sc.integrate();
                     if(enable_hard) xsection_HnonC6 = hc.integrate();
                     print_to_terminal( "gu > suLsuL*(+X)", xsection_tree6, xsection_virt6, xsection_SC6, xsection_HnonC6 );
                     xsec_to_json(j, "gu->suLsuL*(+X)", xsection_tree6, xsection_virt6, xsection_SC6, xsection_HnonC6);
                  }

                  xsection_tree_total = xsection_tree1 + xsection_tree2 + xsection_tree3 + xsection_tree4;
                  xsection_virt_total = xsection_virt1 + xsection_virt2 + xsection_virt3;
                  xsection_SC_total = xsection_SC1 + xsection_SC2 + xsection_SC3 + xsection_SC4 + xsection_SC5;
                  xsection_HnonC_total = xsection_HnonC1 + xsection_HnonC2 + xsection_HnonC3
                          + xsection_HnonC4 + xsection_HnonC5;
                  print_to_terminal( "total", xsection_tree_total, xsection_virt_total, xsection_SC_total, xsection_HnonC_total);
                  break;
               }
               default:
                  cout << "NLO process not implemented\n";
            }
         }
         break;
         case Model::Sgluons:
            Sgluons sgluons(sgluon_params);
            switch(channel) {
               case Channel::pp_OO:
               {
                  const double m1 = pt.get<double>("masses.sgluon");
                  if (subprocess == "") {
                     std::vector<std::array<int, 3>> flav {{2, -2, 2}};
                     XSection_Tree tree(
                        parameters, m1, m1,
                        std::bind(&Sgluons::matrixSgluonsTree_qqbar_OO, sgluons, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     /*
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMVirt_uubar_suLsuLdagger, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMSoft_uubar_suLsuLdaggerg, mrssm, _1, _2, _3, _4, _5),
                        dS, dC,
                        flav,
                        {
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_uubar_suLsuLdagger, mrssm, _1, _2)},
                           std::pair<SplittingKernel, std::function<double(double, double)>>{SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_uubar_suLsuLdagger, mrssm, _1, _2)}
                        }

                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_uubar_suLsuLdaggerg, mrssm, _1, _2),
                        dS, dC,
                        flav
                     );
                     if(enable_born) xsection_tree1 = tree.integrate();
                     if(enable_virt) xsection_virt1 = virt.integrate();
                     if(enable_sc) xsection_SC1 = sc.integrate();
                     if(enable_hard) xsection_HnonC1 = hc.integrate();
                     print_to_terminal( "uubar > suLsuL*", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
                     xsec_to_json(j, "uubar->suLsuL*", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
                     */
               }
            }
         }
      }
      j["technical parameters"] = {
         {"dS", dS},
         {"dC", dC},
         {"WidthOverMass", pt.get<double>("technical parameters.WidthOverMass")},
         {"eta_sign", eta_sign},
         {"delta", delta}
      };
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
