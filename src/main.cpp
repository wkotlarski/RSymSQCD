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
using namespace std::placeholders;

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
      ("precision-virt", po::value<int>() -> default_value(4), "")
      // gu_suLsuLdaggeru with SC precision 5 for BMP2 gives p-value 1
      ("precision-sc",   po::value<int>() -> default_value(7), "")
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
      std::cerr << "Error while trying to parse " << card << " (" << e.what() << ")\n";
      return 1;
   }

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
      sgluon_params.mt = pt.get<double>("masses.top");
      model = Model::Sgluons;
   }
   else {
      std::cerr << "\nError: Model " << pt.get<string>("process.model") << " not implemented\n\n";
      return 1;
   }

   enum class Channel {
       pp_OO,
       pp_suLsuR,
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

   nlohmann::json j;
   j["process"] = pt.get<string>("process.process");
   j["sqrt(S)"] = pt.get<double>("collider setup.sqrt_S");
   j["mu_r"] = pt.get<double>("collider setup.mu_r");
   j["mu_f"] = pt.get<double>("collider setup.mu_f");
   j["pdf"] = pt.get<string>("collider setup.pdf");
   switch (model) {
      case Model::MRSSM:
         j["masses"] = {
            {"gluino", pt.get<double>("masses.gluino")},
            {"pseudoscalar sgluon", pt.get<double>("masses.pseudoscalar_sgluon")},
            {"top", pt.get<double>("masses.top")},
            {"squarks", pt.get<double>("masses.squarks")}
         };
         break;
   }

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

   switch (order) {
      case Order::LO:
      {
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
                     auto xsection_current = tree.integrate();
                     xsec_to_json(j, "uu->suLsuR", xsection_current);
                     print_to_terminal("uu > suLsuR", xsection_current);
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
                  std::array<double, 3> result {};
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
                     print_to_terminal("qqbar -> sqsq*", chan_res);
                     xsec_to_json(j, "qqbar->sqsq*", chan_res);
                     result += chan_res;
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
                     print_to_terminal("qq'bar -> sqsq'*", chan_res);
                     xsec_to_json(j, "qq'bar->sqsq'*", chan_res);
                     result += chan_res;
                  }
                  {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {1, 2, 3, 4, 5}) {
                        // 4 (squark flavours) * 2 (pp symmetry) * 2 (L and R squarks)
                        flav.push_back({i, -i, 2*4*2});
                     }
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMTree_ddbar_suLsuLdagger, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     auto chan_res = tree.integrate();
                     print_to_terminal("qqbar -> sq'sq'*", chan_res);
                     xsec_to_json(j, "ddbar->sq'sq'*", chan_res);
                     result += chan_res;
                  }
                  {
                     // 5 squark flavours * (L + R)
                     std::vector<std::array<int, 3>> flav {{21, 21, 2*5}};
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMTree_GG_suLsuLdagger, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     auto chan_res = tree.integrate();
                     print_to_terminal( "gg -> suLsuL*", chan_res);
                     xsec_to_json(j, "gg->sqsq*", chan_res);
                     result += chan_res;
                  }
                  print_to_terminal("total", result);
                  break;
               }
               case Channel::pp_glglbar:
               {
                  const double m = pt.get<double>("masses.gluino");
                  std::array<double, 3> result {};
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
                     result += chan_res;
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
                     result += chan_res;
                  }
                  print_to_terminal("total", result);
               }
               default:
                  break;
            } // end of process block
         break;
      } // end of model block
      break;
      } // end of LO block
      case Order::NLO: {
      const double dS = pt.get<double>("technical parameters.dS", 1e-5);
      // for the matrix elements that are regular in the limit dS -> 0 because the phase space parametrization
      // fails if we are exactly on the threshold
      constexpr double dS0 = 1e-10;
      const double dC = pt.get<double>("technical parameters.dC", 1e-6);
      cout << "\nINFO: Using phase space slicing parameters δS=" << std::scientific << std::setprecision(1) << dS
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
                  array<double, 3> born_xsec_total {};
                  array<double, 3> soft_xsec_total {};
                  array<double, 3> virt_xsec_total {};
                  array<double, 3> hard_xsec_total {};
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
                        {{
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_uu_suLsuR, mrssm, _1, _2)},
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_uu_suLsuR, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_uu_suLsuRg, mrssm, _1, _2),
                        dS, dC,
                        flav,
                        hard_precision, hard_verbosity
                     );
                     array<double, 3> born_xsec_current {};
                     array<double, 3> soft_xsec_current {};
                     array<double, 3> virt_xsec_current {};
                     array<double, 3> hard_xsec_current {};
                     if(enable_born) born_xsec_current = tree.integrate();
                     if(enable_virt) virt_xsec_current = virt.integrate();
                     if(enable_sc) soft_xsec_current = sc.integrate();
                     if(enable_hard) hard_xsec_current = hc.integrate();
                     print_to_terminal("uu > suLsuR(+X)", born_xsec_current, virt_xsec_current, soft_xsec_current, hard_xsec_current);
                     xsec_to_json(j, "uu->suLsuR(+X)", born_xsec_current, virt_xsec_current, soft_xsec_current, hard_xsec_current);
                     born_xsec_total += born_xsec_current;
                     virt_xsec_total += virt_xsec_current;
                     soft_xsec_total += soft_xsec_current;
                     hard_xsec_total += hard_xsec_current;
                  }

                  // gu > suL suR ubar process
                  if (subprocess == "gu_suLsuRubar" || subprocess == "" ) {
                     const std::vector<std::array<int, 3>> flav {{21, 2, 2}};
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::nullopt,
                        dS0, dC,
                        flav,
                        {{
                           {SplittingKernel::Pqg, std::bind(&MRSSM::sigmaMRSSMTree_uu_suLsuR, mrssm, _1, _2)},
                           {SplittingKernel::Pgq, std::nullopt}
                        }},
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
                     array<double, 3> soft_xsec_current {};
                     array<double, 3> hard_xsec_current {};
                     if(enable_sc) soft_xsec_current = sc.integrate();
                     if(enable_hard) hard_xsec_current = hc.integrate();
                     print_to_terminal("gu > suLsuR(+X)", {}, {}, soft_xsec_current, hard_xsec_current);
                     xsec_to_json(j, "gu->suLsuR(+X)", {}, {}, soft_xsec_current, hard_xsec_current);
                     soft_xsec_total += soft_xsec_current;
                     hard_xsec_total += hard_xsec_current;
                  }

                  print_to_terminal("sum", born_xsec_total, virt_xsec_total, soft_xsec_total, hard_xsec_total);
                  break;
               }
               case Channel::pp_sqLsqR:
               {
                  // qq > sqL sqR (+g) process
                  const double m1 = pt.get<double>("masses.squarks");
                  const double m2 = pt.get<double>("masses.squarks");
                  std::array<double, 3> total_xsec_tree {};
                  std::array<double, 3> total_xsec_virt {};
                  std::array<double, 3> total_xsec_soft {};
                  std::array<double, 3> total_xsec_hard {};
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
                        {{
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_uu_suLsuR, mrssm, _1, _2)},
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_uu_suLsuR, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );

                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_uu_suLsuRg, mrssm, _1, _2),
                        dS, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     std::array<double, 3> current_tree {};
                     std::array<double, 3> current_virt {};
                     std::array<double, 3> current_soft {};
                     std::array<double, 3> current_hard {};
                     if(enable_born) current_tree = tree.integrate();
                     if(enable_virt) current_virt = virt.integrate();
                     if(enable_sc) current_soft = sc.integrate();
                     if(enable_hard) current_hard = hc.integrate();
                     print_to_terminal("qq > sqLsqR(+X)", current_tree, current_virt, current_soft, current_hard);
                     xsec_to_json(j, "qq->sqLsqR(+X)", current_tree, current_virt, current_soft, current_hard);
                     total_xsec_tree += current_tree;
                     total_xsec_virt += current_virt;
                     total_xsec_soft += current_soft;
                     total_xsec_hard += current_hard;
                  }
                  // gu > suL suR ubar process
                  if (subprocess == "gu_suLsuRubar" || subprocess == "" ) {
                     std::vector<std::array<int, 3>> flav {};
                     for(int el : {1, 2, 3, 4, 5}) flav.push_back({21, el, 5});
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::nullopt,
                        dS0, dC,
                        flav,
                        {{
                           {SplittingKernel::Pqg, std::bind(&MRSSM::sigmaMRSSMTree_uu_suLsuR, mrssm, _1, _2)},
                           {SplittingKernel::Pgq, std::nullopt}
                        }},
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
                     std::array<double, 3> current_soft {};
                     std::array<double, 3> current_hard {};
                     if(enable_sc) current_soft = sc.integrate();
                     if(enable_hard) current_hard = hc.integrate();
                     print_to_terminal( "gu > suLsuR(+X)", {}, {}, current_soft, current_hard);
                     xsec_to_json(j, "gu->suLsuR(+X)", {}, {}, current_soft, current_hard);
                     total_xsec_soft += current_soft;
                     total_xsec_hard += current_hard;
                  }

                  print_to_terminal("sum", total_xsec_tree, total_xsec_virt, total_xsec_soft, total_xsec_hard);
                  break;
               }
               case Channel::pp_sqLsqR_w_cc:
               {
                  // qq > sqL sqR (+g) process
                  const double m1 = pt.get<double>("masses.squarks");
                  const double m2 = pt.get<double>("masses.squarks");
                  std::array<double, 3> total_xsec_tree {};
                  std::array<double, 3> total_xsec_virt {};
                  std::array<double, 3> total_xsec_soft {};
                  std::array<double, 3> total_xsec_hard {};
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
                        {{
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_uu_suLsuR, mrssm, _1, _2)},
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_uu_suLsuR, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );

                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_uu_suLsuRg, mrssm, _1, _2),
                        dS, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     std::array<double, 3> current_tree {};
                     std::array<double, 3> current_virt {};
                     std::array<double, 3> current_soft {};
                     std::array<double, 3> current_hard {};
                     if(enable_born) current_tree = tree.integrate();
                     if(enable_virt) current_virt = virt.integrate();
                     if(enable_sc) current_soft = sc.integrate();
                     if(enable_hard) current_hard = hc.integrate();
                     print_to_terminal( "qq > sqLsqR(+X)", current_tree, current_virt, current_soft, current_hard);
                     xsec_to_json(j, "qq->sqLsqR(+X)", current_tree, current_virt, current_soft, current_hard);
                     total_xsec_tree += current_tree;
                     total_xsec_virt += current_virt;
                     total_xsec_soft += current_soft;
                     total_xsec_hard += current_hard;
                  }
                  // gu > suL suR ubar process
                  if (subprocess == "gu_suLsuRubar" || subprocess == "" ) {
                     std::vector<std::array<int, 3>> flav {};
                     for(int el : {1, -1, 2, -2, 3, -3, 4, -4, 5, -5}) flav.push_back({21, el, 2*5});
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::nullopt,
                        dS0, dC,
                        flav,
                        {{
                           {SplittingKernel::Pqg, std::bind(&MRSSM::sigmaMRSSMTree_uu_suLsuR, mrssm, _1, _2)},
                           {SplittingKernel::Pgq, std::nullopt}
                        }},
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
                     std::array<double, 3> current_soft {};
                     std::array<double, 3> current_hard {};
                     if(enable_sc) current_soft = sc.integrate();
                     if(enable_hard) current_hard = hc.integrate();
                     print_to_terminal( "gu > suLsuR(+X)", {}, {}, current_soft, current_hard);
                     xsec_to_json(j, "gu->suLsuR(+X)", {}, {}, current_soft, current_hard);
                     total_xsec_soft += current_soft;
                     total_xsec_hard += current_hard;
                  }

                  print_to_terminal("sum", total_xsec_tree, total_xsec_virt, total_xsec_soft, total_xsec_hard);
                  break;
               }
               case Channel::pp_suLsuLdagger:
               {
                  const double m1 = pt.get<double>("masses.squarks");
                  const double m2 = pt.get<double>("masses.squarks");
                  std::array<double, 3> total_xsec_tree {};
                  std::array<double, 3> total_xsec_virt {};
                  std::array<double, 3> total_xsec_soft {};
                  std::array<double, 3> total_xsec_hard {};
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
                        {{
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_uubar_suLsuLdagger, mrssm, _1, _2)},
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_uubar_suLsuLdagger, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_uubar_suLsuLdaggerg, mrssm, _1, _2),
                        dS, dC,
                        flav,
                        hard_precision, hard_verbosity
                     );
                     std::array<double, 3> current_tree {};
                     std::array<double, 3> current_virt {};
                     std::array<double, 3> current_soft {};
                     std::array<double, 3> current_hard {};
                     if(enable_born) current_tree = tree.integrate();
                     if(enable_virt) current_virt = virt.integrate();
                     if(enable_sc) current_soft = sc.integrate();
                     if(enable_hard) current_hard = hc.integrate();
                     print_to_terminal( "uubar -> suLsuL*(+X)", current_tree, current_virt, current_soft, current_hard);
                     xsec_to_json(j, "uubar->suLsuL*", current_tree, current_virt, current_soft, current_hard);
                     total_xsec_tree += current_tree;
                     total_xsec_virt += current_virt;
                     total_xsec_soft += current_soft;
                     total_xsec_hard += current_hard;
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
                        {{
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_ddbar_suLsuLdagger, mrssm, _1, _2)},
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_ddbar_suLsuLdagger, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_ddbar_suLsuLdaggerg, mrssm, _1, _2),
                        dS, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     std::array<double, 3> current_tree {};
                     std::array<double, 3> current_virt {};
                     std::array<double, 3> current_soft {};
                     std::array<double, 3> current_hard {};
                     if(enable_born) current_tree = tree.integrate();
                     if(enable_virt) current_virt = virt.integrate();
                     if(enable_sc) current_soft = sc.integrate();
                     if(enable_hard) current_hard = hc.integrate();
                     print_to_terminal("ddbar -> suLsuL*(+X)", current_tree, current_virt, current_soft, current_hard);
                     xsec_to_json(j, "ddbar->suLsuL*", current_tree, current_virt, current_soft, current_hard);
                     total_xsec_tree += current_tree;
                     total_xsec_virt += current_virt;
                     total_xsec_soft += current_soft;
                     total_xsec_hard += current_hard;
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
                        {{
                           {SplittingKernel::Pgg, std::bind(&MRSSM::sigmaMRSSMTree_gg_suLsuLdagger, mrssm, _1, _2)},
                           {SplittingKernel::Pgg, std::bind(&MRSSM::sigmaMRSSMTree_gg_suLsuLdagger, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_gg_suLsuLdaggerg, mrssm, _1, _2),
                        dS, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     std::array<double, 3> current_tree {};
                     std::array<double, 3> current_virt {};
                     std::array<double, 3> current_soft {};
                     std::array<double, 3> current_hard {};
                     if(enable_born) current_tree = tree.integrate();
                     if(enable_virt) current_virt = virt.integrate();
                     if(enable_sc) current_soft = sc.integrate();
                     if(enable_hard) current_hard = hc.integrate();
                     print_to_terminal( "gg -> suLsuL*(+X)", current_tree, current_virt, current_soft, current_hard);
                     xsec_to_json(j, "gg->suLsuL*", current_tree, current_virt, current_soft, current_hard);
                     total_xsec_tree += current_tree;
                     total_xsec_virt += current_virt;
                     total_xsec_soft += current_soft;
                     total_xsec_hard += current_hard;
                  }

                  if( subprocess == "") {
                     std::vector<std::array<int, 3>> flav {};
                     for( int el : { 1, -1, 3, -3, 4, -4, 5, -5}) flav.push_back({21, el, 2});
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::nullopt,
                        dS0, dC, flav,
                        {{
                           {SplittingKernel::Pqg, std::bind(&MRSSM::sigmaMRSSMTree_ddbar_suLsuLdagger, mrssm, _1, _2)},
                           {SplittingKernel::Pgq, std::bind(&MRSSM::sigmaMRSSMTree_gg_suLsuLdagger, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_gd_suLsuLdaggerd, mrssm, _1, _2),
                        dS0, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     std::array<double, 3> current_soft {};
                     std::array<double, 3> current_hard {};
                     if(enable_sc) current_soft = sc.integrate();
                     if(enable_hard) current_hard = hc.integrate();
                     print_to_terminal( "gq -> suLsuL*(+X)", {}, {}, current_soft, current_hard);
                     xsec_to_json(j, "gq->suLsuL*(+X)", {}, {}, current_soft, current_hard);
                     total_xsec_soft += current_soft;
                     total_xsec_hard += current_hard;
                  }

                  // g u > suL suLdagger
                  if( subprocess == "gu_suLsuLdaggeru" || subprocess == "" ) {
                     std::vector<std::array<int, 3>> flav {};
                     for( int el : { 2, -2 }) flav.push_back({21, el, 2});
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::nullopt,
                        dS0, dC, flav,
                        {{
                           {SplittingKernel::Pqg, std::bind(&MRSSM::sigmaMRSSMTree_uubar_suLsuLdagger, mrssm, _1, _2)},
                           {SplittingKernel::Pgq, std::bind(&MRSSM::sigmaMRSSMTree_gg_suLsuLdagger, mrssm, _1, _2)}
                        }},
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
                     std::array<double, 3> current_soft {};
                     std::array<double, 3> current_hard {};
                     if(enable_sc) current_soft = sc.integrate();
                     if(enable_hard) current_hard = hc.integrate();
                     print_to_terminal( "gu -> suLsuL*(+X)", {}, {}, current_soft, current_hard);
                     xsec_to_json(j, "gu->suLsuL*(+X)", {}, {}, current_soft, current_hard);
                     total_xsec_soft += current_soft;
                     total_xsec_hard += current_hard;
                  }

                  print_to_terminal("sum", total_xsec_tree, total_xsec_virt, total_xsec_soft, total_xsec_hard);
                  break;
               }
               case Channel::pp_sqsqdagger:
               {
                  const double m1 = pt.get<double>("masses.squarks");
                  const double m2 = pt.get<double>("masses.squarks");
                  std::array<double, 3> total_xsec_tree {};
                  std::array<double, 3> total_xsec_virt {};
                  std::array<double, 3> total_xsec_soft {};
                  std::array<double, 3> total_xsec_hard {};
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
                        {{
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_uubar_suLsuLdagger, mrssm, _1, _2)},
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_uubar_suLsuLdagger, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_uubar_suLsuLdaggerg, mrssm, _1, _2),
                        dS, dC,
                        flav,
                        hard_precision, hard_verbosity
                     );
                     std::array<double, 3> current_tree {};
                     std::array<double, 3> current_virt {};
                     std::array<double, 3> current_soft {};
                     std::array<double, 3> current_hard {};
                     if(enable_born) current_tree = tree.integrate();
                     if(enable_virt) current_virt = virt.integrate();
                     if(enable_sc) current_soft = sc.integrate();
                     if(enable_hard) current_hard = hc.integrate();
                     print_to_terminal("qqbar > sqsq*", current_tree, current_virt, current_soft, current_hard);
                     xsec_to_json(j, "qqbar->sqsq*", current_tree, current_virt, current_soft, current_hard);
                     total_xsec_tree += current_tree;
                     total_xsec_virt += current_virt;
                     total_xsec_soft += current_soft;
                     total_xsec_hard += current_hard;
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
                        {{
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_uubar_suLsuLdagger, mrssm, _1, _2)},
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_uubar_suLsuLdagger, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_uubar_suLsuLdaggerg, mrssm, _1, _2),
                        dS, dC,
                        flav,
                        hard_precision, hard_verbosity
                     );
                     std::array<double, 3> current_tree {};
                     std::array<double, 3> current_virt {};
                     std::array<double, 3> current_soft {};
                     std::array<double, 3> current_hard {};
                     if(enable_born) current_tree = tree.integrate();
                     if(enable_virt) current_virt = virt.integrate();
                     if(enable_sc) current_soft = sc.integrate();
                     if(enable_hard) current_hard = hc.integrate();
                     print_to_terminal("qqbar > sqsq*", current_tree, current_virt, current_soft, current_hard);
                     xsec_to_json(j, "qqbar->sqsq*", current_tree, current_virt, current_soft, current_hard);
                     total_xsec_tree += current_tree;
                     total_xsec_virt += current_virt;
                     total_xsec_soft += current_soft;
                     total_xsec_hard += current_hard;
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
                        {{
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_ddbar_suLsuLdagger, mrssm, _1, _2)},
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaMRSSMTree_ddbar_suLsuLdagger, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_ddbar_suLsuLdaggerg, mrssm, _1, _2),
                        dS, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     std::array<double, 3> current_tree {};
                     std::array<double, 3> current_virt {};
                     std::array<double, 3> current_soft {};
                     std::array<double, 3> current_hard {};
                     if(enable_born) current_tree = tree.integrate();
                     if(enable_virt) current_virt = virt.integrate();
                     if(enable_sc) current_soft = sc.integrate();
                     if(enable_hard) current_hard = hc.integrate();
                     print_to_terminal("ddbar->suLsuL*", current_tree, current_virt, current_soft, current_hard);
                     xsec_to_json(j, "ddbar->suLsuL*", current_tree, current_virt, current_soft, current_hard);
                     total_xsec_tree += current_tree;
                     total_xsec_virt += current_virt;
                     total_xsec_soft += current_soft;
                     total_xsec_hard += current_hard;
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
                        {{
                           {SplittingKernel::Pgg, std::bind(&MRSSM::sigmaMRSSMTree_gg_suLsuLdagger, mrssm, _1, _2)},
                           {SplittingKernel::Pgg, std::bind(&MRSSM::sigmaMRSSMTree_gg_suLsuLdagger, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_gg_suLsuLdaggerg, mrssm, _1, _2),
                        dS, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     std::array<double, 3> current_tree {};
                     std::array<double, 3> current_virt {};
                     std::array<double, 3> current_soft {};
                     std::array<double, 3> current_hard {};
                     if(enable_born) current_tree = tree.integrate();
                     if(enable_virt) current_virt = virt.integrate();
                     if(enable_sc) current_soft = sc.integrate();
                     if(enable_hard) current_hard = hc.integrate();
                     print_to_terminal("gg > suLsuL*", current_tree, current_virt, current_soft, current_hard);
                     xsec_to_json(j, "gg->suLsuL*", current_tree, current_virt, current_soft, current_hard);
                     total_xsec_tree += current_tree;
                     total_xsec_virt += current_virt;
                     total_xsec_soft += current_soft;
                     total_xsec_hard += current_hard;
                  }

                  if( subprocess == "") {
                     std::vector<std::array<int, 3>> flav {};
                     for( int el : { 1, -1, 3, -3, 4, -4, 5, -5}) flav.push_back({21, el, 2});
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::nullopt,
                        dS0, dC, flav,
                        {{
                           {SplittingKernel::Pqg, std::bind(&MRSSM::sigmaMRSSMTree_ddbar_suLsuLdagger, mrssm, _1, _2)},
                           {SplittingKernel::Pgq, std::bind(&MRSSM::sigmaMRSSMTree_gg_suLsuLdagger, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_gd_suLsuLdaggerd, mrssm, _1, _2),
                        dS0, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     std::array<double, 3> current_soft {};
                     std::array<double, 3> current_hard {};
                     if(enable_sc) current_soft = sc.integrate();
                     if(enable_hard) current_hard = hc.integrate();
                     print_to_terminal("gq > suLsuL*(+X)", {}, {}, current_soft, current_hard);
                     xsec_to_json(j, "gq->suLsuL*(+X)", {}, {}, current_soft, current_hard);
                     total_xsec_soft += current_soft;
                     total_xsec_hard += current_hard;
                  }

                  // g u > suL suLdagger
                  if( subprocess == "gu_suLsuLdaggeru" || subprocess == "" ) {
                     std::vector<std::array<int, 3>> flav {};
                     for( int el : { 2, -2 }) flav.push_back({21, el, 2});
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::nullopt,
                        dS0, dC, flav,
                        {{
                           {SplittingKernel::Pqg, std::bind(&MRSSM::sigmaMRSSMTree_ddbar_suLsuLdagger, mrssm, _1, _2)},
                           {SplittingKernel::Pgq, std::bind(&MRSSM::sigmaMRSSMTree_gg_suLsuLdagger, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMHard_gu_suLsuLdaggeru, mrssm, _1, _2),
                        dS0, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     std::array<double, 3> current_soft {};
                     std::array<double, 3> current_hard {};
                     if(enable_sc) current_soft = sc.integrate();
                     if(enable_hard) current_hard = hc.integrate();
                     print_to_terminal("gu > suLsuL*(+X)", {}, {}, current_soft, current_hard);
                     xsec_to_json(j, "gu->suLsuL*(+X)", {}, {}, current_soft, current_hard);
                     total_xsec_soft += current_soft;
                     total_xsec_hard += current_hard;
                  }

                  print_to_terminal("sum", total_xsec_tree, total_xsec_virt, total_xsec_soft, total_xsec_hard);
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
                  const double m1 = pt.get<double>("masses.sgluons");
                  std::array<double, 3> total_xsec_tree {};
                  std::array<double, 3> total_xsec_virt {};
                  std::array<double, 3> total_xsec_soft {};
                  std::array<double, 3> total_xsec_hard {};
                  {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {1, 2, 3, 4, 5}) {
                      flav.push_back({i, -i, 2});
                     }
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
                     */
                     XSection_SC sc(
                        parameters, m1, m1,
                        std::bind(&Sgluons::sgluons_qqbar_OOg_soft, sgluons, _1, _2, _3, _4, _5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pqq, std::bind(&Sgluons::sigmaSgluonsTree_qqbar_OO, sgluons, _1, _2)},
                           {SplittingKernel::Pqq, std::bind(&Sgluons::sigmaSgluonsTree_qqbar_OO, sgluons, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m1,
                        std::bind(&Sgluons::sgluons_qqbar_OOg_hard, sgluons, _1, _2),
                        dS, dC,
                        flav,
                        hard_precision, hard_verbosity
                     );
                     std::array<double, 3> current_tree {};
                     std::array<double, 3> current_virt {};
                     std::array<double, 3> current_soft {};
                     std::array<double, 3> current_hard {};
                     if(enable_born) current_tree = tree.integrate();
                     // if(enable_virt) current_virt = virt.integrate();
                     if(enable_sc) current_soft = sc.integrate();
                     if(enable_hard) current_hard = hc.integrate();
                     print_to_terminal("qqbar -> OO", current_tree, current_virt, current_soft, current_hard);
                     xsec_to_json(j, "qqbar->OO", current_tree, current_virt, current_soft, current_hard);
                     total_xsec_tree += current_tree;
                     total_xsec_virt += current_virt;
                     total_xsec_soft += current_soft;
                     total_xsec_hard += current_hard;
                  }
                  {
                     std::vector<std::array<int, 3>> flav {{21, 21, 1}};
                     XSection_Tree tree(
                        parameters, m1, m1,
                        std::bind(&Sgluons::matrixSgluonsTree_gg_OO, sgluons, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     /*
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixMRSSMVirt_uubar_suLsuLdagger, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav
                     );
                     */
                     XSection_SC sc(
                        parameters, m1, m1,
                        std::bind(&Sgluons::sgluons_gg_OOg_soft, sgluons, _1, _2, _3, _4, _5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pgg, std::bind(&Sgluons::sigmaSgluonsTree_gg_OO, sgluons, _1, _2)},
                           {SplittingKernel::Pgg, std::bind(&Sgluons::sigmaSgluonsTree_gg_OO, sgluons, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m1,
                        std::bind(&Sgluons::sgluons_gg_OOg_hard, sgluons, _1, _2),
                        dS, dC,
                        flav,
                        hard_precision, hard_verbosity
                     );
                     std::array<double, 3> current_tree {};
                     std::array<double, 3> current_virt {};
                     std::array<double, 3> current_soft {};
                     std::array<double, 3> current_hard {};
                     if(enable_born) current_tree = tree.integrate();
                     // if(enable_virt) current_virt = virt.integrate();
                     if(enable_sc) current_soft = sc.integrate();
                     if(enable_hard) current_hard = hc.integrate();
                     print_to_terminal("gg -> OO", current_tree, current_virt, current_soft, current_hard);
                     xsec_to_json(j, "gg->OO", current_tree, current_virt, current_soft, current_hard);
                     total_xsec_tree += current_tree;
                     total_xsec_virt += current_virt;
                     total_xsec_soft += current_soft;
                     total_xsec_hard += current_hard;
               }
               print_to_terminal("sum", total_xsec_tree, total_xsec_virt, total_xsec_soft, total_xsec_hard);
            }
            break;
         }
      }
      j["technical parameters"] = {
         {"dS", dS},
         {"dC", dC},
         {"WidthOverMass", pt.get<double>("technical parameters.WidthOverMass")},
         {"eta_sign", eta_sign},
         {"delta", delta}
      };
      break;
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
           + to_string(muR) + "_"
           + to_string(muF) + "_"
           + pt.get<string>("collider setup.pdf")
           + ".json";
   std::ofstream o(json_outputfile_name);
   o << std::setw(3) << j << std::endl;

   return 0;
}
