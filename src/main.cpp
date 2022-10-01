#include "version.hpp"
#include "XSections/XSection_Tree.hpp"
#include "XSections/XSection_Virt.hpp"
#include "XSections/XSection_SC.hpp"
#include "XSections/XSection_HnonC.hpp"
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
      sgluon_params.mO = pt.get<double>("masses.sgluon");
      sgluon_params.mt = pt.get<double>("masses.top");
      model = Model::Sgluons;
   }
   else {
      std::cerr << "\nError: Model " << pt.get<string>("process.model") << " not implemented\n\n";
      return 1;
   }

   enum class Process {
       pp_OO,
       pp_suLsuR,
       pp_sqLsqR,
       pp_sqLsqR_w_cc,
       pp_suLsuLdagger,
       pp_sqsqdagger,
       pp_glglbar
   };

   Process channel;
   if (pt.get<string>("process.process") == "pp_suLsuR") {
      channel = Process::pp_suLsuR;
   }
   else if (pt.get<string>("process.process") == "pp_sqLsqR") {
      channel = Process::pp_sqLsqR;
   }
   else if (pt.get<string>("process.process") == "pp_sqLsqR+cc") {
      channel = Process::pp_sqLsqR_w_cc;
   }
   else if (pt.get<string>("process.process") == "pp_suLsuLdagger") {
      channel = Process::pp_suLsuLdagger;
   }
   else if (pt.get<string>("process.process") == "pp_sqsqdagger") {
      channel = Process::pp_sqsqdagger;
   }
   else if (pt.get<string>("process.process") == "pp_OO") {
      channel = Process::pp_OO;
   }
   else if (pt.get<string>("process.process") == "pp_glglbar") {
      channel = Process::pp_glglbar;
   }
   else {
      cerr << "Error: Process not implemented\n";
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
      case Model::Sgluons:
         j["masses"] = {
            {"sgluon", pt.get<double>("masses.sgluon")},
            {"top",    pt.get<double>("masses.top")}
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

   std::vector<ChannelResult> allChannels;

   switch (order) {
      case Order::LO:
      {
      switch(model) {
         case Model::MRSSM:
            const MRSSM mrssm(mrssm_params);
            switch(channel) {
               case Process::pp_suLsuR:
               {
                  const double m1 = pt.get<double>("masses.squarks");
                  const double m2 = pt.get<double>("masses.squarks");
                  {
                     std::vector<std::array<int, 3>> flav {{2,2,1}};
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, _1, _2, _3),
                        flav,
                        born_precision, born_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "uu->suLsuR";
                     if(enable_born) chan.b = tree.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }
                  break;
               }
               case Process::pp_suLsuLdagger:
               {
                  const double m1 = pt.get<double>("masses.squarks");
                  const double m2 = pt.get<double>("masses.squarks");
                  {
                     std::vector<std::array<int, 3>> flav {{2, -2, 2}};
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixTree_uubar_suLsuLdagger, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "uubar->suLsuL*";
                     if(enable_born) chan.b = tree.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }
                  {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {1, 3, 4, 5}) {
                        flav.push_back({i,-i, 2});
                     }
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixTree_ddbar_suLsuLdagger, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "ddbar->suLsuL*";
                     if(enable_born) chan.b = tree.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }
                  {
                     std::vector<std::array<int, 3>> flav {{21, 21, 1}};
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixTree_GG_suLsuLdagger, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "gg->suLsuL*";
                     if(enable_born) chan.b = tree.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }
                  break;
               }
               case Process::pp_sqLsqR_w_cc:
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
                        std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "qq->sqLsqR+cc";
                     if(enable_born) chan.b = tree.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }
                  break;
               }
               case Process::pp_sqLsqR:
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
                        std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "qq->sqLsqR";
                     if(enable_born) chan.b = tree.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }
                  break;
               }
               case Process::pp_sqsqdagger:
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
                        std::bind(&MRSSM::matrixTree_uubar_suLsuLdagger, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "qqbar->sqsq*";
                     if(enable_born) chan.b = tree.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
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
                        std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "qq'bar->sqsq'*";
                     if(enable_born) chan.b = tree.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }
                  {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {1, 2, 3, 4, 5}) {
                        // 4 (squark flavours) * 2 (pp symmetry) * 2 (L and R squarks)
                        flav.push_back({i, -i, 2*4*2});
                     }
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixTree_ddbar_suLsuLdagger, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "ddbar->sq'sq'*";
                     if(enable_born) chan.b = tree.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }
                  {
                     // 5 squark flavours * (L + R)
                     std::vector<std::array<int, 3>> flav {{21, 21, 2*5}};
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixTree_GG_suLsuLdagger, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "gg->sqsq*";
                     if(enable_born) chan.b = tree.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }
                  break;
               }
               case Process::pp_glglbar:
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
                        std::bind(&MRSSM::matrixTree_uubar_glglbar, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "qqbar->gluglubar";
                     if(enable_born) chan.b = tree.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }
                  {
                     std::vector<std::array<int, 3>> flav {{21, 21, 1}};
                     XSection_Tree tree(
                        parameters, m, m,
                        std::bind(&MRSSM::matrixTree_gg_glglbar, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "gg->gluglubar";
                     if(enable_born) chan.b = tree.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }
                  break;
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
      constexpr double dS0 = 1e-9;
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
            const MRSSM mrssm(mrssm_params);
            switch(channel) {
               case Process::pp_suLsuR:
               {
                  const double m1 = pt.get<double>("masses.squarks");
                  const double m2 = pt.get<double>("masses.squarks");
                  // uu > suL suR (+g) process
                  if( subprocess == "" ) {
                     const std::vector<std::array<int, 3>> flav {{2, 2, 1}};
                     XSection_Tree tree(
                        parameters, m1, m2,
                        // same ME as in the MSSM
                        std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, _1, _2, _3),
                        flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixVirt_uu_suLsuR, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixSoft_uu_suLsuRg, mrssm, _1, _2, _3, _4, _5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaTree_uu_suLsuR, mrssm, _1, _2)},
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaTree_uu_suLsuR, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_uu_suLsuRg, mrssm, _1, _2),
                        dS, dC,
                        flav,
                        hard_precision, hard_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "uu->suLsuR(+g)";
                     if(enable_born) chan.b = tree.integrate();
                     if(enable_virt) chan.v = virt.integrate();
                     if(enable_sc) chan.s = sc.integrate();
                     if(enable_hard) chan.h = hc.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }

                  // gu > suL suR ubar process
                  if (subprocess == "gu_suLsuRubar" || subprocess == "" ) {
                     const std::vector<std::array<int, 3>> flav {{21, 2, 2}};
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::nullopt,
                        0., dC,
                        flav,
                        {{
                           {SplittingKernel::Pqg, std::bind(&MRSSM::sigmaTree_uu_suLsuR, mrssm, _1, _2)},
                           {SplittingKernel::Pgq, std::nullopt}
                        }},
                        sc_precision, sc_verbosity
                     );
                     auto f =
                        mrssm_params.MassGlu < mrssm_params.MassSq
                           ? &MRSSM::matrixHard_gu_suLsuRubar
                              : (mrssm_params.WidthGlu < 0 ? &MRSSM::matrixHard_gu_suLsuRubar_DR : &MRSSM::matrixHard_gu_suLsuRubar_DS);
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(f, mrssm, _1, _2),
                        dS0, dC,
                        flav,
                        hard_precision, hard_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "gu->suLsuRubar";
                     if(enable_sc) chan.s = sc.integrate();
                     if(enable_hard) chan.h = hc.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }
                  break;
               }
               case Process::pp_sqLsqR:
               {
                  // qq > sqL sqR (+g) process
                  const double m1 = pt.get<double>("masses.squarks");
                  const double m2 = pt.get<double>("masses.squarks");
                  {
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
                        std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixVirt_uu_suLsuR, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixSoft_uu_suLsuRg, mrssm, _1, _2, _3, _4, _5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaTree_uu_suLsuR, mrssm, _1, _2)},
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaTree_uu_suLsuR, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );

                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_uu_suLsuRg, mrssm, _1, _2),
                        dS, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "qq->sqLsqR(+g)";
                     if(enable_born) chan.b = tree.integrate();
                     if(enable_virt) chan.v = virt.integrate();
                     if(enable_sc) chan.s = sc.integrate();
                     if(enable_hard) chan.h = hc.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }
                  // gq > sqL sqR qbar process
                  {
                     std::vector<std::array<int, 3>> flav {};
                     for(int el : {1, 2, 3, 4, 5}) flav.push_back({21, el, 5});
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::nullopt,
                        0., dC,
                        flav,
                        {{
                           {SplittingKernel::Pqg, std::bind(&MRSSM::sigmaTree_uu_suLsuR, mrssm, _1, _2)},
                           {SplittingKernel::Pgq, std::nullopt}
                        }},
                        sc_precision, sc_verbosity
                     );
                     auto f =
                        mrssm_params.MassGlu < mrssm_params.MassSq
                           ? &MRSSM::matrixHard_gu_suLsuRubar
                              : (mrssm_params.WidthGlu < 0 ? &MRSSM::matrixHard_gu_suLsuRubar_DR : &MRSSM::matrixHard_gu_suLsuRubar_DS);
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(f, mrssm, _1, _2),
                        dS0, dC,
                        flav,
                        hard_precision, hard_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "gq->sqLsqRqbar";
                     if(enable_sc) chan.s = sc.integrate();
                     if(enable_hard) chan.h = hc.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }
                  break;
               }
               case Process::pp_sqLsqR_w_cc:
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
                        std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixVirt_uu_suLsuR, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixSoft_uu_suLsuRg, mrssm, _1, _2, _3, _4, _5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaTree_uu_suLsuR, mrssm, _1, _2)},
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaTree_uu_suLsuR, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );

                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_uu_suLsuRg, mrssm, _1, _2),
                        dS, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "qq->sqLsqR(+g)+cc";
                     if(enable_born) chan.b = tree.integrate();
                     if(enable_virt) chan.v = virt.integrate();
                     if(enable_sc) chan.s = sc.integrate();
                     if(enable_hard) chan.h = hc.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }
                  // gq > sqL sqR qbar process
                  {
                     std::vector<std::array<int, 3>> flav {};
                     for(int el : {1, -1, 2, -2, 3, -3, 4, -4, 5, -5}) flav.push_back({21, el, 2*5});
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::nullopt,
                        0., dC,
                        flav,
                        {{
                           {SplittingKernel::Pqg, std::bind(&MRSSM::sigmaTree_uu_suLsuR, mrssm, _1, _2)},
                           {SplittingKernel::Pgq, std::nullopt}
                        }},
                        sc_precision, sc_verbosity
                     );
                     auto f =
                        mrssm_params.MassGlu < mrssm_params.MassSq
                           ? &MRSSM::matrixHard_gu_suLsuRubar
                              : (mrssm_params.WidthGlu < 0 ? &MRSSM::matrixHard_gu_suLsuRubar_DR : &MRSSM::matrixHard_gu_suLsuRubar_DS);
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(f, mrssm, _1, _2),
                        dS0, dC,
                        flav,
                        hard_precision, hard_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "gq->sqLsqRqbar+cc";
                     if(enable_sc) chan.s = sc.integrate();
                     if(enable_hard) chan.h = hc.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }
                  break;
               }
               case Process::pp_suLsuLdagger:
               {
                  const double m1 = pt.get<double>("masses.squarks");
                  const double m2 = pt.get<double>("masses.squarks");
                  if (subprocess == "" || subprocess == "uubar_suLsuLdagger") {
                     std::vector<std::array<int, 3>> flav {{2, -2, 2}};
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixTree_uubar_suLsuLdagger, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixVirt_uubar_suLsuLdagger, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixSoft_uubar_suLsuLdaggerg, mrssm, _1, _2, _3, _4, _5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaTree_uubar_suLsuLdagger, mrssm, _1, _2)},
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaTree_uubar_suLsuLdagger, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_uubar_suLsuLdaggerg, mrssm, _1, _2),
                        dS, dC,
                        flav,
                        hard_precision, hard_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "uubar->suLsuL*(+g)";
                     if(enable_born) chan.b = tree.integrate();
                     if(enable_virt) chan.v = virt.integrate();
                     if(enable_sc) chan.s = sc.integrate();
                     if(enable_hard) chan.h = hc.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }

                  if( subprocess == "") {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {1, 3, 4, 5}) {
                        flav.push_back({i,-i, 2});
                     }
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixTree_ddbar_suLsuLdagger, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixVirt_ddbar_suLsuLdagger, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixSoft_ddbar_suLsuLdaggerg, mrssm, _1, _2, _3, _4, _5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaTree_ddbar_suLsuLdagger, mrssm, _1, _2)},
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaTree_ddbar_suLsuLdagger, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_ddbar_suLsuLdaggerg, mrssm, _1, _2),
                        dS, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "ddbar->suLsuL*(+g)";
                     if(enable_born) chan.b = tree.integrate();
                     if(enable_virt) chan.v = virt.integrate();
                     if(enable_sc) chan.s = sc.integrate();
                     if(enable_hard) chan.h = hc.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }

                  {
                     std::vector<std::array<int, 3>> flav {{21, 21, 1}};
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixTree_GG_suLsuLdagger, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixVirt_GG_suLsuLdagger, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixSoft_gg_suLsuLdaggerg, mrssm, _1, _2, _3, _4, _5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pgg, std::bind(&MRSSM::sigmaTree_gg_suLsuLdagger, mrssm, _1, _2)},
                           {SplittingKernel::Pgg, std::bind(&MRSSM::sigmaTree_gg_suLsuLdagger, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_gg_suLsuLdaggerg, mrssm, _1, _2),
                        dS, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "gg->suLsuL*(+g)";
                     if(enable_born) chan.b = tree.integrate();
                     if(enable_virt) chan.v = virt.integrate();
                     if(enable_sc) chan.s = sc.integrate();
                     if(enable_hard) chan.h = hc.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }

                  {
                     std::vector<std::array<int, 3>> flav {};
                     for( int el : { 1, -1, 3, -3, 4, -4, 5, -5}) flav.push_back({21, el, 2});
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::nullopt,
                        0., dC, flav,
                        {{
                           {SplittingKernel::Pqg, std::bind(&MRSSM::sigmaTree_ddbar_suLsuLdagger, mrssm, _1, _2)},
                           {SplittingKernel::Pgq, std::bind(&MRSSM::sigmaTree_gg_suLsuLdagger, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_gd_suLsuLdaggerd, mrssm, _1, _2),
                        dS0, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "gq->suLsuL*(+X)";
                     if(enable_sc) chan.s = sc.integrate();
                     if(enable_hard) chan.h = hc.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }

                  // g u > suL suLdagger
                  if( subprocess == "gu_suLsuLdaggeru" || subprocess == "" ) {
                     std::vector<std::array<int, 3>> flav {};
                     for( int el : { 2, -2 }) flav.push_back({21, el, 2});
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::nullopt,
                        0., dC, flav,
                        {{
                           {SplittingKernel::Pqg, std::bind(&MRSSM::sigmaTree_uubar_suLsuLdagger, mrssm, _1, _2)},
                           {SplittingKernel::Pgq, std::bind(&MRSSM::sigmaTree_gg_suLsuLdagger, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     auto f =
                        mrssm_params.MassGlu < mrssm_params.MassSq
                           ? &MRSSM::matrixHard_gu_suLsuLdaggeru
                              : (mrssm_params.WidthGlu < 0 ? &MRSSM::matrixHard_gu_suLsuLdaggeru_DR : &MRSSM::matrixHard_gu_suLsuLdaggeru_DS);
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(f, mrssm, _1, _2),
                        dS0, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "gu->suLsuL*(+X)";
                     if(enable_sc) chan.s = sc.integrate();
                     if(enable_hard) chan.h = hc.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }
                  break;
               }
               case Process::pp_sqsqdagger:
               {
                  const double m1 = pt.get<double>("masses.squarks");
                  const double m2 = pt.get<double>("masses.squarks");
                  {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {1, 2, 3, 4, 5}) {
                        flav.push_back({i, -i, 4});
                     }
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixTree_uubar_suLsuLdagger, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixVirt_uubar_suLsuLdagger, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixSoft_uubar_suLsuLdaggerg, mrssm, _1, _2, _3, _4, _5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaTree_uubar_suLsuLdagger, mrssm, _1, _2)},
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaTree_uubar_suLsuLdagger, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_uubar_suLsuLdaggerg, mrssm, _1, _2),
                        dS, dC,
                        flav,
                        hard_precision, hard_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "qqbar->sqsq*(+g)";
                     if(enable_born) chan.b = tree.integrate();
                     if(enable_virt) chan.v = virt.integrate();
                     if(enable_sc) chan.s = sc.integrate();
                     if(enable_hard) chan.h = hc.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
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
                        std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixVirt_uu_suLsuR, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixSoft_uu_suLsuRg, mrssm, _1, _2, _3, _4, _5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaTree_uu_suLsuR, mrssm, _1, _2)},
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaTree_uu_suLsuR, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_uu_suLsuRg, mrssm, _1, _2),
                        dS, dC,
                        flav,
                        hard_precision, hard_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "qq'bar->sqsq'*(+g)";
                     if(enable_born) chan.b = tree.integrate();
                     if(enable_virt) chan.v = virt.integrate();
                     if(enable_sc) chan.s = sc.integrate();
                     if(enable_hard) chan.h = hc.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }

                  {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {1, 2, 3, 4, 5}) {
                        // 4 (squark flavours) * 2 (pp symmetry) * 2 (L and R squarks)
                        flav.push_back({i,-i, 2*4*2});
                     }
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixTree_ddbar_suLsuLdagger, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixVirt_ddbar_suLsuLdagger, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixSoft_ddbar_suLsuLdaggerg, mrssm, _1, _2, _3, _4, _5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaTree_ddbar_suLsuLdagger, mrssm, _1, _2)},
                           {SplittingKernel::Pqq, std::bind(&MRSSM::sigmaTree_ddbar_suLsuLdagger, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_ddbar_suLsuLdaggerg, mrssm, _1, _2),
                        dS, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "qqbar->sq'sq'*(+g)";
                     if(enable_born) chan.b = tree.integrate();
                     if(enable_virt) chan.v = virt.integrate();
                     if(enable_sc) chan.s = sc.integrate();
                     if(enable_hard) chan.h = hc.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }

                  {
                     // 5 squark flavours * L and R
                     std::vector<std::array<int, 3>> flav {{21, 21, 2*5}};
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixTree_GG_suLsuLdagger, mrssm, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixVirt_GG_suLsuLdagger, mrssm, _1, _2, _3, _4, _5, _6, _7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixSoft_gg_suLsuLdaggerg, mrssm, _1, _2, _3, _4, _5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pgg, std::bind(&MRSSM::sigmaTree_gg_suLsuLdagger, mrssm, _1, _2)},
                           {SplittingKernel::Pgg, std::bind(&MRSSM::sigmaTree_gg_suLsuLdagger, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_gg_suLsuLdaggerg, mrssm, _1, _2),
                        dS, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "gg->sqsq*(+g)";
                     if(enable_born) chan.b = tree.integrate();
                     if(enable_virt) chan.v = virt.integrate();
                     if(enable_sc) chan.s = sc.integrate();
                     if(enable_hard) chan.h = hc.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }

                  {
                     std::vector<std::array<int, 3>> flav {};
                     for(int el : { 1, -1, 2, -2, 3, -3, 4, -4, 5, -5}) flav.push_back({21, el, 2*2*4 /* initial state permutation x L+R squarks x 4-squark flavours */});
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::nullopt,
                        dS0, dC, flav,
                        {{
                           {SplittingKernel::Pqg, std::bind(&MRSSM::sigmaTree_ddbar_suLsuLdagger, mrssm, _1, _2)},
                           {SplittingKernel::Pgq, std::bind(&MRSSM::sigmaTree_gg_suLsuLdagger, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_gd_suLsuLdaggerd, mrssm, _1, _2),
                        dS0, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "gq->sq'sq'*(+X)";
                     if(enable_sc) chan.s = sc.integrate();
                     if(enable_hard) chan.h = hc.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }

                  // g q > sq sqdagger
                  {
                     std::vector<std::array<int, 3>> flav {};
                     for( int el : { 1, -1, 2, -2, 3, -3, 4, -4, 5, -5 }) flav.push_back({21, el, 2*2});
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::nullopt,
                        dS0, dC, flav,
                        {{
                           {SplittingKernel::Pqg, std::bind(&MRSSM::sigmaTree_uubar_suLsuLdagger, mrssm, _1, _2)},
                           {SplittingKernel::Pgq, std::bind(&MRSSM::sigmaTree_gg_suLsuLdagger, mrssm, _1, _2)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     auto f =
                        mrssm_params.MassGlu < mrssm_params.MassSq
                           ? &MRSSM::matrixHard_gu_suLsuLdaggeru
                              : (mrssm_params.WidthGlu < 0 ? &MRSSM::matrixHard_gu_suLsuLdaggeru_DR : &MRSSM::matrixHard_gu_suLsuLdaggeru_DS);
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(f, mrssm, _1, _2),
                        dS0, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "gq->sqsq*(+X)";
                     if(enable_sc) chan.s = sc.integrate();
                     if(enable_hard) chan.h = hc.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }
                  {
                     std::vector<std::array<int, 3>> flav {};
                     for( int el : { 1, -1, 2, -2, 3, -3, 4, -4, 5, -5 }) flav.push_back({21, el, 2*4});
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::nullopt,
                        dS0, dC, flav,
                        {{
                           {SplittingKernel::Pqg, std::bind(&MRSSM::sigmaTree_uu_suLsuR, mrssm, _1, _2)},
                           {SplittingKernel::Pgq, std::nullopt}
                        }},
                        sc_precision, sc_verbosity
                     );
                     auto f =
                        mrssm_params.MassGlu < mrssm_params.MassSq
                           ? &MRSSM::matrixHard_gu_suLsuRubar
                              : (mrssm_params.WidthGlu < 0 ? &MRSSM::matrixHard_gu_suLsuRubar_DR : &MRSSM::matrixHard_gu_suLsuRubar_DS);
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(f, mrssm, _1, _2),
                        dS0, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "gq->sqsq'*(+X)+h.c.";
                     if(enable_sc) chan.s = sc.integrate();
                     if(enable_hard) chan.h = hc.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }
                  break;
               }
            } //channels
         break;
         } // MRSSM
         case Model::Sgluons:
            Sgluons sgluons(sgluon_params);
            switch(channel) {
               case Process::pp_OO:
               {
                  const double m1 = pt.get<double>("masses.sgluon");
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
                        std::bind(&MRSSM::matrixVirt_uubar_suLsuLdagger, mrssm, _1, _2, _3, _4, _5, _6, _7),
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
                     ChannelResult chan;
                     chan.channel_name = "qqbar->OO";
                     // if(enable_born) chan.b = tree.integrate();
                     // if(enable_virt) chan.v = virt.integrate();
                     if(enable_sc) chan.s = sc.integrate();
                     if(enable_hard) chan.h = hc.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }
                  {
                     const std::vector<std::array<int, 3>> flav {{21, 21, 1}};
                     XSection_Tree tree(
                        parameters, m1, m1,
                        std::bind(&Sgluons::matrixSgluonsTree_gg_OO, sgluons, _1, _2, _3), flav,
                        born_precision, born_verbosity
                     );
                     /*
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixVirt_uubar_suLsuLdagger, mrssm, _1, _2, _3, _4, _5, _6, _7),
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
                     ChannelResult chan;
                     chan.channel_name = "gg->OO";
                     // if(enable_born) chan.b = tree.integrate();
                     // if(enable_virt) chan.v = virt.integrate();
                     if(enable_sc) chan.s = sc.integrate();
                     if(enable_hard) chan.h = hc.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
               }
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

   switch (order) {
      case Order::LO:
      {
         std::array<double, 3> tot_b {};
         for (ChannelResult const& ch : allChannels) {
            xsec_to_json(j, ch);
            tot_b += ch.b;
         }
         ChannelResult total {"total", tot_b};
         print_to_terminal(total);
         break;
      }
      case Order::NLO:
      {
         std::array<double, 3> tot_b {};
         std::array<double, 3> tot_v {};
         std::array<double, 3> tot_s {};
         std::array<double, 3> tot_h {};
         for (ChannelResult const& ch : allChannels) {
            xsec_to_json(j, ch);
            tot_b += ch.b;
            tot_v += ch.v.value_or(std::array<double, 3>{0., 0., 0.});
            tot_s += ch.s.value_or(std::array<double, 3>{0., 0., 0.});
            tot_h += ch.h.value_or(std::array<double, 3>{0., 0., 0.});
         }
         ChannelResult total {"total", tot_b, tot_v, tot_s, tot_h};
         print_to_terminal(total);
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
   string json_outputfile_name =
      vm.count("json-outputfile-name")
         ? vm["json-outputfile-name"].as<string>()
         : pt.get<string>("process.process") + "_"
           + to_string(pt.get<double>("masses.squarks")) + "_"
           + to_string(pt.get<double>("masses.gluino")) + "_"
           + to_string(pt.get<double>("masses.pseudoscalar_sgluon")) + "_"
           + to_string(pt.get<double>("collider setup.sqrt_S")) + "_"
           + to_string(muR) + "_"
           + to_string(muF) + "_"
           + pt.get<string>("collider setup.pdf");
   int file_index = 0;
   while (std::ifstream(json_outputfile_name + ".json")) {
      file_index += 1;
      std::cout
         << "File " << json_outputfile_name + ".json already exists. "
         << "Trying " << json_outputfile_name + "_" + std::to_string(file_index)  + ".json instead.\n";
      json_outputfile_name += "_" + std::to_string(file_index);
   }
   std::ofstream o(vm.count("json-outputfile-name") ? vm["json-outputfile-name"].as<string>() : json_outputfile_name + ".json");
   o << std::setw(3) << j << std::endl;

   return 0;
}
