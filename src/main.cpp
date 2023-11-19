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
#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/stopwatch.h"

#include <chrono>
#include <filesystem>
#include <functional>
#include <sstream>

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
      ("card,f", po::value<string>()->required(), "path to a run card")
      ("subprocess", po::value<string>() -> default_value(""), "")
      ("log-level", po::value<string>() -> default_value("info"), "")
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

   spdlog::set_pattern("[%H:%M:%S] [%^%l%$] : %v");
   auto console = spdlog::stdout_color_mt("console");
   if (vm["log-level"].as<std::string>() == "trace") {
      console->set_level(spdlog::level::trace);
   }
   else if (vm["log-level"].as<std::string>() == "info") {
      console->set_level(spdlog::level::info);
   }
   else if (vm["log-level"].as<std::string>() == "warn") {
      console->set_level(spdlog::level::warn);
   }
   else if (vm["log-level"].as<std::string>() == "error") {
      console->set_level(spdlog::level::err);
   }
   else if (vm["log-level"].as<std::string>() == "critical") {
      console->set_level(spdlog::level::critical);
   }
   else if (vm["log-level"].as<std::string>() == "off") {
      console->set_level(spdlog::level::off);
   }
   else {
      spdlog::get("console")->error("Unknown log-level {}", vm["log-level"].as<std::string>());
      return 1;
   }
   spdlog::get("console")->info("Log level set to {}", vm["log-level"].as<std::string>());

   if (vm["precision-virt"].as<int>() > 4) {
      spdlog::get("console")->warn("Virtual integration with such a high precision is known to use very large amounts of RAM. In case of problems consider reducing precision");
   }

   bool enable_born = vm["enable-born"].as<bool>();
   bool enable_virt = vm["enable-virt"].as<bool>();
   bool enable_sc = vm["enable-sc"].as<bool>();
   bool enable_hard = vm["enable-hard"].as<bool>();

   string card = vm["card"].as<string>();
   string subprocess = vm["subprocess"].as<string>();

   boost::property_tree::ptree pt;
   try {
      if (card == "-") {
         boost::property_tree::ini_parser::read_ini(std::cin, pt);
      }
      else {
         boost::property_tree::ini_parser::read_ini(card, pt);
      }
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

   MSSMParameters mssm_params;
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
      mssm_params.MassTop = pt.get<double>("masses.top");
      mssm_params.MassGlu = pt.get<double>("masses.gluino");
      mssm_params.MassSq = pt.get<double>("masses.squarks");
      mssm_params.eta_sign = pt.get<int>("technical parameters.eta_sign", -1);
      mssm_params.delta = pt.get<double>("technical parameters.delta", 0.);
      mssm_params.WidthGlu = pt.get<double>("technical parameters.WidthOverMass", -1.) * mrssm_params.MassGlu;
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
       pp_glglbar,
       pp_sqgl_w_cc
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
   else if (pt.get<string>("process.process") == "pp_sqgl+cc") {
      channel = Process::pp_sqgl_w_cc;
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
   LHAPDF::setVerbosity(0);
   std::unique_ptr<LHAPDF::PDF> pdf(LHAPDF::mkPDF(pt.get<std::string>("collider setup.pdf")));
   stringstream ss;
   ss << "Using " << pdf->set().name() << " PDF set, member #" << pdf->memberID()
      << ", version " << pdf->dataversion() << "; LHAPDF ID = " << pdf->lhapdfID();
   console->info(ss.str());

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

   spdlog::stopwatch sw;
   switch (order) {
      case Order::LO:
      {
      switch(model) {
         case Model::MRSSM:
         {
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
                        std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4),
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
                        std::bind(&MRSSM::matrixTree_uubar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
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
                        std::bind(&MRSSM::matrixTree_ddbar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
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
                        std::bind(&MRSSM::matrixTree_gg_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
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
                        std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
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
                        std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
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
                        std::bind(&MRSSM::matrixTree_uubar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
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
                        std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
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
                        std::bind(&MRSSM::matrixTree_ddbar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
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
                        std::bind(&MRSSM::matrixTree_gg_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
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
                        std::bind(&MRSSM::matrixTree_uubar_glglbar, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
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
                        std::bind(&MRSSM::matrixTree_gg_glglbar, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
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
               case Process::pp_sqgl_w_cc:
               {
                  const double m1 = pt.get<double>("masses.squarks");
                  const double m2 = pt.get<double>("masses.gluino");
                  std::array<double, 3> result {};
                  {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {1,2,3,4,5}) {
                        flav.push_back({i, 21, 2});
                     }
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixTree_ug_suLglbar, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
                        born_precision, born_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "qg->sqLglbar";
                     if(enable_born) chan.b = tree.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }
                  {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {1,2,3,4,5}) {
                        flav.push_back({i, 21, 2});
                     }
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixTree_ug_suLglbar, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
                        born_precision, born_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "qg->sqRgl";
                     if(enable_born) chan.b = tree.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }
                  {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {-1,-2,-3,-4,-5}) {
                        flav.push_back({i, 21, 2});
                     }
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixTree_ug_suLglbar, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
                        born_precision, born_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "qbarg->sqLdaggergl";
                     if(enable_born) chan.b = tree.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }
                  {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {-1, -2, -3, -4, -5}) {
                        flav.push_back({i, 21, 2});
                     }
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixTree_ug_suLglbar, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
                        born_precision, born_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "qbarg->sqRdaggerglbar";
                     if(enable_born) chan.b = tree.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }
               }
               default:
                  break;
            } // end of process block
         }
         case Model::MSSM:
            const MSSM mssm(mssm_params);
            switch(channel) {
               case Process::pp_suLsuR:
                  {
                  }
               break;
            }
      } // end of model block
      break;
      } // end of LO block
      case Order::NLO: {
      const double dS = pt.get<double>("technical parameters.dS", 1e-5);
      // for the matrix elements that are regular in the limit dS -> 0 because the phase space parametrization
      // fails if we are exactly on the threshold
      constexpr double dS0 = 1e-9;
      const double dC = pt.get<double>("technical parameters.dC", 1e-6);
      spdlog::get("console")->info("Using phase space slicing parameters δS={0} and δC={1}", dS, dC);
      if (dC > dS) {
         spdlog::get("console")->warn("Warning: δC should be always << than δS");
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
                        std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4),
                        flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixVirt_uu_suLsuR, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixSoft_uu_suLsuRg_finite, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pqq, std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
                           {SplittingKernel::Pqq, std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_uu_suLsuRg, mrssm, std::placeholders::_1, std::placeholders::_2),
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
                           {SplittingKernel::Pqg, std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
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
                        std::bind(f, mrssm, std::placeholders::_1, std::placeholders::_2),
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
                        std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixVirt_uu_suLsuR, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixSoft_uu_suLsuRg_finite, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pqq, std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
                           {SplittingKernel::Pqq, std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)}
                        }},
                        sc_precision, sc_verbosity
                     );

                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_uu_suLsuRg, mrssm, std::placeholders::_1, std::placeholders::_2),
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
                           {SplittingKernel::Pqg, std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
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
                        std::bind(f, mrssm, std::placeholders::_1, std::placeholders::_2),
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
                        std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixVirt_uu_suLsuR, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixSoft_uu_suLsuRg_finite, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pqq, std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
                           {SplittingKernel::Pqq, std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)}
                        }},
                        sc_precision, sc_verbosity
                     );

                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_uu_suLsuRg, mrssm, std::placeholders::_1, std::placeholders::_2),
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
                           {SplittingKernel::Pqg, std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
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
                        std::bind(f, mrssm, std::placeholders::_1, std::placeholders::_2),
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
                        std::bind(&MRSSM::matrixTree_uubar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixVirt_uubar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixSoft_uubar_suLsuLdaggerg_finite, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pqq, std::bind(&MRSSM::matrixTree_uubar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
                           {SplittingKernel::Pqq, std::bind(&MRSSM::matrixTree_uubar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_uubar_suLsuLdaggerg, mrssm, std::placeholders::_1, std::placeholders::_2),
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

                  if (subprocess == "" || subprocess == "ddbar_suLsuLdagger") {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {1, 3, 4, 5}) {
                        flav.push_back({i, -i, 2});
                     }
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixTree_ddbar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixVirt_ddbar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixSoft_ddbar_suLsuLdaggerg_finite, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pqq, std::bind(&MRSSM::matrixTree_ddbar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
                           {SplittingKernel::Pqq, std::bind(&MRSSM::matrixTree_ddbar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_ddbar_suLsuLdaggerg, mrssm, std::placeholders::_1, std::placeholders::_2),
                        dS, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "qqbar->suLsuL*(+g), q!=u";
                     if(enable_born) chan.b = tree.integrate();
                     if(enable_virt) chan.v = virt.integrate();
                     if(enable_sc) chan.s = sc.integrate();
                     if(enable_hard) chan.h = hc.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }

                  if (subprocess == "" || subprocess == "gg_suLsuLdagger") {
                     std::vector<std::array<int, 3>> flav {{21, 21, 1}};
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixTree_gg_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixVirt_gg_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixSoft_gg_suLsuLdaggerg_finite, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pgg, std::bind(&MRSSM::matrixTree_gg_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
                           {SplittingKernel::Pgg, std::bind(&MRSSM::matrixTree_gg_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_gg_suLsuLdaggerg, mrssm, std::placeholders::_1, std::placeholders::_2),
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

                  if (subprocess == "" || subprocess == "gq_suLsuLdaggerq") {
                     std::vector<std::array<int, 3>> flav {};
                     for( int el : { 1, 3, 4, 5}) flav.push_back({21, el, 2});
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::nullopt,
                        0., dC, flav,
                        {{
                           {SplittingKernel::Pqg, std::bind(&MRSSM::matrixTree_ddbar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
                           {SplittingKernel::Pgq, std::bind(&MRSSM::matrixTree_gg_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_gd_suLsuLdaggerd, mrssm, std::placeholders::_1, std::placeholders::_2),
                        dS0, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "gq->suLsuL*q";
                     if(enable_sc) chan.s = sc.integrate();
                     if(enable_hard) chan.h = hc.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }

                  if (subprocess == "" || subprocess == "gqbar_suLsuLdaggerqbar") {
                     std::vector<std::array<int, 3>> flav {};
                     for( int el : { -1, -3, -4, -5}) flav.push_back({21, el, 2});
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::nullopt,
                        0., dC, flav,
                        {{
                           {SplittingKernel::Pqg, std::bind(&MRSSM::matrixTree_ddbar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
                           {SplittingKernel::Pgq, std::bind(&MRSSM::matrixTree_gg_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_gdbar_suLsuLdaggerdbar, mrssm, std::placeholders::_1, std::placeholders::_2),
                        dS0, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "gqbar->suLsuL*qbar";
                     if(enable_sc) chan.s = sc.integrate();
                     if(enable_hard) chan.h = hc.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }

                  // g u > suL suLdagger
                  if( subprocess == "gu_suLsuLdaggeru" || subprocess == "" ) {
                     std::vector<std::array<int, 3>> flav {};
                     for( int el : { 2 }) flav.push_back({21, el, 2});
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::nullopt,
                        0., dC, flav,
                        {{
                           {SplittingKernel::Pqg, std::bind(&MRSSM::matrixTree_uubar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
                           {SplittingKernel::Pgq, std::bind(&MRSSM::matrixTree_gg_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     auto f =
                        mrssm_params.MassGlu < mrssm_params.MassSq
                           ? &MRSSM::matrixHard_gu_suLsuLdaggeru
                              : (mrssm_params.WidthGlu < 0 ? &MRSSM::matrixHard_gu_suLsuLdaggeru_DR : &MRSSM::matrixHard_gu_suLsuLdaggeru_DS);
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(f, mrssm, std::placeholders::_1, std::placeholders::_2),
                        dS0, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "gu->suLsuL*u";
                     if(enable_sc) chan.s = sc.integrate();
                     if(enable_hard) chan.h = hc.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }

                  // g ubar > suL suLdagger
                  if( subprocess == "gubar_suLsuLdaggerubar" || subprocess == "" ) {
                     std::vector<std::array<int, 3>> flav {};
                     for( int el : { -2 }) flav.push_back({21, el, 2});
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::nullopt,
                        0., dC, flav,
                        {{
                           {SplittingKernel::Pqg, std::bind(&MRSSM::matrixTree_uubar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
                           {SplittingKernel::Pgq, std::bind(&MRSSM::matrixTree_gg_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     auto f =
                        mrssm_params.MassGlu < mrssm_params.MassSq
                           ? &MRSSM::matrixHard_gubar_suLsuLdaggerubar
                              : (mrssm_params.WidthGlu < 0 ? &MRSSM::matrixHard_gubar_suLsuLdaggerubar_DR : &MRSSM::matrixHard_gu_suLsuLdaggeru_DS);
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(f, mrssm, std::placeholders::_1, std::placeholders::_2),
                        dS0, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "gubar->suLsuL*ubar";
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
                        std::bind(&MRSSM::matrixTree_uubar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixVirt_uubar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixSoft_uubar_suLsuLdaggerg_finite, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pqq, std::bind(&MRSSM::matrixTree_uubar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
                           {SplittingKernel::Pqq, std::bind(&MRSSM::matrixTree_uubar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_uubar_suLsuLdaggerg, mrssm, std::placeholders::_1, std::placeholders::_2),
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
                        std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixVirt_uu_suLsuR, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixSoft_uu_suLsuRg_finite, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pqq, std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
                           {SplittingKernel::Pqq, std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_uu_suLsuRg, mrssm, std::placeholders::_1, std::placeholders::_2),
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
                        std::bind(&MRSSM::matrixTree_ddbar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixVirt_ddbar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixSoft_ddbar_suLsuLdaggerg_finite, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pqq, std::bind(&MRSSM::matrixTree_ddbar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
                           {SplittingKernel::Pqq, std::bind(&MRSSM::matrixTree_ddbar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_ddbar_suLsuLdaggerg, mrssm, std::placeholders::_1, std::placeholders::_2),
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
                        std::bind(&MRSSM::matrixTree_gg_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixVirt_gg_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixSoft_gg_suLsuLdaggerg_finite, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pgg, std::bind(&MRSSM::matrixTree_gg_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
                           {SplittingKernel::Pgg, std::bind(&MRSSM::matrixTree_gg_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_gg_suLsuLdaggerg, mrssm, std::placeholders::_1, std::placeholders::_2),
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
                           {SplittingKernel::Pqg, std::bind(&MRSSM::matrixTree_ddbar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
                           {SplittingKernel::Pgq, std::bind(&MRSSM::matrixTree_gg_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_gd_suLsuLdaggerd, mrssm, std::placeholders::_1, std::placeholders::_2),
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
                           {SplittingKernel::Pqg, std::bind(&MRSSM::matrixTree_uubar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
                           {SplittingKernel::Pgq, std::bind(&MRSSM::matrixTree_gg_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     auto f =
                        mrssm_params.MassGlu < mrssm_params.MassSq
                           ? &MRSSM::matrixHard_gu_suLsuLdaggeru
                              : (mrssm_params.WidthGlu < 0 ? &MRSSM::matrixHard_gu_suLsuLdaggeru_DR : &MRSSM::matrixHard_gu_suLsuLdaggeru_DS);
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(f, mrssm, std::placeholders::_1, std::placeholders::_2),
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
                           {SplittingKernel::Pqg, std::bind(&MRSSM::matrixTree_uu_suLsuR, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
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
                        std::bind(f, mrssm, std::placeholders::_1, std::placeholders::_2),
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
         case Model::MSSM:
         {
            MSSM mssm(mssm_params);
            switch (channel) {
               case Process::pp_suLsuLdagger:
               {
                  /*
                  const double m1 = pt.get<double>("masses.squarks");
                  const double m2 = pt.get<double>("masses.squarks");
                  if (subprocess == "" || subprocess == "uubar_suLsuLdagger") {
                     std::vector<std::array<int, 3>> flav {{2, -2, 2}};
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MSSM::matrixTree_uubar_suLsuLdagger, mssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MSSM::matrixVirt_uubar_suLsuLdagger, mssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MSSM::matrixSoft_uubar_suLsuLdaggerg_finite, mssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pqq, std::bind(&MSSM::matrixTree_uubar_suLsuLdagger, mssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
                           {SplittingKernel::Pqq, std::bind(&MSSM::matrixTree_uubar_suLsuLdagger, mssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MSSM::matrixHard_uubar_suLsuLdaggerg, mssm, std::placeholders::_1, std::placeholders::_2),
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

                  if (subprocess == "ddbar_suLsuLdagger" || subprocess == "") {
                     std::vector<std::array<int, 3>> flav {};
                     for (int i : {1, 3, 4, 5}) {
                        flav.push_back({i, -i, 2});
                     }
                     XSection_Tree tree(
                        parameters, m1, m2,
                        std::bind(&MSSM::matrixTree_ddbar_suLsuLdagger, mssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MSSM::matrixVirt_ddbar_suLsuLdagger, mssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MSSM::matrixSoft_ddbar_suLsuLdaggerg_finite, mssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pqq, std::bind(&MSSM::matrixTree_ddbar_suLsuLdagger, mssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
                           {SplittingKernel::Pqq, std::bind(&MSSM::matrixTree_ddbar_suLsuLdagger, mssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MSSM::matrixHard_ddbar_suLsuLdaggerg, mssm, std::placeholders::_1, std::placeholders::_2),
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
                        std::bind(&MRSSM::matrixTree_gg_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
                        born_precision, born_verbosity
                     );
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixVirt_gg_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7),
                        flav,
                        virt_precision, virt_verbosity
                     );
                     XSection_SC sc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixSoft_gg_suLsuLdaggerg_finite, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pgg, std::bind(&MRSSM::matrixTree_gg_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
                           {SplittingKernel::Pgg, std::bind(&MRSSM::matrixTree_gg_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_gg_suLsuLdaggerg, mrssm, std::placeholders::_1, std::placeholders::_2),
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
                           {SplittingKernel::Pqg, std::bind(&MRSSM::matrixTree_ddbar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
                           {SplittingKernel::Pgq, std::bind(&MRSSM::matrixTree_gg_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixHard_gd_suLsuLdaggerd, mrssm, std::placeholders::_1, std::placeholders::_2),
                        dS0, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "gq->suLsuL*q";
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
                           {SplittingKernel::Pqg, std::bind(&MRSSM::matrixTree_uubar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
                           {SplittingKernel::Pgq, std::bind(&MRSSM::matrixTree_gg_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     auto f =
                        mrssm_params.MassGlu < mrssm_params.MassSq
                           ? &MRSSM::matrixHard_gu_suLsuLdaggeru
                              : (mrssm_params.WidthGlu < 0 ? &MRSSM::matrixHard_gu_suLsuLdaggeru_DR : &MRSSM::matrixHard_gu_suLsuLdaggeru_DS);
                     XSection_HnonC hc(
                        parameters, m1, m2,
                        std::bind(f, mrssm, std::placeholders::_1, std::placeholders::_2),
                        dS0, dC, flav,
                        hard_precision, hard_verbosity
                     );
                     ChannelResult chan;
                     chan.channel_name = "gu->suLsuL*u";
                     if(enable_sc) chan.s = sc.integrate();
                     if(enable_hard) chan.h = hc.integrate();
                     print_to_terminal(chan);
                     allChannels.push_back(std::move(chan));
                  }
                  */
                  break;
               }
            } //channels
            break;
         } // MSSM
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
                        std::bind(&Sgluons::matrixSgluonsTree_qqbar_OO, sgluons, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
                        born_precision, born_verbosity
                     );
                     /*
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixVirt_uubar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7),
                        flav
                     );
                     */
                     XSection_SC sc(
                        parameters, m1, m1,
                        std::bind(&Sgluons::sgluons_qqbar_OOg_soft, sgluons, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pqq, std::bind(&Sgluons::matrixSgluonsTree_qqbar_OO, sgluons, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
                           {SplittingKernel::Pqq, std::bind(&Sgluons::matrixSgluonsTree_qqbar_OO, sgluons, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m1,
                        std::bind(&Sgluons::sgluons_qqbar_OOg_hard, sgluons, std::placeholders::_1, std::placeholders::_2),
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
                        std::bind(&Sgluons::matrixSgluonsTree_gg_OO, sgluons, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4), flav,
                        born_precision, born_verbosity
                     );
                     /*
                     XSection_Virt virt(
                        parameters, m1, m2,
                        std::bind(&MRSSM::matrixVirt_uubar_suLsuLdagger, mrssm, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7),
                        flav
                     );
                     */
                     XSection_SC sc(
                        parameters, m1, m1,
                        std::bind(&Sgluons::sgluons_gg_OOg_soft, sgluons, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5),
                        dS, dC,
                        flav,
                        {{
                           {SplittingKernel::Pgg, std::bind(&Sgluons::matrixSgluonsTree_gg_OO, sgluons, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
                           {SplittingKernel::Pgg, std::bind(&Sgluons::matrixSgluonsTree_gg_OO, sgluons, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)}
                        }},
                        sc_precision, sc_verbosity
                     );
                     XSection_HnonC hc(
                        parameters, m1, m1,
                        std::bind(&Sgluons::sgluons_gg_OOg_hard, sgluons, std::placeholders::_1, std::placeholders::_2),
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

   std::cout << std::endl;

   console->info("Elapsed {:.3}", sw);

   // print out time statistics
   auto end = chrono::steady_clock::now();
   std::stringstream s;
   if (end - start > 1h) {
      s << chrono::duration_cast<chrono::hours>(end-start).count() << " hour(s), ";
   }
   if (end - start > 1min) {
      s << chrono::duration_cast<chrono::minutes>(end-start).count() %  60 << " minute(s) and ";
   }
   s << chrono::duration_cast<chrono::seconds>(end-start).count() % 60 << " second(s)";
   spdlog::get("console")->info("Calculation ended after {}", s.str());

   // write results to JSON file
   std::filesystem::path path(card);
   path.replace_extension();
   string json_outputfile_name =
      vm.count("json-outputfile-name")
         ? vm["json-outputfile-name"].as<string>()
         :  path.filename().string();
   int file_index = 0;
   std::string file_suffix = "";
   while (std::ifstream(json_outputfile_name + file_suffix + ".json")) {
      file_index += 1;
      file_suffix = "_" + std::to_string(file_index);
   }
   const std::string fileName = json_outputfile_name + file_suffix + ".json";
   std::ofstream o(fileName);
   o << std::setw(3) << j << std::endl;
   console->info("JSON output written to {}", fileName);
   o.close();

   return 0;
}
