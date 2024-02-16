#include "utils.hpp"
#include "mathematica_wrapper.hpp"

#include "nlohmann/json.hpp"
#include "spdlog/fmt/fmt.h"

#include <boost/math/special_functions/pow.hpp>

#include <cmath>
#include <iostream>
#include <iomanip>

using boost::math::pow;

void print_to_terminal(
   std::string_view str,
   std::array<double, 3> const& tree, std::array<double, 3> const& virt, std::array<double, 3> const& soft, std:: array<double,3> const& hard)
{
   static constexpr int width = 19;
   static constexpr int precision = 5;
   static constexpr std::string_view line = "-------------------------------------------------------------------\n";

   std::cout << "\nResults for subprocess " << str << '\n';
   std::cout << line;
   fmt::print("          Born (B): {: .5e} +/- {:.1e} fb (p-value = {:.1e})\n", tree.at(0), tree.at(1), tree.at(2));
   fmt::print("       Virtual (V): {: .5e} +/- {:.1e} fb (p-value = {:.1e})\n", virt.at(0), virt.at(1), virt.at(2));
   fmt::print("       Real (S+HC): {: .5e} +/- {:.1e} fb (p-value = {:.1e})\n", soft.at(0), soft.at(1), soft.at(2));
   fmt::print("        Real (HnC): {: .5e} +/- {:.1e} fb (p-value = {:.1e})\n", hard.at(0), hard.at(1), hard.at(2));
   std::cout << line;
   fmt::print(
      "sum (B+V+S+HC+HnC): {: .5e} +/- {:.1e} fb (p-value = {:.1e})\n",
      tree.at(0) + virt.at(0) + hard.at(0) + soft.at(0),
      std::sqrt(pow<2>(tree.at(1)) + pow<2>(virt.at(1)) + pow<2>(hard.at(1)) + pow<2>(soft.at(1))),
      std::max({tree.at(2), virt.at(2), soft.at(2), hard.at(2)})
   );
}

void print_to_terminal(std::string_view str, std::array<double, 3> const& tree) {
   std::cout << "\nResults for subprocess " << str << '\n';
   std::cout << std::scientific;
   //print out LO run statistics
   std::cout << "---------------------------------------------------------------\n";
   std::cout << std::setprecision(5);
   std::cout << std::setw(12) << "born:" << std::setw(13) << tree.at(0)
         << " +/- " << std::setprecision(1) << tree.at(1)
         << " fb ( p-value = " << std::setw(8) << tree.at(2) << " )\n";
}

void xsec_to_json(nlohmann::json& j, std::string const& str, std::array<double, 3> const& tree) {
   j["cross sections"][str] = {
      {"tree", {{"res", tree.at(0)}, {"err", tree.at(1)}, {"p-val", tree.at(2)}}}
   };
}

void xsec_to_json(
      nlohmann::json& j,
      std::string const& str,
      std::array<double, 3> const& tree, std::array<double, 3> const& virt, std::array<double, 3> const& soft, std::array<double, 3> const& hard) {

   xsec_to_json(j, str, tree);
   j["cross sections"][str] = {
      {"tree", {{"res", tree.at(0)}, {"err", tree.at(1)}, {"p-val", tree.at(2)}}},
      {"virtual", {{"res", virt.at(0)}, {"err", virt.at(1)}, {"p-val", virt.at(2)}}},
      {"SC", {{"res", soft.at(0)}, {"err", soft.at(1)}, {"p-val", soft.at(2)}}},
      {"HnonC", {{"res", hard.at(0)}, {"err", hard.at(1)}, {"p-val", hard.at(2)}}}
   };
}


void print_to_terminal(ChannelResult const& c) {
   if (c.v || c.s || c.h) {
      print_to_terminal(c.channel_name, c.b, c.v.value_or(std::array<double, 3>{0., 0., 0.}), c.s.value_or(std::array<double, 3>{0., 0., 0.}), c.h.value_or(std::array<double, 3>{0., 0., 0.}));
   }
   else {
      print_to_terminal(c.channel_name, c.b);
   }
}

void xsec_to_json(
      nlohmann::json& j,
      ChannelResult const& c) {
   if (c.v || c.s || c.h) {
      xsec_to_json(j, c.channel_name, c.b, c.v.value_or(std::array<double, 3>{0., 0., 0.}), c.s.value_or(std::array<double, 3>{0., 0., 0.}), c.h.value_or(std::array<double, 3>{0., 0., 0.}));
   }
   else {
      xsec_to_json(j, c.channel_name, c.b);
   }
}
