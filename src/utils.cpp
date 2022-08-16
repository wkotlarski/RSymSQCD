#include "utils.hpp"
#include "mathematica_wrapper.hpp"

#include "nlohmann/json.hpp"

#include <cmath>
#include <iostream>
#include <iomanip>

void print_to_terminal(
   std::string_view str,
   std::array<double, 3> const& tree, std::array<double, 3> const& virt, std::array<double, 3> const& soft, std:: array<double,3> const& hard)
{
   std::cout << "\nResults for subprocess " << str << '\n';
   std::cout << std::scientific;
   // print out LO run statistics
   std::cout << "---------------------------------------------------------------\n";
   std::cout << std::setprecision(5);
   std::cout << std::setw(12) << "tree:" << std::setw(13) << tree.at(0)
         << " +/- " << std::setprecision(1) << tree.at(1)
         << " fb ( p-value = " << std::setw(8) << tree.at(2) << " )\n";
   std::cout << std::setprecision(5);
   std::cout << std::setw(12) << "virtual:" << std::setw(13) << virt.at(0) << " +/- "
           << std::setprecision(1) << virt.at(1) << " fb ( p-value = "
           << std::setw(8) << virt.at(2) << " )\n";

      std::cout << std::setprecision(5);
      std::cout << std::setw(12) << "real (soft):" << std::setw(13) << soft.at(0) << " +/- " << std::setprecision(1) << soft.at(1)
           << " fb ( p-value = " << std::setw(8) << soft.at(2) << " )\n";
      std::cout << std::setprecision(5);
      std::cout << std::setw(12) << "real (hard):" << std::setw(13) << hard.at(0) << " +/- " << std::setprecision(1) << hard.at(1)
           << " fb ( p-value = " << std::setw(8) << hard.at(2) << " )\n";
      std::cout << "---------------------------------------------------------------\n";
      std::cout << std::setprecision(5);
      std::cout << std::setw(12) << "sum:" << std::setw(13)
           << tree.at(0) + virt.at(0) + hard.at(0) + soft.at(0)
           << " +/- " << std::setprecision(1) << std::sqrt(pow(tree.at(1),2) +
           Sqr(virt.at(1)) + Sqr(hard.at(1)) +
           Sqr(soft.at(1))) << " fb\n";
}

void print_to_terminal(std::string_view str, std::array<double, 3> const& tree) {
   std::cout << "\nResults for subprocess " << str << '\n';
   std::cout << std::scientific;
   //print out LO run statistics
   std::cout << "---------------------------------------------------------------\n";
   std::cout << std::setprecision(5);
   std::cout << std::setw(12) << "tree:" << std::setw(13) << tree.at(0)
         << " +/- " << std::setprecision(1) << tree.at(1)
         << " fb ( p-value = " << std::setw(8) << tree.at(2) << " )\n";
}

void xsec_to_json(nlohmann::json& j, std::string_view str, std::array<double, 3> const& tree) {
   j["cross sections"][str] = {
      {"tree", {{"res", tree.at(0)}, {"err", tree.at(1)}, {"p-val", tree.at(2)}}}
   };
}

void xsec_to_json(
      nlohmann::json& j,
      std::string_view str,
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
