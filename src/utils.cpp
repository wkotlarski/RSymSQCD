#include "utils.hpp"
#include "mathematica_wrapper.hpp"

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



