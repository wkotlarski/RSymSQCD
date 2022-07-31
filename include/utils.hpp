#include <array>
#include <string_view>

void print_to_terminal(
   std::string_view,
   std::array<double, 3> const&, std::array<double, 3> const&, std::array<double, 3> const&, std::array<double, 3> const&);

void print_to_terminal(std::string_view, std::array<double, 3> const&);
