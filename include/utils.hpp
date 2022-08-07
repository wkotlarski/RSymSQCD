#include <nlohmann/json_fwd.hpp>

#include <array>
#include <cmath>
#include <string_view>

struct ChannelResult {
   std::string channel_name;
   std::array<double, 3> b {0., 0., 0.};
   std::array<double, 3> v {0., 0., 0.};
   std::array<double, 3> s {0., 0., 0.};
   std::array<double, 3> h {0., 0., 0.};
};

void print_to_terminal(
   std::string_view,
   std::array<double, 3> const&, std::array<double, 3> const&, std::array<double, 3> const&, std::array<double, 3> const&);
void print_to_terminal(ChannelResult const&);

void print_to_terminal(std::string_view, std::array<double, 3> const&);

void xsec_to_json(nlohmann::json&, std::string_view, std::array<double, 3> const&);
void xsec_to_json(
   nlohmann::json&,
   std::string_view,
   std::array<double, 3> const&, std::array<double, 3> const&, std::array<double, 3> const&, std::array<double, 3> const&
);
void xsec_to_json(nlohmann::json&, ChannelResult const&);

inline std::array<double, 3> operator+(std::array<double, 3> const& x, std::array<double, 3> const& y) {
   return {x.at(0) + y.at(0), std::hypot(x.at(1), y.at(1)), std::max(x.at(2), y.at(2))};
}

inline std::array<double, 3>& operator+=(std::array<double, 3>& x, std::array<double, 3> const& y) {
   x.at(0) += y.at(0);
   x.at(1) = std::hypot(x.at(1), y.at(1));
   x.at(2) = std::max(x.at(2), y.at(2));
   return x;
}
