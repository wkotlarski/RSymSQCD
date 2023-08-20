#include <nlohmann/json_fwd.hpp>

#include <array>
#include <cmath>
#include <optional>
#include <string>

struct ChannelResult {
   std::string channel_name;
   std::array<double, 3> b {0., 0., 0.};
   std::optional<std::array<double, 3>> v;
   std::optional<std::array<double, 3>> s;
   std::optional<std::array<double, 3>> h;
};

void print_to_terminal(ChannelResult const&);
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
