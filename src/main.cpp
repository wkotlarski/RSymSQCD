#include <iostream>

#include "XSection_SC.h"
#include "XSection_HnonC.h"

using namespace std;

int main(int argc, char* argv[]) {
  XSection_HnonC hc;
  hc.show_settings();
  cout << hc.integrate() << endl;
  return 1;
}
