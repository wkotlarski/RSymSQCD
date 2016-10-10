#ifndef PROCESS_H_
#define PROCESS_H_

#include <string>

class Process {
    private:
        
              
    public: 
        Process(std::string); 
        double (*matrixelementTree)(double, double, double, double, double, double); // matrix element squared 
        double (*matrixelementVirt)(double, double, double, double, double, double,
                                    double, double, double, double, double, double); // Re[M^1L M^(B star)]
        double f1,f2; // flavours of initial partons
        double k;     // 1/k = average over initial state colors and helicities
        double h;     // h = sum over initial and final state helicities of fermions (_hel = 0 in FormCalc)
};
#endif /* PROCESS_H_ */
