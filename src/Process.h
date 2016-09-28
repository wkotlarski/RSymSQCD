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
};
#endif /* PROCESS_H_ */
