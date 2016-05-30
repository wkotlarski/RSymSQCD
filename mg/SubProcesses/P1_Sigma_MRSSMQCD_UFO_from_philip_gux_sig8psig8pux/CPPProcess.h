//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.3.3, 2015-10-25
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MG5_Sigma_MRSSMQCD_UFO_from_philip_gux_sig8psig8pux_H
#define MG5_Sigma_MRSSMQCD_UFO_from_philip_gux_sig8psig8pux_H

#include <complex> 
#include <vector> 

#include "Parameters_MRSSMQCD_UFO_from_philip.h"

using namespace std; 

//==========================================================================
// A class for calculating the matrix elements for
// Process: g u~ > sig8p sig8p u~ NP=1 WEIGHTED=3 @1
// Process: g c~ > sig8p sig8p c~ NP=1 WEIGHTED=3 @1
// Process: g d~ > sig8p sig8p d~ NP=1 WEIGHTED=3 @1
// Process: g s~ > sig8p sig8p s~ NP=1 WEIGHTED=3 @1
// Process: g b~ > sig8p sig8p b~ NP=1 WEIGHTED=3 @1
//--------------------------------------------------------------------------

class CPPProcess
{
  public:

    // Constructor.
    CPPProcess() {}

    // Initialize process.
    virtual void initProc(string param_card_name); 

    // Calculate flavour-independent parts of cross section.
    virtual void sigmaKin(); 

    // Evaluate sigmaHat(sHat).
    virtual double sigmaHat(); 

    // Info on the subprocess.
    virtual string name() const {return "g u~ > sig8p sig8p u~ (MRSSMQCD_UFO_from_philip)";}

    virtual int code() const {return 1;}

    const vector<double> & getMasses() const {return mME;}

    // Get and set momenta for matrix element evaluation
    vector < double * > getMomenta(){return p;}
    void setMomenta(vector < double * > & momenta){p = momenta;}
    void setInitial(int inid1, int inid2){id1 = inid1; id2 = inid2;}

    // Get matrix element vector
    const double * getMatrixElements() const {return matrix_element;}

    // Constants for array limits
    static const int ninitial = 2; 
    static const int nexternal = 5; 
    static const int nprocesses = 2; 

  private:

    // Private functions to calculate the matrix element for all subprocesses
    // Calculate wavefunctions
    void calculate_wavefunctions(const int perm[], const int hel[]); 
    static const int nwavefuncs = 13; 
    std::complex<double> w[nwavefuncs][18]; 
    static const int namplitudes = 7; 
    std::complex<double> amp[namplitudes]; 
    double matrix_1_gux_sig8psig8pux(); 

    // Store the matrix element value from sigmaKin
    double matrix_element[nprocesses]; 

    // Color flows, used when selecting color
    double * jamp2[nprocesses]; 

    // Pointer to the model parameters
    Parameters_MRSSMQCD_UFO_from_philip * pars; 

    // vector with external particle masses
    vector<double> mME; 

    // vector with momenta (to be changed each event)
    vector < double * > p; 
    // Initial particle ids
    int id1, id2; 

}; 


#endif  // MG5_Sigma_MRSSMQCD_UFO_from_philip_gux_sig8psig8pux_H
