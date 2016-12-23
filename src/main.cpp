#include <chrono>

#include "XSection_Tree.hpp"
#include "XSection_Virt.hpp"
#include "XSection_SC.hpp"
#include "XSection_HnonC.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

using namespace std;

array<double,3> add(array<double,3> x, array<double,3> y) {
   return array<double,3> { x.at(0) + y.at(0), sqrt( pow(x.at(1), 2) + pow(y.at(1), 2) ),
     x.at(2) + y.at(2) };
}

// why do I have to write this?
// why isn't init() enough?
std::array< std::array<double, 2>, 6 > XSection::squark_mass; 
Process *XSection::processID; 
double XSection::muF;
double XSection::muR;
double XSection::mu_r;
double XSection::mu_f;
double XSection::m1;
double XSection::m2;
double XSection::prec_virt;
double XSection::prec_sc;
double XSection::prec_hnc;
double XSection::S;
double XSection::S_sqrt;
double XSection::dC;
double XSection::dS;
boost::property_tree::ptree XSection::pt;
const LHAPDF::PDF* XSection::pdf;

void print( string str, array<double,3> tree, array<double,3> virt, array<double,3> soft, array<double,3> hard) {
   cout << "\nResults for subprocess " << str << '\n';
      cout << scientific;
   //print out LO run statistics
   cout << "---------------------------------------------------------------" << endl;
   cout << setprecision(5);
   cout << setw(12) << "tree:" << setw(13) << tree.at(0) 
         << " +/- " << setprecision(1) << tree.at(1)
         << " fb ( p-value = " << setw(8) << tree.at(2) << " )\n";
   cout << setprecision(5);
   cout << setw(12) << "virtual:" << setw(13) << virt.at(0) << " +/- " 
           << setprecision(1) << virt.at(1) << " fb ( p-value = " 
           << setw(8) << virt.at(2) << " )\n";

      cout << setprecision(5);
      cout << setw(12) << "real (soft):" << setw(13) << soft.at(0) << " +/- " << setprecision(1) << soft.at(1)
           << " fb ( p-value = " << setw(8) << soft.at(2) << " )\n"; 
      cout << setprecision(5);
      cout << setw(12) << "real (hard):" << setw(13) << hard.at(0) << " +/- " << setprecision(1) << hard.at(1)
           << " fb ( p-value = " << setw(8) << hard.at(2) << " )\n";
      cout << "---------------------------------------------------------------" << endl;
      cout << setprecision(5);
      cout << setw(12) << "sum:" << setw(13)
           << tree.at(0) + virt.at(0) + hard.at(0) + soft.at(0) 
           << " +/- " << setprecision(1) << sqrt(pow(tree.at(1),2) + 
           pow(virt.at(1),2) + pow(hard.at(1),2) +
           pow(soft.at(1),2)) << " fb" << endl;   
}

int main(int argc, char* argv[]) {
   
/* invoke programm like 
   "./RSymSQCD MSSM pp_suLsuR NLO 1 1 1" for  NLO calculation  
       with last three numbers giving the desired accuracy of the virtual,
       soft-collinear and hard-noncollinear part, respectively 
   "./RSymSQCD MSSM pp_suLsuR LO" for  LO calculation  */
   
   cout << "\nPlease do take care of used pdf-set. For LO/NLO calculations "
        << "LO/NLO pdf's are NOT used automatically, but need to be "
        << "specified in config.ini!\n" << endl;
        
   boost::property_tree::ptree pt;
   boost::property_tree::ini_parser::read_ini("config.ini", pt);   
   
   // local arrays are not aumatically initialized to 0
   // need to use {}
   array<double,3> temp {}, xsection_tree {}, xsection_virt {}, xsection_SC {}, xsection_HnonC {},
           xsection_tree1 {}, xsection_virt1 {}, xsection_SC1 {}, xsection_HnonC1 {},
           xsection_tree2 {}, xsection_virt2 {}, xsection_SC2 {}, xsection_HnonC2 {},
           xsection_tree3 {}, xsection_virt3 {}, xsection_SC3 {}, xsection_HnonC3 {};
   enum Model {
       MRSSM,
       MSSM, 
       Simplified,
       no_model
   };
   enum Channel {
       pp_OO,
       pp_OsOs,
       pp_suLsuR,  
       pp_suLsuL,
       pp_suLsdR,
       pp_suLsdL,
       pp_suLsuLdagger,
       pp_suLsuRdagger,
       pp_suLsdLdagger,
       pp_suLsdRdagger,
       no_channel
   };  
   
   Model model = no_model;
   Channel channel = no_channel; 
   
   if ( string(argv[1]) == "MRSSM" ) {
      model = MRSSM;
   } else if ( string(argv[1]) == "MSSM" ) {
      model = MSSM;
   } else if ( string(argv[1]) == "Simplified" ) {
      model = Simplified;
   } else {
	   cout << "\n Model not implemented! \n\n";
   }
   
   if ( string(argv[2]) == "pp_OsOs" ) {
      channel = pp_OsOs;  
   } else if ( string(argv[2]) == "pp_suLsuR" ) {
      channel = pp_suLsuR;
   } else if ( string(argv[2]) == "pp_suLsuL" ) {
      channel = pp_suLsuL;   
   } else if ( string(argv[2]) == "pp_suLsdR" ) {
      channel = pp_suLsdR;
   } else if ( string(argv[2]) == "pp_suLsdL" ) {
      channel = pp_suLsdL;
   } else if ( string(argv[2]) == "pp_suLsuLdagger" ) {
      channel = pp_suLsuLdagger;
   } else if ( string(argv[2]) == "pp_suLsuRdagger" ) {
      channel = pp_suLsuRdagger;
   } else if ( string(argv[2]) == "pp_suLsdLdagger" ) {
      channel = pp_suLsdLdagger;
   } else if ( string(argv[2]) == "pp_suLsdRdagger" ) {
      channel = pp_suLsdRdagger;
         } else if ( string(argv[2]) == "pp_OO" ) {
      channel = pp_OO;
   } else {
	   cout << "\n Process not implemented! \n\n";
   } 
     
   auto start = chrono::steady_clock::now();

   if( string( argv[3] ) == "LO" ) {
	  switch(model) {
         case MRSSM:
            switch(channel) {
               case pp_OsOs:
                  {
                  Process process1("sgluons-gg_OO", pt);  
                  XSection::init( &process1, pt, 1, 1, 1 );
                  XSection_Tree tree;
                  temp = tree.integrate();     
                  Process process2("sgluons-qqbar_OO", pt);
                  XSection::init( &process2, pt, 1, 1, 1 );
                  xsection_tree = add(tree.integrate(), temp);   
                  break;
			      }                  
               case pp_suLsuR: 
                  {                                                     // checked with MadGraph and Philip
                  Process process1("MRSSM,uu_suLsuR", pt);
                  XSection::init( &process1, pt, 1, 1, 1 );
                  XSection_Tree tree;
                  xsection_tree = tree.integrate();
                  break;      	
			      }	
			   case pp_suLsdR: 
                  {  
				                                                        // checked with MadGraph and Philip 
                  Process process1("MRSSM,ud_suLsdR", pt);
                  XSection::init( &process1, pt, 1, 1, 1 );
                  XSection_Tree tree;
                  xsection_tree = tree.integrate();
                  break;      	
			      }	

            case pp_suLsuLdagger:
            {
               Process process1("MRSSM,GG_suLsuLdagger", pt);
               XSection::init( &process1, pt, 1, 1, 1 );
               XSection_Tree tree;
               xsection_tree = tree.integrate();
               Process process2("MRSSM,uubar_suLsuLdagger", pt);
               XSection::init( &process2, pt, 1, 1, 1 );
               xsection_tree = add(xsection_tree, tree.integrate());
               Process process3("MRSSM,ddbar_suLsuLdagger", pt);
               XSection::init( &process3, pt, 1, 1, 1 );
               xsection_tree = add(xsection_tree, tree.integrate());
               break; 
            }
			   default:
			      {
			      xsection_tree = {0,0,0};
			      break;
			      }		
            }
            break;
         case MSSM:
            switch(channel) {
			   case pp_suLsuR:  
			      {                                                     // checked with MadGraph 
                  Process process1("MSSM,uu_suLsuR", pt);
                  XSection::init( &process1, pt, 1, 1, 1 );
                  XSection_Tree tree;
                  xsection_tree = tree.integrate();
                  break;
		          }
		       case pp_suLsuL:
		          {
				  Process process1("MSSM,uu_suLsuL", pt);
                  XSection::init( &process1, pt, 1, 1, 1 );
                  XSection_Tree tree;
                  xsection_tree = tree.integrate();
                  break;	  
				  }
			   case pp_suLsdR: 
                  {                                                     // checked with MadGraph 
                  Process process1("MSSM,ud_suLsdR", pt);
                  XSection::init( &process1, pt, 1, 1, 1 );
                  XSection_Tree tree;
                  xsection_tree = add(tree.integrate(), tree.integrate()); // twice as there is ud and du initial state
                  break;      	
			      }
			   case pp_suLsdL: 
                  {                                                     // checked with MadGraph 
                  Process process1("MSSM,ud_suLsdL", pt);
                  XSection::init( &process1, pt, 1, 1, 1 );
                  XSection_Tree tree;
                  xsection_tree = add(tree.integrate(), tree.integrate()); // twice as there is ud and du initial state
                  break;      	
			      }
			   case pp_suLsuLdagger:
			      { 
				  // todo
                  break;
			      }
			   case pp_suLsuRdagger:
			      { 
				  // todo
                  break;
			      }
			   case pp_suLsdLdagger:
			      { 
				  // todo
                  break;
			      }
			   case pp_suLsdRdagger:
			      { 
				  // todo
                  break;
			      }
			   default:
			      {
			      xsection_tree = {0,0,0};
			      break;
			      }		       
            }
         break;								
      }
   }   		 

   else if( string( argv[3] ) == "NLO" ) {
	  switch(model) {
         case MRSSM:
            switch(channel) {
               case pp_OsOs:
                  {
                  // todo   
                  break;
			      }                  
               case pp_suLsuR: 
                  {     
                  Process process1("MRSSM,uu_suLsuR", pt);
	              XSection::init( &process1, pt, pow(10, -atoi(argv[4])), pow(10, -atoi(argv[5])), pow(10, -atoi(argv[6])) );
                  XSection_Tree tree;
                  xsection_tree1 = tree.integrate();      
                  XSection_Virt virt;
                  xsection_virt1 = virt.integrate();
                  XSection_SC sc;
                  xsection_SC1 = sc.integrate();
                  XSection_HnonC hc;
                  xsection_HnonC1 = hc.integrate();
                  print( "uu > suLsuR", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
                  break;      	
			      }	
			   case pp_suLsdR:  // result doubled up, as there is ud and du initial state
                  {                                                     
                  Process process1("MRSSM,ud_suLsdR", pt);
	              XSection::init( &process1, pt, pow(10, -atoi(argv[4])), pow(10, -atoi(argv[5])), pow(10, -atoi(argv[6])) );
                  XSection_Tree tree;
                  xsection_tree = add(tree.integrate(),tree.integrate());      

                  XSection_Virt virt;
                  xsection_virt = virt.integrate();
                  //cout << "Virtual correction = " << xsection_virt.at(0) << endl;
                  XSection_SC sc;
                  xsection_SC = add(sc.integrate(),sc.integrate());
        

          XSection_HnonC hc;
          //xsection_HnonC = hc.integrate();  
      } 
               case pp_suLsuLdagger:
               {
         
         //Process process1("MRSSM,uubar_suLsuLdagger", pt);
         //XSection::init( &process1, pt, pow(10, -atoi(argv[4])), pow(10, -atoi(argv[5])), pow(10, -atoi(argv[6])) );
         
         //XSection_Tree tree;
         //xsection_tree1 = tree.integrate();      
         //XSection_Virt virt;
         //xsection_virt1 = virt.integrate();
         XSection_SC sc;
         //xsection_SC1 = sc.integrate();
         XSection_HnonC hc;
         //xsection_HnonC1 = hc.integrate();
         //print( "uubar > suLsuL*", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
         
         Process process2("MRSSM,ddbar_suLsuLdagger", pt);
         XSection::init( &process2, pt, pow(10, -atoi(argv[4])), pow(10, -atoi(argv[5])), pow(10, -atoi(argv[6])) );  
         
         //xsection_tree2 = tree.integrate();
         //xsection_virt2 = virt.integrate();
         xsection_SC2 = sc.integrate();
         xsection_HnonC2 = hc.integrate();
         print( "ddbar > suLsuL*", xsection_tree2, xsection_virt2, xsection_SC2, xsection_HnonC2);
         
         Process process3("MRSSM,GG_suLsuLdagger", pt);
         XSection::init( &process3, pt, pow(10, -atoi(argv[4])), pow(10, -atoi(argv[5])), pow(10, -atoi(argv[6])) );
      //xsection_tree3 = tree.integrate();
         //xsection_virt3 = virt.integrate();
         //xsection_SC3 = sc.integrate();
         //xsection_HnonC3 = hc.integrate();
         print( "gg > suLsuL*", xsection_tree3, xsection_virt3, xsection_SC3, xsection_HnonC3);
      } 
               default:
         cout << "NLO process not implemented\n";
            }
     case Simplified:
        switch(channel) {
           case pp_OO:
           {
              {     
                  Process process1("Simplified,uubar_OO", pt);
                  XSection::init( &process1, pt, pow(10, -atoi(argv[4])), pow(10, -atoi(argv[5])), pow(10, -atoi(argv[6])) );
                  XSection_Tree tree;
                  //xsection_tree1 = tree.integrate();      
                  XSection_Virt virt;
                  //xsection_virt1 = virt.integrate();
                  XSection_SC sc;
                  xsection_SC1 = sc.integrate();
                  XSection_HnonC hc;
                  xsection_HnonC1 = hc.integrate();
                  print( "qqbar > OO", xsection_tree1, xsection_virt1, xsection_SC1, xsection_HnonC1);
                  break;      	
           }
        }
        }
     }
//                  XSection_HnonC hc;
//                  xsection_HnonC = add(hc.integrate(),hc.integrate());
//                  break;      	
//			      }	
//			   case pp_suLsuLdagger:
//			      { 
//				  // todo
//                  break;
//			      }
//			   case pp_suLsdLdagger:
//			      { 
//				  // todo
//                  break;
//			      }
//			   default:
//			      {
//			      xsection_tree = {0,0,0};
//			      break;
//			      }		
//            }
//            break;
//         case MSSM:
//            switch(channel) {
//			   case pp_suLsuR:  
//			      {                                                      
//                  Process process1("MSSM,uu_suLsuR", pt);
//	              XSection::init( &process1, pt, pow(10, -atoi(argv[4])), pow(10, -atoi(argv[5])), pow(10, -atoi(argv[6])) );
//                  XSection_Tree tree;
//                  xsection_tree = tree.integrate();      
//
//                  XSection_Virt virt;
//                  xsection_virt = virt.integrate();
//                  //cout << "Virtual correction = " << xsection_virt.at(0) << endl;
//                  XSection_SC sc;
//                  xsection_SC = sc.integrate();
//        
//                  XSection_HnonC hc;
//                  xsection_HnonC = hc.integrate();
//                  break;
//		          }
//		       case pp_suLsuL:
//		          {
//				  // todo
//                  break;	  
//				  }
//			   case pp_suLsdR: // result doubled up, as there is ud and du initial state
//                  {          	                                             
//                  Process process1("MSSM,ud_suLsdR", pt);
//	              XSection::init( &process1, pt, pow(10, -atoi(argv[4])), pow(10, -atoi(argv[5])), pow(10, -atoi(argv[6])) );
//                  XSection_Tree tree;
//                  xsection_tree = add(tree.integrate(),tree.integrate());      
//
//                  XSection_Virt virt;
//                  xsection_virt = virt.integrate();
//                  //cout << "Virtual correction = " << xsection_virt.at(0) << endl;
//                  XSection_SC sc;
//                  xsection_SC = add(sc.integrate(),sc.integrate());
//        
//                  XSection_HnonC hc;
//                  xsection_HnonC = add(hc.integrate(),hc.integrate());
//                  break;      	
//			      }
//			   case pp_suLsdL: 
//                  {                                                      
//                  // todo
//                  break;      	
//			      }
//			   case pp_suLsuLdagger:
//			      { 
//				  // todo
//                  break;
//			      }
//			   case pp_suLsuRdagger:
//			      { 
//				  // todo
//                  break;
//			      }
//			   case pp_suLsdLdagger:
//			      { 
//				  // todo
//                  break;
//			      }
//			   case pp_suLsdRdagger:
//			      { 
//				  // todo
//                  break;
//			      }
//			   default:
//			      {
//			      xsection_tree = {0,0,0};
//			      break;
//			      }		       
//            }
//         break;								

 //     }
   } else {
      cout << "Third command line argument must be 'LO' or 'NLO'." << endl;
      return 0;
   }


   auto end = chrono::steady_clock::now();
   
   // print out total run time
   cout << "\nRun summary\n";
   cout << "Time: " << chrono::duration_cast<chrono::minutes>(end-start).count()
        << " minutes\n";
   
   cout << scientific;
   //print out LO run statistics
   cout << "---------------------------------------------------------------" << endl;
   cout << setprecision(5);
   cout << setw(12) << "tree:" << setw(13) << xsection_tree.at(0) 
         << " +/- " << setprecision(1) << xsection_tree.at(1)
         << " fb ( p-value = " << setw(8) << xsection_tree.at(2) << " )\n";
   
   return 0;
}
