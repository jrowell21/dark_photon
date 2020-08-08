#include "mass_calibration.cpp"
//#include "omega.cpp"
//#include "fit_res.cpp"
//#include "sig_sub.cpp"
//#include "plot_chi2s.cpp"
//#include "resids.cpp"

int main (int argc, char *argv[]) {
  
  if (argc > 1){
    return 0;
  }
  //resids();
  //plot_chi2s();
  mass_calibration("2018");
  //mass_calibration("2018");	
  //omega();
  //sig_sub();
  //fit_res()
  return 0;
}
