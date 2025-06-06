#include <vector>
#include <complex>


void initial(int N05, double schrittweite, std::vector<std::complex<double> > trotter_ze[5], std::vector<std::complex<double> > trotter_J_ij_1[18], std::vector<std::complex<double> > trotter_J_ij_2[18], std::vector<double> J_ij[11], std::vector<double> h_i[4]);

void initial_imag(int N05, double schrittweite, std::vector<std::complex<double> > trotter_ze[5], std::vector<std::complex<double> > trotter_J_ij_1[18], std::vector<std::complex<double> > trotter_J_ij_2[18], std::vector<double> J_ij[11], std::vector<double> h_i[4]);
