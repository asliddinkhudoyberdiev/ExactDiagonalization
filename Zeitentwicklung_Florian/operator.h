#include <vector>
#include <complex>
#include <omp.h>

//Methoden, die Operatoren auf den Zustand psi wirken lassen

//beliebiger spin05 Operator (psi = op*psi)
void spin05(std::vector<std::complex<double> >& psi, int d05, int i, std::complex<double> a00, std::complex<double> a01, std::complex<double> a10, std::complex<double> a11);

//beliebiger spin05 Operator (psi = op*psi)
void spin05_tmp(std::vector<std::complex<double> >& psi,std::vector<std::complex<double> >& psi_tmp, int d05, int i, std::complex<double> a00, std::complex<double> a01, std::complex<double> a10, std::complex<double> a11, int flag);

//beliebiger spin05spin05 Operator (i != j   !!!!) (psi = op*psi)
void spin05spin05(std::vector<std::complex<double> >& psi, int d05, int i, int j, std::complex<double> a00, std::complex<double> a01, std::complex<double> a02, std::complex<double> a03, std::complex<double> a10, std::complex<double> a11, std::complex<double> a12, std::complex<double> a13, std::complex<double> a20, std::complex<double> a21, std::complex<double> a22, std::complex<double> a23, std::complex<double> a30, std::complex<double> a31, std::complex<double> a32, std::complex<double> a33);


//Operatoren, welche psi_tmp nutzen, um Observablen auszurechnen (flag entscheidet ob += oder =)

// psi_tmp (+)= k*s_+*psi
void s_plus(std::vector<std::complex<double> >& psi, std::vector<std::complex<double> >& psi_tmp, int d05, int i, std::complex<double> k, int flag);

// psi_tmp (+)= k*s_-*psi
void s_minus(std::vector<std::complex<double> >& psi, std::vector<std::complex<double> >& psi_tmp, int d05, int i, std::complex<double> k, int flag);

void const_eins(std::vector<std::complex<double> >& psi, std::vector<std::complex<double> >& psi_tmp, int d05, int i, double k, int flag);

// psi_tmp (+)= k*s_x*psi
void s_x(std::vector<std::complex<double> >& psi, std::vector<std::complex<double> >& psi_tmp, int d05, int i, double k, int flag);


// psi_tmp (+)= k*s_y*psi
void s_y(std::vector<std::complex<double> >& psi, std::vector<std::complex<double> >& psi_tmp, int d05, int i, double k, int flag);


// psi_tmp (+)= k*s_z*psi
void s_z(std::vector<std::complex<double> >& psi, std::vector<std::complex<double> >& psi_tmp, int d05, int i, double k, int flag);


// psi_tmp (+)= k*s_x*s_x*psi
void s_x_s_x(std::vector<std::complex<double> >& psi, std::vector<std::complex<double> >& psi_tmp, int d05, int i, int j, double k, int flag);


// psi_tmp (+)= k*s_x*s_y*psi
void s_x_s_y(std::vector<std::complex<double> >& psi, std::vector<std::complex<double> >& psi_tmp, int d05, int i, int j, double k, int flag);


// psi_tmp (+)= k*s_x*s_z*psi
void s_x_s_z(std::vector<std::complex<double> >& psi, std::vector<std::complex<double> >& psi_tmp, int d05, int i, int j, double k, int flag);


// psi_tmp (+)= k*s_y*s_y*psi
void s_y_s_y(std::vector<std::complex<double> >& psi, std::vector<std::complex<double> >& psi_tmp, int d05, int i, int j, double k, int flag);


// psi_tmp (+)= k*s_y*s_z*psi
void s_y_s_z(std::vector<std::complex<double> >& psi, std::vector<std::complex<double> >& psi_tmp, int d05, int i, int j, double k, int flag);


// psi_tmp (+)= k*s_z*s_z*psi
void s_z_s_z(std::vector<std::complex<double> >& psi, std::vector<std::complex<double> >& psi_tmp, int d05, int i, int j, double k, int flag);


// psi_tmp = H*psi
void H(std::vector<std::complex<double> >& psi, std::vector<std::complex<double> >& psi_tmp, int d05, std::vector<double> J_ij[11], std::vector<double> h_i[4]);


// psi_tmp = H*psi
void H_mit_const(std::vector<std::complex<double> >& psi, std::vector<std::complex<double> >& psi_tmp, int d05, std::vector<double> J_ij[11], std::vector<double> h_i[4], double C);
