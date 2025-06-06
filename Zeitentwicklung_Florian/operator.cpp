#include "operator.h"

using namespace std;


void spin05(vector<complex<double> >& psi, int d05, int i, complex<double> a00, complex<double> a01, complex<double> a10, complex<double> a11){
	
	#pragma omp parallel 
	{
		#pragma omp for 
		for (int z = 0; z<(d05>>1); z++){							//geht alle "Paare" durch
		
			int dummy = (1<<i) - 1;
			int z0 = dummy & z;
			int z1 = z - z0;
			int index0 = z0 + (z1<<1);
			int index1 = index0 + (1<<i);
			
			complex<double> tmp = a00*psi[index0] + a01*psi[index1];
			psi[index1] = a10*psi[index0] + a11*psi[index1];
			psi[index0] = tmp;
		}
	}
}


void spin05_tmp(vector<complex<double> >& psi, vector<complex<double> >& psi_tmp, int d05, int i, complex<double> a00, complex<double> a01, complex<double> a10, complex<double> a11, int flag){
	
	
	if (flag == 0){
		
		#pragma omp parallel 
		{
			#pragma omp for 
			for (int z = 0; z<(d05>>1); z++){							//geht alle "Paare" durch
		
				int dummy = (1<<i) - 1;
				int z0 = dummy & z;
				int z1 = z - z0;
				int index0 = z0 + (z1<<1);
				int index1 = index0 + (1<<i);
				
				complex<double> tmp = a00*psi[index0] + a01*psi[index1];
				psi_tmp[index1] += a10*psi[index0] + a11*psi[index1];
				psi_tmp[index0] += tmp;
			}
		}
	}
	else{
		
		#pragma omp parallel 
		{
			#pragma omp for 
			for (int z = 0; z<(d05>>1); z++){							//geht alle "Paare" durch
		
				int dummy = (1<<i) - 1;
				int z0 = dummy & z;
				int z1 = z - z0;
				int index0 = z0 + (z1<<1);
				int index1 = index0 + (1<<i);
			
				complex<double> tmp = a00*psi[index0] + a01*psi[index1];
				psi_tmp[index1] = a10*psi[index0] + a11*psi[index1];
				psi_tmp[index0] = tmp;
			}
		}
	}
}


void spin05spin05(vector<complex<double> >& psi, int d05, int i, int j, complex<double> a00, complex<double> a01, complex<double> a02, complex<double> a03, complex<double> a10, complex<double> a11, complex<double> a12, complex<double> a13, complex<double> a20, complex<double> a21, complex<double> a22, complex<double> a23, complex<double> a30, complex<double> a31, complex<double> a32, complex<double> a33){
	
	//i<j erzwingen
	if (j < i){
		int k2 = j;
		j = i;
		i = k2;
	}
	
	#pragma omp parallel 
	{
		#pragma omp for 
		for (int z = 0; z<(d05>>2); z++){							//geht alle "Quadrupel" durch
		
			int dummy1 = (1<<i) - 1;
			int z0 = dummy1 & z;
			int z12 = z - z0;
			int dummy2 = (1<<(j-1)) - (1<<i);						// = ((1<<(j-1)) - 1) - ((1<<i) - 1);
			int z1 = dummy2 & z12;
			int z2 = z - z0 - z1;
			
			int index0 = z0 + (z1<<1) + (z2<<2);
			int index1 = index0 + (1<<i);
			int index2 = index0 + (1<<j);
			int index3 = index1 + (1<<j);
			
			complex<double> tmp0 = a00*psi[index0] + a01*psi[index1] + a02*psi[index2] + a03*psi[index3];
			complex<double> tmp1 = a10*psi[index0] + a11*psi[index1] + a12*psi[index2] + a13*psi[index3];
			complex<double> tmp2 = a20*psi[index0] + a21*psi[index1] + a22*psi[index2] + a23*psi[index3];
			psi[index3] = a30*psi[index0] + a31*psi[index1] + a32*psi[index2] + a33*psi[index3];
			psi[index0] = tmp0;
			psi[index1] = tmp1;
			psi[index2] = tmp2;
		}
	}
}


void s_plus(vector<complex<double> >& psi, vector<complex<double> >& psi_tmp, int d05, int i, complex<double> k, int flag){
	
	if (flag == 0){
		
		#pragma omp parallel 
		{
			#pragma omp for 
			for (int z = 0; z<(d05>>1); z++){							//geht alle "Paare" durch
		
				int dummy = (1<<i) - 1;
				int z0 = dummy & z;
				int z1 = z - z0;
				int index0 = z0 + (z1<<1);
				int index1 = index0 + (1<<i);
			
				psi_tmp[index0] += k*1.0*psi[index1];
				psi_tmp[index1] += 0;
			}
		}
	}
	else{
		
		#pragma omp parallel 
		{
			#pragma omp for 
			for (int z = 0; z<(d05>>1); z++){							//geht alle "Paare" durch
		
				int dummy = (1<<i) - 1;
				int z0 = dummy & z;
				int z1 = z - z0;
				int index0 = z0 + (z1<<1);
				int index1 = index0 + (1<<i);
			
				psi_tmp[index0] = k*1.0*psi[index1];
				psi_tmp[index1] = 0;
			}
		}
	}
}


void s_minus(vector<complex<double> >& psi, vector<complex<double> >& psi_tmp, int d05, int i, complex<double> k, int flag){
	
	if (flag == 0){
		
		#pragma omp parallel 
		{
			#pragma omp for 
			for (int z = 0; z<(d05>>1); z++){							//geht alle "Paare" durch
		
				int dummy = (1<<i) - 1;
				int z0 = dummy & z;
				int z1 = z - z0;
				int index0 = z0 + (z1<<1);
				int index1 = index0 + (1<<i);
			
				psi_tmp[index0] += 0;
				psi_tmp[index1] += k*1.0*psi[index0];
			}
		}
	}
	else{
		
		#pragma omp parallel 
		{
			#pragma omp for 
			for (int z = 0; z<(d05>>1); z++){							//geht alle "Paare" durch
		
				int dummy = (1<<i) - 1;
				int z0 = dummy & z;
				int z1 = z - z0;
				int index0 = z0 + (z1<<1);
				int index1 = index0 + (1<<i);
			
				psi_tmp[index0] = 0;
				psi_tmp[index1] = k*1.0*psi[index0];
			}
		}
	}
}


void const_eins(vector<complex<double> >& psi, vector<complex<double> >& psi_tmp, int d05, int i, double k, int flag){
	
	if (flag == 0){
		
		#pragma omp parallel 
		{
			#pragma omp for 
			for (int z = 0; z<(d05>>1); z++){							//geht alle "Paare" durch
		
				int dummy = (1<<i) - 1;
				int z0 = dummy & z;
				int z1 = z - z0;
				int index0 = z0 + (z1<<1);
				int index1 = index0 + (1<<i);
			
				psi_tmp[index0] += k*psi[index0];
				psi_tmp[index1] += k*psi[index1];
			}
		}
	}
	else{
		
		#pragma omp parallel 
		{
			#pragma omp for 
			for (int z = 0; z<(d05>>1); z++){							//geht alle "Paare" durch
		
				int dummy = (1<<i) - 1;
				int z0 = dummy & z;
				int z1 = z - z0;
				int index0 = z0 + (z1<<1);
				int index1 = index0 + (1<<i);
			
				psi_tmp[index0] = k*psi[index0];
				psi_tmp[index1] = k*psi[index1];
			}
		}
	}
}


void s_x(vector<complex<double> >& psi, vector<complex<double> >& psi_tmp, int d05, int i, double k, int flag){
	
	if (flag == 0){
		
		#pragma omp parallel 
		{
			#pragma omp for 
			for (int z = 0; z<(d05>>1); z++){							//geht alle "Paare" durch
		
				int dummy = (1<<i) - 1;
				int z0 = dummy & z;
				int z1 = z - z0;
				int index0 = z0 + (z1<<1);
				int index1 = index0 + (1<<i);
			
				psi_tmp[index0] += k*0.5*psi[index1];
				psi_tmp[index1] += k*0.5*psi[index0];
			}
		}
	}
	else{
		
		#pragma omp parallel 
		{
			#pragma omp for 
			for (int z = 0; z<(d05>>1); z++){							//geht alle "Paare" durch
		
				int dummy = (1<<i) - 1;
				int z0 = dummy & z;
				int z1 = z - z0;
				int index0 = z0 + (z1<<1);
				int index1 = index0 + (1<<i);
			
				psi_tmp[index0] = k*0.5*psi[index1];
				psi_tmp[index1] = k*0.5*psi[index0];
			}
		}
	}
}


void s_y(vector<complex<double> >& psi, vector<complex<double> >& psi_tmp, int d05, int i, double k, int flag){
	
	if (flag == 0){
	
		#pragma omp parallel 
		{
			#pragma omp for 
			for (int z = 0; z<(d05>>1); z++){							//geht alle "Paare" durch
		
				int dummy = (1<<i) - 1;
				int z0 = dummy & z;
				int z1 = z - z0;
				int index0 = z0 + (z1<<1);
				int index1 = index0 + (1<<i);
			
				psi_tmp[index0] += (complex<double> (0,-0.5))*k*psi[index1];
				psi_tmp[index1] += (complex<double> (0,0.5))*k*psi[index0];
			}
		}
	}
	else{
		
		#pragma omp parallel 
		{
			#pragma omp for 
			for (int z = 0; z<(d05>>1); z++){							//geht alle "Paare" durch
		
				int dummy = (1<<i) - 1;
				int z0 = dummy & z;
				int z1 = z - z0;
				int index0 = z0 + (z1<<1);
				int index1 = index0 + (1<<i);
			
				psi_tmp[index0] = (complex<double> (0,-0.5))*k*psi[index1];
				psi_tmp[index1] = (complex<double> (0,0.5))*k*psi[index0];
			}
		}
	}
}


void s_z(vector<complex<double> >& psi, vector<complex<double> >& psi_tmp, int d05, int i, double k, int flag){
	
	if (flag == 0){
		
		#pragma omp parallel 
		{
			#pragma omp for 
			for (int z = 0; z<(d05>>1); z++){							//geht alle "Paare" durch
		
				int dummy = (1<<i) - 1;
				int z0 = dummy & z;
				int z1 = z - z0;
				int index0 = z0 + (z1<<1);
				int index1 = index0 + (1<<i);
			
				psi_tmp[index0] += 0.5*k*psi[index0];
				psi_tmp[index1] += -0.5*k*psi[index1];
			}
		}
	}
	else{
		
		#pragma omp parallel 
		{
			#pragma omp for 
			for (int z = 0; z<(d05>>1); z++){							//geht alle "Paare" durch
		
				int dummy = (1<<i) - 1;
				int z0 = dummy & z;
				int z1 = z - z0;
				int index0 = z0 + (z1<<1);
				int index1 = index0 + (1<<i);
			
				psi_tmp[index0] = 0.5*k*psi[index0];
				psi_tmp[index1] = -0.5*k*psi[index1];
			}
		}
	}
}


//Im Folgenden: s_i*s_j mit i<j !!! 
void s_x_s_x(vector<complex<double> >& psi, vector<complex<double> >& psi_tmp, int d05, int i, int j, double k, int flag){
	
	//i<j erzwingen
	if (j < i){
		int k2 = j;
		j = i;
		i = k2;
	}
	
	if (flag == 0){
		
		#pragma omp parallel 
		{
			#pragma omp for 
			for (int z = 0; z<(d05>>2); z++){							//geht alle "Quadrupel" durch
		
				int dummy1 = (1<<i) - 1;
				int z0 = dummy1 & z;
				int z12 = z - z0;
				int dummy2 = (1<<(j-1)) - (1<<i);						// = ((1<<(j-1)) - 1) - ((1<<i) - 1);
				int z1 = dummy2 & z12;
				int z2 = z - z0 - z1;
			
				int index0 = z0 + (z1<<1) + (z2<<2);
				int index1 = index0 + (1<<i);
				int index2 = index0 + (1<<j);
				int index3 = index1 + (1<<j);
	
				psi_tmp[index0] += k*0.25*psi[index3];
				psi_tmp[index1] += k*0.25*psi[index2];
				psi_tmp[index2] += k*0.25*psi[index1];
				psi_tmp[index3] += k*0.25*psi[index0];
			}
		}
	}
	else{
		
		#pragma omp parallel 
		{
			#pragma omp for 
			for (int z = 0; z<(d05>>2); z++){							//geht alle "Quadrupel" durch
		
				int dummy1 = (1<<i) - 1;
				int z0 = dummy1 & z;
				int z12 = z - z0;
				int dummy2 = (1<<(j-1)) - (1<<i);						// = ((1<<(j-1)) - 1) - ((1<<i) - 1);
				int z1 = dummy2 & z12;
				int z2 = z - z0 - z1;
			
				int index0 = z0 + (z1<<1) + (z2<<2);
				int index1 = index0 + (1<<i);
				int index2 = index0 + (1<<j);
				int index3 = index1 + (1<<j);
	
				psi_tmp[index0] = k*0.25*psi[index3];
				psi_tmp[index1] = k*0.25*psi[index2];
				psi_tmp[index2] = k*0.25*psi[index1];
				psi_tmp[index3] = k*0.25*psi[index0];
			}
		}
	}
}


void s_x_s_y(vector<complex<double> >& psi, vector<complex<double> >& psi_tmp, int d05, int i, int j, double k, int flag){
	
	if (i < j){
		
		if (flag == 0){
		
			#pragma omp parallel 
			{
				#pragma omp for 
				for (int z = 0; z<(d05>>2); z++){							//geht alle "Quadrupel" durch
		
					int dummy1 = (1<<i) - 1;
					int z0 = dummy1 & z;
					int z12 = z - z0;
					int dummy2 = (1<<(j-1)) - (1<<i);						// = ((1<<(j-1)) - 1) - ((1<<i) - 1);
					int z1 = dummy2 & z12;
					int z2 = z - z0 - z1;
			
					int index0 = z0 + (z1<<1) + (z2<<2);
					int index1 = index0 + (1<<i);
					int index2 = index0 + (1<<j);
					int index3 = index1 + (1<<j);
	
					psi_tmp[index0] += (complex<double> (0,-0.25))*k*psi[index3];
					psi_tmp[index1] += (complex<double> (0,-0.25))*k*psi[index2];
					psi_tmp[index2] += (complex<double> (0,0.25))*k*psi[index1];
					psi_tmp[index3] += (complex<double> (0,0.25))*k*psi[index0];
				}
			}
		}
		else{
		
			#pragma omp parallel 
			{
				#pragma omp for 
				for (int z = 0; z<(d05>>2); z++){							//geht alle "Quadrupel" durch
		
					int dummy1 = (1<<i) - 1;
					int z0 = dummy1 & z;
					int z12 = z - z0;
					int dummy2 = (1<<(j-1)) - (1<<i);						// = ((1<<(j-1)) - 1) - ((1<<i) - 1);
					int z1 = dummy2 & z12;
					int z2 = z - z0 - z1;
				
					int index0 = z0 + (z1<<1) + (z2<<2);
					int index1 = index0 + (1<<i);
					int index2 = index0 + (1<<j);
					int index3 = index1 + (1<<j);
	
					psi_tmp[index0] = (complex<double> (0,-0.25))*k*psi[index3];
					psi_tmp[index1] = (complex<double> (0,-0.25))*k*psi[index2];
					psi_tmp[index2] = (complex<double> (0,0.25))*k*psi[index1];
					psi_tmp[index3] = (complex<double> (0,0.25))*k*psi[index0];
				}
			}
		}
	}
	
	if (j < i){
		
		int k2 = j;
		j = i;
		i = k2;
		
		if (flag == 0){
		
			#pragma omp parallel 
			{
				#pragma omp for 
				for (int z = 0; z<(d05>>2); z++){							//geht alle "Quadrupel" durch
		
					int dummy1 = (1<<i) - 1;
					int z0 = dummy1 & z;
					int z12 = z - z0;
					int dummy2 = (1<<(j-1)) - (1<<i);						// = ((1<<(j-1)) - 1) - ((1<<i) - 1);
					int z1 = dummy2 & z12;
					int z2 = z - z0 - z1;
			
					int index0 = z0 + (z1<<1) + (z2<<2);
					int index1 = index0 + (1<<i);
					int index2 = index0 + (1<<j);
					int index3 = index1 + (1<<j);
	
					psi_tmp[index0] += (complex<double> (0,-0.25))*k*psi[index3];
					psi_tmp[index1] += (complex<double> (0,0.25))*k*psi[index2];
					psi_tmp[index2] += (complex<double> (0,-0.25))*k*psi[index1];
					psi_tmp[index3] += (complex<double> (0,0.25))*k*psi[index0];
				}
			}
		}
		else{
		
			#pragma omp parallel 
			{
				#pragma omp for 
				for (int z = 0; z<(d05>>2); z++){							//geht alle "Quadrupel" durch
		
					int dummy1 = (1<<i) - 1;
					int z0 = dummy1 & z;
					int z12 = z - z0;
					int dummy2 = (1<<(j-1)) - (1<<i);						// = ((1<<(j-1)) - 1) - ((1<<i) - 1);
					int z1 = dummy2 & z12;
					int z2 = z - z0 - z1;
			
					int index0 = z0 + (z1<<1) + (z2<<2);
					int index1 = index0 + (1<<i);
					int index2 = index0 + (1<<j);
					int index3 = index1 + (1<<j);
	
					psi_tmp[index0] = (complex<double> (0,-0.25))*k*psi[index3];
					psi_tmp[index1] = (complex<double> (0,0.25))*k*psi[index2];
					psi_tmp[index2] = (complex<double> (0,-0.25))*k*psi[index1];
					psi_tmp[index3] = (complex<double> (0,0.25))*k*psi[index0];
				}
			}
		}
	
	}	
}


void s_x_s_z(vector<complex<double> >& psi, vector<complex<double> >& psi_tmp, int d05, int i, int j, double k, int flag){
	
	if (i < j){
	
		if (flag == 0){
		
			#pragma omp parallel 
			{
				#pragma omp for 
				for (int z = 0; z<(d05>>2); z++){							//geht alle "Quadrupel" durch
		
					int dummy1 = (1<<i) - 1;
					int z0 = dummy1 & z;
					int z12 = z - z0;
					int dummy2 = (1<<(j-1)) - (1<<i);						// = ((1<<(j-1)) - 1) - ((1<<i) - 1);
					int z1 = dummy2 & z12;
					int z2 = z - z0 - z1;
			
					int index0 = z0 + (z1<<1) + (z2<<2);
					int index1 = index0 + (1<<i);
					int index2 = index0 + (1<<j);
					int index3 = index1 + (1<<j);
	
					psi_tmp[index0] += 0.25*k*psi[index1];
					psi_tmp[index1] += 0.25*k*psi[index0];
					psi_tmp[index2] += -0.25*k*psi[index3];
					psi_tmp[index3] += -0.25*k*psi[index2];
				}
			}
		}
		else{
		
			#pragma omp parallel 
			{
				#pragma omp for 
				for (int z = 0; z<(d05>>2); z++){							//geht alle "Quadrupel" durch
		
					int dummy1 = (1<<i) - 1;
					int z0 = dummy1 & z;
					int z12 = z - z0;
					int dummy2 = (1<<(j-1)) - (1<<i);						// = ((1<<(j-1)) - 1) - ((1<<i) - 1);
					int z1 = dummy2 & z12;
					int z2 = z - z0 - z1;
			
					int index0 = z0 + (z1<<1) + (z2<<2);
					int index1 = index0 + (1<<i);
					int index2 = index0 + (1<<j);
					int index3 = index1 + (1<<j);
	
					psi_tmp[index0] = 0.25*k*psi[index1];
					psi_tmp[index1] = 0.25*k*psi[index0];
					psi_tmp[index2] = -0.25*k*psi[index3];
					psi_tmp[index3] = -0.25*k*psi[index2];
				}
			}
		}
	}
	
	if (j < i){
		
		int k2 = j;
		j = i;
		i = k2;
		
		if (flag == 0){
		
			#pragma omp parallel 
			{
				#pragma omp for 
				for (int z = 0; z<(d05>>2); z++){							//geht alle "Quadrupel" durch
		
					int dummy1 = (1<<i) - 1;
					int z0 = dummy1 & z;
					int z12 = z - z0;
					int dummy2 = (1<<(j-1)) - (1<<i);						// = ((1<<(j-1)) - 1) - ((1<<i) - 1);
					int z1 = dummy2 & z12;
					int z2 = z - z0 - z1;
			
					int index0 = z0 + (z1<<1) + (z2<<2);
					int index1 = index0 + (1<<i);
					int index2 = index0 + (1<<j);
					int index3 = index1 + (1<<j);
	
					psi_tmp[index0] += 0.25*k*psi[index2];
					psi_tmp[index1] += -0.25*k*psi[index3];
					psi_tmp[index2] += 0.25*k*psi[index0];
					psi_tmp[index3] += -0.25*k*psi[index1];
				}
			}
		}
		else{
		
			#pragma omp parallel 
			{
				#pragma omp for 
				for (int z = 0; z<(d05>>2); z++){							//geht alle "Quadrupel" durch
		
					int dummy1 = (1<<i) - 1;
					int z0 = dummy1 & z;
					int z12 = z - z0;
					int dummy2 = (1<<(j-1)) - (1<<i);						// = ((1<<(j-1)) - 1) - ((1<<i) - 1);
					int z1 = dummy2 & z12;
					int z2 = z - z0 - z1;
			
					int index0 = z0 + (z1<<1) + (z2<<2);
					int index1 = index0 + (1<<i);
					int index2 = index0 + (1<<j);
					int index3 = index1 + (1<<j);
	
					psi_tmp[index0] = 0.25*k*psi[index2];
					psi_tmp[index1] = -0.25*k*psi[index3];
					psi_tmp[index2] = 0.25*k*psi[index0];
					psi_tmp[index3] = -0.25*k*psi[index1];
				}
			}
		}
	}
}



void s_y_s_y(vector<complex<double> >& psi, vector<complex<double> >& psi_tmp, int d05, int i, int j, double k, int flag){
	
	//i<j erzwingen
	if (j < i){
		int k2 = j;
		j = i;
		i = k2;
	}
	
	
	if (flag == 0){
		
		#pragma omp parallel 
		{
			#pragma omp for 
			for (int z = 0; z<(d05>>2); z++){							//geht alle "Quadrupel" durch
		
				int dummy1 = (1<<i) - 1;
				int z0 = dummy1 & z;
				int z12 = z - z0;
				int dummy2 = (1<<(j-1)) - (1<<i);						// = ((1<<(j-1)) - 1) - ((1<<i) - 1);
				int z1 = dummy2 & z12;
				int z2 = z - z0 - z1;
			
				int index0 = z0 + (z1<<1) + (z2<<2);
				int index1 = index0 + (1<<i);
				int index2 = index0 + (1<<j);
				int index3 = index1 + (1<<j);
	
				psi_tmp[index0] += (-1.0)*k*0.25*psi[index3];
				psi_tmp[index1] += k*0.25*psi[index2];
				psi_tmp[index2] += k*0.25*psi[index1];
				psi_tmp[index3] += (-1.0)*k*0.25*psi[index0];
			}
		}
	}
	else{
		
		#pragma omp parallel 
		{
			#pragma omp for 
			for (int z = 0; z<(d05>>2); z++){							//geht alle "Quadrupel" durch
		
				int dummy1 = (1<<i) - 1;
				int z0 = dummy1 & z;
				int z12 = z - z0;
				int dummy2 = (1<<(j-1)) - (1<<i);						// = ((1<<(j-1)) - 1) - ((1<<i) - 1);
				int z1 = dummy2 & z12;
				int z2 = z - z0 - z1;
			
				int index0 = z0 + (z1<<1) + (z2<<2);
				int index1 = index0 + (1<<i);
				int index2 = index0 + (1<<j);
				int index3 = index1 + (1<<j);
	
				psi_tmp[index0] = (-1.0)*k*0.25*psi[index3];
				psi_tmp[index1] = k*0.25*psi[index2];
				psi_tmp[index2] = k*0.25*psi[index1];
				psi_tmp[index3] = (-1.0)*k*0.25*psi[index0];
			}
		}
	}
}


void s_y_s_z(vector<complex<double> >& psi, vector<complex<double> >& psi_tmp, int d05, int i, int j, double k, int flag){
	
	if (i < j){
		
		if (flag == 0){
		
			#pragma omp parallel 
			{
				#pragma omp for 
				for (int z = 0; z<(d05>>2); z++){							//geht alle "Quadrupel" durch
		
					int dummy1 = (1<<i) - 1;
					int z0 = dummy1 & z;
					int z12 = z - z0;
					int dummy2 = (1<<(j-1)) - (1<<i);						// = ((1<<(j-1)) - 1) - ((1<<i) - 1);
					int z1 = dummy2 & z12;
					int z2 = z - z0 - z1;
			
					int index0 = z0 + (z1<<1) + (z2<<2);
					int index1 = index0 + (1<<i);
					int index2 = index0 + (1<<j);
					int index3 = index1 + (1<<j);
	
					psi_tmp[index0] += (complex<double> (0,-0.25))*k*psi[index1];
					psi_tmp[index1] += (complex<double> (0,0.25))*k*psi[index0];
					psi_tmp[index2] += (complex<double> (0,0.25))*k*psi[index3];
					psi_tmp[index3] += (complex<double> (0,-0.25))*k*psi[index2];
				}
			}
		}
		else{
		
			#pragma omp parallel 
			{
				#pragma omp for 
				for (int z = 0; z<(d05>>2); z++){							//geht alle "Quadrupel" durch
		
					int dummy1 = (1<<i) - 1;
					int z0 = dummy1 & z;
					int z12 = z - z0;
					int dummy2 = (1<<(j-1)) - (1<<i);						// = ((1<<(j-1)) - 1) - ((1<<i) - 1);
					int z1 = dummy2 & z12;
					int z2 = z - z0 - z1;
			
					int index0 = z0 + (z1<<1) + (z2<<2);
					int index1 = index0 + (1<<i);
					int index2 = index0 + (1<<j);
					int index3 = index1 + (1<<j);
	
					psi_tmp[index0] = (complex<double> (0,-0.25))*k*psi[index1];
					psi_tmp[index1] = (complex<double> (0,0.25))*k*psi[index0];
					psi_tmp[index2] = (complex<double> (0,0.25))*k*psi[index3];
					psi_tmp[index3] = (complex<double> (0,-0.25))*k*psi[index2];
				}
			}
		}
	}
	
	if (j < i){
		
		int k2 = j;
		j = i;
		i = k2;
		
			if (flag == 0){
		
			#pragma omp parallel 
			{
				#pragma omp for 
				for (int z = 0; z<(d05>>2); z++){							//geht alle "Quadrupel" durch
		
					int dummy1 = (1<<i) - 1;
					int z0 = dummy1 & z;
					int z12 = z - z0;
					int dummy2 = (1<<(j-1)) - (1<<i);						// = ((1<<(j-1)) - 1) - ((1<<i) - 1);
					int z1 = dummy2 & z12;
					int z2 = z - z0 - z1;
			
					int index0 = z0 + (z1<<1) + (z2<<2);
					int index1 = index0 + (1<<i);
					int index2 = index0 + (1<<j);
					int index3 = index1 + (1<<j);
	
					psi_tmp[index0] += (complex<double> (0,-0.25))*k*psi[index2];
					psi_tmp[index1] += (complex<double> (0,0.25))*k*psi[index3];
					psi_tmp[index2] += (complex<double> (0,0.25))*k*psi[index0];
					psi_tmp[index3] += (complex<double> (0,-0.25))*k*psi[index1];
				}
			}
		}
		else{
		
			#pragma omp parallel 
			{
				#pragma omp for 
				for (int z = 0; z<(d05>>2); z++){							//geht alle "Quadrupel" durch
		
					int dummy1 = (1<<i) - 1;
					int z0 = dummy1 & z;
					int z12 = z - z0;
					int dummy2 = (1<<(j-1)) - (1<<i);						// = ((1<<(j-1)) - 1) - ((1<<i) - 1);
					int z1 = dummy2 & z12;
					int z2 = z - z0 - z1;
			
					int index0 = z0 + (z1<<1) + (z2<<2);
					int index1 = index0 + (1<<i);
					int index2 = index0 + (1<<j);
					int index3 = index1 + (1<<j);
	
					psi_tmp[index0] = (complex<double> (0,-0.25))*k*psi[index2];
					psi_tmp[index1] = (complex<double> (0,0.25))*k*psi[index3];
					psi_tmp[index2] = (complex<double> (0,0.25))*k*psi[index0];
					psi_tmp[index3] = (complex<double> (0,-0.25))*k*psi[index1];
				}
			}
		}
	}
}




void s_z_s_z(vector<complex<double> >& psi, vector<complex<double> >& psi_tmp, int d05, int i, int j, double k, int flag){
	
	//i<j erzwingen
	if (j < i){
		int k2 = j;
		j = i;
		i = k2;
	}
	
	if (flag == 0){
		
		#pragma omp parallel 
		{
			#pragma omp for 
			for (int z = 0; z<(d05>>2); z++){							//geht alle "Quadrupel" durch
		
				int dummy1 = (1<<i) - 1;
				int z0 = dummy1 & z;
				int z12 = z - z0;
				int dummy2 = (1<<(j-1)) - (1<<i);						// = ((1<<(j-1)) - 1) - ((1<<i) - 1);
				int z1 = dummy2 & z12;
				int z2 = z - z0 - z1;
			
				int index0 = z0 + (z1<<1) + (z2<<2);
				int index1 = index0 + (1<<i);
				int index2 = index0 + (1<<j);
				int index3 = index1 + (1<<j);
	
				psi_tmp[index0] += k*0.25*psi[index0];
				psi_tmp[index1] += (-1.0)*k*0.25*psi[index1];
				psi_tmp[index2] += (-1.0)*k*0.25*psi[index2];
				psi_tmp[index3] += k*0.25*psi[index3];
			}
		}
	}
	else{
		
		#pragma omp parallel 
		{
			#pragma omp for 
			for (int z = 0; z<(d05>>2); z++){							//geht alle "Quadrupel" durch
		
				int dummy1 = (1<<i) - 1;
				int z0 = dummy1 & z;
				int z12 = z - z0;
				int dummy2 = (1<<(j-1)) - (1<<i);						// = ((1<<(j-1)) - 1) - ((1<<i) - 1);
				int z1 = dummy2 & z12;
				int z2 = z - z0 - z1;
			
				int index0 = z0 + (z1<<1) + (z2<<2);
				int index1 = index0 + (1<<i);
				int index2 = index0 + (1<<j);
				int index3 = index1 + (1<<j);
	
				psi_tmp[index0] = k*0.25*psi[index0];
				psi_tmp[index1] = (-1.0)*k*0.25*psi[index1];
				psi_tmp[index2] = (-1.0)*k*0.25*psi[index2];
				psi_tmp[index3] = k*0.25*psi[index3];
			}
		}
	}
}


//Hamiltonoperator H = - sum(s_i J s_j) - sum(h_i s_i) + C
void H(vector<complex<double> >& psi, std::vector<std::complex<double> >& psi_tmp, int d05, vector<double> J_ij[11], vector<double> h_i[4]){
	
	for (int z = 0; z < d05; z++){
		psi_tmp[z] = complex<double> (0,0);
	}
	
	for(int z = 0; z < h_i[0].size(); z++){
		
		if(h_i[1][z] != 0){
			s_x(psi, psi_tmp, d05,(int) h_i[0][z],(-1.0)*h_i[1][z],0);		//psi_tmp += s_i_x*h_i_x*psi		
		}
		if(h_i[2][z] != 0){
			s_y(psi, psi_tmp, d05,(int) h_i[0][z],(-1.0)*h_i[2][z],0);		//psi_tmp += s_i_y*h_i_y*psi		
		}
		if(h_i[3][z] != 0){
			s_z(psi, psi_tmp, d05,(int) h_i[0][z],(-1.0)*h_i[3][z],0);		//psi_tmp += s_i_z*h_i_z*psi		
		}	
	}
		
		
	for(int z = 0; z < J_ij[0].size(); z++){
	
		int i = (int) J_ij[0][z];
		int j = (int) J_ij[1][z];
		
		//Die Folgenden J_ab Konstanten beziehen sich auf s_i^a*s_j^b mit i<j !!
		double Jxx = J_ij[2][z];
		double Jyy = J_ij[6][z];
		double Jzz = J_ij[10][z];
		double Jxy = 0;
		double Jyx = 0;
		double Jxz = 0;
		double Jzx = 0;
		double Jyz = 0;
		double Jzy = 0;
		
		if (i < j){
			Jxy = J_ij[3][z];
			Jyx = J_ij[5][z];
			Jxz = J_ij[4][z];
			Jzx = J_ij[8][z];
			Jyz = J_ij[7][z];
			Jzy = J_ij[9][z];
		}
		else{
			Jxy = J_ij[5][z];
			Jyx = J_ij[3][z];
			Jxz = J_ij[8][z];
			Jzx = J_ij[4][z];
			Jyz = J_ij[9][z];
			Jzy = J_ij[7][z];
		}
		
		if(Jxx != 0){
			s_x_s_x(psi, psi_tmp, d05, i, j, (-1.0)*Jxx,0);
		}
		if(Jxy != 0){
			s_x_s_y(psi, psi_tmp, d05, i, j, (-1.0)*Jxy,0);
		}
		if(Jyx != 0){
			s_x_s_y(psi, psi_tmp, d05, j, i, (-1.0)*Jyx,0);
		}
		if(Jxz != 0){
			s_x_s_z(psi, psi_tmp, d05, i, j, (-1.0)*Jxz,0);
		}
		if(Jzx != 0){
			s_x_s_z(psi, psi_tmp, d05, j, i, (-1.0)*Jzx,0);
		}
		if(Jyy != 0){
			s_y_s_y(psi, psi_tmp, d05, i, j, (-1.0)*Jyy,0);
		}
		if(Jyz != 0){
			s_y_s_z(psi, psi_tmp, d05, i, j, (-1.0)*Jyz,0);
		}
		if(Jzy != 0){
			s_y_s_z(psi, psi_tmp, d05, j, i, (-1.0)*Jzy,0);
		}
		if(Jzz != 0){
			s_z_s_z(psi, psi_tmp, d05, i, j, (-1.0)*Jzz,0);
		}
	}
}


//Hamiltonoperator H = - sum(s_i J s_j) - sum(h_i s_i) + C
void H_mit_const(vector<complex<double> >& psi, std::vector<std::complex<double> >& psi_tmp, int d05, vector<double> J_ij[11], vector<double> h_i[4], double C){
	
	for (int z = 0; z < d05; z++){
		psi_tmp[z] = complex<double> (0,0);
	}
	
	for(int z = 0; z < h_i[0].size(); z++){
		
		if(h_i[1][z] != 0){
			s_x(psi, psi_tmp, d05,(int) h_i[0][z],(-1.0)*h_i[1][z],0);		//psi_tmp += s_i_x*h_i_x*psi		
		}
		if(h_i[2][z] != 0){
			s_y(psi, psi_tmp, d05,(int) h_i[0][z],(-1.0)*h_i[2][z],0);		//psi_tmp += s_i_y*h_i_y*psi		
		}
		if(h_i[3][z] != 0){
			s_z(psi, psi_tmp, d05,(int) h_i[0][z],(-1.0)*h_i[3][z],0);		//psi_tmp += s_i_z*h_i_z*psi		
		}	
	}
		
		
	for(int z = 0; z < J_ij[0].size(); z++){
	
		int i = (int) J_ij[0][z];
		int j = (int) J_ij[1][z];
		
		//Die Folgenden J_ab Konstanten beziehen sich auf s_i^a*s_j^b mit i<j !!
		double Jxx = J_ij[2][z];
		double Jyy = J_ij[6][z];
		double Jzz = J_ij[10][z];
		double Jxy = 0;
		double Jyx = 0;
		double Jxz = 0;
		double Jzx = 0;
		double Jyz = 0;
		double Jzy = 0;
		
		if (i < j){
			Jxy = J_ij[3][z];
			Jyx = J_ij[5][z];
			Jxz = J_ij[4][z];
			Jzx = J_ij[8][z];
			Jyz = J_ij[7][z];
			Jzy = J_ij[9][z];
		}
		else{
			Jxy = J_ij[5][z];
			Jyx = J_ij[3][z];
			Jxz = J_ij[8][z];
			Jzx = J_ij[4][z];
			Jyz = J_ij[9][z];
			Jzy = J_ij[7][z];
		}
		
		if(Jxx != 0){
			s_x_s_x(psi, psi_tmp, d05, i, j, (-1.0)*Jxx,0);
		}
		if(Jxy != 0){
			s_x_s_y(psi, psi_tmp, d05, i, j, (-1.0)*Jxy,0);
		}
		if(Jyx != 0){
			s_x_s_y(psi, psi_tmp, d05, j, i, (-1.0)*Jyx,0);
		}
		if(Jxz != 0){
			s_x_s_z(psi, psi_tmp, d05, i, j, (-1.0)*Jxz,0);
		}
		if(Jzx != 0){
			s_x_s_z(psi, psi_tmp, d05, j, i, (-1.0)*Jzx,0);
		}
		if(Jyy != 0){
			s_y_s_y(psi, psi_tmp, d05, i, j, (-1.0)*Jyy,0);
		}
		if(Jyz != 0){
			s_y_s_z(psi, psi_tmp, d05, i, j, (-1.0)*Jyz,0);
		}
		if(Jzy != 0){
			s_y_s_z(psi, psi_tmp, d05, j, i, (-1.0)*Jzy,0);
		}
		if(Jzz != 0){
			s_z_s_z(psi, psi_tmp, d05, i, j, (-1.0)*Jzz,0);
		}
	}
	
	
	if(C != 0){
		const_eins(psi, psi_tmp, d05,0 ,C ,0);		
	}
	
	
	
}