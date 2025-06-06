#include "initial.h"

using namespace std;


void initial(int N05, double schrittweite, vector<complex<double> > trotter_ze[5], vector<complex<double> > trotter_J_ij_1[18], vector<complex<double> > trotter_J_ij_2[18], vector<double> J_ij[11], vector<double> h_i[4]){
	
	
	//Suzuki-Trotter Faktoren der Zeeman-Terme bauen
	for(int z = 0; z < h_i[0].size(); z++){
		
		double h = sqrt(h_i[1][z]*h_i[1][z] + h_i[2][z]*h_i[2][z] + h_i[3][z]*h_i[3][z]);
		
		if (h != 0){
			complex<double> a00( cos(schrittweite*h/2.0), sin(schrittweite*h/2.0)*h_i[3][z]/h );
			complex<double> a01( sin(schrittweite*h/2.0)*h_i[2][z]/h , sin(schrittweite*h/2.0)*h_i[1][z]/h );
			complex<double> a10( (-1)*sin(schrittweite*h/2.0)*h_i[2][z]/h , sin(schrittweite*h/2.0)*h_i[1][z]/h );
			complex<double> a11( cos(schrittweite*h/2.0), (-1)*sin(schrittweite*h/2.0)*h_i[3][z]/h );
			
			trotter_ze[0].push_back(a00);
			trotter_ze[1].push_back(a01);
			trotter_ze[2].push_back(a10);
			trotter_ze[3].push_back(a11);
			trotter_ze[4].push_back(complex<double> (h_i[0][z], 0));
		}
	}
	
	
	//Suzuki-Trotter Faktoren der J_ij Terme (i != j) bauen
	vector<complex<double> > trotter_xx_yy_zz[18];			
	vector<complex<double> > trotter_xy_yx[18];				
	vector<complex<double> > trotter_xz[18];					
	vector<complex<double> > trotter_zx[18];					
	vector<complex<double> > trotter_yz[18];					
	vector<complex<double> > trotter_zy[18];					
	
	
	for(int z = 0; z < J_ij[0].size(); z++){
		
		complex<double> a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33;
		int i = (int) J_ij[0][z];
		int j = (int) J_ij[1][z];
		
		//Die Folgenden J_ab Konstanten beziehen sich auf s_i^a*s_j^b mit i>j !!
		double Jxx = J_ij[2][z];
		double Jyy = J_ij[6][z];
		double Jzz = J_ij[10][z];
		double Jxy = 0;
		double Jyx = 0;
		double Jxz = 0;
		double Jzx = 0;
		double Jyz = 0;
		double Jzy = 0;
		
		if (j < i){
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
		
		
		if (Jxx != 0 || Jyy != 0 || Jzz != 0){
					
			double a = Jzz/4.0;
			double b = (Jxx - Jyy)/4.0;
			double c = (Jxx + Jyy)/4.0;
					
			a00 = complex<double> (cos(a*schrittweite)*cos(b*schrittweite), sin(a*schrittweite)*cos(b*schrittweite));
			a01 = complex<double> (0, 0);
			a02 = complex<double> (0, 0);
			a03 = complex<double> ((-1.0)*sin(a*schrittweite)*sin(b*schrittweite), cos(a*schrittweite)*sin(b*schrittweite));
					
			a10 = complex<double> (0, 0);
			a11 = complex<double> (cos(a*schrittweite)*cos(c*schrittweite), (-1.0)*sin(a*schrittweite)*cos(c*schrittweite));
			a12 = complex<double> (sin(a*schrittweite)*sin(c*schrittweite), cos(a*schrittweite)*sin(c*schrittweite));
			a13 = complex<double> (0, 0);
					
			a20 = complex<double> (0, 0);
			a21 = a12;
			a22 = a11;
			a23 = complex<double> (0, 0);
					
			a30 = a03;
			a31 = complex<double> (0, 0);
			a32 = complex<double> (0, 0);
			a33 = a00;
					
			trotter_xx_yy_zz[0].push_back(a00);
			trotter_xx_yy_zz[1].push_back(a01);
			trotter_xx_yy_zz[2].push_back(a02);
			trotter_xx_yy_zz[3].push_back(a03);
			trotter_xx_yy_zz[4].push_back(a10);
			trotter_xx_yy_zz[5].push_back(a11);
			trotter_xx_yy_zz[6].push_back(a12);
			trotter_xx_yy_zz[7].push_back(a13);
			trotter_xx_yy_zz[8].push_back(a20);
			trotter_xx_yy_zz[9].push_back(a21);
			trotter_xx_yy_zz[10].push_back(a22);
			trotter_xx_yy_zz[11].push_back(a23);
			trotter_xx_yy_zz[12].push_back(a30);
			trotter_xx_yy_zz[13].push_back(a31);
			trotter_xx_yy_zz[14].push_back(a32);
			trotter_xx_yy_zz[15].push_back(a33);
			trotter_xx_yy_zz[16].push_back(complex<double> (i, 0));
			trotter_xx_yy_zz[17].push_back(complex<double> (j, 0));
		}
		
		
		if (Jxy != 0 || Jyx != 0){
					
			a00 = complex<double> (cos((Jxy + Jyx)*schrittweite/4.0), 0);
			a01 = complex<double> (0, 0);
			a02 = complex<double> (0, 0);
			a03 = complex<double> (sin((Jxy + Jyx)*schrittweite/4.0), 0);
					
			a10 = complex<double> (0, 0);
			a11 = complex<double> (cos((Jxy - Jyx)*schrittweite/4.0), 0);
			a12 = complex<double> (-1.0*sin((Jxy - Jyx)*schrittweite/4.0), 0);
			a13 = complex<double> (0, 0);
					
			a20 = complex<double> (0, 0);
			a21 = complex<double> (sin((Jxy - Jyx)*schrittweite/4.0), 0);
			a22 = complex<double> (cos((Jxy - Jyx)*schrittweite/4.0), 0);
			a23 = complex<double> (0, 0);
					
			a30 = complex<double> (-1.0*sin((Jxy + Jyx)*schrittweite/4.0), 0);
			a31 = complex<double> (0, 0);
			a32 = complex<double> (0, 0);
			a33 = complex<double> (cos((Jxy + Jyx)*schrittweite/4.0), 0);
					
			trotter_xy_yx[0].push_back(a00);
			trotter_xy_yx[1].push_back(a01);
			trotter_xy_yx[2].push_back(a02);
			trotter_xy_yx[3].push_back(a03);
			trotter_xy_yx[4].push_back(a10);
			trotter_xy_yx[5].push_back(a11);
			trotter_xy_yx[6].push_back(a12);
			trotter_xy_yx[7].push_back(a13);
			trotter_xy_yx[8].push_back(a20);
			trotter_xy_yx[9].push_back(a21);
			trotter_xy_yx[10].push_back(a22);
			trotter_xy_yx[11].push_back(a23);
			trotter_xy_yx[12].push_back(a30);
			trotter_xy_yx[13].push_back(a31);
			trotter_xy_yx[14].push_back(a32);
			trotter_xy_yx[15].push_back(a33);
			trotter_xy_yx[16].push_back(complex<double> (i, 0));
			trotter_xy_yx[17].push_back(complex<double> (j, 0));
		}
		
		
		if (Jxz != 0){
					
			a00 = complex<double> (cos(Jxz*schrittweite/4.0), 0);
			a01 = complex<double> (0, 0);
			a02 = complex<double> (0, sin(Jxz*schrittweite/4.0));
			a03 = complex<double> (0, 0);
					
			a10 = complex<double> (0, 0);
			a11 = complex<double> (cos(Jxz*schrittweite/4.0), 0);
			a12 = complex<double> (0, 0);
			a13 = complex<double> (0, -1.0*sin(Jxz*schrittweite/4.0));
					
			a20 = complex<double> (0, sin(Jxz*schrittweite/4.0));
			a21 = complex<double> (0, 0);
			a22 = complex<double> (cos(Jxz*schrittweite/4.0), 0);
			a23 = complex<double> (0, 0);
					
			a30 = complex<double> (0, 0);
			a31 = complex<double> (0, -1.0*sin(Jxz*schrittweite/4.0));
			a32 = complex<double> (0, 0);
			a33 = complex<double> (cos(Jxz*schrittweite/4.0), 0);
					
			trotter_xz[0].push_back(a00);
			trotter_xz[1].push_back(a01);
			trotter_xz[2].push_back(a02);
			trotter_xz[3].push_back(a03);
			trotter_xz[4].push_back(a10);
			trotter_xz[5].push_back(a11);
			trotter_xz[6].push_back(a12);
			trotter_xz[7].push_back(a13);
			trotter_xz[8].push_back(a20);
			trotter_xz[9].push_back(a21);
			trotter_xz[10].push_back(a22);
			trotter_xz[11].push_back(a23);
			trotter_xz[12].push_back(a30);
			trotter_xz[13].push_back(a31);
			trotter_xz[14].push_back(a32);
			trotter_xz[15].push_back(a33);
			trotter_xz[16].push_back(complex<double> (i, 0));
			trotter_xz[17].push_back(complex<double> (j, 0));
		}
				
				
		if (Jzx != 0){
					
			a00 = complex<double> (cos(Jzx*schrittweite/4.0), 0);
			a01 = complex<double> (0, sin(Jzx*schrittweite/4.0));
			a02 = complex<double> (0, 0);
			a03 = complex<double> (0, 0);
					
			a10 = complex<double> (0, sin(Jzx*schrittweite/4.0));
			a11 = complex<double> (cos(Jzx*schrittweite/4.0), 0);
			a12 = complex<double> (0, 0);
			a13 = complex<double> (0, 0);
					
			a20 = complex<double> (0, 0);
			a21 = complex<double> (0, 0);
			a22 = complex<double> (cos(Jzx*schrittweite/4.0), 0);
			a23 = complex<double> (0, -1.0*sin(Jzx*schrittweite/4.0));
					
			a30 = complex<double> (0, 0);
			a31 = complex<double> (0, 0);
			a32 = complex<double> (0, -1.0*sin(Jzx*schrittweite/4.0));
			a33 = complex<double> (cos(Jzx*schrittweite/4.0), 0);
					
			trotter_zx[0].push_back(a00);
			trotter_zx[1].push_back(a01);
			trotter_zx[2].push_back(a02);
			trotter_zx[3].push_back(a03);
			trotter_zx[4].push_back(a10);
			trotter_zx[5].push_back(a11);
			trotter_zx[6].push_back(a12);
			trotter_zx[7].push_back(a13);
			trotter_zx[8].push_back(a20);
			trotter_zx[9].push_back(a21);
			trotter_zx[10].push_back(a22);
			trotter_zx[11].push_back(a23);
			trotter_zx[12].push_back(a30);
			trotter_zx[13].push_back(a31);
			trotter_zx[14].push_back(a32);
			trotter_zx[15].push_back(a33);
			trotter_zx[16].push_back(complex<double> (i, 0));
			trotter_zx[17].push_back(complex<double> (j, 0));
		}
				
				
		if (Jyz != 0){
					
			a00 = complex<double> (cos(Jyz*schrittweite/4.0), 0);
			a01 = complex<double> (0, 0);
			a02 = complex<double> (sin(Jyz*schrittweite/4.0), 0);
			a03 = complex<double> (0, 0);
					
			a10 = complex<double> (0, 0);
			a11 = complex<double> (cos(Jyz*schrittweite/4.0), 0);
			a12 = complex<double> (0, 0);
			a13 = complex<double> (-1.0*sin(Jyz*schrittweite/4.0), 0);
					
			a20 = complex<double> (-1.0*sin(Jyz*schrittweite/4.0), 0);
			a21 = complex<double> (0, 0);
			a22 = complex<double> (cos(Jyz*schrittweite/4.0), 0);
			a23 = complex<double> (0, 0);
					
			a30 = complex<double> (0, 0);
			a31 = complex<double> (sin(Jyz*schrittweite/4.0), 0);
			a32 = complex<double> (0, 0);
			a33 = complex<double> (cos(Jyz*schrittweite/4.0), 0);
					
			trotter_yz[0].push_back(a00);
			trotter_yz[1].push_back(a01);
			trotter_yz[2].push_back(a02);
			trotter_yz[3].push_back(a03);
			trotter_yz[4].push_back(a10);
			trotter_yz[5].push_back(a11);
			trotter_yz[6].push_back(a12);
			trotter_yz[7].push_back(a13);
			trotter_yz[8].push_back(a20);
			trotter_yz[9].push_back(a21);
			trotter_yz[10].push_back(a22);
			trotter_yz[11].push_back(a23);
			trotter_yz[12].push_back(a30);
			trotter_yz[13].push_back(a31);
			trotter_yz[14].push_back(a32);
			trotter_yz[15].push_back(a33);
			trotter_yz[16].push_back(complex<double> (i, 0));
			trotter_yz[17].push_back(complex<double> (j, 0));
		}
				
				
		if (Jzy != 0){
				
			a00 = complex<double> (cos(Jzy*schrittweite/4.0), 0);
			a01 = complex<double> (sin(Jzy*schrittweite/4.0), 0);
			a02 = complex<double> (0, 0);
			a03 = complex<double> (0, 0);
					
			a10 = complex<double> (-1.0*sin(Jzy*schrittweite/4.0), 0);
			a11 = complex<double> (cos(Jzy*schrittweite/4.0), 0);
			a12 = complex<double> (0, 0);
			a13 = complex<double> (0, 0);
					
			a20 = complex<double> (0, 0);
			a21 = complex<double> (0, 0);
			a22 = complex<double> (cos(Jzy*schrittweite/4.0), 0);
			a23 = complex<double> (-1.0*sin(Jzy*schrittweite/4.0), 0);
					
			a30 = complex<double> (0, 0);
			a31 = complex<double> (0, 0);
			a32 = complex<double> (sin(Jzy*schrittweite/4.0), 0);
			a33 = complex<double> (cos(Jzy*schrittweite/4.0), 0);
					
			trotter_zy[0].push_back(a00);
			trotter_zy[1].push_back(a01);
			trotter_zy[2].push_back(a02);
			trotter_zy[3].push_back(a03);
			trotter_zy[4].push_back(a10);
			trotter_zy[5].push_back(a11);
			trotter_zy[6].push_back(a12);
			trotter_zy[7].push_back(a13);
			trotter_zy[8].push_back(a20);
			trotter_zy[9].push_back(a21);
			trotter_zy[10].push_back(a22);
			trotter_zy[11].push_back(a23);
			trotter_zy[12].push_back(a30);
			trotter_zy[13].push_back(a31);
			trotter_zy[14].push_back(a32);
			trotter_zy[15].push_back(a33);
			trotter_zy[16].push_back(complex<double> (i, 0));
			trotter_zy[17].push_back(complex<double> (j, 0));
		}
		
		
		//Zusammenfassung
		if (Jxx != 0 || Jyy != 0 || Jzz != 0 || Jxy != 0 || Jyx != 0 || Jxz != 0 || Jzx != 0 || Jyz != 0 || Jzy != 0){
					
			trotter_J_ij_1[0].push_back(complex<double> (1,0));
			trotter_J_ij_1[1].push_back(complex<double> (0,0));
			trotter_J_ij_1[2].push_back(complex<double> (0,0));
			trotter_J_ij_1[3].push_back(complex<double> (0,0));
			trotter_J_ij_1[4].push_back(complex<double> (0,0));
			trotter_J_ij_1[5].push_back(complex<double> (1,0));
			trotter_J_ij_1[6].push_back(complex<double> (0,0));
			trotter_J_ij_1[7].push_back(complex<double> (0,0));
			trotter_J_ij_1[8].push_back(complex<double> (0,0));
			trotter_J_ij_1[9].push_back(complex<double> (0,0));
			trotter_J_ij_1[10].push_back(complex<double> (1,0));
			trotter_J_ij_1[11].push_back(complex<double> (0,0));
			trotter_J_ij_1[12].push_back(complex<double> (0,0));
			trotter_J_ij_1[13].push_back(complex<double> (0,0));
			trotter_J_ij_1[14].push_back(complex<double> (0,0));
			trotter_J_ij_1[15].push_back(complex<double> (1,0));
			trotter_J_ij_1[16].push_back(complex<double> (i,0));
			trotter_J_ij_1[17].push_back(complex<double> (j,0));
					
			trotter_J_ij_2[0].push_back(complex<double> (1,0));
			trotter_J_ij_2[1].push_back(complex<double> (0,0));
			trotter_J_ij_2[2].push_back(complex<double> (0,0));
			trotter_J_ij_2[3].push_back(complex<double> (0,0));
			trotter_J_ij_2[4].push_back(complex<double> (0,0));
			trotter_J_ij_2[5].push_back(complex<double> (1,0));
			trotter_J_ij_2[6].push_back(complex<double> (0,0));
			trotter_J_ij_2[7].push_back(complex<double> (0,0));
			trotter_J_ij_2[8].push_back(complex<double> (0,0));
			trotter_J_ij_2[9].push_back(complex<double> (0,0));
			trotter_J_ij_2[10].push_back(complex<double> (1,0));
			trotter_J_ij_2[11].push_back(complex<double> (0,0));
			trotter_J_ij_2[12].push_back(complex<double> (0,0));
			trotter_J_ij_2[13].push_back(complex<double> (0,0));
			trotter_J_ij_2[14].push_back(complex<double> (0,0));
			trotter_J_ij_2[15].push_back(complex<double> (1,0));
			trotter_J_ij_2[16].push_back(complex<double> (i,0));
			trotter_J_ij_2[17].push_back(complex<double> (j,0));
					
			int index = trotter_J_ij_1[0].size() - 1;
					
			if (Jxx != 0 || Jyy != 0 || Jzz != 0){
						
				vector<complex<double> > tmp;
				for(int z3 = 0; z3 < 16; z3++){
					tmp.push_back(complex<double> (0,0));
				}
				
				int index2 = trotter_xx_yy_zz[0].size() - 1;
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						for(int z3 = 0; z3 < 4; z3++){
							tmp[spalte + zeile*4] += trotter_xx_yy_zz[z3 + zeile*4][index2]*trotter_J_ij_1[spalte + z3*4][index];
						}	
					}
				}
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						trotter_J_ij_1[spalte + zeile*4][index] = tmp[spalte + zeile*4];
					}
				}
			}
					
					
			if (Jxy != 0 || Jyx != 0){
						
				vector<complex<double> > tmp;
				for(int z3 = 0; z3 < 16; z3++){
					tmp.push_back(complex<double> (0,0));
				}
				
				int index2 = trotter_xy_yx[0].size() - 1;
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						for(int z3 = 0; z3 < 4; z3++){
							tmp[spalte + zeile*4] += trotter_xy_yx[z3 + zeile*4][index2]*trotter_J_ij_1[spalte + z3*4][index];
						}	
					}
				}
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						trotter_J_ij_1[spalte + zeile*4][index] = tmp[spalte + zeile*4];
					}
				}
			}
					
					
			if (Jxz != 0){
						
				vector<complex<double> > tmp;
				for(int z3 = 0; z3 < 16; z3++){
					tmp.push_back(complex<double> (0,0));
				}
				
				int index2 = trotter_xz[0].size() - 1;
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						for(int z3 = 0; z3 < 4; z3++){
							tmp[spalte + zeile*4] += trotter_xz[z3 + zeile*4][index2]*trotter_J_ij_1[spalte + z3*4][index];
						}	
					}
				}
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						trotter_J_ij_1[spalte + zeile*4][index] = tmp[spalte + zeile*4];
					}
				}
			}
					
					
			if (Jzx != 0){
						
				vector<complex<double> > tmp;
				for(int z3 = 0; z3 < 16; z3++){
					tmp.push_back(complex<double> (0,0));
				}
				
				int index2 = trotter_zx[0].size() - 1;
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						for(int z3 = 0; z3 < 4; z3++){
							tmp[spalte + zeile*4] += trotter_zx[z3 + zeile*4][index2]*trotter_J_ij_1[spalte + z3*4][index];
						}	
					}
				}
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						trotter_J_ij_1[spalte + zeile*4][index] = tmp[spalte + zeile*4];
					}
				}
			}
					
					
			if (Jyz != 0){
						
				vector<complex<double> > tmp;
				for(int z3 = 0; z3 < 16; z3++){
					tmp.push_back(complex<double> (0,0));
				}
				
				int index2 = trotter_yz[0].size() - 1;
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						for(int z3 = 0; z3 < 4; z3++){
							tmp[spalte + zeile*4] += trotter_yz[z3 + zeile*4][index2]*trotter_J_ij_1[spalte + z3*4][index];
						}	
					}
				}
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						trotter_J_ij_1[spalte + zeile*4][index] = tmp[spalte + zeile*4];
					}
				}
			}
					
					
			if (Jzy != 0){
						
				vector<complex<double> > tmp;
				for(int z3 = 0; z3 < 16; z3++){
					tmp.push_back(complex<double> (0,0));
				}
				
				int index2 = trotter_zy[0].size() - 1;
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						for(int z3 = 0; z3 < 4; z3++){
							tmp[spalte + zeile*4] += trotter_zy[z3 + zeile*4][index2]*trotter_J_ij_1[spalte + z3*4][index];
						}	
					}
				}
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						trotter_J_ij_1[spalte + zeile*4][index] = tmp[spalte + zeile*4];
					}
				}
			}
					
					
					
			//Andersherum (für höhere Ordnung)	
			if (Jzy != 0){
						
				vector<complex<double> > tmp;
				for(int z3 = 0; z3 < 16; z3++){
					tmp.push_back(complex<double> (0,0));
				}
				
				int index2 = trotter_zy[0].size() - 1;
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						for(int z3 = 0; z3 < 4; z3++){
							tmp[spalte + zeile*4] += trotter_zy[z3 + zeile*4][index2]*trotter_J_ij_2[spalte + z3*4][index];
						}	
					}
				}
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						trotter_J_ij_2[spalte + zeile*4][index] = tmp[spalte + zeile*4];
					}
				}
			}
					
					
			if (Jyz != 0){
						
				vector<complex<double> > tmp;
				for(int z3 = 0; z3 < 16; z3++){
					tmp.push_back(complex<double> (0,0));
				}
				
				int index2 = trotter_yz[0].size() - 1;
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						for(int z3 = 0; z3 < 4; z3++){
							tmp[spalte + zeile*4] += trotter_yz[z3 + zeile*4][index2]*trotter_J_ij_2[spalte + z3*4][index];
						}	
					}
				}
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						trotter_J_ij_2[spalte + zeile*4][index] = tmp[spalte + zeile*4];
					}
				}
			}
					
					
			if (Jzx != 0){
						
				vector<complex<double> > tmp;
				for(int z3 = 0; z3 < 16; z3++){
					tmp.push_back(complex<double> (0,0));
				}
				
				int index2 = trotter_zx[0].size() - 1;
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						for(int z3 = 0; z3 < 4; z3++){
							tmp[spalte + zeile*4] += trotter_zx[z3 + zeile*4][index2]*trotter_J_ij_2[spalte + z3*4][index];
						}	
					}
				}
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						trotter_J_ij_2[spalte + zeile*4][index] = tmp[spalte + zeile*4];
					}
				}
			}
					
					
			if (Jxz != 0){
						
				vector<complex<double> > tmp;
				for(int z3 = 0; z3 < 16; z3++){
					tmp.push_back(complex<double> (0,0));
				}
				
				int index2 = trotter_xz[0].size() - 1;
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						for(int z3 = 0; z3 < 4; z3++){
							tmp[spalte + zeile*4] += trotter_xz[z3 + zeile*4][index2]*trotter_J_ij_2[spalte + z3*4][index];
						}	
					}
				}
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						trotter_J_ij_2[spalte + zeile*4][index] = tmp[spalte + zeile*4];
					}
				}
			}
					
					
			if (Jxy != 0 || Jyx != 0){
						
				vector<complex<double> > tmp;
				for(int z3 = 0; z3 < 16; z3++){
					tmp.push_back(complex<double> (0,0));
				}
				
				int index2 = trotter_xy_yx[0].size() - 1;
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						for(int z3 = 0; z3 < 4; z3++){
							tmp[spalte + zeile*4] += trotter_xy_yx[z3 + zeile*4][index2]*trotter_J_ij_2[spalte + z3*4][index];
						}	
					}
				}
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						trotter_J_ij_2[spalte + zeile*4][index] = tmp[spalte + zeile*4];
					}
				}
			}
					
					
			if (Jxx != 0 || Jyy != 0 || Jzz != 0){
						
				vector<complex<double> > tmp;
				for(int z3 = 0; z3 < 16; z3++){
					tmp.push_back(complex<double> (0,0));
				}
				
				int index2 = trotter_xx_yy_zz[0].size() - 1;
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						for(int z3 = 0; z3 < 4; z3++){
							tmp[spalte + zeile*4] += trotter_xx_yy_zz[z3 + zeile*4][index2]*trotter_J_ij_2[spalte + z3*4][index];
						}	
					}
				}
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						trotter_J_ij_2[spalte + zeile*4][index] = tmp[spalte + zeile*4];
					}
				}
			}
		}
	}
}
















void initial_imag(int N05, double schrittweite, vector<complex<double> > trotter_ze[5], vector<complex<double> > trotter_J_ij_1[18], vector<complex<double> > trotter_J_ij_2[18], vector<double> J_ij[11], vector<double> h_i[4]){
	
	
	//Suzuki-Trotter Faktoren der Zeeman-Terme bauen
	for(int z = 0; z < h_i[0].size(); z++){
		
		double h = sqrt(h_i[1][z]*h_i[1][z] + h_i[2][z]*h_i[2][z] + h_i[3][z]*h_i[3][z]);
		
		if (h != 0){
			complex<double> a00(exp(schrittweite*(-1.0)*h/2.0)*((-1.0)*h_i[3][z] + h)/(2.0*h) +  exp(schrittweite*h/2.0)*(h_i[3][z] + h)/(2.0*h), 0);
			complex<double> a01(exp(schrittweite*(-1.0)*h/2.0)*((-1.0)*h_i[1][z])/(2.0*h) + exp(schrittweite*h/2.0)*(h_i[1][z])/(2.0*h), exp(schrittweite*(-1.0)*h/2.0)*(h_i[2][z])/(2.0*h) + exp(schrittweite*h/2.0)*((-1.0)*h_i[2][z])/(2.0*h));
			complex<double> a10(exp(schrittweite*(-1.0)*h/2.0)*((-1.0)*h_i[1][z])/(2.0*h) + exp(schrittweite*h/2.0)*(h_i[1][z])/(2.0*h), exp(schrittweite*(-1.0)*h/2.0)*((-1.0)*h_i[2][z])/(2.0*h) + exp(schrittweite*h/2.0)*(h_i[2][z])/(2.0*h));
			complex<double> a11(exp(schrittweite*(-1.0)*h/2.0)*(h_i[3][z] + h)/(2.0*h) +  exp(schrittweite*h/2.0)*((-1.0)*h_i[3][z] + h)/(2.0*h), 0);
			
			trotter_ze[0].push_back(a00);
			trotter_ze[1].push_back(a01);
			trotter_ze[2].push_back(a10);
			trotter_ze[3].push_back(a11);
			trotter_ze[4].push_back(complex<double> (h_i[0][z], 0));
		}
	}
	
	
	//Suzuki-Trotter Faktoren der J_ij Terme (i != j) bauen
	vector<complex<double> > trotter_xx_yy_zz[18];			
	vector<complex<double> > trotter_xy_yx[18];				
	vector<complex<double> > trotter_xz[18];					
	vector<complex<double> > trotter_zx[18];					
	vector<complex<double> > trotter_yz[18];					
	vector<complex<double> > trotter_zy[18];					
	
	
	for(int z = 0; z < J_ij[0].size(); z++){
		
		complex<double> a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33;
		int i = (int) J_ij[0][z];
		int j = (int) J_ij[1][z];
		
		//Die Folgenden J_ab Konstanten beziehen sich auf s_i^a*s_j^b mit i>j !!
		double Jxx = J_ij[2][z];
		double Jyy = J_ij[6][z];
		double Jzz = J_ij[10][z];
		double Jxy = 0;
		double Jyx = 0;
		double Jxz = 0;
		double Jzx = 0;
		double Jyz = 0;
		double Jzy = 0;
		
		if (j < i){
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
		
		
		if (Jxx != 0 || Jyy != 0 || Jzz != 0){
							
			a00 = complex<double> (0.5*exp(0.25*schrittweite*(Jxx - Jyy + Jzz)) + 0.5*exp(0.25*schrittweite*(-Jxx + Jyy + Jzz)), 0);
			a01 = complex<double> (0, 0);
			a02 = complex<double> (0, 0);
			a03 = complex<double> (0.5*exp(0.25*schrittweite*(Jxx - Jyy + Jzz)) - 0.5*exp(0.25*schrittweite*(-Jxx + Jyy + Jzz)), 0);
					
			a10 = complex<double> (0, 0);
			a11 = complex<double> (0.5*exp(0.25*schrittweite*(Jxx + Jyy - Jzz)) + 0.5*exp(0.25*schrittweite*(-Jxx - Jyy - Jzz)), 0);
			a12 = complex<double> (0.5*exp(0.25*schrittweite*(Jxx + Jyy - Jzz)) - 0.5*exp(0.25*schrittweite*(-Jxx - Jyy - Jzz)), 0);
			a13 = complex<double> (0, 0);
					
			a20 = complex<double> (0, 0);
			a21 = complex<double> (0.5*exp(0.25*schrittweite*(Jxx + Jyy - Jzz)) - 0.5*exp(0.25*schrittweite*(-Jxx - Jyy - Jzz)), 0);
			a22 = complex<double> (0.5*exp(0.25*schrittweite*(Jxx + Jyy - Jzz)) + 0.5*exp(0.25*schrittweite*(-Jxx - Jyy - Jzz)), 0);
			a23 = complex<double> (0, 0);
					
			a30 = complex<double> (0.5*exp(0.25*schrittweite*(Jxx - Jyy + Jzz)) - 0.5*exp(0.25*schrittweite*(-Jxx + Jyy + Jzz)), 0);
			a31 = complex<double> (0, 0);
			a32 = complex<double> (0, 0);
			a33 = complex<double> (0.5*exp(0.25*schrittweite*(Jxx - Jyy + Jzz)) + 0.5*exp(0.25*schrittweite*(-Jxx + Jyy + Jzz)), 0);
					
			trotter_xx_yy_zz[0].push_back(a00);
			trotter_xx_yy_zz[1].push_back(a01);
			trotter_xx_yy_zz[2].push_back(a02);
			trotter_xx_yy_zz[3].push_back(a03);
			trotter_xx_yy_zz[4].push_back(a10);
			trotter_xx_yy_zz[5].push_back(a11);
			trotter_xx_yy_zz[6].push_back(a12);
			trotter_xx_yy_zz[7].push_back(a13);
			trotter_xx_yy_zz[8].push_back(a20);
			trotter_xx_yy_zz[9].push_back(a21);
			trotter_xx_yy_zz[10].push_back(a22);
			trotter_xx_yy_zz[11].push_back(a23);
			trotter_xx_yy_zz[12].push_back(a30);
			trotter_xx_yy_zz[13].push_back(a31);
			trotter_xx_yy_zz[14].push_back(a32);
			trotter_xx_yy_zz[15].push_back(a33);
			trotter_xx_yy_zz[16].push_back(complex<double> (i, 0));
			trotter_xx_yy_zz[17].push_back(complex<double> (j, 0));
		}
		
		
		if (Jxy != 0 || Jyx != 0){
					
			a00 = complex<double> (cosh((Jxy + Jyx)*schrittweite/4.0), 0);
			a01 = complex<double> (0, 0);
			a02 = complex<double> (0, 0);
			a03 = complex<double> (0, -1.0*sinh((Jxy + Jyx)*schrittweite/4.0));
					
			a10 = complex<double> (0, 0);
			a11 = complex<double> (cosh((Jxy - Jyx)*schrittweite/4.0), 0);
			a12 = complex<double> (0, sinh((Jxy - Jyx)*schrittweite/4.0));
			a13 = complex<double> (0, 0);
					
			a20 = complex<double> (0, 0);
			a21 = complex<double> (0, -1.0*sinh((Jxy - Jyx)*schrittweite/4.0));
			a22 = complex<double> (cosh((Jxy - Jyx)*schrittweite/4.0), 0);
			a23 = complex<double> (0, 0);
					
			a30 = complex<double> (0, sinh((Jxy + Jyx)*schrittweite/4.0));
			a31 = complex<double> (0, 0);
			a32 = complex<double> (0, 0);
			a33 = complex<double> (cosh((Jxy + Jyx)*schrittweite/4.0), 0);
					
			trotter_xy_yx[0].push_back(a00);
			trotter_xy_yx[1].push_back(a01);
			trotter_xy_yx[2].push_back(a02);
			trotter_xy_yx[3].push_back(a03);
			trotter_xy_yx[4].push_back(a10);
			trotter_xy_yx[5].push_back(a11);
			trotter_xy_yx[6].push_back(a12);
			trotter_xy_yx[7].push_back(a13);
			trotter_xy_yx[8].push_back(a20);
			trotter_xy_yx[9].push_back(a21);
			trotter_xy_yx[10].push_back(a22);
			trotter_xy_yx[11].push_back(a23);
			trotter_xy_yx[12].push_back(a30);
			trotter_xy_yx[13].push_back(a31);
			trotter_xy_yx[14].push_back(a32);
			trotter_xy_yx[15].push_back(a33);
			trotter_xy_yx[16].push_back(complex<double> (i, 0));
			trotter_xy_yx[17].push_back(complex<double> (j, 0));
		}
		
		
		if (Jxz != 0){
					
			a00 = complex<double> (0.5*exp((-1.0)*Jxz*schrittweite/4.0) + 0.5*exp(Jxz*schrittweite/4.0), 0);
			a01 = complex<double> (0, 0);
			a02 = complex<double> (-0.5*exp((-1.0)*Jxz*schrittweite/4.0) + 0.5*exp(Jxz*schrittweite/4.0), 0);
			a03 = complex<double> (0, 0);
					
			a10 = complex<double> (0, 0);
			a11 = complex<double> (0.5*exp((-1.0)*Jxz*schrittweite/4.0) + 0.5*exp(Jxz*schrittweite/4.0), 0);
			a12 = complex<double> (0, 0);
			a13 = complex<double> (0.5*exp((-1.0)*Jxz*schrittweite/4.0) - 0.5*exp(Jxz*schrittweite/4.0), 0);
					
			a20 = complex<double> (-0.5*exp((-1.0)*Jxz*schrittweite/4.0) + 0.5*exp(Jxz*schrittweite/4.0), 0);
			a21 = complex<double> (0, 0);
			a22 = complex<double> (0.5*exp((-1.0)*Jxz*schrittweite/4.0) + 0.5*exp(Jxz*schrittweite/4.0), 0);
			a23 = complex<double> (0, 0);
					
			a30 = complex<double> (0, 0);
			a31 = complex<double> (0.5*exp((-1.0)*Jxz*schrittweite/4.0) - 0.5*exp(Jxz*schrittweite/4.0), 0);
			a32 = complex<double> (0, 0);
			a33 = complex<double> (0.5*exp((-1.0)*Jxz*schrittweite/4.0) + 0.5*exp(Jxz*schrittweite/4.0), 0);
					
			trotter_xz[0].push_back(a00);
			trotter_xz[1].push_back(a01);
			trotter_xz[2].push_back(a02);
			trotter_xz[3].push_back(a03);
			trotter_xz[4].push_back(a10);
			trotter_xz[5].push_back(a11);
			trotter_xz[6].push_back(a12);
			trotter_xz[7].push_back(a13);
			trotter_xz[8].push_back(a20);
			trotter_xz[9].push_back(a21);
			trotter_xz[10].push_back(a22);
			trotter_xz[11].push_back(a23);
			trotter_xz[12].push_back(a30);
			trotter_xz[13].push_back(a31);
			trotter_xz[14].push_back(a32);
			trotter_xz[15].push_back(a33);
			trotter_xz[16].push_back(complex<double> (i, 0));
			trotter_xz[17].push_back(complex<double> (j, 0));
		}
				
				
		if (Jzx != 0){
					
			a00 = complex<double> (0.5*exp((-1.0)*Jzx*schrittweite/4.0) + 0.5*exp(Jzx*schrittweite/4.0), 0);
			a01 = complex<double> (-0.5*exp((-1.0)*Jzx*schrittweite/4.0) + 0.5*exp(Jzx*schrittweite/4.0), 0);
			a02 = complex<double> (0, 0);
			a03 = complex<double> (0, 0);
					
			a10 = complex<double> (-0.5*exp((-1.0)*Jzx*schrittweite/4.0) + 0.5*exp(Jzx*schrittweite/4.0), 0);
			a11 = complex<double> (0.5*exp((-1.0)*Jzx*schrittweite/4.0) + 0.5*exp(Jzx*schrittweite/4.0), 0);
			a12 = complex<double> (0, 0);
			a13 = complex<double> (0, 0);
					
			a20 = complex<double> (0, 0);
			a21 = complex<double> (0, 0);
			a22 = complex<double> (0.5*exp((-1.0)*Jzx*schrittweite/4.0) + 0.5*exp(Jzx*schrittweite/4.0), 0);
			a23 = complex<double> (0.5*exp((-1.0)*Jzx*schrittweite/4.0) - 0.5*exp(Jzx*schrittweite/4.0), 0);
					
			a30 = complex<double> (0, 0);
			a31 = complex<double> (0, 0);
			a32 = complex<double> (0.5*exp((-1.0)*Jzx*schrittweite/4.0) - 0.5*exp(Jzx*schrittweite/4.0), 0);
			a33 = complex<double> (0.5*exp((-1.0)*Jzx*schrittweite/4.0) + 0.5*exp(Jzx*schrittweite/4.0), 0);
					
			trotter_zx[0].push_back(a00);
			trotter_zx[1].push_back(a01);
			trotter_zx[2].push_back(a02);
			trotter_zx[3].push_back(a03);
			trotter_zx[4].push_back(a10);
			trotter_zx[5].push_back(a11);
			trotter_zx[6].push_back(a12);
			trotter_zx[7].push_back(a13);
			trotter_zx[8].push_back(a20);
			trotter_zx[9].push_back(a21);
			trotter_zx[10].push_back(a22);
			trotter_zx[11].push_back(a23);
			trotter_zx[12].push_back(a30);
			trotter_zx[13].push_back(a31);
			trotter_zx[14].push_back(a32);
			trotter_zx[15].push_back(a33);
			trotter_zx[16].push_back(complex<double> (i, 0));
			trotter_zx[17].push_back(complex<double> (j, 0));
		}
				
				
		if (Jyz != 0){
					
			a00 = complex<double> (cosh(Jyz*schrittweite/4.0), 0);
			a01 = complex<double> (0, 0);
			a02 = complex<double> (0, -sinh(Jyz*schrittweite/4.0));
			a03 = complex<double> (0, 0);
					
			a10 = complex<double> (0, 0);
			a11 = complex<double> (cosh(Jyz*schrittweite/4.0), 0);
			a12 = complex<double> (0, 0);
			a13 = complex<double> (0, sinh(Jyz*schrittweite/4.0));
					
			a20 = complex<double> (0, sinh(Jyz*schrittweite/4.0));
			a21 = complex<double> (0, 0);
			a22 = complex<double> (cosh(Jyz*schrittweite/4.0), 0);
			a23 = complex<double> (0, 0);
					
			a30 = complex<double> (0, 0);
			a31 = complex<double> (0, -1.0*sinh(Jyz*schrittweite/4.0));
			a32 = complex<double> (0, 0);
			a33 = complex<double> (cosh(Jyz*schrittweite/4.0), 0);
					
			trotter_yz[0].push_back(a00);
			trotter_yz[1].push_back(a01);
			trotter_yz[2].push_back(a02);
			trotter_yz[3].push_back(a03);
			trotter_yz[4].push_back(a10);
			trotter_yz[5].push_back(a11);
			trotter_yz[6].push_back(a12);
			trotter_yz[7].push_back(a13);
			trotter_yz[8].push_back(a20);
			trotter_yz[9].push_back(a21);
			trotter_yz[10].push_back(a22);
			trotter_yz[11].push_back(a23);
			trotter_yz[12].push_back(a30);
			trotter_yz[13].push_back(a31);
			trotter_yz[14].push_back(a32);
			trotter_yz[15].push_back(a33);
			trotter_yz[16].push_back(complex<double> (i, 0));
			trotter_yz[17].push_back(complex<double> (j, 0));
		}
				
				
		if (Jzy != 0){
				
			a00 = complex<double> (cosh(Jzy*schrittweite/4.0), 0);
			a01 = complex<double> (0, -1.0*sinh(Jzy*schrittweite/4.0));
			a02 = complex<double> (0, 0);
			a03 = complex<double> (0, 0);
					
			a10 = complex<double> (0, sinh(Jzy*schrittweite/4.0));
			a11 = complex<double> (cosh(Jzy*schrittweite/4.0), 0);
			a12 = complex<double> (0, 0);
			a13 = complex<double> (0, 0);
					
			a20 = complex<double> (0, 0);
			a21 = complex<double> (0, 0);
			a22 = complex<double> (cosh(Jzy*schrittweite/4.0), 0);
			a23 = complex<double> (0, sinh(Jzy*schrittweite/4.0));
					
			a30 = complex<double> (0, 0);
			a31 = complex<double> (0, 0);
			a32 = complex<double> (0, -1.0*sinh(Jzy*schrittweite/4.0));
			a33 = complex<double> (cosh(Jzy*schrittweite/4.0), 0);
					
			trotter_zy[0].push_back(a00);
			trotter_zy[1].push_back(a01);
			trotter_zy[2].push_back(a02);
			trotter_zy[3].push_back(a03);
			trotter_zy[4].push_back(a10);
			trotter_zy[5].push_back(a11);
			trotter_zy[6].push_back(a12);
			trotter_zy[7].push_back(a13);
			trotter_zy[8].push_back(a20);
			trotter_zy[9].push_back(a21);
			trotter_zy[10].push_back(a22);
			trotter_zy[11].push_back(a23);
			trotter_zy[12].push_back(a30);
			trotter_zy[13].push_back(a31);
			trotter_zy[14].push_back(a32);
			trotter_zy[15].push_back(a33);
			trotter_zy[16].push_back(complex<double> (i, 0));
			trotter_zy[17].push_back(complex<double> (j, 0));
		}
		
		
		//Zusammenfassung
		if (Jxx != 0 || Jyy != 0 || Jzz != 0 || Jxy != 0 || Jyx != 0 || Jxz != 0 || Jzx != 0 || Jyz != 0 || Jzy != 0){
					
			trotter_J_ij_1[0].push_back(complex<double> (1,0));
			trotter_J_ij_1[1].push_back(complex<double> (0,0));
			trotter_J_ij_1[2].push_back(complex<double> (0,0));
			trotter_J_ij_1[3].push_back(complex<double> (0,0));
			trotter_J_ij_1[4].push_back(complex<double> (0,0));
			trotter_J_ij_1[5].push_back(complex<double> (1,0));
			trotter_J_ij_1[6].push_back(complex<double> (0,0));
			trotter_J_ij_1[7].push_back(complex<double> (0,0));
			trotter_J_ij_1[8].push_back(complex<double> (0,0));
			trotter_J_ij_1[9].push_back(complex<double> (0,0));
			trotter_J_ij_1[10].push_back(complex<double> (1,0));
			trotter_J_ij_1[11].push_back(complex<double> (0,0));
			trotter_J_ij_1[12].push_back(complex<double> (0,0));
			trotter_J_ij_1[13].push_back(complex<double> (0,0));
			trotter_J_ij_1[14].push_back(complex<double> (0,0));
			trotter_J_ij_1[15].push_back(complex<double> (1,0));
			trotter_J_ij_1[16].push_back(complex<double> (i,0));
			trotter_J_ij_1[17].push_back(complex<double> (j,0));
					
			trotter_J_ij_2[0].push_back(complex<double> (1,0));
			trotter_J_ij_2[1].push_back(complex<double> (0,0));
			trotter_J_ij_2[2].push_back(complex<double> (0,0));
			trotter_J_ij_2[3].push_back(complex<double> (0,0));
			trotter_J_ij_2[4].push_back(complex<double> (0,0));
			trotter_J_ij_2[5].push_back(complex<double> (1,0));
			trotter_J_ij_2[6].push_back(complex<double> (0,0));
			trotter_J_ij_2[7].push_back(complex<double> (0,0));
			trotter_J_ij_2[8].push_back(complex<double> (0,0));
			trotter_J_ij_2[9].push_back(complex<double> (0,0));
			trotter_J_ij_2[10].push_back(complex<double> (1,0));
			trotter_J_ij_2[11].push_back(complex<double> (0,0));
			trotter_J_ij_2[12].push_back(complex<double> (0,0));
			trotter_J_ij_2[13].push_back(complex<double> (0,0));
			trotter_J_ij_2[14].push_back(complex<double> (0,0));
			trotter_J_ij_2[15].push_back(complex<double> (1,0));
			trotter_J_ij_2[16].push_back(complex<double> (i,0));
			trotter_J_ij_2[17].push_back(complex<double> (j,0));
					
			int index = trotter_J_ij_1[0].size() - 1;
					
			if (Jxx != 0 || Jyy != 0 || Jzz != 0){
						
				vector<complex<double> > tmp;
				for(int z3 = 0; z3 < 16; z3++){
					tmp.push_back(complex<double> (0,0));
				}
				
				int index2 = trotter_xx_yy_zz[0].size() - 1;
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						for(int z3 = 0; z3 < 4; z3++){
							tmp[spalte + zeile*4] += trotter_xx_yy_zz[z3 + zeile*4][index2]*trotter_J_ij_1[spalte + z3*4][index];
						}	
					}
				}
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						trotter_J_ij_1[spalte + zeile*4][index] = tmp[spalte + zeile*4];
					}
				}
			}
					
					
			if (Jxy != 0 || Jyx != 0){
						
				vector<complex<double> > tmp;
				for(int z3 = 0; z3 < 16; z3++){
					tmp.push_back(complex<double> (0,0));
				}
				
				int index2 = trotter_xy_yx[0].size() - 1;
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						for(int z3 = 0; z3 < 4; z3++){
							tmp[spalte + zeile*4] += trotter_xy_yx[z3 + zeile*4][index2]*trotter_J_ij_1[spalte + z3*4][index];
						}	
					}
				}
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						trotter_J_ij_1[spalte + zeile*4][index] = tmp[spalte + zeile*4];
					}
				}
			}
					
					
			if (Jxz != 0){
						
				vector<complex<double> > tmp;
				for(int z3 = 0; z3 < 16; z3++){
					tmp.push_back(complex<double> (0,0));
				}
				
				int index2 = trotter_xz[0].size() - 1;
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						for(int z3 = 0; z3 < 4; z3++){
							tmp[spalte + zeile*4] += trotter_xz[z3 + zeile*4][index2]*trotter_J_ij_1[spalte + z3*4][index];
						}	
					}
				}
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						trotter_J_ij_1[spalte + zeile*4][index] = tmp[spalte + zeile*4];
					}
				}
			}
					
					
			if (Jzx != 0){
						
				vector<complex<double> > tmp;
				for(int z3 = 0; z3 < 16; z3++){
					tmp.push_back(complex<double> (0,0));
				}
				
				int index2 = trotter_zx[0].size() - 1;
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						for(int z3 = 0; z3 < 4; z3++){
							tmp[spalte + zeile*4] += trotter_zx[z3 + zeile*4][index2]*trotter_J_ij_1[spalte + z3*4][index];
						}	
					}
				}
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						trotter_J_ij_1[spalte + zeile*4][index] = tmp[spalte + zeile*4];
					}
				}
			}
					
					
			if (Jyz != 0){
						
				vector<complex<double> > tmp;
				for(int z3 = 0; z3 < 16; z3++){
					tmp.push_back(complex<double> (0,0));
				}
				
				int index2 = trotter_yz[0].size() - 1;
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						for(int z3 = 0; z3 < 4; z3++){
							tmp[spalte + zeile*4] += trotter_yz[z3 + zeile*4][index2]*trotter_J_ij_1[spalte + z3*4][index];
						}	
					}
				}
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						trotter_J_ij_1[spalte + zeile*4][index] = tmp[spalte + zeile*4];
					}
				}
			}
					
					
			if (Jzy != 0){
						
				vector<complex<double> > tmp;
				for(int z3 = 0; z3 < 16; z3++){
					tmp.push_back(complex<double> (0,0));
				}
				
				int index2 = trotter_zy[0].size() - 1;
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						for(int z3 = 0; z3 < 4; z3++){
							tmp[spalte + zeile*4] += trotter_zy[z3 + zeile*4][index2]*trotter_J_ij_1[spalte + z3*4][index];
						}	
					}
				}
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						trotter_J_ij_1[spalte + zeile*4][index] = tmp[spalte + zeile*4];
					}
				}
			}
					
					
					
			//Andersherum (für höhere Ordnung)	
			if (Jzy != 0){
						
				vector<complex<double> > tmp;
				for(int z3 = 0; z3 < 16; z3++){
					tmp.push_back(complex<double> (0,0));
				}
				
				int index2 = trotter_zy[0].size() - 1;
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						for(int z3 = 0; z3 < 4; z3++){
							tmp[spalte + zeile*4] += trotter_zy[z3 + zeile*4][index2]*trotter_J_ij_2[spalte + z3*4][index];
						}	
					}
				}
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						trotter_J_ij_2[spalte + zeile*4][index] = tmp[spalte + zeile*4];
					}
				}
			}
					
					
			if (Jyz != 0){
						
				vector<complex<double> > tmp;
				for(int z3 = 0; z3 < 16; z3++){
					tmp.push_back(complex<double> (0,0));
				}
				
				int index2 = trotter_yz[0].size() - 1;
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						for(int z3 = 0; z3 < 4; z3++){
							tmp[spalte + zeile*4] += trotter_yz[z3 + zeile*4][index2]*trotter_J_ij_2[spalte + z3*4][index];
						}	
					}
				}
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						trotter_J_ij_2[spalte + zeile*4][index] = tmp[spalte + zeile*4];
					}
				}
			}
					
					
			if (Jzx != 0){
						
				vector<complex<double> > tmp;
				for(int z3 = 0; z3 < 16; z3++){
					tmp.push_back(complex<double> (0,0));
				}
				
				int index2 = trotter_zx[0].size() - 1;
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						for(int z3 = 0; z3 < 4; z3++){
							tmp[spalte + zeile*4] += trotter_zx[z3 + zeile*4][index2]*trotter_J_ij_2[spalte + z3*4][index];
						}	
					}
				}
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						trotter_J_ij_2[spalte + zeile*4][index] = tmp[spalte + zeile*4];
					}
				}
			}
					
					
			if (Jxz != 0){
						
				vector<complex<double> > tmp;
				for(int z3 = 0; z3 < 16; z3++){
					tmp.push_back(complex<double> (0,0));
				}
				
				int index2 = trotter_xz[0].size() - 1;
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						for(int z3 = 0; z3 < 4; z3++){
							tmp[spalte + zeile*4] += trotter_xz[z3 + zeile*4][index2]*trotter_J_ij_2[spalte + z3*4][index];
						}	
					}
				}
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						trotter_J_ij_2[spalte + zeile*4][index] = tmp[spalte + zeile*4];
					}
				}
			}
					
					
			if (Jxy != 0 || Jyx != 0){
						
				vector<complex<double> > tmp;
				for(int z3 = 0; z3 < 16; z3++){
					tmp.push_back(complex<double> (0,0));
				}
				
				int index2 = trotter_xy_yx[0].size() - 1;
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						for(int z3 = 0; z3 < 4; z3++){
							tmp[spalte + zeile*4] += trotter_xy_yx[z3 + zeile*4][index2]*trotter_J_ij_2[spalte + z3*4][index];
						}	
					}
				}
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						trotter_J_ij_2[spalte + zeile*4][index] = tmp[spalte + zeile*4];
					}
				}
			}
					
					
			if (Jxx != 0 || Jyy != 0 || Jzz != 0){
						
				vector<complex<double> > tmp;
				for(int z3 = 0; z3 < 16; z3++){
					tmp.push_back(complex<double> (0,0));
				}
				
				int index2 = trotter_xx_yy_zz[0].size() - 1;
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						for(int z3 = 0; z3 < 4; z3++){
							tmp[spalte + zeile*4] += trotter_xx_yy_zz[z3 + zeile*4][index2]*trotter_J_ij_2[spalte + z3*4][index];
						}	
					}
				}
						
				for(int zeile = 0; zeile < 4; zeile++){
					for(int spalte = 0; spalte < 4; spalte++){
						trotter_J_ij_2[spalte + zeile*4][index] = tmp[spalte + zeile*4];
					}
				}
			}
		}
	}
}




