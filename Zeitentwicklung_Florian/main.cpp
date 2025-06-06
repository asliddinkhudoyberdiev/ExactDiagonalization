#include <iostream>
#include <ctime>	
#include "operator.h"
#include "initial.h"
#include "einlesen.h"
#include <fstream>
#include <sstream>

using namespace std;
vector<complex<double> > psi;				//Zustandsvektor Psi
vector<complex<double> > psi_tmp;			//temporäre Ablage
vector<complex<double> > psi_tmp2;          //temporäre Ablage
vector<complex<double> > psi_tmp3;          //temporäre Ablage
vector<complex<double> > zz;	            //temporäre Ablage

int thread_number = 8;                      // Anzahl Threads
int N05 = 16;                               // Anzahl Spins
double gesamtZeit = 50;		                // Gesamtzeit t, die simuliert werden soll
double schrittweite;	                    // (in main festgelegt)
int schrittzahl;                            // (in main festgelegt)
int d05 = 1<<N05;                           // Dimension des Hilbertraums



//Zeitentwicklung zweiter Ordnung (zeitunabhängiger Hamiltonoperator)
void entwicklung2(vector<double> J_ij[11], vector<double> h_i[4]){
	
    
    // Für Zeitentwicklung benötigte Suzuki-Trotter Faktoren, die aus J_ij, h_i generiert werden (nicht anfassen ;))
	vector<complex<double> > trotter_ze[5];
	vector<complex<double> > trotter_J_ij_1[18];
	vector<complex<double> > trotter_J_ij_2[18];
	initial(N05, schrittweite/2.0, trotter_ze, trotter_J_ij_1, trotter_J_ij_2, J_ij, h_i);
	int faktoren_1 = trotter_ze[0].size();
	int faktoren_2 = trotter_J_ij_1[0].size();
	
	
    // Generierung einer Ausgabedatei, wo z.B. während der Zeitentwicklung zu jedem Zeitschritt Observablen abgespeichert werden können
	FILE *outputFile;
	char name[64];
	sprintf (name,"output_time_evolution/Test_mit_N%d_Spins.dat",N05);
	outputFile = fopen(name, "w");
    
    
	double t = 0;
	int fortschritt = 0;
    
    
    // Beginn der Zeitentwicklung
	for (int k = 0; k < schrittzahl; k++){
		
        
        //Fortschrittsanzeige in 5% Schritten
		cout.precision(2);								
		if ((100*k/schrittzahl) % 5 == 0 && 100*k/schrittzahl != fortschritt){
			fortschritt = 100*k/schrittzahl;			
			cout << "Fortschritt: " << (100*k/schrittzahl) << " %" << "\n";
		}
		
    
		//Ausgabe von Observablen in Textdatei alle dt=0.05
		if (k % ((int)(0.01/schrittweite)) == 0){
			
			//Energie-Erwartungswert
            complex<double> energie (0,0);
			H(psi, psi_tmp, d05, J_ij, h_i);
			for (int z = 0; z < d05; z++){
				energie += conj(psi[z])*psi_tmp[z];
			}
            
         
            // Erwartungswerte von s_x aller Spins im Gitter
			vector<complex<double> > sxj;	
			for (int z = 0; z < N05; z++){
				sxj.push_back(complex<double> (0,0));
			}
			
			for (int z = 0; z < N05; z++){
				
				s_x(psi, psi_tmp, d05,z, 1.0, 1);
				
				for (int z2 = 0; z2 < d05; z2++){
					sxj[z] += conj(psi[z2])*psi_tmp[z2];	
				}
			}
			
			// Erwartungswerte von s_y aller Spins im Gitter
			vector<complex<double> > syj;	
			for (int z = 0; z < N05; z++){
				syj.push_back(complex<double> (0,0));
			}
			
			for (int z = 0; z < N05; z++){
				
				s_y(psi, psi_tmp, d05,z, 1.0, 1);
				
				for (int z2 = 0; z2 < d05; z2++){
					syj[z] += conj(psi[z2])*psi_tmp[z2];	
				}
			}
			
			// Erwartungswerte von s_z aller Spins im Gitter
			vector<complex<double> > szj;	
			for (int z = 0; z < N05; z++){
				szj.push_back(complex<double> (0,0));
			}
			
			for (int z = 0; z < N05; z++){
				
				s_z(psi, psi_tmp, d05,z, 1.0, 1);
				
				for (int z2 = 0; z2 < d05; z2++){
					szj[z] += conj(psi[z2])*psi_tmp[z2];	
				}
			}
			
			// Magnetisierung des Untergitters für 4x4 gerade
			complex<double> szH_even (0,0);
			for (int z = 0; z < 3; z=z+2){
				szH_even += szj[z];
				szH_even += szj[z+5];
				szH_even += szj[z+8];
				szH_even += szj[z+13];
			}
			
			// Magnetisierung des Untergitters für 4x4 ungerade
			complex<double> szH_odd (0,0);
			for (int z = 1; z < 4; z=z+2){
				szH_odd += szj[z];
				szH_odd += szj[z+3];
				szH_odd += szj[z+8];
				szH_odd += szj[z+11];
			}
		
			
			// Schreibe in Ausgabedatei (Zeit, Energie, Erwartungswert s_0^z)
			fprintf(outputFile, "%lf\t %.17e\t %.17e\t %.17e\t %.17e\t \n", t, real(energie), real(szj[0]), real(szH_even), real(szH_odd));
			
		}
        
		t += schrittweite;
		
		
		//Zustandsvektor updaten (übrigens zweite Ordnung weil die Trotter Faktoren in jedem Zeitschritt einmal sorum und dann nochmal mit inverser Reihenfolge angewandt werden)
		for(int z = 0; z<faktoren_1; z++){
			spin05(psi, d05, (int) real(trotter_ze[4][z]), trotter_ze[0][z], trotter_ze[1][z], trotter_ze[2][z], trotter_ze[3][z]);
		}
		for(int z = 0; z<faktoren_2; z++){
			spin05spin05(psi, d05, (int) real(trotter_J_ij_1[16][z]), (int) real(trotter_J_ij_1[17][z]), trotter_J_ij_1[0][z], trotter_J_ij_1[1][z], trotter_J_ij_1[2][z], trotter_J_ij_1[3][z], trotter_J_ij_1[4][z], trotter_J_ij_1[5][z], trotter_J_ij_1[6][z], trotter_J_ij_1[7][z], trotter_J_ij_1[8][z], trotter_J_ij_1[9][z], trotter_J_ij_1[10][z], trotter_J_ij_1[11][z], trotter_J_ij_1[12][z], trotter_J_ij_1[13][z], trotter_J_ij_1[14][z], trotter_J_ij_1[15][z]);
		}
		for(int z = faktoren_2 -1; z>(-1); z--){
			spin05spin05(psi, d05, (int) real(trotter_J_ij_2[16][z]), (int) real(trotter_J_ij_2[17][z]), trotter_J_ij_2[0][z], trotter_J_ij_2[1][z], trotter_J_ij_2[2][z], trotter_J_ij_2[3][z], trotter_J_ij_2[4][z], trotter_J_ij_2[5][z], trotter_J_ij_2[6][z], trotter_J_ij_2[7][z], trotter_J_ij_2[8][z], trotter_J_ij_2[9][z], trotter_J_ij_2[10][z], trotter_J_ij_2[11][z], trotter_J_ij_2[12][z], trotter_J_ij_2[13][z], trotter_J_ij_2[14][z], trotter_J_ij_2[15][z]);
		}
		for(int z = faktoren_1 -1; z>(-1); z--){
			spin05(psi, d05, (int) real(trotter_ze[4][z]), trotter_ze[0][z], trotter_ze[1][z], trotter_ze[2][z], trotter_ze[3][z]);
		}
	} 
	fclose(outputFile);
}



//Entwicklung zweiter Ordnung (zeitabhängig, nur ein Zeitschritt)
void entwicklung2_t(){
	
	vector<complex<double> > trotter_ze[5];					// = Suzuki-Trotter Zeeman-Faktoren
	vector<complex<double> > trotter_J_ij_1[18];			// = trotter_xx_yy_zz * trotter_xy_yx * trotter_xz * ... * trotter_zy
	vector<complex<double> > trotter_J_ij_2[18];			// = trotter_zy * ... * trotter_xz * trotter_xy_yx * trotter_xx_yy_zz
	
	vector<double> h_i[4];								
	vector<double> J_ij[11];
	
	FILE *outputFile;
	char name[64];
	sprintf (name,"output_time_evolution/Test_mit_N%d_Spins.dat",N05);
	outputFile = fopen(name, "w");
	double t = 0;
	
	int fortschritt = 0;									//Fortschrittsanzeige in 5% Schritten							 
	for (int k = 0; k < schrittzahl; k++){
		
		//Variablen zurücksetzen, um sie anschließend für neues H festzulegen
		for (int z = 0; z<18; z++){
			if (z < 4){
				h_i[z].clear();
			}
			if (z < 5){
				trotter_ze[z].clear();
			}
			if (z < 11){
				J_ij[z].clear();
			}
			if (z < 18){
				trotter_J_ij_1[z].clear();
				trotter_J_ij_2[z].clear();
			}
		}
			
			
		//Neue Wechselwirkungen 
		
		double J = -1.0;
		double h = sin(t);		//(z.B. zeitabhängiges h)
	
		//Zeilen (offene Randbedingungen) 
		for(int z = 0; z < 3; z++){
		
			//Zeilen
			J_ij[0].push_back(z);               // Spinindex1
			J_ij[1].push_back((z+1));           // Spinindex2
			J_ij[2].push_back(J);            // J_xx
			J_ij[3].push_back(0);               // J_xy
			J_ij[4].push_back(0);               // J_xz
			J_ij[5].push_back(0);               // J_yx
			J_ij[6].push_back(J);            // J_yy
			J_ij[7].push_back(0);               // J_yz
			J_ij[8].push_back(0);               // J_zx
			J_ij[9].push_back(0);               // J_zy
			J_ij[10].push_back(J);           // J_zz
		
			J_ij[0].push_back(z+4);             // Spinindex1
			J_ij[1].push_back((z+4+1));         // Spinindex2
			J_ij[2].push_back(J);            // J_xx
			J_ij[3].push_back(0);               // J_xy
			J_ij[4].push_back(0);               // J_xz
			J_ij[5].push_back(0);               // J_yx
			J_ij[6].push_back(J);            // J_yy
			J_ij[7].push_back(0);               // J_yz
			J_ij[8].push_back(0);               // J_zx
			J_ij[9].push_back(0);               // J_zy
			J_ij[10].push_back(J);           // J_zz
		
			J_ij[0].push_back(z+8);             // Spinindex1
			J_ij[1].push_back((z+8+1));         // Spinindex2
			J_ij[2].push_back(J);            // J_xx
			J_ij[3].push_back(0);               // J_xy
			J_ij[4].push_back(0);               // J_xz
			J_ij[5].push_back(0);               // J_yx
			J_ij[6].push_back(J);            // J_yy
			J_ij[7].push_back(0);               // J_yz
			J_ij[8].push_back(0);               // J_zx
			J_ij[9].push_back(0);               // J_zy
			J_ij[10].push_back(J);           // J_zz
		
			J_ij[0].push_back(z+12);            // Spinindex1
			J_ij[1].push_back((z+12+1));        // Spinindex2
			J_ij[2].push_back(J);            // J_xx
			J_ij[3].push_back(0);               // J_xy
			J_ij[4].push_back(0);               // J_xz
			J_ij[5].push_back(0);               // J_yx
			J_ij[6].push_back(J);            // J_yy
			J_ij[7].push_back(0);               // J_yz
			J_ij[8].push_back(0);               // J_zx
			J_ij[9].push_back(0);               // J_zy
			J_ij[10].push_back(J);           // J_zz
		
		}
	
		//Spalten (offene Randbedingungen) 
		for(int z = 0; z < 4; z++){
		
			//Zeilen
			J_ij[0].push_back(z);               // Spinindex1
			J_ij[1].push_back((z+4));           // Spinindex2
			J_ij[2].push_back(J);            // J_xx
			J_ij[3].push_back(0);               // J_xy
			J_ij[4].push_back(0);               // J_xz
			J_ij[5].push_back(0);               // J_yx
			J_ij[6].push_back(J);            // J_yy
			J_ij[7].push_back(0);               // J_yz
			J_ij[8].push_back(0);               // J_zx
			J_ij[9].push_back(0);               // J_zy
			J_ij[10].push_back(J);           // J_zz
		
			J_ij[0].push_back(z+4);             // Spinindex1
			J_ij[1].push_back((z+8));          // Spinindex2
			J_ij[2].push_back(J);            // J_xx
			J_ij[3].push_back(0);               // J_xy
			J_ij[4].push_back(0);               // J_xz
			J_ij[5].push_back(0);               // J_yx
			J_ij[6].push_back(J);            // J_yy
			J_ij[7].push_back(0);               // J_yz
			J_ij[8].push_back(0);               // J_zx
			J_ij[9].push_back(0);               // J_zy
			J_ij[10].push_back(J);           // J_zz
		
			J_ij[0].push_back(z+8);             // Spinindex1
			J_ij[1].push_back((z+12));          // Spinindex2
			J_ij[2].push_back(J);            // J_xx
			J_ij[3].push_back(0);               // J_xy
			J_ij[4].push_back(0);               // J_xz
			J_ij[5].push_back(0);               // J_yx
			J_ij[6].push_back(J);            // J_yy
			J_ij[7].push_back(0);               // J_yz
			J_ij[8].push_back(0);               // J_zx
			J_ij[9].push_back(0);               // J_zy
			J_ij[10].push_back(J);           // J_zz
	
		}
	
		//Randbedingungen schließen
	
		//Zeilen
		J_ij[0].push_back(0);               // Spinindex1
		J_ij[1].push_back(3);               // Spinindex2
		J_ij[2].push_back(J);            // J_xx
		J_ij[3].push_back(0);               // J_xy
		J_ij[4].push_back(0);               // J_xz
		J_ij[5].push_back(0);               // J_yx
		J_ij[6].push_back(J);            // J_yy
		J_ij[7].push_back(0);               // J_yz
		J_ij[8].push_back(0);               // J_zx
		J_ij[9].push_back(0);               // J_zy
		J_ij[10].push_back(J);           // J_zz
	
		J_ij[0].push_back(4);               // Spinindex1
		J_ij[1].push_back(7);               // Spinindex2
		J_ij[2].push_back(J);            // J_xx
		J_ij[3].push_back(0);               // J_xy
		J_ij[4].push_back(0);               // J_xz
		J_ij[5].push_back(0);               // J_yx
		J_ij[6].push_back(J);            // J_yy
		J_ij[7].push_back(0);               // J_yz
		J_ij[8].push_back(0);               // J_zx
		J_ij[9].push_back(0);               // J_zy
		J_ij[10].push_back(J);           // J_zz
	
		J_ij[0].push_back(8);               // Spinindex1
		J_ij[1].push_back(11);               // Spinindex2
		J_ij[2].push_back(J);            // J_xx
		J_ij[3].push_back(0);               // J_xy
		J_ij[4].push_back(0);               // J_xz
		J_ij[5].push_back(0);               // J_yx
		J_ij[6].push_back(J);            // J_yy
		J_ij[7].push_back(0);               // J_yz
		J_ij[8].push_back(0);               // J_zx
		J_ij[9].push_back(0);               // J_zy
		J_ij[10].push_back(J);           // J_zz
	
		J_ij[0].push_back(12);               // Spinindex1
		J_ij[1].push_back(15);               // Spinindex2
		J_ij[2].push_back(J);            // J_xx
		J_ij[3].push_back(0);               // J_xy
		J_ij[4].push_back(0);               // J_xz
		J_ij[5].push_back(0);               // J_yx
		J_ij[6].push_back(J);            // J_yy
		J_ij[7].push_back(0);               // J_yz
		J_ij[8].push_back(0);               // J_zx
		J_ij[9].push_back(0);               // J_zy
		J_ij[10].push_back(J);           // J_zz
		
	
		//Spalten
		J_ij[0].push_back(0);               // Spinindex1
		J_ij[1].push_back(12);               // Spinindex2
		J_ij[2].push_back(J);            // J_xx
		J_ij[3].push_back(0);               // J_xy
		J_ij[4].push_back(0);               // J_xz
		J_ij[5].push_back(0);               // J_yx
		J_ij[6].push_back(J);            // J_yy
		J_ij[7].push_back(0);               // J_yz
		J_ij[8].push_back(0);               // J_zx
		J_ij[9].push_back(0);               // J_zy
		J_ij[10].push_back(J);           // J_zz
	
		J_ij[0].push_back(1);               // Spinindex1
		J_ij[1].push_back(13);               // Spinindex2
		J_ij[2].push_back(J);            // J_xx
		J_ij[3].push_back(0);               // J_xy
		J_ij[4].push_back(0);               // J_xz
		J_ij[5].push_back(0);               // J_yx
		J_ij[6].push_back(J);            // J_yy
		J_ij[7].push_back(0);               // J_yz
		J_ij[8].push_back(0);               // J_zx
		J_ij[9].push_back(0);               // J_zy
		J_ij[10].push_back(J);           // J_zz
	
		J_ij[0].push_back(2);               // Spinindex1
		J_ij[1].push_back(14);               // Spinindex2
		J_ij[2].push_back(J);            // J_xx
		J_ij[3].push_back(0);               // J_xy
		J_ij[4].push_back(0);               // J_xz
		J_ij[5].push_back(0);               // J_yx
		J_ij[6].push_back(J);            // J_yy
		J_ij[7].push_back(0);               // J_yz
		J_ij[8].push_back(0);               // J_zx
		J_ij[9].push_back(0);               // J_zy
		J_ij[10].push_back(J);           // J_zz
	
		J_ij[0].push_back(3);               // Spinindex1
		J_ij[1].push_back(15);               // Spinindex2
		J_ij[2].push_back(J);            // J_xx
		J_ij[3].push_back(0);               // J_xy
		J_ij[4].push_back(0);               // J_xz
		J_ij[5].push_back(0);               // J_yx
		J_ij[6].push_back(J);            // J_yy
		J_ij[7].push_back(0);               // J_yz
		J_ij[8].push_back(0);               // J_zx
		J_ij[9].push_back(0);               // J_zy
		J_ij[10].push_back(J);           // J_zz
		
		
	
		//Neues externes Magnetfeld für alle Spins
		for(int z = 0; z<N05; z++){
	
			h_i[0].push_back(z);                // Spinindex
			h_i[1].push_back(0);                // x-Feld
			h_i[2].push_back(0);                // y-Feld
			h_i[3].push_back(h);              // z-Feld
		
		} 
		
		
		//initialisiere neue Trotter-Faktoren		
		initial(N05, schrittweite/2.0, trotter_ze, trotter_J_ij_1, trotter_J_ij_2, J_ij, h_i);
		
		
		int faktoren_1 = trotter_ze[0].size();
		int faktoren_2 = trotter_J_ij_1[0].size();
	
		cout.precision(3);	
		if ((100*k/schrittzahl) % 5 == 0 && 100*k/schrittzahl != fortschritt){
			fortschritt = 100*k/schrittzahl;			
			cout << "Fortschritt: " << (100*k/schrittzahl) << " %" << "\n";
		}
		
		
		//Ausgabe von Observablen in Textdatei alle dt=0.05 
		if (k % ((int)(0.05/schrittweite)) == 0){
			
			//Energie-Erwartungswert
            complex<double> energie (0,0);
			H(psi, psi_tmp, d05, J_ij, h_i);
			for (int z = 0; z < d05; z++){
				energie += conj(psi[z])*psi_tmp[z];
			}
            
            
            //s_z Erwartungswert vom ersten Spin
            complex<double> sz (0,0);
            int spinindex = 0;
            double konstante = 1.0;
            
            // Hier wird in psi_tmp = konstante * s_0^x * psi gespeichert (letzte Zahl = 1 sagt,
            // dass vorheriges psi_tmp überschrieben werden soll, ansonsten psi_tmp += konstante * s_0^x * psi)
            s_z(psi, psi_tmp, d05, spinindex, konstante, 1);
            
            // In psi_tmp steht jetzt also s^x * psi, für den Erwartungswert muss jetzt noch Skalarprodukt <psi, psi_tmp> ausgerechnet werden
            for (int z = 0; z < d05; z++){
                sz += conj(psi[z])*psi_tmp[z];
            }
        
			// Schreibe in Ausgabedatei (Zeit, Energie, Erwartungswert s_0^z)
			fprintf(outputFile, "%lf\t %.17e\t %.17e\t \n", t, real(energie), real(sz));
		}
		
		
		//Zustandsvektor updaten
		for(int z = 0; z<faktoren_1; z++){
			spin05(psi, d05, (int) real(trotter_ze[4][z]), trotter_ze[0][z], trotter_ze[1][z], trotter_ze[2][z], trotter_ze[3][z]);
		}
		for(int z = 0; z<faktoren_2; z++){
			spin05spin05(psi, d05, (int) real(trotter_J_ij_1[16][z]), (int) real(trotter_J_ij_1[17][z]), trotter_J_ij_1[0][z], trotter_J_ij_1[1][z], trotter_J_ij_1[2][z], trotter_J_ij_1[3][z], trotter_J_ij_1[4][z], trotter_J_ij_1[5][z], trotter_J_ij_1[6][z], trotter_J_ij_1[7][z], trotter_J_ij_1[8][z], trotter_J_ij_1[9][z], trotter_J_ij_1[10][z], trotter_J_ij_1[11][z], trotter_J_ij_1[12][z], trotter_J_ij_1[13][z], trotter_J_ij_1[14][z], trotter_J_ij_1[15][z]);
		}
		for(int z = faktoren_2 -1; z>(-1); z--){
			spin05spin05(psi, d05, (int) real(trotter_J_ij_2[16][z]), (int) real(trotter_J_ij_2[17][z]), trotter_J_ij_2[0][z], trotter_J_ij_2[1][z], trotter_J_ij_2[2][z], trotter_J_ij_2[3][z], trotter_J_ij_2[4][z], trotter_J_ij_2[5][z], trotter_J_ij_2[6][z], trotter_J_ij_2[7][z], trotter_J_ij_2[8][z], trotter_J_ij_2[9][z], trotter_J_ij_2[10][z], trotter_J_ij_2[11][z], trotter_J_ij_2[12][z], trotter_J_ij_2[13][z], trotter_J_ij_2[14][z], trotter_J_ij_2[15][z]);
		}
		for(int z = faktoren_1 -1; z>(-1); z--){
			spin05(psi, d05, (int) real(trotter_ze[4][z]), trotter_ze[0][z], trotter_ze[1][z], trotter_ze[2][z], trotter_ze[3][z]);
		}
		
		t += schrittweite;
		
	} 
	
	fclose(outputFile);
}



//Initialisiere Anfangszustand und Hamiltonoperator
void initialize(vector<double> J_ij[11], vector<double> h_i[4]){
	
    // Zur Sicherheit alle Variablen noch einmal komplett zurücksetzen
	psi.clear();
	psi_tmp.clear();
	psi_tmp2.clear();
	psi_tmp3.clear();
	zz.clear();
    
	for (int z = 0; z<11; z++){
		J_ij[z].clear();
		if (z < 4){
			h_i[z].clear();
		}
	}
    
    // Zustandsvektor erstmal komplett mit 0en füllen damit er die richtige Größe hat (Hilbertraumdimension)
	for (int z = 0; z < d05; z++){
		psi.push_back(complex<double> (0,0));
		psi_tmp.push_back(complex<double> (0,0));
		psi_tmp2.push_back(complex<double> (0,0));
		psi_tmp3.push_back(complex<double> (0,0));
		zz.push_back(complex<double> (0,0));
	}

	// Beispiel (antiferro) Heisenberg-Wechselwirkungen auf einem 4x4-Gitter von Spins (Vorzeichenkonvention: H = -J*s*s ... - h*s
	
	double J1 = -0.9;
	double J2 = -1.0;
	double h = 0.0;
	
	//Zeilen (offene Randbedingungen) 
	for(int z = 0; z < 3; z++){
		
		//Zeilen
        J_ij[0].push_back(z);               // Spinindex1
        J_ij[1].push_back((z+1));           // Spinindex2
		J_ij[2].push_back(J1);            // J_xx
		J_ij[3].push_back(0);               // J_xy
		J_ij[4].push_back(0);               // J_xz
		J_ij[5].push_back(0);               // J_yx
		J_ij[6].push_back(J1);            // J_yy
		J_ij[7].push_back(0);               // J_yz
		J_ij[8].push_back(0);               // J_zx
		J_ij[9].push_back(0);               // J_zy
		J_ij[10].push_back(J2);           // J_zz
		
		J_ij[0].push_back(z+4);             // Spinindex1
        J_ij[1].push_back((z+4+1));         // Spinindex2
		J_ij[2].push_back(J1);            // J_xx
		J_ij[3].push_back(0);               // J_xy
		J_ij[4].push_back(0);               // J_xz
		J_ij[5].push_back(0);               // J_yx
		J_ij[6].push_back(J1);            // J_yy
		J_ij[7].push_back(0);               // J_yz
		J_ij[8].push_back(0);               // J_zx
		J_ij[9].push_back(0);               // J_zy
		J_ij[10].push_back(J2);           // J_zz
		
		J_ij[0].push_back(z+8);             // Spinindex1
        J_ij[1].push_back((z+8+1));         // Spinindex2
		J_ij[2].push_back(J1);            // J_xx
		J_ij[3].push_back(0);               // J_xy
		J_ij[4].push_back(0);               // J_xz
		J_ij[5].push_back(0);               // J_yx
		J_ij[6].push_back(J1);            // J_yy
		J_ij[7].push_back(0);               // J_yz
		J_ij[8].push_back(0);               // J_zx
		J_ij[9].push_back(0);               // J_zy
		J_ij[10].push_back(J2);           // J_zz
		
		J_ij[0].push_back(z+12);            // Spinindex1
        J_ij[1].push_back((z+12+1));        // Spinindex2
		J_ij[2].push_back(J1);            // J_xx
		J_ij[3].push_back(0);               // J_xy
		J_ij[4].push_back(0);               // J_xz
		J_ij[5].push_back(0);               // J_yx
		J_ij[6].push_back(J1);            // J_yy
		J_ij[7].push_back(0);               // J_yz
		J_ij[8].push_back(0);               // J_zx
		J_ij[9].push_back(0);               // J_zy
		J_ij[10].push_back(J2);           // J_zz
		
	}
	
	//Spalten (offene Randbedingungen) 
	for(int z = 0; z < 4; z++){
		
		//Zeilen
        J_ij[0].push_back(z);               // Spinindex1
        J_ij[1].push_back((z+4));           // Spinindex2
		J_ij[2].push_back(J1);            // J_xx
		J_ij[3].push_back(0);               // J_xy
		J_ij[4].push_back(0);               // J_xz
		J_ij[5].push_back(0);               // J_yx
		J_ij[6].push_back(J1);            // J_yy
		J_ij[7].push_back(0);               // J_yz
		J_ij[8].push_back(0);               // J_zx
		J_ij[9].push_back(0);               // J_zy
		J_ij[10].push_back(J2);           // J_zz
		
		J_ij[0].push_back(z+4);             // Spinindex1
        J_ij[1].push_back((z+8));          // Spinindex2
		J_ij[2].push_back(J1);            // J_xx
		J_ij[3].push_back(0);               // J_xy
		J_ij[4].push_back(0);               // J_xz
		J_ij[5].push_back(0);               // J_yx
		J_ij[6].push_back(J1);            // J_yy
		J_ij[7].push_back(0);               // J_yz
		J_ij[8].push_back(0);               // J_zx
		J_ij[9].push_back(0);               // J_zy
		J_ij[10].push_back(J2);           // J_zz
		
		J_ij[0].push_back(z+8);             // Spinindex1
        J_ij[1].push_back((z+12));          // Spinindex2
		J_ij[2].push_back(J1);            // J_xx
		J_ij[3].push_back(0);               // J_xy
		J_ij[4].push_back(0);               // J_xz
		J_ij[5].push_back(0);               // J_yx
		J_ij[6].push_back(J1);            // J_yy
		J_ij[7].push_back(0);               // J_yz
		J_ij[8].push_back(0);               // J_zx
		J_ij[9].push_back(0);               // J_zy
		J_ij[10].push_back(J2);           // J_zz
	
	}
	
	
	//Randbedingungen schließen
	
	//Zeilen
    J_ij[0].push_back(0);               // Spinindex1
    J_ij[1].push_back(3);               // Spinindex2
	J_ij[2].push_back(J1);            // J_xx
	J_ij[3].push_back(0);               // J_xy
	J_ij[4].push_back(0);               // J_xz
	J_ij[5].push_back(0);               // J_yx
	J_ij[6].push_back(J1);            // J_yy
	J_ij[7].push_back(0);               // J_yz
	J_ij[8].push_back(0);               // J_zx
	J_ij[9].push_back(0);               // J_zy
	J_ij[10].push_back(J2);           // J_zz
	
	J_ij[0].push_back(4);               // Spinindex1
    J_ij[1].push_back(7);               // Spinindex2
	J_ij[2].push_back(J1);            // J_xx
	J_ij[3].push_back(0);               // J_xy
	J_ij[4].push_back(0);               // J_xz
	J_ij[5].push_back(0);               // J_yx
	J_ij[6].push_back(J1);            // J_yy
	J_ij[7].push_back(0);               // J_yz
	J_ij[8].push_back(0);               // J_zx
	J_ij[9].push_back(0);               // J_zy
	J_ij[10].push_back(J2);           // J_zz
	
	J_ij[0].push_back(8);               // Spinindex1
    J_ij[1].push_back(11);               // Spinindex2
	J_ij[2].push_back(J1);            // J_xx
	J_ij[3].push_back(0);               // J_xy
	J_ij[4].push_back(0);               // J_xz
	J_ij[5].push_back(0);               // J_yx
	J_ij[6].push_back(J1);            // J_yy
	J_ij[7].push_back(0);               // J_yz
	J_ij[8].push_back(0);               // J_zx
	J_ij[9].push_back(0);               // J_zy
	J_ij[10].push_back(J2);           // J_zz
	
	J_ij[0].push_back(12);               // Spinindex1
    J_ij[1].push_back(15);               // Spinindex2
	J_ij[2].push_back(J1);            // J_xx
	J_ij[3].push_back(0);               // J_xy
	J_ij[4].push_back(0);               // J_xz
	J_ij[5].push_back(0);               // J_yx
	J_ij[6].push_back(J1);            // J_yy
	J_ij[7].push_back(0);               // J_yz
	J_ij[8].push_back(0);               // J_zx
	J_ij[9].push_back(0);               // J_zy
	J_ij[10].push_back(J2);           // J_zz
	
	
	//Spalten
    J_ij[0].push_back(0);               // Spinindex1
    J_ij[1].push_back(12);               // Spinindex2
	J_ij[2].push_back(J1);            // J_xx
	J_ij[3].push_back(0);               // J_xy
	J_ij[4].push_back(0);               // J_xz
	J_ij[5].push_back(0);               // J_yx
	J_ij[6].push_back(J1);            // J_yy
	J_ij[7].push_back(0);               // J_yz
	J_ij[8].push_back(0);               // J_zx
	J_ij[9].push_back(0);               // J_zy
	J_ij[10].push_back(J2);           // J_zz
	
	J_ij[0].push_back(1);               // Spinindex1
    J_ij[1].push_back(13);               // Spinindex2
	J_ij[2].push_back(J1);            // J_xx
	J_ij[3].push_back(0);               // J_xy
	J_ij[4].push_back(0);               // J_xz
	J_ij[5].push_back(0);               // J_yx
	J_ij[6].push_back(J1);            // J_yy
	J_ij[7].push_back(0);               // J_yz
	J_ij[8].push_back(0);               // J_zx
	J_ij[9].push_back(0);               // J_zy
	J_ij[10].push_back(J2);           // J_zz
	
	J_ij[0].push_back(2);               // Spinindex1
    J_ij[1].push_back(14);               // Spinindex2
	J_ij[2].push_back(J1);            // J_xx
	J_ij[3].push_back(0);               // J_xy
	J_ij[4].push_back(0);               // J_xz
	J_ij[5].push_back(0);               // J_yx
	J_ij[6].push_back(J1);            // J_yy
	J_ij[7].push_back(0);               // J_yz
	J_ij[8].push_back(0);               // J_zx
	J_ij[9].push_back(0);               // J_zy
	J_ij[10].push_back(J2);           // J_zz
	
	J_ij[0].push_back(3);               // Spinindex1
    J_ij[1].push_back(15);               // Spinindex2
	J_ij[2].push_back(J1);            // J_xx
	J_ij[3].push_back(0);               // J_xy
	J_ij[4].push_back(0);               // J_xz
	J_ij[5].push_back(0);               // J_yx
	J_ij[6].push_back(J1);            // J_yy
	J_ij[7].push_back(0);               // J_yz
	J_ij[8].push_back(0);               // J_zx
	J_ij[9].push_back(0);               // J_zy
	J_ij[10].push_back(J2);           // J_zz
		

	
	
    // Magnetfelder
	for(int z = 0; z<N05; z++){
	
		h_i[0].push_back(z);                // Spinindex
		h_i[1].push_back(0);               // x-Feld
		h_i[2].push_back(h);               // y-Feld
		h_i[3].push_back(0);              // z-Feld
		
	} 
	
	// Lese externen Koeffizienten für den Anfangszustand ein 
    std::ifstream infile("/home/asliddin/CPP/Exact_DIagonalization/Mit_Florian/output/coefficients_state_0.txt");
	
    if (!infile.is_open()) {
        std::cout << "Fehler beim Öffnen der Datei" << std::endl;
    }
    std::vector<std::complex<double>> coef;
    std::string line;
    double real, imag;
    while (std::getline(infile, line)) {
        // Kommentarzeilen oder leere Zeilen überspringen
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        if (iss >> real >> imag) {
            coef.emplace_back(real, imag);
        } else {
            std::cerr << "Fehler beim Parsen der Zeile: " << line << std::endl;
        }
    }
	infile.close();
	
	
	for (int i = 0; i < d05; i++) {
		psi[i] = coef[i];
	}
	
	//Anfangszustand festlegen (Beispiel: Spins sind alternierend)
	//psi[23130] = 1.0; 
	
	



	/*// Anfangszustand, in diesem Beispiel ein Produktzustand wo alle Spin-Erwartungswerte in einem Intervall delta*N05 in der xy-Ebene verteilt sind
	double delta = 1.0*M_PI/(N05);
	for(int z = 0; z < d05; z++){

		double phase = 0;
		double koeff = 1.0;
		for(int z2 = 0; z2 < N05; z2++){	
			int dummy = 1<<z2;
			if ((dummy & z) != 0){
				phase += z2*delta;
				koeff *= sqrt(0.5);
			}
			else{
				koeff *= sqrt(0.5);
			}
		}		
		psi[z] = complex<double> (cos(phase)*koeff,sin(phase)*koeff);
	}*/
	
	
	
	
}



// main des Programms, also das was wirklich ausgeführt wird und die vorherigen Funktionen aufgerufen werden
int main(void){
	
	omp_set_num_threads(thread_number);

	vector<double> h_i[4];								       // Liste für lokale Magnetfelder (Spinindex, x-Feld, y-Feld, z-Feld)
	vector<double> J_ij[11];                                   // Wechselwirkungen zwischen Spins (Spindex1, Spinindex2, J-Matrix)
	
	schrittweite = 0.01;
    schrittzahl = ((int) (gesamtZeit/schrittweite)) + 1;
	
	initialize(J_ij, h_i);                                     // In der Funktion wird der Hamiltonoperator festgelegt (für den Fall H zeitunabhängig)
															   // und auch der Anfangszustand wird festgelegt, deshalb muss die Funktion auf jeden Fall ausgeführt werden
	
	entwicklung2(J_ij, h_i);                                   // Zeitentwicklung (H nicht zeitabhängig)
	//entwicklung2_t();										   // Zeitentwicklung (H zeitabhängig, in der Funktion festlegen)
	
	psi.clear();                                               // Rechnung fertig, Arbeitsspeicher wieder löschen
	psi_tmp.clear();
	psi_tmp2.clear();
	psi_tmp3.clear();
	zz.clear();
		
	return 0;
}
