//                               
//  _____ _     _    __    _____ 
// |     |_|___|_|  |  |  | __  |
// | | | | |   | |  |  |__| __ -|
// |_|_|_|_|_|_|_|  |_____|_____|                                 
//

// MIT License
//
// Copyright (c) 2021 J-Y Delenne, V. Richefeu, F. Radjai
// Jean-Yves Delenne <jean-yves.delenne@inrae.fr> - INRAE, University of Montpellier
// Vincent Richefeu <vincent.richefeu@3sr-grenoble.fr> - University Grenoble Alpes
// Farhang Radjai <franck.radjai@umontpellier.fr> - CNRS, University of Montpellier
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <iostream>
#include <fstream>

using namespace std;

// *** Parameters are setted has const values ********************************
// Fluid domain
const int lx=200;
const int ly=20;
// Boundary conditions
const double rhoin=1.0003;
const double rhoout=1.0;
// Relaxation time
const double tau=0.55;
// Output results
const int period=10000;
// ***************************************************************************

// D2Q9 
const double w[9]={4./9.,1./9.,1./9.,1./9.,1./9.,1./36.,1./36.,1./36.,1./36.};
const int cx[9]={0,1,0,-1, 0,1,-1,-1, 1};
const int cy[9]={0,0,1, 0,-1,1, 1,-1,-1};
double f[lx][ly][9],fout[lx][ly][9];
unsigned long int t;

// *** See full definitions below *******************************************
void init_Densities();
void collisions_BGK();
void streaming_BasicSwap();
void boundaries_yBounceBack();
void boundaries_xPoiseuillePressure();
void check_Densities();
void output_Velocity(int);
void output_pressure_velocity_fields();

// *** Main Program *********************************************************
int main(int argc, char *argv[]) {
	init_Densities();
	for (t=0;;t++) {
		collisions_BGK();
		streaming_BasicSwap();
		boundaries_yBounceBack();
		boundaries_xPoiseuillePressure();
		if (t%period==0) {
			cout << "Iteration " << t/period << endl << flush;
			// Output profile for x = lx/2
			output_Velocity(lx/2);
			output_pressure_velocity_fields();
		}
	}
	return 0;
}

// *** Begining of LBM functions ********************************************

// Initialise densities
void init_Densities() {
	for (int x=0; x<lx;x++)
		for (int y=0; y<ly;y++)
			for (int i=0; i<9;i++) f[x][y][i]=w[i];				  
}

// Compute collisions according to BGK
void collisions_BGK() {
	for (int x=0; x<lx;x++)
		for (int y=0; y<ly;y++) {
		double rho=0.;
		double ux=0., uy=0.;

		for (int i=0; i<9;i++) rho+=f[x][y][i];
		for (int i=0; i<9;i++) { 
			ux+=f[x][y][i]*cx[i];
			uy+=f[x][y][i]*cy[i];
		}
		ux/=rho;
		uy/=rho;
		double u2=ux*ux+uy*uy;
		for (int i=0; i<9;i++) { 
			double cu=cx[i]*ux+cy[i]*uy;
			fout[x][y][i]=f[x][y][i]+(w[i]*rho*(1.+3*cu+4.5*cu*cu-1.5*u2)-f[x][y][i])/tau;
		}
	}
}

// Streaming
void streaming_BasicSwap() {
	for (int x=0; x<lx;x++)
		for (int y=0; y<ly;y++)
			for (int i=0; i<9;i++)
				if(x+cx[i]>=0 && y+cy[i]>=0 && x+cx[i]<lx && y+cy[i]<ly)
					f[x+cx[i]][y+cy[i]][i]=fout[x][y][i];
}

// Boundary conditions
void boundaries_yBounceBack() {
	for (int x=0;x<lx;x++) {
	// Bounce back at y = 0
			f[x][0][5]=f[x][0][7];
			f[x][0][2]=f[x][0][4];
			f[x][0][6]=f[x][0][8];
			f[x][0][3]=f[x][0][1];
	// Bounce back at y = ly - 1 
			f[x][ly-1][7]=f[x][ly-1][5];
			f[x][ly-1][4]=f[x][ly-1][2];
			f[x][ly-1][8]=f[x][ly-1][6];
			f[x][ly-1][3]=f[x][ly-1][1];
	}
}

void boundaries_xPoiseuillePressure() {
	// Pressure boundary at x = 0
	for (int y=1; y<ly-1;y++) {
		double qin=rhoin-(f[0][y][0]+f[0][y][2]+f[0][y][4]+2.*(f[0][y][3]+f[0][y][6]+f[0][y][7]));
		f[0][y][1]=f[0][y][3]+qin*2./3.;
		f[0][y][5]=f[0][y][7]+qin/6.-(f[0][y][2]-f[0][y][4])/2.;
		f[0][y][8]=f[0][y][6]+qin/6.+(f[0][y][2]-f[0][y][4])/2.;
	}
	// Pressure boundary at x = lx - 1
	for (int y=1; y<ly-1;y++) {
		double qout=f[lx-1][y][0]+f[lx-1][y][2]+f[lx-1][y][4]+2.*(f[lx-1][y][1]+f[lx-1][y][8]+f[lx-1][y][5])-rhoout;
		f[lx-1][y][3]=f[lx-1][y][1]-qout*2./3.;
		f[lx-1][y][7]=f[lx-1][y][5]-qout/6.+(f[lx-1][y][2]-f[lx-1][y][4])/2.;
		f[lx-1][y][6]=f[lx-1][y][8]-qout/6.-(f[lx-1][y][2]-f[lx-1][y][4])/2.;
	}	
}

// *** End of LBM functions *************************************************

// Outputs: velocity profile and pressure and velocity fields ***************

// Output velocity profile
void output_Velocity(int xpos) {
	char fname[30];
	sprintf(fname,"output_U_Profile_%03ld.txt",t/period);
	ofstream file(fname);
	file << "# y ux" << endl;
	for (int y=0;y<ly;y++) {
		double ux,rho;
		ux=0.; rho=0.;
		for (int i=0; i<9;i++) rho += f[xpos][y][i];
		for (int i=0; i<9;i++) ux += f[xpos][y][i]*cx[i];
		ux/=rho;
		file << y << " " << ux;
	}
}

// Save pressure and velocity fields in paraview format
void output_pressure_velocity_fields() {
	int x,y,i;
    double resolxy;
	double P,ux,uy;
	char nomfic[25];
	FILE * sortie;
	
	sprintf(nomfic,"output_P_U_Fields_%.3i.vtk",(int) t/period);
	resolxy=1./lx; 
	sortie = fopen(nomfic, "w");
	fprintf(sortie,"# vtk DataFile Version 2.0\n");
	fprintf(sortie,"Sortie domaine LB+LINK t %li\n",t);
	fprintf(sortie,"ASCII\n");
	fprintf(sortie,"DATASET RECTILINEAR_GRID\n");
	fprintf(sortie,"DIMENSIONS %d %d 1\n",lx,ly);
	fprintf(sortie,"X_COORDINATES %d float\n",lx);
	for(i=0;i<=lx-1;i++) fprintf(sortie,"%e ",(float)i * resolxy);
	fprintf(sortie,"\n");
	fprintf(sortie,"Y_COORDINATES %d float\n",ly);
	for(i=0;i<=ly-1;i++) fprintf(sortie,"%e ",(float)i * resolxy);
	fprintf(sortie,"\n");	
	fprintf(sortie,"Z_COORDINATES 1 float\n");
	fprintf(sortie,"0\n");
	fprintf(sortie,"POINT_DATA %d\n",lx*ly);
	fprintf(sortie,"SCALARS Pression float 1\n");
	fprintf(sortie,"LOOKUP_TABLE default\n");
	for (y=0; y<ly;y++) {
		for (x=0; x<lx;x++) {
			P=0.;
			for (i=0; i<9; i++) P += f[x][y][i];
			P = (1./3.)*(P-1.);	
			fprintf(sortie,"%.4e\n",P);
		}	
	}
	fprintf(sortie,"VECTORS VecVelocity float\n");
	for (y=0; y<ly;y++) {
		for (x=0; x<lx;x++) {
			ux=0.;
			uy=0.;
			for (i=0; i<9; i++) {
				ux+= f[x][y][i]*cx[i];
				uy+= f[x][y][i]*cy[i];
			}
			fprintf(sortie,"%.4e %.4e 0.\n",ux,uy);
		}
	}
	fclose(sortie);
}

