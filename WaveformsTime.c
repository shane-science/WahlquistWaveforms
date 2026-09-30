// ***************************************************************************
//                                                                            
//  PROGRAM :  WaveformsTime                                                  
//  AUTHOR  :  Shane L. Larson                                                
//  VERSION :  3.1 (first github repository version)
//  DATE    :  20 SEPT 2026
//  ORIGINAL:  15 March 2000
//
//  This program is a derivative of my original 'Waveform.p' program.  I need 
//    some more practice writing C code, so I decided to do this rather than  
//    modify the waveform code.  What this code does is uses the Wahlquist    
//    formalism for computing waveforms from arbitrarily oriented binaries.   
//                                                                            
//  Reference Paper: WAHLQUIST, GRG _19_, 1101 [1987]                         
//                                                                            
//  NB: This waveform code is for binaries which are not evolving appreciably 
//      due to the emission of GW (or otherwise).  The orbital parameters are 
//      assumed to be fixed over the course of an orbit.                      
//                                                                            
//  HISTORY                                                                   
//  ========================================================================  
//  1.0 ALPHA       : First RK4, CodeWarrior Pascal compiler (v 4)            
//  2.0             : Rewrite in C to output locus of amplitudes for SGR A*
//  2.1             : Update to be used for binary simulations
//  3.0             : Modified to give const time output                      
//  3.1             : 30 Sept 2026 -- Moved to Github Repository
//
// ***************************************************************************


//  ============ INCLUDE LIBRARIES ============ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

//  ============ CONSTANTS DEFINED ============ 

#define MSun            1.989E30                   //  solar mass in kg   
#define AU              1.495978706910e11          //  m in an AU         
#define lyr             9.46052840488E15           //  m in a lyr         
#define pc              3.08567818585E16           //  m in a pc          
#define YR              31558149.7635456           // seconds per sidereal year    

#define PI              3.14159265358979323846264338327950288     //  Why ask pi?        
#define degrees         1.7453292519943295769236907684886e-2      //  radians per degree
#define radians         57.295779513082320876798154814105         //  degrees per radian

#define G               6.6726E-11                 //  Newton G, SI       
#define C               299792458                  //  Speed of light, SI 

#define kEPS            1.0e-8                     //  tiny number 
#define bEPS            1.0e-6                     //  bisection tolerance 

#define MAX_LENGTH      200                        // max string length


//  ============ FUNCTION PROTOTYPES ============ 

int GetData(char baseName[], double *a, double *ec, double *inc, double *m1, double *m2, double *R, double *thN,
            double *thP, double *phi, double *th1, double *tsim, long int *Nt, long int *fftFLAG);

double OrbTime(double theta, double Ei, double P);

double GetEccentricAnomaly(double Time, double PsiLo, double P, double ecc);

double GetTrueAnomaly(double Psi, double ecc);

void ErrorExit(char routine[], char errorMsg[]);



//  ========================================================================= 



// ***************************************************************************
//                                                                            
//  FUNCTION:  main                                                           
//  TYPE    :  int                                                            
//  RETURNS :  0 (void)                                                       
//                                                                            
//  Program takes the input user parameters, computes the amplitude for       
//     both gravitaitonal wave polarizations for one orbital period, and      
//     and dumps the waveforms to file.                                       
//                                                                            
//                                                                            
//  VARIABLES                                                                 
//  ------------------------------------------------------------------------- 
//  fhandle       : (FILE) handle for output file                             
//  basename[]    : (char) User input filename
//  a             : (double) Semi-major axis [-> meters]
//  e             : (double) Orbital eccentricity                             
//  inc           : (double) Orbital inclination  [-> radians]                
//  m1            : (double) mass 1 of binary [-> kg]                         
//  m2            : (double) mass 2 of binary [-> kg]                         
//  D             : (double) distance to binary [-> meters]                   
//  thn           : (double) anomaly at line of nodes [-> radians]            
//  thp           : (double) anomaly at periapse [-> radians]                 
//  Phi           : (double) orientation of line of nodes on sky [-> radians] 
//  r             : (double) Orbital radius from Shape Equation               
//  x,y           : (double) cartesian coordinate values in orbital plane     
//  deltaTh       : (double) Step value in anomaly, computed from nPts        
//  Period        : (double) Computed orbital period, Kepler III              
//  Root          : (double) eccentricity fraction                            
//  MA            : (double) Mean Anomaly                                     
//  Psi           : (double) Eccentric Anomaly                                
//  Theta         : (double) True Anomaly                                     
//  hp            : (double) h+ polarization amplitude                        
//  hx            : (double) hx polarization amplitude                        
//  Time          : (double) elapsed time since periapse (from mean anomaly)  
//  pMx, cMx      : (double) |h+| and |hx| maximum amplitudes                 
//  Ho            : (double) Wahlquist scaling amplitdue                      
//  A0, A1, A2    : (double) Wahlquist expansion coefficients {Ai}            
//  B0, B1, B2    : (double) Wahlquist expansion coefficients {Bi}            
//  Tmp1, Tmp2, Tmp3 : (double) Temp vars                                     
//                                                                            
// ***************************************************************************

int main(void)
{
	//  ============= VARIABLE DECLARATIONS ============= 
	
	FILE *fhandle, *ffthandleP, *ffthandleX;
	char basename[200],tmpName[200];
	
	double a, e, inc, m1, m2, D, thn, thp, Phi, th1;
	
	long int jj, fftFLAG;              //  looping index 
	
	double r, x, y;
	
	double Period, Root, Psi, Theta;
	
	double hp, hx, Time;
	double pMx, cMx;
	
	double Ho, A0, A1, A2, B0, B1, B2;
	double Tmp1, Tmp2, Tmp3;
	
	long int nT;
	double dT, t1, Tsim;
	
	
	//  =============== START ROUTINE HERE =============== 
	
	printf("Welcome to WAHLQUIST WAVEFORM-TIME DELUXE.\n\n");
	
	//  -------- Initialize variables -------- 
	
	GetData(basename, &a, &e, &inc, &m1, &m2, &D, &thn, &thp, &Phi, &th1, &Tsim, &nT, &fftFLAG);
	
	
	//  ----------- setup Orbital Period and spacing of time points -------------- 
	
	Period = sqrt(pow(2.0*PI,2.0)*pow(a,3.0)/(G*(m1+m2)));  //  orbital period, Kepler III 
	
	
    //  --- first, find time associated with th1 and th2 --- 
    
    t1 = OrbTime(th1, e, Period);  //  time at starting angle, offset for Tperiapse = 0.0 
	
    dT = Tsim/nT;   //  spacing of time interval points 
    
    //  --- Initial value of Eccentric Anomaly, associated with starting point --- 
    Psi = th1 - fabs(th1);
    Psi = GetEccentricAnomaly(t1, Psi, Period, e);
    
	
    printf("True Anomaly = %lf\n",th1/degrees);
    printf("Eccentric Anomaly = %lf\n",Psi/degrees);
    printf("time 1 = %lf\n",t1);
    printf("dT = %lf s\n",dT);
    printf("Period = %lf yrs = %lf sec\n",Period/YR,Period);
    
	
	//  -------- OUTPUT FILE HANDLING -------- 
	
	// printf("\nEnter output file name: ");       //  get filename from user for output file
	// gets(basename);
    
	tmpName[0] = '\0';
	strcat(tmpName,basename);
	strcat(tmpName,".tplot");                // from basename, make a *.tplot file == CSV file for plotting time-domain data
	fhandle = fopen(tmpName,"w");            // Open file for WRITE
	
	
	if (fftFLAG == 1)                         //  If user asked for FFT time series, open file 
	{
		tmpName[0] = '\0';
		strcat(tmpName,basename);
		strcat(tmpName,"_P.time");          // from basename, make a + polarization *.time file == formated file for my FFT codes to read in
		ffthandleP = fopen(tmpName,"w");    //  + Polarization file
        fprintf(ffthandleP,"0\n%ld\n%lf\n",nT,dT); // write fft parameters == nDAT and deltaT
		
		tmpName[0] = '\0';
		strcat(tmpName,basename);
		strcat(tmpName,"_X.time");          // from basename, make a x polarization *.time file == formated file for my FFT codes to read in
		ffthandleX = fopen(tmpName,"w");   //  x Polarization file
		fprintf(ffthandleX,"0\n%ld\n%lf\n",nT,dT); // write fft parameters == nDAT and deltaT
		
		fflush(ffthandleP);
		fflush(ffthandleX);
	}
	
	
	
    //  -------- Output file headers -------- 
	
	fprintf(fhandle,"# Gravitational waveform data file\n");
	fprintf(fhandle,"# Code: WaveformsTime.c v. 4.0 alpha (sll - 06 Aug 2003 build)\n");
	fprintf(fhandle,"# Based on WAHLQUIST, GRG _19_, 1101 [1987]\n#\n");
	
	fprintf(fhandle,"# Binary Parameters used in this file:\n");
	fprintf(fhandle,"# ------------------------------------\n");
	
	fprintf(fhandle,"# Semi-major axis a: %f AU\n",a/AU);
	fprintf(fhandle,"# Mass 1: %f MSun\n",m1/MSun);
	fprintf(fhandle,"# Mass 2: %f MSun\n",m2/MSun);
	fprintf(fhandle,"# Eccentricity e: %f \n",e);
	fprintf(fhandle,"# Distance to binary: %f pc\n",D/pc);
	fprintf(fhandle,"# Inclination i: %f deg\n",inc/degrees);
	fprintf(fhandle,"# Line of nodes on sky phi: %f deg\n",Phi/degrees);
	fprintf(fhandle,"# Anomaly at line of nodes thetaN: %f deg\n",thn/degrees);
	fprintf(fhandle,"# Anomaly at periapse thetaP: %f deg\n",thp/degrees);
	fprintf(fhandle,"# Computed orbital period P: %f s\n",Period);
	fprintf(fhandle,"# Number Time Points: %ld \n",nT);
	
    fprintf(fhandle,"#\n#\n#\n# IDX,TIME(s),TRUE ANOMALY(d),rOrb(m),xOrb(m),yOrb(m),hp,hx\n");
	fflush(fhandle);
	
	
	//  -------- MAIN WAVEFORM COMPUTATION LOOP -------- 
	
	Time = 0.0;                          //  Initial time value          
	
	pMx = 0.0;                          //  Initialize max h+ value     
	cMx = 0.0;                          //  Initialize max hx value     
	
	for (jj = 0; jj < nT; jj++)    //  Generate full waveform, loop over theta 
	{
		
        //  -- get new eccentric anomaly values for current time 
        Psi = GetEccentricAnomaly(Time + t1, Psi, Period, e);
		
		
        //  -- get new true anomaly values for current eccentric anomaly 
        Theta = GetTrueAnomaly(Psi, e);
		
        //  -- compute new orbital radii for given anomaly value 
        Root = sqrt((1.0 - e)/(1.0 + e));
        Psi = 2.0*atan(Root*tan(Theta/2.0));  //  Eccentric Anomaly 
		
        r = a*(1.0 - (e*cos(Psi)));      //  shape equation 
        x = r*cos(Theta);                //  x coord in orbital plane 
        y = r*sin(Theta);                //  y coord in orbital plane 
        
		
        //  ============ GET WAVEFORMS ============ 
        
        //  ----- Calculate Waveforms: Wahlquist Formulae ----- 
        //  ----- Wahlquist, GRG _19_, 1101 [1987]        -----
        //  ----- Equations (30), (31), (32)              -----
		
		//  ---- Scaling Amplitude ---- 
        Tmp1 = 4.0*pow(G,2.0)*m1*m2;
		
        Tmp3 = pow(C,4.0)*a*(1.0 - pow(e,2.0))*D;
        Ho = Tmp1/Tmp3;
        
		//  ---- Ai coefficients  ---- 
        A0 = (-1.0/2.0)*(1.0 + pow(cos(inc),2.0))*cos(2.0*(Theta - thn));
		
        Tmp1 = (1.0/4.0)*pow(sin(inc),2.0)*cos(Theta - thp);
        Tmp2 = (1.0/8.0)*(1.0 + pow(cos(inc),2.0));
        Tmp3 = (5.0*cos(Theta - (2.0*thn) + thp)) + cos((3.0*Theta) - (2.0*thn) - thp);
        A1 = Tmp1 - (Tmp2*Tmp3);
		
        A2 = (1.0/4.0)*((pow(sin(inc),2.0)) - ((1.0 + pow(cos(inc),2.0))*cos(2.0*(thn - thp))));
		
		
		//  ---- Bi coefficients  ---- 
        B0 = (-1.0)*cos(inc)*sin(2.0*(Theta - thn));
		
        Tmp1 = (5.0*sin(Theta - (2.0*thn) + thp)) + sin((3.0*Theta) - (2.0*thn) - thp);
        B1 = (-1.0/4.0)*cos(inc)*Tmp1;
		
        B2 = (1.0/2.0)*cos(inc)*sin(2.0*(thn - thp));
		
		
		//  ---- Polarizatin Amplitudes  ---- 
        Tmp1 = (A0 + (e*A1) + (pow(e,2.0)*A2));
        Tmp2 = (B0 + (e*B1) + (pow(e,2.0)*B2));
		
        hp = Ho*((Tmp1*cos(2*Phi)) - (Tmp2*sin(2.0*Phi)));
		
        hx = Ho*((Tmp1*sin(2*Phi)) + (Tmp2*cos(2.0*Phi)));
        
        
		//  ---- WRITE to File  ---- 
        fprintf(fhandle,"%ld,%f,%f,%e,%e,%e,%e,%e\n",jj+1,Time,Theta*180.0/PI,r,x,y,hp,hx);
        fflush(fhandle);
        
        if (fftFLAG == 1)
        {
            // each file has pairs of numbers on sequential lines: REAL part, then IMAG part of the
            // time domain data
            
			fprintf(ffthandleP,"%e\n0.0\n",hp);  //  The zero is because this is real data, makes
			fprintf(ffthandleX,"%e\n0.0\n",hx);  //  suitable input for my FFT program            
			
			fflush(ffthandleP);
			fflush(ffthandleX);
        }
		
		//  ---- Look for max amplitudes  ---- 
        if (fabs(hx) > cMx)      //  ---------- search for hx max -------- 
			cMx = fabs(hx);      //  Stores the value of the max amp
		
        
        if (fabs(hp) > pMx)      //  ---------- search for h+ max -------- 
			pMx = fabs(hp);      //  Stores the value of the max amp
		
		
		//  ---- increment time  ---- 
        Time += dT;
        
		
	} //  <--- for jj LOOP over theta <--- 
	
	
	
	//  ------ OUTPUT Final Summary TO stdout ------ 
	
	printf("\n\n|hx| max = %e\n",cMx);
	printf("|h+| max = %e\n",pMx);
	
	fclose(fhandle);
	
	printf("\n\nALL DONE!  Whoo hoo!\n");
	
	return 0;
}





// ***************************************************************************
//                                                                            
//  FUNCTION:  GetData                                                        
//  TYPE    :  int                                                            
//  RETURNS :  0 (void)                                                       
//  MODIFIED:  15 February 2002                                               
//                                                                            
//  'GetData' obtains the orbital parameters, and the units used.             
//                                                                            
//  FILE FORMAT FOR "binary.dat" ---------------------------                  
//    row 1 == basename to build filenames around
//    row 2 == semi major axis in AU
//    row 3 == eccentricity
//    row 4 == inclination in DEGREES
//    row 5 == mass1 in SOLAR MASSES
//    row 6 == mass2 in SOLAR MASSES
//    row 7 == distance to source in PARSECS
//    row 8 == anomaly of node in DEGREES
//    row 9 == anomaly of periapse in DEGREES
//    row 10 == orientation of nodes in DEGREES
//    row 11 == Initial Anomaly in DEGREES
//    row 12 == length of simulation in YRS
//    row 13 == number of time points
//    row 14 == 0 = make orbit only, 1 = Make FFT time series file
//
//  PARAMETERS                                                                
//  ------------------------------------------------------------------------- 
//  basename : (string) user input basename to build filenames
//  a        : (double) user input Semi-major axis
//  e        : (double) user input Orbital eccentricity
//  inc      : (double) user input Orbital inclination
//  m1       : (double) user input mass 1 of binary
//  m2       : (double) user input mass 2 of binary
//  R        : (double) user input distance to binary
//  thN      : (double) user input anomaly at line of nodes
//  thP      : (double) user input anomaly at periapse
//  phi      : (double) user input orientation of line of nodes on sky
//  th1      : (double) user input inital anomaly
//  tsim     : (double) user input simulation length in YRs
//  Nt       : (long int) user requested number of time points
//
// ***************************************************************************

int GetData(char baseName[], double *a, double *ec, double *inc, double *m1, double *m2, double *R, double *thN,
            double *thP, double *phi, double *th1, double *tsim, long int *Nt, long int *fftFLAG)
{
	//  =============== VARIABLES =============== 
	
	FILE *phandle;
	char tmpRead[MAX_LENGTH + 1];
	
	//  =============== START ROUTINE HERE =============== 
	
    phandle = fopen("binary.dat","r");   //  open data file with binary parameters 
    
	// ------ Parameter Inquisition ------ 
	
	fgets(tmpRead,MAX_LENGTH + 1,phandle);
	sscanf(tmpRead,"basename  = %s",baseName);
	
	fgets(tmpRead,MAX_LENGTH + 1,phandle);
	sscanf(tmpRead,"a (AU)    = %lf",a);
	
	fgets(tmpRead,MAX_LENGTH + 1,phandle);
	sscanf(tmpRead,"ecc       = %lf",ec);
	
	fgets(tmpRead,MAX_LENGTH + 1,phandle);
	sscanf(tmpRead,"inc (deg) = %lf",inc);
	
	fgets(tmpRead,MAX_LENGTH + 1,phandle);
	sscanf(tmpRead,"m1 (msun) = %lf",m1);
	
	fgets(tmpRead,MAX_LENGTH + 1,phandle);
	sscanf(tmpRead,"m2 (msun) = %lf",m2);
	
	fgets(tmpRead,MAX_LENGTH + 1,phandle);
	sscanf(tmpRead,"R (pc)    = %lf",R);
	
	fgets(tmpRead,MAX_LENGTH + 1,phandle);
	sscanf(tmpRead,"thN (deg) = %lf",thN);
	
	fgets(tmpRead,MAX_LENGTH + 1,phandle);
	sscanf(tmpRead,"thP (deg) = %lf",thP);
	
	fgets(tmpRead,MAX_LENGTH + 1,phandle);
	sscanf(tmpRead,"phi (deg) = %lf",phi);
	
	fgets(tmpRead,MAX_LENGTH + 1,phandle);
	sscanf(tmpRead,"th1 (deg) = %lf",th1);
	
	fgets(tmpRead,MAX_LENGTH + 1,phandle);
	sscanf(tmpRead,"tsim (yr) = %lf",tsim);
	
	fgets(tmpRead,MAX_LENGTH + 1,phandle);
	sscanf(tmpRead,"Npts      = %ld",Nt);
	
	fgets(tmpRead,MAX_LENGTH + 1,phandle);
	sscanf(tmpRead,"fftFlag   = %ld",fftFLAG);
	
	fclose(phandle);
	
	
	//  ------ Printout Inputs ------ 
	
	printf("Binary SemiMajor (AU)  =%lf\n", *a);
	printf("Binary eccentricity    =%lf\n", *ec);
	printf("Binary inclination (d) = %lf\n", *inc);
	printf("Mass 1 (solar masses)  = %lf\n", *m1);
	printf("Mass 2 (solar masses)  = %lf\n", *m2);
	printf("Distance (parsec)      = %lf\n", *R);
	printf("Anomaly at nodes (d)   = %lf\n", *thN);
	printf("Anomaly at perictr (d) = %lf\n", *thP);
	printf("Orientation@nodes (d)  = %lf\n", *phi);
	printf("Initial anomaly (d)    = %lf\n", *th1);
	printf("Total time (yr)        = %lf\n", *tsim);
	printf("Number of time pts     = %ld\n", *Nt);
	printf("fftFLAG (1 = FFTfile)  = %ld\n\n", *fftFLAG);
	
	fflush(stdout);
	
	//  ------ Scale values to normal SI units ------ 
	
	*a *= AU;
	*inc *= degrees;
	*m1 *= MSun;
	*m2 *= MSun;
	*R *= pc;
	*thN *= degrees;
	*thP *= degrees;
	*phi *= degrees;
	*th1 *= degrees;
	*tsim *= YR;
	
	return 0;
	
} //  FUNCTION: GetData 



// ***************************************************************************
//                                                                            
//  FUNCTION:  OrbTime                                                        
//  TYPE    :  int                                                            
//  RETURNS :  t(theta)                                                       
//  MODIFIED:  9 October 2002                                                 
//                                                                            
//  'OrbTime' takes current theta value and computes the time in seconds      
//     from periapse.                                                         
//                                                                            
//  PARAMETERS                                                                
//  ------------------------------------------------------------------------- 
//                                                                            
// ***************************************************************************

double OrbTime(double theta, double Ei, double P)
{
	//  =============== VARIABLES =============== 
	
	double root, psi, trm2, timeVal;
	
	//  =============== START ROUTINE HERE =============== 
	
	root = sqrt((1.0 - Ei)/(1.0 + Ei));
	
	psi = 2.0*atan(root*tan(theta/2.0));  //  Eccentric Anomaly 
	
	trm2 = Ei*sqrt(1.0 - Ei*Ei)*sin(theta)/(1.0 + Ei*cos(theta));
	
	timeVal = (P/(2.0*PI))*(psi - trm2);
	
	return timeVal;
	
}



// ***************************************************************************
//                                                                            
//  FUNCTION:  GetEccentricAnomaly                                            
//  TYPE    :  double                                                         
//  RETURNS :  psi(t)                                                         
//  MODIFIED:  30 July 2003                                                   
//                                                                            
//  'GetEccentricAnomaly' is a bisection routine which solves the Kepler      
//     Equation for Psi given any time t.                                     
//                                                                            
//  NB: THIS ROUTINE ONLY WORKS FOR ECCENTRICITY <= 1.0                       
//                                                                            
//  PARAMETERS                                                                
//  ------------------------------------------------------------------------- 
//                                                                            
// ***************************************************************************

double GetEccentricAnomaly(double Time, double PsiLo, double P, double ecc)
{
	//  =============== VARIABLES =============== 
	
	int tFlag, N;
	
	double dPsi, PsiUp, pL, pU, p1, pVal;   //  Psi variables 
	
	double tEstimate, u;   //  Time variables 
	
	//  =============== START ROUTINE HERE =============== 
	
	
	//  ----- Search for Psi Upper Bracket for Bisection ----- 
	
	tFlag = 0;     //  Set while flag to make sure we do loop 1 
	
	dPsi = 1.0;    //  delta Psi, to search for upper bracket for bisection 
	PsiUp = PsiLo; //  Start search at lower bracket 
	
	while (tFlag != 1)
	{
		PsiUp += dPsi;
		tEstimate = (P/(2.0*PI))*(PsiUp - ecc*sin(PsiUp));
		
		if (tEstimate > Time)       //  If current bracket is greater than time we seek 
			tFlag = 1;               //  kick us out of loop and head for bisection      
		else                        //  else                                            
			PsiLo = PsiUp;           //  Hi bracket still too low so make it low bracket 
	}
	
	
	//  ----------- MAIN BISECTION CODE --------------- 
	
	pL = PsiLo;
	pU = PsiUp;
	
	N = 0;          //  initialize bisection counter   
	
	tFlag = 0;      //  set flag to make sure we bisect at least once 
	
	while (tFlag == 0)
	{
		p1 = pL + (pU - pL)/2.0;                 //  start at middle of bracket interval 
		
		u = Time - (P/(2.0*PI))*(p1 - ecc*sin(p1));  //  current time - time @ psi value 
		
		if (u > 0.0) pL = p1;                      //  if under t, move lower bracket up 
		
		if (u < 0.0) pU = p1;                     //  if over t, move upper bracket down 
		
		if ( fabs(u) <= bEPS ) tFlag = 1;     //  if reached tolerance, kick out of bisect 
		
		if (N > 100)                      //  Terminate for a bisection error          
		{
			printf("N = %d    Time = %g      u = %g\n",N,Time,fabs(u));
            ErrorExit("GetEccentricAnomaly()","Bisections exceed 10000");
		}
		
		N++;                                //  Increment bisection counter   
		
	} //  <--- end WHILE <--- 
	
	pVal = p1;                           //  set return value to bisected value 
	
	return pVal;
	
}



// ***************************************************************************
//                                                                            
//  FUNCTION:  GetTrueAnomaly                                                 
//  TYPE    :  double                                                         
//  RETURNS :  theta(t)                                                       
//  MODIFIED:  30 July 2003                                                   
//                                                                            
//  'GetTrueAnomaly' is an inversion routine to get the true anomaly theta    
//     from the eccentric anomaly Psi.                                        
//                                                                            
//  PARAMETERS                                                                
//  ------------------------------------------------------------------------- 
//                                                                            
// ***************************************************************************

double GetTrueAnomaly(double Psi, double ecc)
{
	//  =============== VARIABLES =============== 
	
	long int nWIND;
	
	double tho, thVal;
	
	//  =============== START ROUTINE HERE =============== 
	
	
	nWIND = (int)floor(Psi/PI);   //  Number of orbits in current anomaly specification 
	
	//  make Theta estimate from Psi 
	
	tho = 2.0*atan(sqrt((1.0 + ecc)/(1.0 - ecc))*tan(Psi/2.0));
	
	
	//  ---- Principle Value           ---------------------------- 
	//  Recover correct value of theta given the input value of Psi 
	
	thVal = tho; 
	
	if (tho > 0.0) thVal += nWIND*PI;
	
	if (tho < 0.0) thVal += (nWIND + 1)*PI;
	
	
	return thVal;
	
}





// ***************************************************************************
//
//  FUNCTION:  ErrorExit
//
//  General routine for error report and exit of program. Takes name of routine,
//  error message, and reports to stdout with the time of the exception, then
//  exits.
//
// ***************************************************************************

void ErrorExit(char routine[], char errorMsg[])
{
    /* =============== VARIABLES =============== */
    
    // --- date & time processing ---
    time_t errorTime;
    struct tm *errorQuery;
    char timeEnded[128];
    
    /* =============== START ROUTINE HERE =============== */

    // --- get the end time of the run ---
    time(&errorTime);                                 // get time from the system
    errorQuery = localtime(&errorTime);               // process time
    strftime(timeEnded,128,"%x - %I:%M%p",errorQuery);  // human readable string
    

    printf("\n\n");
    printf("============== ERROR & EXIT ===============\n\n");
    printf("ROUTINE   : %s\n",routine);
    printf("PROBLEM   : %s\n\n",errorMsg);
    printf("EXIT TIME : %s\n\n",timeEnded);
    printf("===========================================\n\n");

    fflush(stdin);
    exit(0);
    
    return;
    
} // FUNCTION: Error Exit
