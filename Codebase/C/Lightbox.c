#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <windows.h>
#include <math.h>


double rand_double(double min, double max);
float get_IR(float theta,float nhat, float uhat);

int main() {

//Timer Start
    LARGE_INTEGER frequency;        
    LARGE_INTEGER t1, t2;           
    double elapsedTime;
    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&t1);

//Output file declaration
FILE *fptr;
fptr = fopen("Path.dat", "w");

/*__________________________________________________________________________________
    INDOOR LIGHTBOX CODE by Mostafa Sedky replicated from Tarek Zohdi
_____________________________________________________________________________________*/

    int tar,i,j,NTargetsX1,NTargetsX2,NTargetsX3,count,NTargets,WFlag=0,TFLAG=0;
    const int RayCountArraySize=13500, TargetArraySize=1000;
    int NVoxelX1,NVoxelX2,NVoxelX3,IRNUM,override;
    int RANDSEED,reflections,I1,I2,I3,IKEY,printscreen;
    float X1,X2,X3,DISTOL=0.00000001;
    float TARDELTAX1,TARDELTAX2,TARDELTAX3,SUM,radiusFactor;
    float TargetX1[TargetArraySize],TargetX2[TargetArraySize],TargetX3[TargetArraySize];
    float TargetRADIUSX1[TargetArraySize],TargetRADIUSX2[TargetArraySize],TargetRADIUSX3[TargetArraySize],CONDITION,TargetABSORB[TargetArraySize];
    float PX1[TargetArraySize],PX2[TargetArraySize],PX3[TargetArraySize];
    float DELTAX1,DELTAX2,DELTAX3;
    float RANDNUMBER,TIME,TIMELIMIT;
    float RLIGHTX1TPDT[RayCountArraySize],RLIGHTX2TPDT[RayCountArraySize],RLIGHTX3TPDT[RayCountArraySize],WAVELENGTH[RayCountArraySize];
    float RLIGHTX1T[RayCountArraySize],RLIGHTX2T[RayCountArraySize],RLIGHTX3T[RayCountArraySize];
    float VLIGHTX1[RayCountArraySize],VLIGHTX2[RayCountArraySize],VLIGHTX3[RayCountArraySize];
    float VLIGHT0,NORMX1,NORMX2,NORMX3,NORMXX1,NORMXX2,NORMXX3,DELTAT;
    float POPOUT,A1,A2,A3,A4,A5,A6,A7,A8,A9;
    float A10,A11,A12,A13,A14,A15,A16,A17,A18;
    int RCOUNT,TargetFlag[RayCountArraySize];
    float RLIGHTS;
    float WALLX1MINUSCUTOFF,WALLX1PLUSCUTOFF;
    float WALLX2MINUSCUTOFF,WALLX2PLUSCUTOFF;
    float WALLX3MINUSCUTOFF,WALLX3PLUSCUTOFF;
    float GRADF[3],XDIST[3],TERM1,TERM2;
    float GNORM,VNORM,VN,THETAI;
    float NTIWVARIABLE[RayCountArraySize],MAGRATW;
    float NTIVARIABLE[RayCountArraySize],MAGRAT;
    float RPER,RPAR,TIMEMULTIPLIER;
    float BIGRPAR,BIGRPER,BIGR;
    float NX1,NX2,NX3,RADIUS;
    float POWERX1[RayCountArraySize],POWERX2[RayCountArraySize],POWERX3[RayCountArraySize];
    float IRRAD[RayCountArraySize],IRRADO[RayCountArraySize],REF,RAYTOL;
    float TOTALFACEPOWERX1MINUS,TOTALFACEPOWERX1PLUS;
    float TOTALFACEPOWERX2MINUS,TOTALFACEPOWERX2PLUS;
    float TOTALFACEPOWERX3MINUS,TOTALFACEPOWERX3PLUS;
    float TOTALPOWER;
    float IRRADOX1MINUS,IRRADOX1PLUS;
    float IRRADOX2MINUS,IRRADOX2PLUS;
    float IRRADOX3MINUS,IRRADOX3PLUS;
    int RAYCOUNT,RACKX1X2,RACKX2X3,RACKX1X3,RACKX1,RACKX2,RACKX3;
    float FORPLUSX3SOURCETUBEX1,FORPLUSX3SOURCETUBEX2;
    float FORMINUSX3SOURCETUBEX1,FORMINUSX3SOURCETUBEX2;
    float FORPLUSX2SOURCETUBEX1,FORPLUSX2SOURCETUBEX3;
    float FORMINUSX2SOURCETUBEX1,FORMINUSX2SOURCETUBEX3;
    float FORPLUSX1SOURCETUBEX2,FORPLUSX1SOURCETUBEX3;
    float FORMINUSX1SOURCETUBEX2,FORMINUSX1SOURCETUBEX3;
    int X1MINUSON,X1PLUSON,X2MINUSON,X2PLUSON,X3MINUSON,X3PLUSON;
    float RAYDENSITY,STEPFACTOR;
    float NTIWX1MINUS,NTIWX1PLUS,NTIWX2MINUS,NTIWX2PLUS;
    float NTIWX3MINUS,NTIWX3PLUS;
    float NTIX1MINUS,NTIX1PLUS,NTIX2MINUS,NTIX2PLUS;
    float NTIX3MINUS,NTIX3PLUS;
    float WAVELENGTHX1MINUS,WAVELENGTHX1PLUS;
    float WAVELENGTHX2MINUS,WAVELENGTHX2PLUS;
    float WAVELENGTHX3MINUS,WAVELENGTHX3PLUS;
    int NRAYSPLUSX1,NRAYSMINUSX1;
    int NRAYSPLUSX2,NRAYSMINUSX2;
    int NRAYSPLUSX3,NRAYSMINUSX3,MOVIEFRAMES;
    float DELTAMOVIEFRAME,MOVIETIME;
    int TotRays;
    float POWERMAGNITUDE,MAGFACEPOWERX1MINUS,MAGFACEPOWERX1PLUS,MAGFACEPOWERX2MINUS,
    MAGFACEPOWERX2PLUS,MAGFACEPOWERX3MINUS,MAGFACEPOWERX3PLUS,SCORE[10],TOTALCOST;
    float IR;
    float RAYFRAC;
    float TIME2;
    

    //Initial parameters (if no optimization needed). Comment the lines below if you'll optimize these.

    //BEAM SPREAD (18 parameters)
       A1=0.5;
       A2=0.5;
       A3=0.5;
       A4=0.5;
       A5=0.5;
       A6=0.5;
       A7=0.5;
       A8=0.5;
       A9=0.5;
       A10=0.5;
       A11=0.5;
       A12=0.5;
       A13=0.5;
       A14=0.5;
       A15=0.5;
       A16=0.5;
       A17=0.5;
       A18=0.5;

   //SOURCETUBE (12 parameters)
       FORPLUSX3SOURCETUBEX1=3;
       FORPLUSX3SOURCETUBEX2=0.5;
       FORMINUSX3SOURCETUBEX1=3;
       FORMINUSX3SOURCETUBEX2=0.5;
       FORPLUSX2SOURCETUBEX1=3;
       FORPLUSX2SOURCETUBEX3=0.5;
       FORMINUSX2SOURCETUBEX1=3;
       FORMINUSX2SOURCETUBEX3=0.5;
       FORPLUSX1SOURCETUBEX2=0.5;
       FORPLUSX1SOURCETUBEX3=0.5;
       FORMINUSX1SOURCETUBEX2=0.5;
       FORMINUSX1SOURCETUBEX3=0.5;

    //POWER DEPENDENCY (6 parameters)

      TOTALFACEPOWERX1MINUS=10000000.0;
      TOTALFACEPOWERX1PLUS=10000000.0;
      TOTALFACEPOWERX2MINUS=10000000.0;
      TOTALFACEPOWERX2PLUS=10000000.0;
      TOTALFACEPOWERX3MINUS=10000000.0;
      TOTALFACEPOWERX3PLUS=10000000.0;



    //Key parameters
    printscreen=1;
    MOVIEFRAMES=100;
    MOVIETIME=0.0;
    STEPFACTOR=2.0;
    RAYDENSITY=2000;
    TIMEMULTIPLIER=5.0;
    VLIGHT0=299792458.0;
    reflections=1;
    radiusFactor=5;
    RACKX1X3=1;
    RACKX2X3=0;
    RACKX1X2=1;
    RACKX1=0;
    RACKX2=0;
    RACKX3=0;



    //Number of rays in each face

    NRAYSPLUSX1=RAYDENSITY*FORPLUSX1SOURCETUBEX2*FORPLUSX1SOURCETUBEX3;
    printf("Number of rays +X1 = %d \n",NRAYSPLUSX1);
    NRAYSMINUSX1=RAYDENSITY*FORMINUSX1SOURCETUBEX2*FORMINUSX1SOURCETUBEX3;
    printf("Number of rays -X1 = %d \n",NRAYSMINUSX1);
    NRAYSPLUSX2=RAYDENSITY*FORPLUSX2SOURCETUBEX1*FORPLUSX2SOURCETUBEX3;
    printf("Number of rays +X2 = %d \n",NRAYSPLUSX2);
    NRAYSMINUSX2=RAYDENSITY*FORMINUSX2SOURCETUBEX1*FORMINUSX2SOURCETUBEX3;
     printf("Number of rays -X2 = %d \n",NRAYSMINUSX2);
    NRAYSPLUSX3=RAYDENSITY*FORPLUSX3SOURCETUBEX1*FORPLUSX3SOURCETUBEX2;
    printf("Number of rays +X3 = %d \n",NRAYSPLUSX3);
    NRAYSMINUSX3=RAYDENSITY*FORMINUSX3SOURCETUBEX1*FORMINUSX3SOURCETUBEX2;
    printf("Number of rays -X3 = %d \n",NRAYSMINUSX3);
    
    TotRays=NRAYSPLUSX1+NRAYSMINUSX1+NRAYSPLUSX2+NRAYSMINUSX2+NRAYSPLUSX3+NRAYSMINUSX3;
    printf("Total Rays: %d \n",TotRays);
    

    //Which faces of the lightbox are on
    X1MINUSON=1;
    X1PLUSON=1;
    X2MINUSON=1;
    X2PLUSON=1;
    X3MINUSON=1;
    X3PLUSON=1;

    TOTALPOWER=TOTALFACEPOWERX1MINUS*X1MINUSON+TOTALFACEPOWERX1PLUS*X1PLUSON+
               TOTALFACEPOWERX2MINUS*X2MINUSON+TOTALFACEPOWERX2PLUS*X2PLUSON+
               TOTALFACEPOWERX3MINUS*X3MINUSON+TOTALFACEPOWERX3PLUS*X3PLUSON;

    POWERMAGNITUDE=10000000.0;
    MAGFACEPOWERX1MINUS=1.0;
    MAGFACEPOWERX1PLUS=1.0;
    MAGFACEPOWERX2MINUS=1.0;
    MAGFACEPOWERX2PLUS=1.0;
    MAGFACEPOWERX3MINUS=1.0;
    MAGFACEPOWERX3PLUS=1.0;


    //Box voxel dimensions and number of targets
    DELTAX1=0.05;
    DELTAX2=0.05;
    DELTAX3=0.05;
    NVoxelX1=20;
    NVoxelX2=20;
    NVoxelX3=20;

    NTargetsX1=9;
    NTargetsX2=3;
    NTargetsX3=3;
    TARDELTAX1=5.0*DELTAX1;
    TARDELTAX2=5.0*DELTAX2;
    TARDELTAX3=5.0*DELTAX3;

    count=0;

printf("Generating targets... \n");

//X2X3 Targets at X1=0

    if (RACKX2X3 == 1){

     I1=0;
     X1=TARDELTAX1*I1;

     for (I2=-NTargetsX2; I2 <= NTargetsX2; I2++){

        X2=TARDELTAX2*I2;

        for (I3=-NTargetsX3;I3 <= NTargetsX3;I3++){
            if(I3||I2 != 0){
            
            X3=TARDELTAX3*I3;
            TargetX1[count]=X1;
            TargetX2[count]=X2;
            TargetX3[count]=X3;
            TargetRADIUSX1[count]=radiusFactor*DELTAX1;
            TargetRADIUSX2[count]=radiusFactor*DELTAX2;
            TargetRADIUSX3[count]=radiusFactor*DELTAX3;
            PX1[count]=6.0;
            PX2[count]=6.0;
            PX3[count]=6.0;
            TargetABSORB[count]=0.0;
            count++;
            
            }
        }
     }
    }


//X1X2 Targets at X3=0

if (RACKX1X2 == 1){

     I3=0;
     X3=TARDELTAX3*I3;

     for (I2= -NTargetsX2; I2 <= NTargetsX2; I2++){

        X2=TARDELTAX2*I2;
        
        for (I1=-NTargetsX1;I1<=NTargetsX1;I1++){
            if(I1||I2 != 0){
            
            X1=TARDELTAX1*I1;
            TargetX1[count]=X1;
            TargetX2[count]=X2;
            TargetX3[count]=X3;
            TargetRADIUSX1[count]=radiusFactor*DELTAX1;
            TargetRADIUSX2[count]=radiusFactor*DELTAX2;
            TargetRADIUSX3[count]=radiusFactor*DELTAX3;
            PX1[count]=6.0;
            PX2[count]=6.0;
            PX3[count]=6.0;
            TargetABSORB[count]=0.0;
            count++;
            
            }
        }
     }
    }

//X1X3 Targets at X2=0

if (RACKX1X3 == 1){

     I2=0;
     X2=TARDELTAX2*I2;

     for (I1=-NTargetsX1; I1 <= NTargetsX1; I1++){

        X1=TARDELTAX1*I1;
        
        for (I3=-NTargetsX3; I3 <= NTargetsX3;I3++){
            if(I1||I3 != 0){
            
            X3=TARDELTAX3*I3;
            TargetX1[count]=X1;
            TargetX2[count]=X2;
            TargetX3[count]=X3;
            TargetRADIUSX1[count]=radiusFactor*DELTAX1;
            TargetRADIUSX2[count]=radiusFactor*DELTAX2;
            TargetRADIUSX3[count]=radiusFactor*DELTAX3;
            PX1[count]=6.0;
            PX2[count]=6.0;
            PX3[count]=6.0;
            TargetABSORB[count]=0.0;
            count++;
            
            }
        }
     }
    }


//X1X3 Targets at X2=0

if (RACKX1X3 == 1){

     I2=0;
     X2=TARDELTAX2*I2;

     for (I1=-NTargetsX1; I1 <= NTargetsX1; I1++){

        X1=TARDELTAX1*I1;
        
        for (I3=-NTargetsX3; I3 <= NTargetsX3;I3++){
            if(I1||I3 != 0){
            
            X3=TARDELTAX3*I3;
            TargetX1[count]=X1;
            TargetX2[count]=X2;
            TargetX3[count]=X3;
            TargetRADIUSX1[count]=radiusFactor*DELTAX1;
            TargetRADIUSX2[count]=radiusFactor*DELTAX2;
            TargetRADIUSX3[count]=radiusFactor*DELTAX3;
            PX1[count]=6.0;
            PX2[count]=6.0;
            PX3[count]=6.0;
            TargetABSORB[count]=0.0;
            count++;
            
            }
        }
     }
    }
    NTargets=count;
printf("%d Targets Generated. \n",NTargets);
printf("Setting wallcutoffs and initial powers... \n");

//Wall cutoffs or "box" limits

    WALLX1MINUSCUTOFF = -3*(NVoxelX1+1)*DELTAX1;
    WALLX1PLUSCUTOFF = 3*(NVoxelX1+1)*DELTAX1;
    WALLX2MINUSCUTOFF = -(NVoxelX2+1)*DELTAX2; 
    WALLX2PLUSCUTOFF = (NVoxelX2+1)*DELTAX2; 
    WALLX3MINUSCUTOFF = -(NVoxelX3+1)*DELTAX3;
    WALLX3PLUSCUTOFF = (NVoxelX3+1)*DELTAX3;

    printf("-X1 cutoff= %lf, +X1 cutoff= %lf, -X2 cutoff = %lf, +X2 cutoff= %lf, -X3 cutoff = %lf, +X3 cutoff= %lf \n",WALLX1MINUSCUTOFF, 
                                                                                                                     WALLX1PLUSCUTOFF,WALLX2MINUSCUTOFF,
                                                                                                                     WALLX2PLUSCUTOFF, WALLX3MINUSCUTOFF,WALLX3PLUSCUTOFF);

//Wavelengths and indices of refraction for each wall

    NTIWX1MINUS=2.0;
    NTIWX1PLUS=2.0;
    NTIWX2MINUS=2.0;
    NTIWX2PLUS=2.0;
    NTIWX3MINUS=2.0;
    NTIWX3PLUS=2.0;
    MAGRATW=1.0;
    NTIX1MINUS=2.0;
    NTIX1PLUS=2.0;
    NTIX2MINUS=2.0;
    NTIX2PLUS=2.0;
    NTIX3MINUS=2.0;
    NTIX3PLUS=2.0;
    MAGRAT=1.0;
    
    WAVELENGTHX1MINUS=10.0;
    WAVELENGTHX1PLUS=4.0;
    WAVELENGTHX2MINUS=1.5;
    WAVELENGTHX2PLUS=3.0;
    WAVELENGTHX3MINUS=9.0;
    WAVELENGTHX3PLUS=2.0;



//Initial power in each face   

    TOTALFACEPOWERX1MINUS=MAGFACEPOWERX1MINUS*POWERMAGNITUDE;
    TOTALFACEPOWERX1PLUS=MAGFACEPOWERX1PLUS*POWERMAGNITUDE;
    TOTALFACEPOWERX2MINUS=MAGFACEPOWERX2MINUS*POWERMAGNITUDE;
    TOTALFACEPOWERX2PLUS=MAGFACEPOWERX2PLUS*POWERMAGNITUDE;
    TOTALFACEPOWERX3MINUS=MAGFACEPOWERX3MINUS*POWERMAGNITUDE;
    TOTALFACEPOWERX3PLUS=MAGFACEPOWERX3PLUS*POWERMAGNITUDE;

//Initial power in each ray

    IRRADOX1MINUS=TOTALFACEPOWERX1MINUS/NRAYSMINUSX1;
    IRRADOX1PLUS=TOTALFACEPOWERX1PLUS/NRAYSPLUSX1;
    IRRADOX2MINUS=TOTALFACEPOWERX2MINUS/NRAYSMINUSX2;
    IRRADOX2PLUS=TOTALFACEPOWERX2PLUS/NRAYSPLUSX2;
    IRRADOX3MINUS=TOTALFACEPOWERX3MINUS/NRAYSMINUSX3;
    IRRADOX3PLUS=TOTALFACEPOWERX3PLUS/NRAYSPLUSX3;


    RAYTOL=0.000001;
    REF=0.0;
    RADIUS=0.0001;

printf("Done. \n");
//Defining time step and time limit for simulation

    DELTAT=STEPFACTOR*(DELTAX1+DELTAX2+DELTAX3)/(3*VLIGHT0);

    TIMELIMIT=TIMEMULTIPLIER*(2*NVoxelX1*DELTAX1+2*NVoxelX2*DELTAX2+2*NVoxelX3*DELTAX3)/(3*VLIGHT0);

    DELTAMOVIEFRAME=TIMELIMIT/MOVIEFRAMES;

//Generating rays
    printf("Generating Rays... \n");
    RAYCOUNT=0;

//Minus x1 wall
if (X1MINUSON==1){
for (int i=1;i<=NRAYSMINUSX1;i++){


    IRRAD[RAYCOUNT]=IRRADOX1MINUS;
    IRRADO[RAYCOUNT]=IRRADOX1MINUS; //initial ray power
    RLIGHTX1T[RAYCOUNT]=WALLX1PLUSCUTOFF-1.5*DELTAX1; //ray x1 
    RLIGHTX2T[RAYCOUNT]=FORMINUSX1SOURCETUBEX2*(rand_double(0,1)-0.5)/0.5; //ray x2 position
    RLIGHTX3T[RAYCOUNT]=FORMINUSX1SOURCETUBEX3*(rand_double(0,1)-0.5)/0.5; //ray x3 position

    // printf("-X1 Wall Ray number %d positions are %lf, %lf, %lf \n",RAYCOUNT,RLIGHTX1T[RAYCOUNT],RLIGHTX2T[RAYCOUNT],RLIGHTX3T[RAYCOUNT]);

    NORMXX1=A1*((rand_double(0,1)-0.5)/0.5)-1;
    NORMXX2=A2*((rand_double(0,1)-0.5)/0.5);
    NORMXX3=A3*((rand_double(0,1)-0.5)/0.5);
    NORMX1=NORMXX1/pow((pow(NORMXX1,2)+pow(NORMXX2,2)+pow(NORMXX3,2)),0.5);
    NORMX2=NORMXX2/pow((pow(NORMXX1,2)+pow(NORMXX2,2)+pow(NORMXX3,2)),0.5);
    NORMX3=NORMXX3/pow((pow(NORMXX1,2)+pow(NORMXX2,2)+pow(NORMXX3,2)),0.5);
    VLIGHTX1[RAYCOUNT]=NORMX1*VLIGHT0;
    VLIGHTX2[RAYCOUNT]=NORMX2*VLIGHT0;
    VLIGHTX3[RAYCOUNT]=NORMX3*VLIGHT0;
    POWERX1[RAYCOUNT]=IRRAD[RAYCOUNT]*NORMX1;
    POWERX2[RAYCOUNT]=IRRAD[RAYCOUNT]*NORMX2;
    POWERX3[RAYCOUNT]=IRRAD[RAYCOUNT]*NORMX3;

    NTIVARIABLE[RAYCOUNT]=NTIX1MINUS;
    NTIWVARIABLE[RAYCOUNT]=NTIWX1MINUS;
    WAVELENGTH[RAYCOUNT]=WAVELENGTHX1MINUS;
    
    RAYCOUNT++;
}
}

//Plus x1 wall

if (X1PLUSON==1){
for (int i=1;i<=NRAYSPLUSX1;i++){

    IRRAD[RAYCOUNT]=IRRADOX1PLUS;
    IRRADO[RAYCOUNT]=IRRADOX1PLUS; //initial ray power
    RLIGHTX1T[RAYCOUNT]=WALLX1MINUSCUTOFF+1.5*DELTAX1; //ray x1 position
    RLIGHTX2T[RAYCOUNT]=FORPLUSX1SOURCETUBEX2*(rand_double(0,1)-0.5)/0.5; //ray x2 position
    RLIGHTX3T[RAYCOUNT]=FORPLUSX1SOURCETUBEX3*(rand_double(0,1)-0.5)/0.5; //ray x3 position

    // printf("+X1 Wall Ray number %d positions are %lf, %lf, %lf \n",RAYCOUNT,RLIGHTX1T[RAYCOUNT],RLIGHTX2T[RAYCOUNT],RLIGHTX3T[RAYCOUNT]);

    NORMXX1=A4*((rand_double(0,1)-0.5)/0.5)+1;
    NORMXX2=A5*((rand_double(0,1)-0.5)/0.5);
    NORMXX3=A6*((rand_double(0,1)-0.5)/0.5);
    NORMX1=NORMXX1/pow((pow(NORMXX1,2)+pow(NORMXX2,2)+pow(NORMXX3,2)),0.5);
    NORMX2=NORMXX2/pow((pow(NORMXX1,2)+pow(NORMXX2,2)+pow(NORMXX3,2)),0.5);
    NORMX3=NORMXX3/pow((pow(NORMXX1,2)+pow(NORMXX2,2)+pow(NORMXX3,2)),0.5);
    VLIGHTX1[RAYCOUNT]=NORMX1*VLIGHT0;
    VLIGHTX2[RAYCOUNT]=NORMX2*VLIGHT0;
    VLIGHTX3[RAYCOUNT]=NORMX3*VLIGHT0;
    POWERX1[RAYCOUNT]=IRRAD[RAYCOUNT]*NORMX1;
    POWERX2[RAYCOUNT]=IRRAD[RAYCOUNT]*NORMX2;
    POWERX3[RAYCOUNT]=IRRAD[RAYCOUNT]*NORMX3;

    NTIVARIABLE[RAYCOUNT]=NTIX1PLUS;
    NTIWVARIABLE[RAYCOUNT]=NTIWX1PLUS;
    WAVELENGTH[RAYCOUNT]=WAVELENGTHX1PLUS;
    
    RAYCOUNT++;
}
}

//Minus x2 wall

if (X2MINUSON==1){
for (int i=1;i<=NRAYSMINUSX2;i++){

    IRRAD[RAYCOUNT]=IRRADOX2MINUS;
    IRRADO[RAYCOUNT]=IRRADOX2MINUS; //initial ray power
    RLIGHTX1T[RAYCOUNT]=FORMINUSX2SOURCETUBEX1*(rand_double(0,1)-0.5)/0.5; //ray x1 position
    RLIGHTX2T[RAYCOUNT]=WALLX2PLUSCUTOFF-1.5*DELTAX2; //ray x2 position
    RLIGHTX3T[RAYCOUNT]= FORMINUSX2SOURCETUBEX3*(rand_double(0,1)-0.5)/0.5; //ray x3 position

    // printf("-X2 Wall Ray number %d positions are %lf, %lf, %lf \n",RAYCOUNT,RLIGHTX1T[RAYCOUNT],RLIGHTX2T[RAYCOUNT],RLIGHTX3T[RAYCOUNT]);

    NORMXX1=A7*((rand_double(0,1)-0.5)/0.5);
    NORMXX2=A8*((rand_double(0,1)-0.5)/0.5)-1;
    NORMXX3=A9*((rand_double(0,1)-0.5)/0.5);
    NORMX1=NORMXX1/pow((pow(NORMXX1,2)+pow(NORMXX2,2)+pow(NORMXX3,2)),0.5);
    NORMX2=NORMXX2/pow((pow(NORMXX1,2)+pow(NORMXX2,2)+pow(NORMXX3,2)),0.5);
    NORMX3=NORMXX3/pow((pow(NORMXX1,2)+pow(NORMXX2,2)+pow(NORMXX3,2)),0.5);
    VLIGHTX1[RAYCOUNT]=NORMX1*VLIGHT0;
    VLIGHTX2[RAYCOUNT]=NORMX2*VLIGHT0;
    VLIGHTX3[RAYCOUNT]=NORMX3*VLIGHT0;
    POWERX1[RAYCOUNT]=IRRAD[RAYCOUNT]*NORMX1;
    POWERX2[RAYCOUNT]=IRRAD[RAYCOUNT]*NORMX2;
    POWERX3[RAYCOUNT]=IRRAD[RAYCOUNT]*NORMX3;

    NTIVARIABLE[RAYCOUNT]=NTIX2MINUS;
    NTIWVARIABLE[RAYCOUNT]=NTIWX2MINUS;
    WAVELENGTH[RAYCOUNT]=WAVELENGTHX2MINUS;
    
    RAYCOUNT++;
}
}

//Plus x2 wall

if (X2PLUSON==1){
for (int i=1;i<=NRAYSPLUSX2;i++){


    IRRAD[RAYCOUNT]=IRRADOX2PLUS;
    IRRADO[RAYCOUNT]=IRRADOX2PLUS; //initial ray power
    RLIGHTX1T[RAYCOUNT]=FORPLUSX2SOURCETUBEX1*(rand_double(0,1)-0.5)/0.5; //ray x1 position
    RLIGHTX2T[RAYCOUNT]=WALLX2MINUSCUTOFF+1.5*DELTAX2; //ray x2 position
    RLIGHTX3T[RAYCOUNT]= FORPLUSX2SOURCETUBEX3*(rand_double(0,1)-0.5)/0.5; //ray x3 position

    // printf("+X2 Wall Ray number %d positions are %lf, %lf, %lf \n",RAYCOUNT,RLIGHTX1T[RAYCOUNT],RLIGHTX2T[RAYCOUNT],RLIGHTX3T[RAYCOUNT]);

    NORMXX1=A10*((rand_double(0,1)-0.5)/0.5);
    NORMXX2=A11*((rand_double(0,1)-0.5)/0.5)+1;
    NORMXX3=A12*((rand_double(0,1)-0.5)/0.5);
    NORMX1=NORMXX1/pow((pow(NORMXX1,2)+pow(NORMXX2,2)+pow(NORMXX3,2)),0.5);
    NORMX2=NORMXX2/pow((pow(NORMXX1,2)+pow(NORMXX2,2)+pow(NORMXX3,2)),0.5);
    NORMX3=NORMXX3/pow((pow(NORMXX1,2)+pow(NORMXX2,2)+pow(NORMXX3,2)),0.5);
    VLIGHTX1[RAYCOUNT]=NORMX1*VLIGHT0;
    VLIGHTX2[RAYCOUNT]=NORMX2*VLIGHT0;
    VLIGHTX3[RAYCOUNT]=NORMX3*VLIGHT0;
    POWERX1[RAYCOUNT]=IRRAD[RAYCOUNT]*NORMX1;
    POWERX2[RAYCOUNT]=IRRAD[RAYCOUNT]*NORMX2;
    POWERX3[RAYCOUNT]=IRRAD[RAYCOUNT]*NORMX3;

    NTIVARIABLE[RAYCOUNT]=NTIX2PLUS;
    NTIWVARIABLE[RAYCOUNT]=NTIWX2PLUS;
    WAVELENGTH[RAYCOUNT]=WAVELENGTHX2PLUS;
    
    RAYCOUNT++;
}
}

//Minus x3 wall

if (X3MINUSON==1){
for (int i=1;i<=NRAYSMINUSX3;i++){

    IRRAD[RAYCOUNT]=IRRADOX3MINUS;
    IRRADO[RAYCOUNT]=IRRADOX3MINUS; //initial ray power
    RLIGHTX1T[RAYCOUNT]=FORMINUSX3SOURCETUBEX1*(rand_double(0,1)-0.5)/0.5; //ray x1 position
    RLIGHTX2T[RAYCOUNT]= FORMINUSX3SOURCETUBEX2*(rand_double(0,1)-0.5)/0.5; //ray x2 position
    RLIGHTX3T[RAYCOUNT]=  WALLX3PLUSCUTOFF-1.5*DELTAX3; //ray x3 position 

    // printf("-X3 Wall Ray number %d positions are %lf, %lf, %lf \n",RAYCOUNT,RLIGHTX1T[RAYCOUNT],RLIGHTX2T[RAYCOUNT],RLIGHTX3T[RAYCOUNT]);


    NORMXX1=A13*((rand_double(0,1)-0.5)/0.5);
    NORMXX2=A14*((rand_double(0,1)-0.5)/0.5);
    NORMXX3=A15*((rand_double(0,1)-0.5)/0.5)-1;
    NORMX1=NORMXX1/pow((pow(NORMXX1,2)+pow(NORMXX2,2)+pow(NORMXX3,2)),0.5);
    NORMX2=NORMXX2/pow((pow(NORMXX1,2)+pow(NORMXX2,2)+pow(NORMXX3,2)),0.5);
    NORMX3=NORMXX3/pow((pow(NORMXX1,2)+pow(NORMXX2,2)+pow(NORMXX3,2)),0.5);
    VLIGHTX1[RAYCOUNT]=NORMX1*VLIGHT0;
    VLIGHTX2[RAYCOUNT]=NORMX2*VLIGHT0;
    VLIGHTX3[RAYCOUNT]=NORMX3*VLIGHT0;
    POWERX1[RAYCOUNT]=IRRAD[RAYCOUNT]*NORMX1;
    POWERX2[RAYCOUNT]=IRRAD[RAYCOUNT]*NORMX2;
    POWERX3[RAYCOUNT]=IRRAD[RAYCOUNT]*NORMX3;

    NTIVARIABLE[RAYCOUNT]=NTIX3MINUS;
    NTIWVARIABLE[RAYCOUNT]=NTIWX3MINUS;
    WAVELENGTH[RAYCOUNT]=WAVELENGTHX3MINUS;
    
    RAYCOUNT++;
}
}

//Plus x3 wall

if (X3PLUSON==1){
for (int i=1;i<=NRAYSPLUSX3;i++){


    IRRAD[RAYCOUNT]=IRRADOX3PLUS;
    IRRADO[RAYCOUNT]=IRRADOX3PLUS; //initial ray power
    RLIGHTX1T[RAYCOUNT]=FORPLUSX3SOURCETUBEX1*(rand_double(0,1)-0.5)/0.5; //ray x1 position
    RLIGHTX2T[RAYCOUNT]= FORPLUSX3SOURCETUBEX2*(rand_double(0,1)-0.5)/0.5; //ray x2 position 
    RLIGHTX3T[RAYCOUNT]= WALLX3MINUSCUTOFF+1.5*DELTAX3; //ray x3 position

    // printf("+X3 Wall Ray number %d positions are %lf, %lf, %lf \n",RAYCOUNT,RLIGHTX1T[RAYCOUNT],RLIGHTX2T[RAYCOUNT],RLIGHTX3T[RAYCOUNT]);


    NORMXX1=A16*((rand_double(0,1)-0.5)/0.5);
    NORMXX2=A17*((rand_double(0,1)-0.5)/0.5);
    NORMXX3=A18*((rand_double(0,1)-0.5)/0.5)+1;
    NORMX1=NORMXX1/pow((pow(NORMXX1,2)+pow(NORMXX2,2)+pow(NORMXX3,2)),0.5);
    NORMX2=NORMXX2/pow((pow(NORMXX1,2)+pow(NORMXX2,2)+pow(NORMXX3,2)),0.5);
    NORMX3=NORMXX3/pow((pow(NORMXX1,2)+pow(NORMXX2,2)+pow(NORMXX3,2)),0.5);
    VLIGHTX1[RAYCOUNT]=NORMX1*VLIGHT0;
    VLIGHTX2[RAYCOUNT]=NORMX2*VLIGHT0;
    VLIGHTX3[RAYCOUNT]=NORMX3*VLIGHT0;
    POWERX1[RAYCOUNT]=IRRAD[RAYCOUNT]*NORMX1;
    POWERX2[RAYCOUNT]=IRRAD[RAYCOUNT]*NORMX2;
    POWERX3[RAYCOUNT]=IRRAD[RAYCOUNT]*NORMX3;

    NTIVARIABLE[RAYCOUNT]=NTIX3PLUS;
    NTIWVARIABLE[RAYCOUNT]=NTIWX3PLUS;
    WAVELENGTH[RAYCOUNT]=WAVELENGTHX3PLUS;
    
    RAYCOUNT++;
}
}
printf("%d Rays Generated. \n",RAYCOUNT);

// for(int k=0;k<RAYCOUNT;k++){

// printf("Ray number %d positions are %lf, %lf, %lf \n",k,RLIGHTX1T[k],RLIGHTX2T[k],RLIGHTX3T[k]);

// }

// Light propogation
    printf("________________________________________________________________ \n");
    printf("Propagating Rays... \n");
    TIME=0;
    TIME2=0;

while(TIME<TIMELIMIT){

printf("TIME = %lf \n", TIME/TIMELIMIT);

    if (reflections==1){
        
        for (i=0; i<RAYCOUNT; i++){

            for (tar=0; tar<NTargets; tar++){
    
            CONDITION=pow((fabs((RLIGHTX1T[i] - TargetX1[tar])/TargetRADIUSX1[tar])),PX1[tar])+
                      pow((fabs((RLIGHTX2T[i] - TargetX2[tar])/TargetRADIUSX2[tar])),PX2[tar]) +
                      pow((fabs((RLIGHTX3T[i] - TargetX3[tar])/TargetRADIUSX3[tar])),PX3[tar]);       
        
        

        if(CONDITION <= 1 || RLIGHTX1T[i]>=WALLX1PLUSCUTOFF || RLIGHTX1T[i] <= WALLX1MINUSCUTOFF || 
                                RLIGHTX2T[i]>=WALLX2PLUSCUTOFF || RLIGHTX2T[i] <= WALLX2MINUSCUTOFF ||   
                                RLIGHTX3T[i]>=WALLX3PLUSCUTOFF || RLIGHTX3T[i] <= WALLX3MINUSCUTOFF){
                            
                            //printf("%lf \n",CONDITION);
                    
                            if(RLIGHTX1T[i]>=WALLX1PLUSCUTOFF){
                            
                            WFlag=1;
                            NX1=-1;
                            NX2=0;
                            NX3=0;
                            VNORM=sqrt(pow(VLIGHTX1[i],2)+pow(VLIGHTX2[i],2)+pow(VLIGHTX3[i],2));
                            VN=VLIGHTX1[i]*NX1+VLIGHTX2[i]*NX2+VLIGHTX3[i]*NX3;
                            break;
                            //printf("Ray number %d hit wall +X1 \n",i);   
                            //printf("Ray number %d positions are %lf, %lf, %lf \n",i,RLIGHTX1T[i],RLIGHTX2T[i],RLIGHTX3T[i]);
     

                            }else if(RLIGHTX1T[i] <= WALLX1MINUSCUTOFF){
                            
                            WFlag=2;
                            NX1=1;
                            NX2=0;
                            NX3=0;
                            VNORM=sqrt(pow(VLIGHTX1[i],2)+pow(VLIGHTX2[i],2)+pow(VLIGHTX3[i],2));
                            VN=VLIGHTX1[i]*NX1+VLIGHTX2[i]*NX2+VLIGHTX3[i]*NX3;
                            break;
                            //printf("Ray number %d hit wall -X1 \n",i); 
                            //printf("Ray number %d positions are %lf, %lf, %lf \n",i,RLIGHTX1T[i],RLIGHTX2T[i],RLIGHTX3T[i]);

                            }else if(RLIGHTX2T[i] >= WALLX2PLUSCUTOFF){
                            
                            WFlag=3;
                            NX1=0;
                            NX2=-1;
                            NX3=0;
                            VNORM=sqrt(pow(VLIGHTX1[i],2)+pow(VLIGHTX2[i],2)+pow(VLIGHTX3[i],2));
                            VN=VLIGHTX1[i]*NX1+VLIGHTX2[i]*NX2+VLIGHTX3[i]*NX3;
                            break;
                            //printf("Ray number %d hit wall +X2 \n",i); 
                            //printf("Ray number %d positions are %lf, %lf, %lf \n",i,RLIGHTX1T[i],RLIGHTX2T[i],RLIGHTX3T[i]);

                            }else if(RLIGHTX2T[i] <= WALLX2MINUSCUTOFF){
                                
                            WFlag=4;
                            NX1=0;
                            NX2=1;
                            NX3=0; 
                            VNORM=sqrt(pow(VLIGHTX1[i],2)+pow(VLIGHTX2[i],2)+pow(VLIGHTX3[i],2));
                            VN=VLIGHTX1[i]*NX1+VLIGHTX2[i]*NX2+VLIGHTX3[i]*NX3;
                            break;
                            //printf("Ray number %d hit wall -X2 \n",i);
                            //printf("Ray number %d positions are %lf, %lf, %lf \n",i,RLIGHTX1T[i],RLIGHTX2T[i],RLIGHTX3T[i]);

                            }else if(RLIGHTX3T[i]>=WALLX3PLUSCUTOFF){
                            
                            WFlag=5;
                            NX1=0;
                            NX2=0;
                            NX3=-1; 
                            VNORM=sqrt(pow(VLIGHTX1[i],2)+pow(VLIGHTX2[i],2)+pow(VLIGHTX3[i],2));
                            VN=VLIGHTX1[i]*NX1+VLIGHTX2[i]*NX2+VLIGHTX3[i]*NX3;
                            break;
                            //printf("Ray number %d hit wall +X3 \n",i);
                            //printf("Ray number %d positions are %lf, %lf, %lf \n",i,RLIGHTX1T[i],RLIGHTX2T[i],RLIGHTX3T[i]);

                            }else if(RLIGHTX3T[i] <= WALLX3MINUSCUTOFF){
                            
                            WFlag=6;
                            NX1=0;
                            NX2=0;
                            NX3=1; 
                            VNORM=sqrt(pow(VLIGHTX1[i],2)+pow(VLIGHTX2[i],2)+pow(VLIGHTX3[i],2));
                            VN=VLIGHTX1[i]*NX1+VLIGHTX2[i]*NX2+VLIGHTX3[i]*NX3;
                            break;
                            //printf("Ray number %d hit wall -X3 \n",i);
                            //printf("Ray number %d positions are %lf, %lf, %lf \n",i,RLIGHTX1T[i],RLIGHTX2T[i],RLIGHTX3T[i]);

                            }else{

                            //Target hit
                            TFLAG=1;
                            TargetFlag[i]=1.0;
                            //printf("Ray number %d hit wall Target number %d \n",i,j);
                            //printf("Ray number %d positions are %lf, %lf, %lf \n",i,RLIGHTX1T[i],RLIGHTX2T[i],RLIGHTX3T[i]);

                            GRADF[1]=0;
                            GRADF[2]=0;
                            GRADF[3]=0;

                            XDIST[1]=RLIGHTX1T[i]-TargetX1[tar];
                            XDIST[2]=RLIGHTX2T[i]-TargetX2[tar];
                            XDIST[3]=RLIGHTX3T[i]-TargetX3[tar];

                                if(fabs(XDIST[1])>DISTOL){

                                TERM1=PX1[tar]*pow(((fabs(XDIST[1]))/TargetRADIUSX1[tar]),(PX1[tar]-1.0));
                                TERM2=XDIST[1]/(TargetRADIUSX1[tar]*fabs(XDIST[1]));
                                GRADF[1]=TERM1*TERM2;
        
                                }
        
                               if(fabs(XDIST[2])>DISTOL){

                                TERM1=PX2[tar]*pow(((fabs(XDIST[1]))/TargetRADIUSX2[tar]),(PX2[tar]-1.0));
                                TERM2=XDIST[1]/(TargetRADIUSX2[tar]*fabs(XDIST[1]));
                                GRADF[2]=TERM1*TERM2;
        
                                }
        
                                if(fabs(XDIST[3])>DISTOL){

                                TERM1=PX3[tar]*pow(((fabs(XDIST[1]))/TargetRADIUSX3[tar]),(PX3[tar]-1.0));
                                TERM2=XDIST[1]/(TargetRADIUSX3[tar]*fabs(XDIST[1]));
                                GRADF[3]=TERM1*TERM2;
        
                                }

                            GNORM=sqrt(pow(GRADF[1],2)+pow(GRADF[2],2)+pow(GRADF[3],2));
                            NX1=GRADF[1]/GNORM;
                            NX2=GRADF[2]/GNORM;
                            NX3=GRADF[3]/GNORM;
                            
                            VN=VLIGHTX1[i]*NX1+VLIGHTX2[i]*NX2+VLIGHTX3[i]*NX3;
                            VNORM=sqrt(pow(VLIGHTX1[i],2)+pow(VLIGHTX2[i],2)+pow(VLIGHTX3[i],2));
                            break;
                            }
                }
            }

            if(TFLAG==1||WFlag>0){

                //Calculating incident theta
                THETAI=acos(fabs(VN)/VNORM);
                

                if (TFLAG==1){
                TargetFlag[i]=1;

                //Updating ray powers after target is hit
                IR=get_IR(THETAI,NTIVARIABLE[i],MAGRAT);
                //printf("Power= %lf \n",IR);

                TargetABSORB[tar]=TargetABSORB[tar]+(IRRAD[i]-IR*IRRAD[i]);
                IRRAD[i]=IR*IRRAD[i];
                
            
                }else{
                
                //Updating ray powers after wall is hit
                IR=get_IR(THETAI,NTIWVARIABLE[i],MAGRATW);
                IRRAD[i]=IR*IRRAD[i];

                }

                //Updating ray powers after target or wall is hit
                POWERX1[i]=IRRAD[i]*VLIGHTX1[i]/VNORM;
                POWERX2[i]=IRRAD[i]*VLIGHTX2[i]/VNORM;
                POWERX3[i]=IRRAD[i]*VLIGHTX3[i]/VNORM;

                //Updating ray velocities after target or wall is hit
                VLIGHTX1[i]=VLIGHTX1[i]-2*VN*NX1;
                VLIGHTX2[i]=VLIGHTX2[i]-2*VN*NX2;
                VLIGHTX3[i]=VLIGHTX3[i]-2*VN*NX3;

                //Updating ray positions after target or wall is hit
                if(IRRAD[i]/IRRADO[i]<=RAYTOL){ //If ray lost its energy
                //printf("Ray deactivated! \n");
                RLIGHTX1T[i]=RLIGHTX1T[i];
                RLIGHTX2T[i]=RLIGHTX2T[i];
                RLIGHTX3T[i]=RLIGHTX3T[i];

                }else{

                RLIGHTX1T[i]=RLIGHTX1T[i]+VLIGHTX1[i]*DELTAT;
                RLIGHTX2T[i]=RLIGHTX2T[i]+VLIGHTX2[i]*DELTAT;
                RLIGHTX3T[i]=RLIGHTX3T[i]+VLIGHTX3[i]*DELTAT;

                }
                
                TFLAG=0;
                WFlag=0;

            }else{
                //Ray in the air not touching walls or targets
                //Only position is updated. Velocity and power are not changed

                RLIGHTX1T[i]=RLIGHTX1T[i]+VLIGHTX1[i]*DELTAT;
                RLIGHTX2T[i]=RLIGHTX2T[i]+VLIGHTX2[i]*DELTAT;
                RLIGHTX3T[i]=RLIGHTX3T[i]+VLIGHTX3[i]*DELTAT;

            }

        }//end of ray count for loop

    }

if (printscreen==1){
       

        fprintf(fptr,"ZONE T=\"Rays\" , SOLUTIONTIME= %lf \n",TIME/TIMELIMIT);
        // printf("TIME = %lf\n",TIME);

        for(i=0;i<RAYCOUNT;i++){

        fprintf(fptr, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",RLIGHTX1T[i],RLIGHTX2T[i],RLIGHTX3T[i]
                                                                                    ,2.0*RADIUS,IRRAD[i]/IRRADO[i],VLIGHTX1[i]
                                                                                    ,VLIGHTX2[i],VLIGHTX3[i],POWERX1[i]/IRRADO[i]
                                                                                    ,POWERX2[i]/IRRADO[i],POWERX3[i]/IRRADO[i],0.0,WAVELENGTH[i]);
        }

        
        fprintf(fptr,"ZONE T=\"TARGETS\", SOLUTIONTIME=%lf \n",TIME/TIMELIMIT);
        
        for(i=0;i<NTargets;i++){
        fprintf(fptr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",TargetX1[i],TargetX2[i],TargetX3[i],2.0*TargetRADIUSX1[i],
                                                                    0.0,0.0,0.0,0.0,0.0,0.0,0.0,TargetABSORB[i]/(TOTALPOWER),0.0);
        }
}


    TIME=TIME+DELTAT;
}

        RCOUNT=0;
        for(i=0;i<RAYCOUNT;i++){

            RCOUNT=RCOUNT+TargetFlag[i];
        }   

        SUM=0;

        for(i=0;i<NTargets;i++){

            SUM=SUM+TargetABSORB[i];
            // printf("%lf \n",TargetABSORB[i]);
        }

RAYFRAC=(float) RCOUNT/RAYCOUNT;
printf("Rays hitting target = %d \n",RCOUNT);
printf("Total rays = %d \n",RAYCOUNT);
printf("Total of power absorbed by plants= %lf \n", SUM);
printf("Total intitial power= %lf \n",TOTALPOWER);
printf("Fraction of rays hitting target = %lf \n",RAYFRAC);
printf("Fraction of power absorbed by plants= %lf \n", SUM/TOTALPOWER);
printf("Power usage effectiveness by pods= %lf \n",TOTALPOWER/SUM);

//Close files
fclose(fptr);

// Stop timer and print the elapsed time in seconds
    QueryPerformanceCounter(&t2);
    elapsedTime = (t2.QuadPart - t1.QuadPart) * 1000 / frequency.QuadPart;
    printf("%f seconds elapsed.\n", elapsedTime/1000);

  return 0;
}


//_________________________________________________________________________________________


double rand_double(double min, double max) {
    return (double)rand() / RAND_MAX * (max - min) + min;
}

float get_IR(float theta,float nhat,float uhat){
    
    float IR;
    float Rper;
    float Rpar;
    float nsq;

    nsq=nhat*nhat;
    
    if((fabs(sin(theta)/nhat))>1){

        Rper=1;
        Rpar=1;
        // IR=0.5*(pow(Rper,2)+pow(Rpar,2));
    }else{

        Rper=(cos(theta)-(1/uhat)*sqrt(pow(nhat,2)-pow(sin(theta),2)))/
             (cos(theta)+(1/uhat)*sqrt(pow(nhat,2)-pow(sin(theta),2)));

        Rpar=(pow(nhat,2)/uhat * cos(theta) - sqrt(pow(nhat,2) - pow(sin(theta),2)))/
             (pow(nhat,2)/uhat * cos(theta) + sqrt(pow(nhat,2) - pow(sin(theta),2))); 

        // IR=0.5*(powf(((nsq*cos(theta)-sqrtf(nsq-powf((sin(theta)),2)))/(nsq*cos(theta)+sqrtf(nsq-powf((sin(theta)),2)))),2)+
        // powf(((cos(theta)-sqrtf(nsq-powf((sin(theta)),2)))/(cos(theta)+sqrtf(nsq-powf((sin(theta)),2)))),2));
    }

    IR=0.5*(pow(Rper,2)+pow(Rpar,2));
    return(IR);
}
