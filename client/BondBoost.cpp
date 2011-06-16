//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "BondBoost.h"
#include <math.h>

//const string Hyperdynamics::NONE = "none";
//const string Hyperdynamics::BOND_BOOST = "bond_boost";
const char Hyperdynamics::NONE[] = "none";
const char Hyperdynamics::BOND_BOOST[] = "bond_boost";


BondBoost::BondBoost(Matter *matt,Parameters *params)
{
    matter = matt;    
    parameters = params;
    nAtoms = matter->numberOfAtoms();
}


BondBoost::~BondBoost()
{
    return;
}


void BondBoost::initial()
{
    long i, j, k = 0, count = 0;
    bool flag  = 1;

    nBBs = 0;
    nReg = 1;
    nBAs = matter->numberOfFreeAtoms();
    nRAs = nAtoms - nBAs; // nRestAtoms
    nTABs = nBAs*(nBAs-1)/2+nBAs*nRAs; // number of Bonds involved with Tagged Atoms
    BAList = new long[nBAs]; // BoostAtomsList
    RAList = new long[nRAs]; // RestAtomsList
    TABAList = new long[2*nTABs]; // CorrespondingAtomsList of TABLList
    TABLList.setZero(nTABs,1);

    printf("BondBoost Used !\n");
    for(i = 0;i < nAtoms; i++){
        if(!matter->getFixed(i)){
            BAList[k]=i;
            k++;
        }
    }

    for(i = 0 ;i < nAtoms; i++){
        flag = 1;
        for (j = 0; j < nBAs; j++){
            if(i == BAList[j]){
                flag = 0;
            }
        }
        if(flag == 1){
            RAList[count] = i;   
            count ++;
        }
    }

    if(count != nRAs){
        printf("Error: nRestAtoms does not equal to counted number!\n");
    }

// TEST PRINT
/*
    for(i=0 ;i< nBAs;i++){
        printf("nBoostAtoms %ld: %ld\n",i,BAList[i]);
    }

    for(i=0; i< nRAs;i++){
        printf("nRestAtoms %ld: %ld\n",i,RAList[i]);
    }
*/ 
}


double BondBoost::boost()
{
    long RMDS = parameters->bondBoostRMDS;
    double dt = parameters->mdTimeStep;
    double boost_dt = 0.0, AVE_Boost_Fact;
    Matrix<double, Eigen::Dynamic, 1> TABL_tmp(nTABs,1);
    bool flag = 0;

    SDtime_B = 0.0;
    SDtime = 0.0;
    SPtime = 0.0;
    if(nReg <= RMDS){
        flag = 0;
    }else{ 
        flag = 1;
    }

    if(flag == 0){
        TABL_tmp=Rmdsteps();
        TABLList = TABLList + (1.0/RMDS)*TABL_tmp;
        nReg ++;
        SDtime += dt;
        SPtime += dt;
        // TEST PRINT
/*  for(long i=nTABs-1;i<nTABs;i++){
              printf("Distance beweem Atoms %ld and %ld is %lf\n",TABAList[2*i],TABAList[2*i+1],TABLList(i,0));
        }

        printf("\n");
*/        // TEST END
    }else{
       //    printf("nreg = %ld; RMDS=%ld\n", nReg,RMDS);
        if( nReg == RMDS+1){
              // printf("First steps = RMDS+1\n");
            nBBs = BondSelect();

           //   for (long i=0;i<nBBs;i++){
           //          printf("Equilibrium Distance between Atoms %ld and %ld is %lf\n",BBAList[2*i],BBAList[2*i+1],EBBLList(i,0));
           //    }
        }

        Epsr_Q = new double[nBBs];
        CBBLList.setZero(nBBs,1);
        boost_dt = Booststeps();
        SDtime += dt;
        SPtime += boost_dt;
        nReg++;
        AVE_Boost_Fact = (SPtime-SDtime+SDtime_B)/SDtime_B;
        //   printf("AVE_BOOST_FACT=%E;nReg = %ld\n",AVE_Boost_Fact,nReg);
        //  printf("nReg=%ld, boost_dt = %lf,SPtime = %lf\n",nReg-1,boost_dt,SPtime);
        //STC: I'm not sure if this is safe, but this should plug the memory
        //     leak.
        delete Epsr_Q;
    }
    return boost_dt;
}

double BondBoost::Booststeps()
{
    long i,j,Mi;
    long AtomI_1,AtomI_2;
    double QRR, PRR, Epsr_MAX, A_EPS_M, Sum_V, Boost_Fact, DVMAX;
    double Dforce, Fact_1, Fact_2, Mforce, dt, bdt, Temp, kb;
    double Ri[3] = {0.0} , R = 0.0;

    AtomMatrix OldForce(nAtoms,3);
    AtomMatrix Free(nAtoms,3);
    AtomMatrix AddForce(nBBs,3);
    AtomMatrix TADF(nAtoms,3);
    AtomMatrix NewForce(nAtoms,3);
 
    OldForce.setZero();
    TADF.setZero();
    Free.setZero();
    AddForce.setZero();
    NewForce.setZero();

    dt = parameters->mdTimeStep;
    Temp = parameters->temperature;
    QRR = parameters->bondBoostQRR;
    PRR = parameters->bondBoostPRR;
    DVMAX = parameters->bondBoostDVMAX;
    bdt = dt;
    kb = 1.0/11604.5;
    Epsr_MAX = 0.0;
    A_EPS_M = 0.0;
    Boost_Fact = 0.0;
    Sum_V = 0.0;
    Dforce = 0.0;
    Mforce = 0.0;
    Fact_1 = 0.0;
    Fact_2 = 0.0;

    for(i=0; i<nBBs; i++){
        AtomI_1 = BBAList[2*i];
        AtomI_2 = BBAList[2*i+1];
        CBBLList(i,0) = matter->distance(AtomI_1,AtomI_2);
    }

    for (i=0; i<nBBs; i++){
        Epsr_Q[i] = (CBBLList(i,0)-EBBLList(i,0))/EBBLList(i,0)/QRR;
        if (abs(Epsr_Q[i]) >= Epsr_MAX ){
            Epsr_MAX = abs(Epsr_Q[i]);
            Mi = i;
        }
    }

    A_EPS_M = (1.0-Epsr_MAX*Epsr_MAX)*(1.0-Epsr_MAX*Epsr_MAX)/(1.0-PRR*PRR*Epsr_MAX*Epsr_MAX);
/*
    for (i=0;i<nBBs;i++){
        printf("Boost::Distance between Atoms %ld and %ld is %lf, EQ= %lf; EPSR_Q=%lf\n",BBAList[2*i],BBAList[2*i+1],CBBLList(i,0),EBBLList(i,0),Epsr_Q[i]);
    }

    printf("Boost::Epsr_MAX= %lf, Atoms %ld and %ld\n",Epsr_MAX,BBAList[2*Mi],BBAList[2*Mi+1]);
    printf("Boost::A_EPS_M= %lf\n",A_EPS_M);
*/
    if(Epsr_MAX < 1.0){
        for(i=0;i<nBBs;i++){
            Sum_V += DVMAX*(1.0-Epsr_Q[i]*Epsr_Q[i])/double(nBBs);
        }
        Boost_Fact = A_EPS_M*Sum_V;
        SDtime_B += dt;

        for(i=0;i<nBBs;i++){
            if( abs(Epsr_Q[i]) < Epsr_MAX ){
                Fact_1 = 2.0*A_EPS_M*DVMAX*Epsr_Q[i]/QRR/EBBLList(i,0)/double(nBBs);
                Dforce = Fact_1;

                AtomI_1 = BBAList[2*i];
                AtomI_2 = BBAList[2*i+1];

                R = matter->distance(AtomI_1,AtomI_2);

                for(j=0;j<3;j++){
                    Ri[j] = matter->pdistance(AtomI_1,AtomI_2,j);
                    AddForce(i,j) = Ri[j]/R*Dforce;
                    TADF(AtomI_1,j) = TADF(AtomI_1,j)+AddForce(i,j);
                    TADF(AtomI_2,j) = TADF(AtomI_2,j)-AddForce(i,j);
                }
              //  printf("%ld Bond:: Rx=%lf;Ry=%lf;Rz=%lf;Rsum=%lf;R=%lf\n",i,Ri[0],Ri[1],Ri[2],sqrt(Ri[0]*Ri[0]+Ri[1]*Ri[1]+Ri[2]*Ri[2]),R);

            }
            else if( abs(Epsr_Q[i]) == Epsr_MAX ){
                Fact_1 = 2.0*A_EPS_M*DVMAX*Epsr_Q[i]/QRR/EBBLList(i,0)/double(nBBs);
                double Fact_tmp1 = (1.0-PRR*PRR*Epsr_Q[i]*Epsr_Q[i]);
                double Fact_tmp2 = (1.0-Epsr_Q[i]*Epsr_Q[i]);
                Fact_2 = 2.0*Fact_tmp2*Epsr_Q[i]*(2.0*Fact_tmp1-PRR*PRR*Fact_tmp2)/QRR/QRR/EBBLList(i,0)/Fact_tmp1/Fact_tmp1; 
                Dforce = Fact_1+Sum_V*Fact_2;

                AtomI_1 = BBAList[2*i];
                AtomI_2 = BBAList[2*i+1];

                Mi = i;
                Mforce = Dforce;
                R = matter->distance(AtomI_1,AtomI_2);

                for(j=0;j<3;j++){
                    Ri[j] = matter->pdistance(AtomI_1,AtomI_2,j);
                    AddForce(i,j) = Ri[j]/R*Dforce;
                    TADF(AtomI_1,j) = TADF(AtomI_1,j)+AddForce(i,j);
                    TADF(AtomI_2,j) = TADF(AtomI_2,j)-AddForce(i,j);
                }
            }
/*
            if(i != Mi){
            printf("DeltaF(Epsr_Q = %lf ), %lf  %lf  %lf; Atom %ld and Atom %ld\n", Epsr_Q[i],AddForce(i,0),AddForce(i,1),AddForce(i,2),BBAList[2*i],BBAList[2*i+1]);
        }else{
        printf("DeltaF(Epsr_MAX = %lf) = %lf, %lf  %lf  %lf; Atom %ld and Atom %ld\n", Epsr_MAX,Mforce,AddForce(i,0),AddForce(i,1),AddForce(i,2),BBAList[2*i],BBAList[2*i+1]);
       }
*/
    }

    OldForce = matter->getForces();
    NewForce = OldForce + TADF;
    Free = matter->getFree();
    NewForce = NewForce.cwise() *  Free;
    matter->setForces(NewForce);
/*       
        printf("TADF::\n");
    for(i=0;i<nAtoms;i++){
        printf("%lf   %lf   %lf\n",TADF(i,0),TADF(i,1),TADF(i,2));
        }
 
        printf("OldForces::\n");
        for(i=0;i<nAtoms;i++){
                printf("%lf   %lf   %lf\n",OldForce(i,0),OldForce(i,1),OldForce(i,2));
        }
  
        printf("NewForces::\n");
        for(i=0;i<nAtoms;i++){
                printf("%lf   %lf   %lf\n",NewForce(i,0),NewForce(i,1),NewForce(i,2));
        }
        OldForce = matter->getForces();
        printf("NewSettedForces::\n");
        for(i=0;i<nAtoms;i++){
                printf("%lf   %lf   %lf\n",OldForce(i,0),OldForce(i,1),OldForce(i,2));
        }
*/

    }else if(Epsr_MAX >= 1.0){
        Boost_Fact = 0.0;
        //Ave_Boost_Fact =
    }
    //printf("boost_fact=%lf\n",Boost_Fact);
    //update bdt;
    bdt = dt*exp(Boost_Fact/kb/Temp);
    //printf("bdt = %E\n",bdt);
    
    return bdt;
}


Matrix<double, Eigen::Dynamic,1> BondBoost::Rmdsteps()
{
    Matrix<double, Eigen::Dynamic, 1> TABL_t(nTABs,1);
    long count = 0, i, j;

    for(i=0;i<nBAs;i++){
        for(j=i+1;j<nBAs;j++){
            //printf("distance = %lf\n", matter->distance(BAList[i],BAList[j]));
            TABL_t(count,0) = matter->distance(BAList[i],BAList[j]);
            TABAList[2*count] = BAList[i];
            TABAList[2*count+1] = BAList[j];
            count ++;
        }
    }

    for(i=0;i<nBAs;i++){
        for(j=0;j<nRAs;j++){
            //printf("distance = %lf\n", matter->distance(BAList[i],RAList[j]));
            TABL_t(count,0) = matter->distance(BAList[i],RAList[j]);
            TABAList[2*count] = BAList[i];
            TABAList[2*count+1] = RAList[j];
            count ++;
        }
    }
//    printf("count=%ld\n",count);
    if(count != nTABs){
        printf("Total Involved Bond number does not equal counted number\n");
    }

    //printf("test here\n");
/*     for(i=count-1;i<count;i++){
         printf("test ::Distance beweem Atoms %ld and %ld is %lf\n",TABAList[2*i],TABAList[2*i+1],TABL_t(i,0));
    }
*/
    return TABL_t;
}

long BondBoost::BondSelect()
{
    long count = 0,i,nBBs_tmp = 0;
    double Qcutoff = parameters->bondBoostQcut;

    for(i=0;i<nTABs;i++){
        if(TABLList(i,0) <= Qcutoff){
            count ++;
        }
    }
    nBBs_tmp = count;
    EBBLList.setZero(nBBs_tmp,1);
    BBAList = new long[2*nBBs_tmp];
    count = 0;
    for(i=0;i<nTABs;i++){
        if(TABLList(i,0) <= Qcutoff){
            EBBLList(count,0) = TABLList(i,0);
            BBAList[2*count] = TABAList[2*i];
            BBAList[2*count+1] = TABAList[2*i+1];
            count ++;
        }
    }
    return nBBs_tmp;
}

