/** Authors: M.A.Pires and  SÍlvio M. Duarte Queirós
*** Main program to the Monte Carlo Simulation of the model presented in paper
*** Optimal diffusion in ecological dynamics  with Allee effect in a metapopulation
***
*** Specifically this code yields the results shown in Fig.1, but with small
*** modifications all the other results presented in our paper can be obtained.
*** Do not hesitate to contact me: piresma@cbpf.br  or pires.ma.fisica@gmail.com
***/


#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "ran2_new.h"


int simu=2;
double fatorIo = 1.0;

#define nsubp 10
#define sizesubpop 10000
#define N sizesubpop*nsubp    //int N=1000*nsubp;
#define L nsubp               //int nn=2;
#define nn 2

#define tmax      3 //number of time steps to simulate over
#define t_steady  1 // NAO USO AQUI: time to reach the steady state
#define cnt       10 // distancia entre medidas

#define lam 1.0     // DO NOT CHANGE HERE

int state[N];
int state2[N];





/***Choosing a second agent Inside Subpopulation ***/
int secondAgentInsideSubpopulation(int i1, long *Seed)
{
    int flag=0, i2, possible_i2=i1;

    //if( local_population[i]>2 ) as the S and I do not move all the subpopulation have more than 2 agents
    while(flag==0)
    {
        possible_i2 = ((int)N*ran2(Seed));
        if(possible_i2!=i1)
        {
            if( abs(state[possible_i2])==abs(state[i1]) )
            {
                flag=1;
                i2=possible_i2;
                //printf("i:%d site:%d i2:%d site2:%d     ",i1, state[possible_i1],i2,state[possible_i2]);
                return(i2);
            }
        }
    }
}




/***Choosing a third agent Inside Subpopulation ***/
int thirdAgentInsideSubpopulation(int i1, int i2, long *Seed)
{
    int flag=0, i3, possible_i3=i1;

    while(flag==0)
    {
        possible_i3 = ((int)N*ran2(Seed));

        if(possible_i3!=i1)
        {
            if(possible_i3!=i2)
            {
                if( abs(state[possible_i3]) == abs(state[i1]) )
                {
                    flag=1;
                    i3=possible_i3;
                    //printf("i:%d site:%d i2:%d site2:%d     ",i1, state[possible_i1],i2,state[possible_i2]);
                    return(i3);
                }
            }
        }
    }
}





//int main(int argc, char *argv[])
int main(int argc, char *argv[])
{

    /**********************************************************/
    long SEED=time(NULL);
    if (SEED > 0) SEED = -SEED; /* Seed passed to ran2() must initially be negative. */
    //Seed: a negative number means start a new sequence of
    //pseudo-random numbers; a positive number means continue
    //with the same sequence.  S is turned positive by ran2().
    /**********************************************************/


    int pop,i,j,jj,label,count;
    int subpop,index,newsubpop;
    int xx,aux,cont,cont2,tot;
    int site,agent,t;
    int node1,node2,node3;
    int S_cont[L+1],I_cont[L+1],vecpop[2*L];
    int vecn[L],vecIo[L],vecSo[L];
    int neig_mat[L][nn+1];



    /*****************************************************************************/
    /*****************************************************************************/
    // get int x and convert it into real/float
    // double D = 0.0; // #mobility parameter
    // int arg1=atof(argv[1]);
    // if()

    int vecD[6]; // $D=\{0.03,0.05,0.07,0.09,0.14,0.19\}$
    vecD[0] = 0.03;
    vecD[1] = 0.05;
    vecD[2] = 0.07;
    vecD[3] = 0.09;
    vecD[4] = 0.14;
    vecD[5] = 0.19;

    int arg1 = atof(argv[1]);
    double D   = vecD[arg1];

    double alp = 0.005*atof(argv[2]);
    int nSourc =    1*atof(argv[3]);
    int sample =    1*atof(argv[4]);

    //ktot = 2; // 2*atof(argv[2]);

    /*****************************************************************************/
    /*****************************************************************************/

    label=0;
    aux  = -L-1;
    for(jj=0; jj<2*L+1; jj++)
    {
        aux += 1;
        if(aux!=0)
        {
            vecpop[label]=aux;
            label++;
        }
    }

    printf("L:%d   vecpop:",L);
    for(jj=0; jj<2*L; jj++) printf("%d  ",vecpop[jj]);
    printf("\n");
    /************************************************************************/




    /************************************************************************/
    double Ioglobal=0;
    int aux2;
    aux2 = (int) (sizesubpop/nSourc);
    //printf("aux2:%d \n",aux2);
    for(jj=0; jj<nSourc; jj++)
    {
        vecn[jj]  = sizesubpop;
        vecIo[jj] = aux2*fatorIo;  // maxIo*(L-0)=(1/L)*L=1    //aux2*fatorIo*maxIo*(L-0);
        vecSo[jj] = vecn[jj]-vecIo[jj];
        Ioglobal += 1.0*vecIo[jj];
    }
    for(jj=nSourc; jj<L; jj++)
    {
        vecn[jj]  = sizesubpop;
        vecIo[jj] = 0;
        vecSo[jj] = vecn[jj]-vecIo[jj];
        printf("%d  vecn:%d  vecSo:%d   vecIo:%d  vecIo/vecn:%1.3lf\n",jj,vecn[jj],vecSo[jj],vecIo[jj],1.0*vecIo[jj]/(1.0*sizesubpop));
        Ioglobal += 1.0*vecIo[jj];
    }
    printf("fatorIo:%1.2lf  Ioglobal%1.3lf", fatorIo, Ioglobal/(1.0*N));
    printf("*** SEED:%ld  simu:%d  L:%d  ttot:%d\n",SEED,simu,L,tmax);
    /************************************************************************/




    /************************************************************************/
    for(jj=0; jj<L; jj++)
    {
        printf("full state:%d  sign:%1.1f  subpop:%d  Io:%d\n",vecpop[jj]  ,copysign(1.0,vecpop[jj]),abs(vecpop[jj]),vecIo[jj]);
        xx=2*L-jj-1;
        printf("full state:%d  sign:%1.1f  subpop:%d  So:%d\n",vecpop[xx],copysign(1.0,vecpop[xx]),abs(vecpop[xx]),vecSo[jj]);
    }



    printf("neig_mat:\n");
    for(i=0; i<L; i++)// condicoes de contorno periodicas em x:
    {
        neig_mat[i][0] = 1+i;
        neig_mat[i][1] = 1+(i+L-1)%L; //im=(i+L-1)%L;
        neig_mat[i][2] = 1+(i+1)%L; //ip=(i+1)%L;
        printf("%d  %d  %d\n",neig_mat[i][0],neig_mat[i][1],neig_mat[i][2]);
    }



    /**********************************************/
    /**********************************************/
    // make data directory
    char dir[500];
    sprintf(dir,"allee-simu%d",simu);
    mkdir(dir, 0755);// 0755 is related to Folder/Directory Permissions.
    //askubuntu.com/questions/638796/what-is-meaning-of-755-permissions-in-samba-share

    FILE *fM;
    char outNameM[500];
    /**********************************************/
    /**********************************************/

    //  int nSourc =    1*atof(argv[3]);
    sprintf(outNameM,"./%s/allee_simu%d_nSourc%d_alp%1.3lf_N%d_L%d_D%1.2lf_tmax%d_samp%d.dat",dir,simu,nSourc,alp,N,L,D,tmax,sample);
    fM=fopen(outNameM,"w");

    /************ Initialization   ***********/

    label=0;
    for(jj=0; jj<L; jj++)
    {
        if(vecIo[jj]!=0)
        {
            for(i=0; i<vecIo[jj]; i++)
            {
                state[label] = vecpop[jj];
                state2[label]= state[label];
                label++;
            }
        }
        if(vecSo[jj]!=0)
        {
            for(i=0; i<vecSo[jj]; i++)
            {
                state[label] = vecpop[2*L-jj-1];
                state2[label]= state[label];
                label++;
            }
        }
    }


//    printf("state:\n");
//    for(j=0; j<N; j++) printf("%d  ",state[j]);
//    printf(" \n");

    for(pop=0; pop<L+1; pop++)
    {
        S_cont[pop]=0;
        I_cont[pop]=0;
    }

    for(i=0; i<N; i++)
    {
        for(pop=1; pop<L+1; pop++)
        {
            if( state[i] == +pop ) S_cont[pop]++;
            else if( state[i] == -pop ) I_cont[pop]++;
        }
    }


    for(pop=1; pop<L+1; pop++) printf("S_cont[%d]:%d   I_cont[%d]:%d\n",pop,S_cont[pop],pop,I_cont[pop]);





    fprintf(fM,"%d  ",0);
    for(pop=1; pop<L+1; pop++) fprintf(fM,"%d %d  ",S_cont[pop],I_cont[pop]);
    fprintf(fM,"\n");
    /************ PRINT AND CREATE NAMEFILE  ***********/

    printf("****SEED:%ld  alp:%1.1lf   fator:%1.2lf  D:%1.2lf    L:%d  ttot:%d\n",SEED,alp,fatorIo,D,L,tmax);


    /************************   Dynamics (1): BEGIN   ***********************/
    count  = 0;
    for(t=1; t<tmax; t++)
    {
        for(node1=0; node1<N; node1++)
        {
            if(  ran2(&SEED) < D ) //ran2(Seed)<D
            {
                if(  copysign(1.0,state[node1])<0 )
                {
                    subpop        = abs(state[node1]);
                    index         =  1 + (int) ( ran2(&SEED)*nn );  //A+(int)(ran2(Seed)*(B-A))-->[A,B-1]//range [1,nn]
                    //(int) uniform(1,nn+1);//range [1,nn]
                    newsubpop     = neig_mat[subpop-1][index];
                    state2[node1] = copysign(1.0,state[node1])*newsubpop;
                    //printf("%d  %d  %d  %d  %d  %d\n",node1,subpop,index,newsubpop,state[node1],state2[node1]);
                }
            }
            else
            {
                if(  copysign(1.0,state[node1])>0 )
                {
                    if(I_cont[abs(state[node1])]>2) node2=secondAgentInsideSubpopulation(node1,&SEED);
                    else node2=node1;

                    if( state[node2] == -state[node1] )
                    {
                        if(I_cont[abs(state[node1])]>2) node3=thirdAgentInsideSubpopulation(node1,node2,&SEED);
                        else node3=node1;

                        //printf("nodes: %d  %d  %d  states:%d  %d  %d \n",node1,node2,node3,state[node1],state[node2],state[node3]);
                        if( state[node3] == -state[node1] ) state2[node1]=-state[node1]; //lam=1
                    }
                }
                else
                {
                    if( ran2(&SEED) < alp ) state2[node1]=-state[node1];
                }
            }

        }//End of MC step over all the N states

        for(pop=0; pop<L+1; pop++)
        {
            S_cont[pop]=0;
            I_cont[pop]=0;
        }

        for(i=0; i<N; i++)
        {
            state[i]=state2[i]; // Parallel/Syncronous Updating

            for(pop=1; pop<L+1; pop++)
            {
                if( state[i] == +pop ) S_cont[pop]++;
                else if( state[i] == -pop ) I_cont[pop]++;
            }
        }

        //for(pop=1; pop<L+1; pop++) printf("S_cont[%d]:%d   I_cont[%d]:%d\n",pop,S_cont[pop],pop,I_cont[pop]);



        count++;
        if(count==cnt)
        {
            count=0;
            fprintf(fM,"%d  ",t);
            for(pop=1; pop<L+1; pop++) fprintf(fM,"%d %d  ",S_cont[pop],I_cont[pop]);
            fprintf(fM,"\n");
        }
    }// END time evolution
    /************************   Dynamics: END   ***********************/


    fclose(fM);

    return(0);
}


