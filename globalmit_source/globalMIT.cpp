/*
  The following code is originally written by Vinh Nguyen (see below),
  but modified by Peter Lo a little in Apr 2015, to allow it to be
  used as subroutine. The original author has approved (by email) that
  the code be modified.

  Use the following to enable the original standalone exe mode.
  #define STANDALONE_EXE
 */


//////////////////////////////////////////////////////////
//                                                      //
//     globalMIT algorithm for DBN Learning		//
//     High order DBN, parallel, exact version		//
//     Author          : Vinh Nguyen                    //
//							//
//	   Multiple time series version			//
//     Project Members :				//
//                                                      //
//     Last Updated    : 13 Aug, 2011                   //
//     Contact info    : vinh.nguyenx@gmail.com         //
//                       vinh.nguyen@monash.edu         //
//							//
//                                                      //
//////////////////////////////////////////////////////////


#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <time.h>
#include <limits>
#include <vector>
#include <string>
#include <cstdlib>
#include <string.h>

#ifndef STANDALONE_EXE
namespace globalMIT {
#endif

#include "combination.h"
using namespace std;
using namespace stdcomb;
template<class BidIt>
void display(BidIt begin,BidIt end)
{
  for (BidIt it=begin;it!=end;++it)
    cout<<*it<<" ";
  cout<<endl;
}
typedef unsigned long long int UL;

//general functions
void help();
void optionInfo(int argc, int* argv[],char *inputFile);

// I/O
int readInputFile(char *inputFile);
int writeResult(char *outputFile);
 
//DBN related
double getScore(int **net);				   //get the MIT score of a DBN
int getNPa(int* &Pa,int **net, int a);	   //get the number of parents of X_i	
	
double conditional_MI_DBN(int **data,int i, int j,int nPa,int *Pa, int n_state);
double conditional_MI_DBN(int **data,int a, char b,int nPa, int* Pa, int n_state);
double conditional_MI_DBN_DBN_HighOrder(int **data,int a, int b,int nPa, int* Pa, int d,int n_state);
double conditional_MI_DBN_DBN_HighOrder_ab(int **data,int a, int b,int nPa, int* Pa, int d,int n_state);

void Contingency(int Mem1,int Mem2,int n_state);
void Contingency_HighOrder(int a,int order_a,int b,int d,int n_state);
void Contingency_HighOrder_ab(int a,int order_a,int b,int d,int n_state);

double Mutu_Info(int **T, int n_state)		;
void ClearT(int n_state);			//clear the share contingency table
	
int compare_Parent_config(int **data,int nPa,int* Pa,int a,int posi, int posj); 
int compare_Parent_config(int **data,int nPa,char* Pa,int a,int posi, int posj);
int compare_Parent_config(int nPa,int* Pa, int* order_Pa,int posi, int posj);
int compare_Parent_config_ab(int nPa,int* Pa, int* order_Pa,int posi, int posj);
	
double myEntropy(int x);			//calculate entropy of a 1-unit-shifted in time vector
double myEntropy_ab(int x);			//calculate entropy of a 1-unit-shifted in time vector

int findPstar(double *g_score,double g);
void updateBestPa(int * &best_Pa,int * Pa,int p);
	
//helper functions
int getPosition(int* powo,int p, int * Pa);
int findLexicalIndex(int n, int p, int * Pa); //Find lexical index of a parent set, variable IDs are from 0->dim-1
double optimized_double_C(UL n,UL k);//nchoosek function double version
unsigned long long nchoosek(int n,int k);//nchoosek function
unsigned long long Factorial(int num);

time_t time_=time(NULL);
void tic(){time_=time(NULL);}
double toc(){return difftime(time(NULL),time_);}
void cleanUp();

//////////////////////Global variables//////////////////////////////////
//Data & Data structure
int      N = 0;              //number of samples
int      Ne= 0;              //number of effective samples
int      dim=0;              //number of variables
int      d=1;                //high order DBN, d degree 
int**    data=NULL;          //data matrix
int      n_state=0;          //number of states, densely mapped [1:n_state], i.e Matlab style
double*  chi=NULL;           //chi values
int**    net=NULL;           //the DBN
int**    net_order=NULL;     //the DBN matrix with order info
int**    T=NULL;             //the shared contingency table
double*  g_score=NULL;       //the g-score
 
int*     dataLength=NULL;    //store the data length
int      Nseries=1;          //number of time series

// I/O & timing
ofstream outLog;             //log report
ofstream timeLog;            //time report
time_t   Start_t=time(NULL);
time_t   End_t  =time(NULL);
double   setupTime=0;

//Parameters
double   alpha=0.999;        //Description: significance level for MI test
int      maxFanIn=20;        //maxFanIn preset to 20, more than enough for most practical purposes
int      maxTime=-1;
int      allowSelfLink=0;    //allow self loop or not?

//debug variables
int      readChi=1;          //read in the chi value 
int      readNet=0;          //read in the initial net

//-------------------------------------------Main program start----------------------
#ifdef STANDALONE_EXE
int main(int argc, char* argv[]) {

  outLog.open("GlobalMIT_log.txt",ios::app);

  outLog<<endl<<"-----------GlobalMIT for High Order DBN structure learning--------------"<<endl;
  cout  <<endl<<"-----------GlobalMIT for High Order DBN structure learning--------------"<<endl;

  if (argc<2){
    help();
    return 0;
  }

  cout<<"Parameters:\n";
  for(int i = 0; i < argc; i++){
    cout   << "argv[" << i << "] = " << argv[i] << endl;
    outLog << "argv[" << i << "] = " << argv[i] << endl;
  }

	
  d=atoi(argv[1]);			//DBN order
  cout<<"DBN max order: " <<d<<endl;

  int i_start=atoi(argv[2])-1;  //starting variable, calculate parent for varialbes [i_start-i_end]
  int i_end  =atoi(argv[3])-1;
    
  char*  inputFile=argv[4];
  char*  outputFile=argv[5];

  alpha=atof(argv[6]);
  allowSelfLink=atoi(argv[7]);



  /*
    d=1;
    alpha=0.95;
    allowSelfLink=0;
    char*  inputFile="c:/New projects/globalMIT_1.1_Developer/myDBNinput.txt";
    char*  outputFile="c:/New projects/GlobalMIT_1.1_Developer/myDBNoutput.txt";
    //char*  inputFile="D:/New project/globalMIT_1.1_Developer/myDBNinput.txt";
    //char*  outputFile="D:/New project/GlobalMIT_1.1_Developer/myDBNoutput.txt";
    */

  char*  statusFile="./ParallelStatusFile.txt";
  //parsing parameters
  //optionInfo(argc,argv,inputFile);       

  //reading input
  if (readInputFile(inputFile)!=0){exit(-1);}

  //int i_end=dim-1; 
  //int i_start=0;
    
  //pre-calculate the g-score
  g_score=new double[dim+1];
  g_score[0]=0;	   //g_score[0]->0 parent, g_score[1]->1 parent...
  g_score[1]=chi[1]; //chi[0]-> 0 Parent, chi[1]-> 1 parent...
  for(int i=2;i<=dim;i++){
    g_score[i]=g_score[i-1]+chi[i];		
  }

  //the shared contingency table
  T=new int* [n_state];
  for (int i=0;i<n_state;i++){
    T[i]=new int[n_state];	
    for (int j=0;j<n_state;j++){
      T[i][j]=0;
    }
  }

  //initialize the DBN
  net=new int* [dim];
  for (int i=0;i<dim;i++){
    net[i]=new int[dim];
    for(int j=0;j<dim;j++){
      net[i][j]=0;
    }
  }

  net_order=new int* [dim];
  for (int i=0;i<dim;i++){
    net_order[i]=new int[dim];
    for(int j=0;j<dim;j++){
      net_order[i][j]=0;
    }
  }

  //////////////////Main program//////////////////////////////////////////////////

  double *best_score_arr=new double [dim];      //%best score for each individual node
  double *HX_arr=new double[dim];		//entropy of each (shifted) variable
	
  cout<<"GlobalMIT Parallel Started, startNode="<<i_start+1<<", endNode="<<i_end+1<<"... "<<endl;
	
  //main loop
  for(int i=i_start;i<=i_end;i++){//loop through the set of selected variables

    //HX_arr[i]=myEntropy(i);	
    HX_arr[i]=myEntropy_ab(i);		
    //investigate all set from 1->P*-1 elements
    int Pstar=findPstar(g_score,2*Ne*HX_arr[i]);	
    cout<<"Node "<<i+1 <<" Pstar= "<< Pstar<<endl;

	
    int *best_Pa=NULL;//empty set
    double best_s_MIT=2*Ne*HX_arr[i]; //score of the empty network
    int	 best_nPa=0;

    //co-ordinate coding: [Pai] -> p*-1 elements
    double * score=new double[dim*d]; //1-d array to store the scores of all parent combination
    double * new_score=NULL; //1-d array to store the scores of all parent combination
    int* ca=NULL;
    int* cb=NULL;

    int n_allowed_parents=0;
    if(forbidden != NULL) {
      cout<<"Parent(s) allowed for node " << i+1 << ":";
      for(int j=0; j<dim; j++) {
	if(! (forbidden[j*dim + i])) {
	  n_allowed_parents++;
	  cout << " " << j+1;
	}
      }
      n_allowed_parents *= d; /* for each delay */
      cout << endl;
    } else {
      n_allowed_parents = dim*d;
    }

    for(int p=1;p<Pstar && p<=n_allowed_parents;p++) { //loop through parent sets of increasing size

      ca=new int[dim*d];//set of all parents
      cb=new int[p];  //subset of parent of cardinality p

      // check forbidden to see what parents are allowed, if applicable.
      if(forbidden != NULL) {
	int k=0;
	for(int j=0; j<dim*d; j++) {
	  if(! (forbidden[(j%dim)*dim + i])) {ca[k]=j; k++;}
	}
	for(int j=0; j<p; j++) {cb[j]=ca[j];}
      } else {
	for(int j=0;j<dim*d;j++)     {ca[j]=j;}
	for(int j=0;j<p;j++)	   {cb[j]=j;}
      }
      //cout<<"P="<<p<<" Allocating "<<nchoosek(dim,p)<< "*"<< sizeof(double)<<" bytes = "<< nchoosek(dim,p)*sizeof(double) <<" bytes = " << double(nchoosek(dim,p))*sizeof(double)/(1024*1024)<<" Mb of Mutual Information Cache."<<endl;
      cout<<"P="<<p<<" Allocating "<< nchoosek(dim,p)*sizeof(double) <<" bytes = " << double(nchoosek(dim,p))*sizeof(double)/(1024*1024)<<" Mb of MI Cache."<<endl;
      if(p>1) {new_score=new double[nchoosek(dim*d,p)];}
      int combi_count=0; //combination count
			
      do  //generate all parent combination and score
	{
	  if(allowSelfLink==0){//check self-link
	    int selfLink=0;
	    for(int j=0;j<p;j++){
	      if(cb[j]%dim==i){
		selfLink=1;
		break;
	      } 
	    }
	    if(selfLink) continue;
	  }
	  //for(int j=0;j<p;j++) cout<<int(cb[j])<<" ";cout<<endl;
	  int pos=findLexicalIndex(dim*d,p,cb);
	  //cout<<"combi_count="<<combi_count<<"Position= "<<pos<<endl;
	  combi_count++;
	  //score this set of parents
	  if (p==1){ //only canculate the score for the 1st level
	    double CMI=conditional_MI_DBN_DBN_HighOrder_ab(data,cb[0], i, 0,cb,d,n_state);
	    double d_MIT=2*Ne*(HX_arr[i]-CMI);
	    double s_MIT=g_score[p]+d_MIT;
	    if(best_s_MIT-s_MIT>1e-12){
	      best_s_MIT=s_MIT;
	      updateBestPa(best_Pa,cb,p);
	      best_nPa=p;
	    }
	    int pos=cb[0];
	    score[pos]=CMI; //store the score
	  }else{
	    double score_i=0;
	    //get from cache
	    //int pos=getPosition(powo,p-1,cb);
	    int pos=findLexicalIndex(dim*d,p-1,cb);
	    score_i=score[pos];
					 
	    //calculate the last score and store
	    double CMI=conditional_MI_DBN_DBN_HighOrder_ab(data,cb[p-1],i, p-1 ,cb,d,n_state);
	    score_i+=CMI;	
					 
	    //pos=getPosition(powo,p,cb);
	    pos=findLexicalIndex(dim*d,p,cb);
	    new_score[pos]=score_i; //store the last calculated score
					 
	    double d_MIT=2*Ne*(HX_arr[i]-score_i);
	    double s_MIT=g_score[p]+d_MIT;
	    if(best_s_MIT-s_MIT>1e-12){
	      best_s_MIT=s_MIT;
	      updateBestPa(best_Pa,cb,p);
	      best_nPa=p;
	    }

	  }

	} //while(false)
      while(next_combination(ca,ca+n_allowed_parents,cb,cb+p));
			
      delete[] ca;
      ca=NULL;
      delete[] cb;
      cb=NULL;
      if(p>1) {
	delete[] score; 
	score=new_score;
      }
    }// of p loop


    cout<<"Complete scoring all parent sets of node " << i+1 <<"!"<<endl;
    cout<<"Best score: "<<best_s_MIT<< " Best Pa=["; for(int k=0;k<best_nPa;k++) {cout<<best_Pa[k]%dim+1<<"("<<floor(double(best_Pa[k])/dim)+1 <<") ";}cout<<"]"<<endl;
    best_score_arr[i]=best_s_MIT;
    for(int k=0;k<best_nPa;k++) {
      net[best_Pa[k]%dim][i]=1;
      net_order[best_Pa[k]%dim][i]= floor(double(best_Pa[k])/dim)+1;
    }
		
    delete[] score;
  }// of i loop


  double best_score=0;
  for(int i=0;i<dim;i++) best_score+=2*Ne*HX_arr[i]-best_score_arr[i];
  cout<<"Final network best score (original MIT) = " << best_score<<endl;

  //double score1=getScore(net_order);
  //cout<<"Best score again= "<<score1<<endl;

  //writeResult(outputFile);
  ofstream  outFile;
  outFile.open(outputFile);
  cout.precision(20);
  //outFile<<fixed <<best_score<<endl;
  for(int i=i_start;i<=i_end;i++){
    //print node name, score, then the parents
    outFile<<i+1<<" " <<2*(Ne)*HX_arr[i]-best_score_arr[i]<<" ";
    for(int j=0;j<dim;j++){
      if (net_order[j][i]>0) outFile<<j+1<<" "<<net_order[j][i]<<" ";
    }
    outFile<<endl;
  }
  outFile.close();

  ofstream statFile;
  statFile.open(statusFile,ios_base::app);
  for(int i=i_start;i<=i_end;i++){
    statFile<<i+1<<" ";
  }
  statFile.close();

  //cleaning up
  cleanUp();
  delete[] best_score_arr;
  delete[] HX_arr;
  delete[] dataLength;
  return 0;

}
#endif
///////////////////////////////////END of MAIN PROGRAM////////////////////////////////////////
int getPosition(int* powo,int p, int * Pa){ //get position in the 1-d array
  int position=0;
  for(int i=0;i<p;i++){
    position+=powo[i]*Pa[i];
  }
  return position;
}

void updateBestPa(int * &best_Pa,int * Pa,int p){ //update the best Parent set
  if (best_Pa!=NULL) delete[] best_Pa;
  best_Pa=new int[p];
  for(int i=0;i<p;i++){
    best_Pa[i]=Pa[i];
  }
}

int findPstar(double *g_score,double g){ //search for the max-fan-in
  int p=0;
  while (g_score[p]<g && p<dim){
    p++;
  }
  return p;
}

//======================Global RSC functions==========================

void optionInfo(int argc, char* argv[],char *inputFile){
  for(int i=3;i<argc;++i){
    if     (!strcmp(argv[i],"-maxFanIn"     )) maxFanIn    = atoi(argv[++i]);
    else if(!strcmp(argv[i],"-maxTime "     )) maxTime     = atoi(argv[++i]);
    else {cout << "Argument not recognized! " << argv[i] << endl; exit(0);}
  }
  cout  <<"---------------"<<endl;outLog<<"---------------"<<endl;
  cout<<"Options:"<<endl;outLog<<"Options:"<<endl;
  if (maxFanIn==1){
    cout  <<"Max-Fan-In= "<<maxFanIn<<endl;
    outLog<<"Max-Fan-In= "<<maxFanIn<<endl;
  }
  cout  <<"Maximum time: "<<maxTime<<endl;
  outLog<<"Maximum time: "<<maxTime<<endl;

  cout  <<"---------------"<<endl;outLog<<"---------------"<<endl;
}

#ifdef STANDALONE_EXE
int readInputFile(char *inputFile){
  //read the data file
  cout  <<"Reading numeric input file: "<< inputFile<<endl;
  outLog<<"Reading numeric input file: "<< inputFile<<endl;
   
  tic();   
  FILE * inFile;
  inFile = fopen (inputFile,"r");

  if (inFile==NULL) {
    cout << "Unable to open Data file";outLog << "Unable to open Data file";
    outLog.close();
    exit(-1); // terminate with error
  }
  fscanf (inFile, "%d %d %d", &N,&dim,&Nseries);
  cout<< "N= "<< N << " Dim= " <<dim<<endl;

  dataLength=new int[Nseries+1];
  dataLength[0]=-1;
  int length;
  for(int i=1;i<=Nseries;i++){
    fscanf (inFile, "%d", &length);
    dataLength[i]=dataLength[i-1]+length;
  }

  data=new int* [N];
  for(int i=0;i<N;i++){
    data[i]=new int[dim];
    for (int j=0;j<dim;j++){
      fscanf (inFile, "%d",&data[i][j]);
      if (data[i][j]>n_state) {n_state=data[i][j];}
      data[i][j]--;  //remapping the data to 1-> n_state-1;
      //cout<<data[i][j]<<" ";
    }
    //cout<<endl;
  }
   
   
  cout<<"Number of time series= "<<Nseries<<endl;
  Ne=N-d*Nseries;
  cout<<"Number of effective samples= "<<Ne<<endl;
  cout<<"N_state= "<<n_state<<endl;
  fclose (inFile);

  if (readChi){
    inFile = fopen ("./myChiValue.txt","r");
    chi=new double[dim+1];
    chi[0]=0;  //chi[0]-> 0 Parent, chi[1]-> 1 parent...
    if(dim<maxFanIn){
      for (int i=1;i<=dim;i++){
	fscanf (inFile, "%lf",&chi[i]);
	cout<<"Chi["<<i<<"]="<<chi[i]<<endl;
      }
    }
    else{
      for (int i=1;i<=maxFanIn;i++){
	fscanf (inFile, "%lf",&chi[i]);
	cout<<"Chi["<<i<<"]="<<chi[i]<<endl;
      }
      for (int i=maxFanIn;i<=dim;i++){
	//fscanf (inFile, "%lf",&chi[i]);
	//cout<<"Chi["<<i<<"]="<<chi[i]<<endl;
	chi[i]=1e99;
      }
    }
    fclose (inFile);
  }

  if (readNet){
    inFile = fopen ("./myNet.txt","r");
    net=new int* [dim];
    for(int i=0;i<dim;i++){
      net[i]=new int[dim];
      for (int j=0;j<dim;j++){
	fscanf (inFile, "%d",&net[i][j]);	
	cout<<net[i][j]<<" ";
      }
      cout<<endl;
    }
  }

  cout  <<"DONE in " << toc() << " seconds."<<endl;
  outLog<<"DONE in " << toc() << " seconds."<<endl;

  return (0);
}
#endif

void cleanUp(){
  if (data!=NULL){
    for (int i=0;i<N;i++){delete data[i];}
    delete[] data;
  }
	
  if (T!=NULL){
    for (int i=0;i<n_state;i++){delete T[i];}
    delete[] T;
  }

  if (net!=NULL){
    for (int i=0;i<dim;i++){delete net[i];}
    delete[] net;
  }
	
  if (net_order!=NULL){
    for (int i=0;i<dim;i++){delete net_order[i];}
    delete[] net_order;
  }

  if (chi!=NULL)             delete[] chi;
}

void help(){
  using namespace std;
  cout << "Usage: globalMIT.exe <maxOrder> <startVariable> <endVariable>  <inputFile> <outputFile> <alpha> [otherOptions]"<<endl;
  cout << "<maxOrder>				: max DBN order\n";
  cout << "<startVariable>		: Learing from this variable\n";
  cout << "<endVariable>			: To this variable \n";
  cout << "<inputFile>            : Ascii file containing the data matrix. The first rows contain the number of rows and collumns\n";
  cout << "<outputFile>           : Ascii file containing the graph\n";
  cout << "<alpha>				: Significance level for the Mutual Information test";
  cout << "-----------------------------------------------------------------------------------------------------\n";
  cout << "options :\n";
  cout << "-maxFanIn              : set the max-fan-in parameter\n";
  cout << "-maxTime               : set max run time\n";
  cout << "-----------------------------------------------------------------------------------------------------\n";
}

double getScore(int **net){
  double score=0;
  for(int i=0;i<dim;i++){
    int* Pa=NULL;
    int nPa=getNPa(Pa,net,i);
		
    double score_i=0;
    cout<<"Node "<<i << " nPa= "<<nPa<<endl;
    for(int j=0;j<nPa;j++){
      score_i+=2*(Ne)*conditional_MI_DBN_DBN_HighOrder(data, Pa[j],i,j,Pa,d,n_state)-chi[j+1];
    }
    cout<<"Node "<< i <<" score " << score_i<<endl;
    score+=score_i;
  }	
  cout<<"Total score: "<<score <<endl;
  return score;
}

int getNPa(int* &Pa,int **net_order, int a){  //get the number of parents of X_a and put the set in Pa
  int nPa=0;
  for(int j=0;j<dim;j++){
    nPa=nPa+net[j][a];
  }

  if(Pa!=NULL)  {delete [] Pa;}

  Pa=new int[nPa];
  int pos=0;
  for(int j=0;j<dim;j++){
    if (net_order[j][a]>0){ //add this parent
      Pa[pos]=j+dim*(net_order[j][a]-1);
      pos++;
    }
  }
  return nPa;
}


//conditional MI between node a-> node b given other parent Pa
double conditional_MI_DBN(int **data,int a, int b,int nPa, int* Pa, int n_state){
  double MI=0;

  if (nPa==0){ //no parent
    Contingency(a,b,n_state);
    return Mutu_Info(T, n_state);
  }
  else {	//with some parents?
    int  * scanned=new int[N];

    for(int i=0;i<N;i++){scanned[i]=0;}

    for(int i=0;i<N-1;i++){ //scan all rows of data
      if(scanned[i]==0){  //a new  combination of Pa found		
	scanned[i]=1;
	double count=1;
	ClearT(n_state);
	T[data[i][a]][data[i+1][b]]++;

	for(int j=i+1;j<N-1;j++){
	  if(scanned[j]==0 && compare_Parent_config(data,nPa,Pa,b,i,j)){
	    scanned[j]=1;	 				
	    T[data[j][a]][data[j+1][b]]++;
	    count++;
	  }
	}
	MI+=(count/(N-1))*Mutu_Info(T,n_state);
      }
    }
    delete[] scanned;	
  }

  return MI;
}

//conditional MI between node a-> node b given other parent Pa, high order DBN, d
double conditional_MI_DBN_DBN_HighOrder_2(int **data,int a, int b,int nPa, int* Pa, int d,int n_state){
  double MI=0;

  int real_a= a%dim;
  int order_a=floor(double(a)/dim)+1;

  if (nPa==0){ //no parent
    Contingency_HighOrder(real_a,order_a,b,d,n_state);
    return Mutu_Info(T, n_state);
  }
  else {	//with some parents?
    int * scanned=new int[N];
    int * real_Pa=new int[nPa];
    int * order_Pa=new int[nPa];
    int ** small_data;

    //find out real node ID from super node, and order
    for (int i=0;i<nPa;i++){
      order_Pa[i]=floor(double(Pa[i])/dim)+1;
      real_Pa[i]=Pa[i]%dim;
    }
	
    //make new data matrix corresponding to the concerning variables	   
    small_data=new int* [N-d];
    for (int i=0;i<N-d;i++) {small_data[i]=new int[nPa+2];}

    for (int i=0;i<nPa;i++){
      for(int j=0;j<N-d;j++){
	small_data[j][i]  =data[j+(d-order_Pa[i])][real_Pa[i]];	
      }
      real_Pa[i]=i;
    }
    for(int j=0;j<N-d;j++){  //copy parent a
      small_data[j][nPa]  =data[j+(d-order_a)][real_a];	
    }    
    for(int j=0;j<N-d;j++){  // copy target node b
      small_data[j][nPa+1]  =data[j+d][b];	
    }    
    a=nPa;  //re-index node for the new data matrix
    b=nPa+1;


    for(int i=0;i<N-d;i++){scanned[i]=0;}

    for(int i=0;i<N-d;i++){ //scan all rows of data
      if(scanned[i]==0){  //a new  combination of Pa found		
	scanned[i]=1;
	double count=1;
	ClearT(n_state);
	T[small_data[i][a]][small_data[i][b]]++;

	for(int j=i+1;j<N-d;j++){
	  if(scanned[j]==0 && compare_Parent_config(small_data,nPa,real_Pa,b,i,j)){
	    scanned[j]=1;	 				
	    T[small_data[j][a]][small_data[j][b]]++;
	    count++;
	  }
	}
	MI+=(count/(N-d))*Mutu_Info(T,n_state);
      }
    }
    delete[] scanned;	
    for (int i=0;i<N-d;i++){delete small_data[i];}
    delete[] small_data;
    delete[] order_Pa;
    delete[] real_Pa;

  }
  return MI;
}

//conditional MI between node a-> node b given other parent Pa, high order DBN, d
//no copy data
double conditional_MI_DBN_DBN_HighOrder(int **data,int a, int b,int nPa, int* Pa, int d,int n_state){
  double MI=0;

  int real_a= a%dim;
  int order_a=floor(double(a)/dim)+1;

  if (nPa==0){ //no parent
    Contingency_HighOrder(real_a,order_a,b,d,n_state);
    return Mutu_Info(T, n_state);
  }
  else {	//with some parents?
    int * scanned=new int[N];
    int * real_Pa=new int[nPa];
    int * order_Pa=new int[nPa];

    //find out real node ID from super node, and order
    for (int i=0;i<nPa;i++){
      order_Pa[i]=floor(double(Pa[i])/dim)+1;
      real_Pa[i]=Pa[i]%dim;
    }

    for(int i=0;i<N-d;i++){scanned[i]=0;}

    for(int i=0;i<N-d;i++){ //scan all rows of data
      if(scanned[i]==0){  //a new  combination of Pa found		
	scanned[i]=1;
	double count=1;
	ClearT(n_state);
	T[data[i+d-order_a][real_a]][data[i+d][b]]++;

	for(int j=i+1;j<N-d;j++){
	  if(scanned[j]==0 && compare_Parent_config(nPa,real_Pa,order_Pa,i,j)){
	    scanned[j]=1;	 				
	    T[data[j+d-order_a][real_a]][data[j+d][b]]++;
	    count++;
	  }
	}
	MI+=(count/(N-d))*Mutu_Info(T,n_state);
      }
    }
    delete[] scanned;	
    delete[] order_Pa;
    delete[] real_Pa;
  }
  return MI;
}

double conditional_MI_DBN_DBN_HighOrder_ab(int **data,int a, int b,int nPa, int* Pa, int d,int n_state){
  double MI=0;

  int real_a= a%dim;
  int order_a=floor(double(a)/dim)+1;

  if (nPa==0){ //no parent
    Contingency_HighOrder_ab(real_a,order_a,b,d,n_state);
    return Mutu_Info(T, n_state);
  }
  else {	//with some parents?
    int * scanned=new int[N];
    int * real_Pa=new int[nPa];
    int * order_Pa=new int[nPa];

    //find out real node ID from super node, and order
    for (int i=0;i<nPa;i++){
      order_Pa[i]=floor(double(Pa[i])/dim)+1;
      real_Pa[i]=Pa[i]%dim;
    }

    for(int i=0;i<N;i++){scanned[i]=1;}
    for(int l=1;l<=Nseries;l++){
      for(int i=dataLength[l-1]+1+d;i<=dataLength[l];i++){
	scanned[i]=0;
      }	
    }

    for(int i=0;i<N;i++){ //scan all rows of data
      if(scanned[i]==0){  //a new  combination of Pa found		
	scanned[i]=1;
	double count=1;
	ClearT(n_state);
	T[data[i-order_a][real_a]][data[i][b]]++;

	for(int j=i+1;j<N;j++){
	  if(scanned[j]==0 && compare_Parent_config_ab(nPa,real_Pa,order_Pa,i,j)){
	    scanned[j]=1;	 				
	    T[data[j-order_a][real_a]][data[j][b]]++;
	    count++;
	  }
	}
	MI+=(count/Ne)*Mutu_Info(T,n_state);
      }
    }
    delete[] scanned;	
    delete[] order_Pa;
    delete[] real_Pa;
  }
  return MI;
}


//compare a parent set configuration of node a at two position in the data
int compare_Parent_config(int nPa,int* Pa, int* Pa_order,int posi, int posj){
  int	isSame=1;
  for (int i=0;i<nPa;i++){ //scan through the list of parents
    if(data[posi+d-Pa_order[i]][Pa[i]]!=data[posj+d-Pa_order[i]][Pa[i]]){//check this parent value at posi & posj
      return 0;
    }
  }
  return isSame;
}

int compare_Parent_config_ab(int nPa,int* Pa, int* Pa_order,int posi, int posj){
  int	isSame=1;
  for (int i=0;i<nPa;i++){ //scan through the list of parents
    if(data[posi-Pa_order[i]][Pa[i]]!=data[posj-Pa_order[i]][Pa[i]]){//check this parent value at posi & posj
      return 0;
    }
  }
  return isSame;
}

//conditional MI between node a-> node b given other parent Pa: char type
double conditional_MI_DBN(int **data,int a, int b,int nPa, char* Pa, int n_state){
  double MI=0;

  if (nPa==0){ //no parent
    Contingency(a,b,n_state);
    return Mutu_Info(T, n_state);
  }
  else {	//with some parents?
    int  * scanned=new int[N];

    for(int i=0;i<N;i++){scanned[i]=0;}

    for(int i=0;i<N-1;i++){ //scan all rows of data
      if(scanned[i]==0){  //a new  combination of Pa found		
	scanned[i]=1;
	double count=1;
	ClearT(n_state);
	T[data[i][a]][data[i+1][b]]++;

	for(int j=i+1;j<N-1;j++){
	  if(scanned[j]==0 && compare_Parent_config(data,nPa,Pa,b,i,j)){
	    scanned[j]=1;	 				
	    T[data[j][a]][data[j+1][b]]++;
	    count++;
	  }
	}
	MI+=(count/(N-1))*Mutu_Info(T,n_state);
      }
    }
    delete[] scanned;	
  }

  return MI;
}

//compare a parent set configuration of node a at two position in the data
int compare_Parent_config(int **small_data,int nPa,int* Pa,int a,int posi, int posj){
  int	isSame=1;
  for (int i=0;i<nPa;i++){ //scan through the list of parents
    if(small_data[posi][Pa[i]]!=small_data[posj][Pa[i]]){//check this parent value at posi & posj
      return 0;
    }
  }
  return isSame;
}

//compare a parent set configuration of node a at two position in the data: char type
int compare_Parent_config(int **data,int nPa,char* Pa,int a,int posi, int posj){
  int	isSame=1;
  for (int i=0;i<nPa;i++){ //scan through the list of parents
    if(data[posi][Pa[i]]!=data[posj][Pa[i]]){//check this parent value at posi & posj
      return 0;
    }
  }
  return isSame;
}



//calculate the unconditional contingency table between node a=> node b (shifted, store to the global var T
void Contingency(int a,int b,int n_state){
  //clear up T
  ClearT(n_state);

  //build table
  for(int i =0;i<N-1;i++){
    T[data[i][a]][data[i+1][b]]++;  //note: b is one unit shifted in time
  }
}

//calculate MI, high order DBN
void Contingency_HighOrder(int a,int order_a,int b,int d,int n_state){
  //clear up T
  ClearT(n_state);

  //build table
  for(int i =d-order_a;i<N-order_a;i++){
    T[data[i][a]][data[i+order_a][b]]++;  //note: b is d units shifted in time
  }
}

void Contingency_HighOrder_ab(int a,int order_a,int b,int d,int n_state){
  //clear up T
  ClearT(n_state);
  int* scanned=new int[N];
  for(int i=0;i<N;i++){scanned[i]=1;}
  for(int l=1;l<=Nseries;l++){
    for(int i=dataLength[l-1]+1+d;i<=dataLength[l];i++){
      scanned[i]=0;
    }	
  }

  //build table
  for(int i=0;i<N;i++){
    if(scanned[i]==0){
      T[data[i-order_a][a]][data[i][b]]++;  //note: b is d units shifted in time
      scanned[i]=1;
    }
  }
  delete[] scanned;	
}



void ClearT(int n_state){
  for(int i=0;i<n_state;i++){
    for(int j=0;j<n_state;j++){
      T[i][j]=0;
    }
  }
}

double Mutu_Info(int **T, int n_state){  //get the mutual information from a contingency table
  double MI=0;
  int *a = new int[n_state];
  int *b = new int[n_state];
  int N=0;

  for(int i=0;i<n_state;i++){ //row sum
    a[i]=0;
    for(int j=0;j<n_state;j++)
      {a[i]+=T[i][j];}
  }

  for(int i=0;i<n_state;i++){ //col sum
    b[i]=0;
    for(int j=0;j<n_state;j++)
      {b[i]+=T[j][i];}
  }

  for(int i=0;i<n_state;i++) {N+=a[i];}
	
  for(int i=0;i<n_state;i++){
    for(int j=0;j<n_state;j++){
      if(T[i][j]>0){
	MI+= T[i][j]*log(double(T[i][j])*N/a[i]/b[j]);
      }
    }
  }
  delete []a;
  delete []b;

  return MI/N;
}

double myEntropy(int x){
  double *H =new double[n_state];
  for(int i=0;i<n_state;i++) {H[i]=0;}
  //entropy of a 1-unit-shifted in time vector
  for(int i=d;i<N;i++){// i run from 1
    H[data[i][x]]++;	
  }
  double e=0;
  for(int i=0;i<n_state;i++) {
    H[i]/=(N-d);
    if (H[i]!=0) {e-=H[i]*log(H[i]);}		
  }
  return e;
}

double myEntropy_ab(int x){
  double *H =new double[n_state];
  for(int i=0;i<n_state;i++){H[i]=0;}
	
  for(int l=1;l<=Nseries;l++){		
    for(int i=dataLength[l-1]+1+d;i<=dataLength[l];i++){
      H[data[i][x]]++;	
    }
  }

  double e=0;
  for(int i=0;i<n_state;i++) {
    H[i]/=Ne;
    if (H[i]!=0) {e-=H[i]*log(H[i]);}		
  }
  return e;
}

int findLexicalIndex(int n, int p, int * Pa){
  if(p==0) return -1;
  if(p==1) return Pa[0];

  int pos=0;
  int last_pos=0;

  for(int i=0;i<p-1;i++){
    if(i==0) {last_pos=0;}
    else{
      last_pos=Pa[i-1]+1;
    }
    for(int j=last_pos;j<Pa[i];j++){
      pos=pos+nchoosek(n-(j+1),p-(i+1));
    }

  }//for i

  pos=pos+Pa[p-1]-Pa[p-2]-1;
  return pos;

}


unsigned long long nchoosek(int n,int k){
  unsigned long long i,temp = 1;
  if(k > (n/2))
    k = n-k;
  for(i = n; i >= (n-k+1); i--){
    temp = temp * i;
  }
  return (temp/Factorial(k));
}



double optimized_double_C(UL n,UL k){
  double answer = 1.0;
  UL i;
  if(k > (n/2))
    k = n-k;
  for(i = 0; i < k; i++){
    answer = answer * ((double)(n-i)/(double)(i+1));
  }
  return answer;
}

unsigned long long Factorial(int num) {
  unsigned long long res = 1;
  while(num > 0){
    res = res * num;
    num = num - 1;
  }
  return res;
} 

#ifndef STANDALONE_EXE
}; // namespace
#endif

/*
  The following wrapper function is intended to allow linking to C
  program. It is copied and modified from main() and readInputFile()
  and the original main() and readInputFile() are commented out
  through preprocessor directives.

  Written by Peter Lo.
 */

#ifndef STANDALONE_EXE
extern "C" {
  struct array2d; // defined outside
  void extract_carray2d_row(array2d* d, int i, int L, int out[]); // defined outside
  double gsl_cdf_chisq_Pinv (double P, double nu); // from the GSL
  double gsl_cdf_chisq_Qinv (double Q, double nu); // from the GSL

  typedef struct s_link links; // defined outside
  links* obtain_links(int n_links);
  void output_one_link(links* Ls, int i, int from, int to, int delay); // defined outside

  int globalMIT_exe_wrapper(array2d* input_data, int n, int Dim, int n_segs, int lens[],
			    int order, int i_start, int i_end, double Alpha, int allow_self_link,
			    char* forbidden, links** out_links) {
    /*
      input_data has n rows, and dim columns. And data(i,j) is (0-based integer) expression of gene j at time i (considering the stacked series).
      There are n_segs series, and the lengths are in lens.
      maximum delay is order.
      i_start and i_end are both 0-based and inclusive (normally set to 0 and dim-1), calculate parent for varialbes [i_start-i_end].

      if forbidden is non-NULL, then if forbidden[i*dim + j] is 1, then the link i->j is forbidden for 0 <= i,j < dim.
     */
    using namespace globalMIT;

    outLog.open("GlobalMIT_log.txt",ios::app);

    outLog<<endl<<"-----------GlobalMIT for High Order DBN structure learning--------------"<<endl;
    cout  <<endl<<"-----------GlobalMIT for High Order DBN structure learning--------------"<<endl;

    d=order;			//DBN order
    cout<<"DBN max order: " <<d<<endl;

    alpha=Alpha;
    allowSelfLink=allow_self_link;

    //reading input, originally use readInputFile(), but copied and modified to here
    // we do not need to read from file here, since it is already read (or otherwise prepared) from the outside
    N = n;
    dim = Dim;
    Nseries = n_segs;
    cout<< "N= "<< N << " Dim= " <<dim<<endl;

    dataLength=new int[Nseries+1];
    dataLength[0]=-1;
    int length;
    for(int i=1;i<=Nseries;i++){
      length = lens[i-1];
      dataLength[i]=dataLength[i-1]+length;
    }

    data=new int* [N];
    n_state = 0;
    for(int i=0;i<N;i++){
      data[i]=new int[dim];
      extract_carray2d_row(input_data, i, dim, data[i]);

      for (int j=0;j<dim;j++){
	if (data[i][j]>n_state) {n_state=data[i][j];}
      }
    }
    n_state++;


    cout<<"Number of time series= "<<Nseries<<endl;
    Ne=N-d*Nseries;
    cout<<"Number of effective samples= "<<Ne<<endl;
    cout<<"N_state= "<<n_state<<endl;

    // Also, we use GSL to get the chi-square quantiles, without the use of external program to generate them in file
    chi=new double[dim+1];
    chi[0]=0;  //chi[0]-> 0 Parent, chi[1]-> 1 parent...
    double df = (n_state-1)*(n_state-1);
    for (int i=1; i<=dim && i<=maxFanIn;i++){
      chi[i] = gsl_cdf_chisq_Qinv(1-alpha, df); // it seems using Qinv can handle larger df.
      //cout << "Q(" << 1-alpha << ", df=" << df <<") = " << gsl_cdf_chisq_Qinv(1-alpha, df) << endl;
      if(chi[i] != chi[i]) { /* nan in case of error */
	chi[i] = chi[i-1]*2;
      }
      cout<<"Chi["<<i<<"]="<<chi[i]<<endl;

      df *= n_state;
    }
    for (int i=maxFanIn;i<=dim;i++){
      chi[i]=1e99;
    }
    //end of getting the input to the needed format.

    //pre-calculate the g-score
    g_score=new double[dim+1];
    g_score[0]=0;	   //g_score[0]->0 parent, g_score[1]->1 parent...
    g_score[1]=chi[1]; //chi[0]-> 0 Parent, chi[1]-> 1 parent...
    for(int i=2;i<=dim;i++){
      g_score[i]=g_score[i-1]+chi[i];		
    }

    //the shared contingency table
    T=new int* [n_state];
    for (int i=0;i<n_state;i++){
      T[i]=new int[n_state];	
      for (int j=0;j<n_state;j++){
	T[i][j]=0;
      }
    }

    //initialize the DBN
    net=new int* [dim];
    for (int i=0;i<dim;i++){
      net[i]=new int[dim];
      for(int j=0;j<dim;j++){
	net[i][j]=0;
      }
    }

    net_order=new int* [dim];
    for (int i=0;i<dim;i++){
      net_order[i]=new int[dim];
      for(int j=0;j<dim;j++){
	net_order[i][j]=0;
      }
    }

    //////////////////Main program//////////////////////////////////////////////////
	
    double *best_score_arr=new double [dim];      //%best score for each individual node
    double *HX_arr=new double[dim];		//entropy of each (shifted) variable
	
    cout<<"GlobalMIT Parallel Started, startNode="<<i_start+1<<", endNode="<<i_end+1<<"... "<<endl;
	
    //main loop
    for(int i=i_start;i<=i_end;i++){//loop through the set of selected variables

      //HX_arr[i]=myEntropy(i);	
      HX_arr[i]=myEntropy_ab(i);		
      //investigate all set from 1->P*-1 elements
      int Pstar=findPstar(g_score,2*Ne*HX_arr[i]);	
      cout<<"Node "<<i+1 <<" Pstar= "<< Pstar<<endl;

	
      int *best_Pa=NULL;//empty set
      double best_s_MIT=2*Ne*HX_arr[i]; //score of the empty network
      int	 best_nPa=0;

      //co-ordinate coding: [Pai] -> p*-1 elements
      double * score=new double[dim*d]; //1-d array to store the scores of all parent combination
      double * new_score=NULL; //1-d array to store the scores of all parent combination
      int* ca=NULL;
      int* cb=NULL;

      for(int p=1;p<Pstar;p++) { //loop through parent sets of increasing size
		
	ca=new int[dim*d];//set of all parents
	cb=new int[p];  //subset of parent of cardinality p
	for(int j=0;j<dim*d;j++)     {ca[j]=j;}
	for(int j=0;j<p;j++)	   {cb[j]=j;}
	//cout<<"P="<<p<<" Allocating "<<nchoosek(dim,p)<< "*"<< sizeof(double)<<" bytes = "<< nchoosek(dim,p)*sizeof(double) <<" bytes = " << double(nchoosek(dim,p))*sizeof(double)/(1024*1024)<<" Mb of Mutual Information Cache."<<endl;
	cout<<"P="<<p<<" Allocating "<< nchoosek(dim,p)*sizeof(double) <<" bytes = " << double(nchoosek(dim,p))*sizeof(double)/(1024*1024)<<" Mb of MI Cache."<<endl;
	if(p>1) {new_score=new double[nchoosek(dim*d,p)];}
	int combi_count=0; //combination count
			
	do  //generate all parent combination and score
	  {
	    if(allowSelfLink==0){//check self-link
	      int selfLink=0;
	      for(int j=0;j<p;j++){
		if(cb[j]%dim==i){
		  selfLink=1;
		  break;
		} 
	      }
	      if(selfLink) continue;
	    }
	    //for(int j=0;j<p;j++) cout<<int(cb[j])<<" ";cout<<endl;
	    int pos=findLexicalIndex(dim*d,p,cb);
	    //cout<<"combi_count="<<combi_count<<"Position= "<<pos<<endl;
	    combi_count++;
	    //score this set of parents
	    if (p==1){ //only canculate the score for the 1st level
	      double CMI=conditional_MI_DBN_DBN_HighOrder_ab(data,cb[0], i, 0,cb,d,n_state);
	      double d_MIT=2*Ne*(HX_arr[i]-CMI);
	      double s_MIT=g_score[p]+d_MIT;
	      if(best_s_MIT-s_MIT>1e-12){
		best_s_MIT=s_MIT;
		updateBestPa(best_Pa,cb,p);
		best_nPa=p;
	      }
	      int pos=cb[0];
	      score[pos]=CMI; //store the score
	    }else{
	      double score_i=0;
	      //get from cache
	      //int pos=getPosition(powo,p-1,cb);
	      int pos=findLexicalIndex(dim*d,p-1,cb);
	      score_i=score[pos];
					 
	      //calculate the last score and store
	      double CMI=conditional_MI_DBN_DBN_HighOrder_ab(data,cb[p-1],i, p-1 ,cb,d,n_state);
	      score_i+=CMI;	
					 
	      //pos=getPosition(powo,p,cb);
	      pos=findLexicalIndex(dim*d,p,cb);
	      new_score[pos]=score_i; //store the last calculated score
					 
	      double d_MIT=2*Ne*(HX_arr[i]-score_i);
	      double s_MIT=g_score[p]+d_MIT;
	      if(best_s_MIT-s_MIT>1e-12){
		best_s_MIT=s_MIT;
		updateBestPa(best_Pa,cb,p);
		best_nPa=p;
	      }

	    }

	  } //while(false)
	while(next_combination(ca,ca+dim*d,cb,cb+p));
			
	delete[] ca;
	ca=NULL;
	delete[] cb;
	cb=NULL;
	if(p>1) {
	  delete[] score; 
	  score=new_score;
	}
      }// of p loop


      cout<<"Complete scoring all parent sets of node " << i+1 <<"!"<<endl;
      cout<<"Best score: "<<best_s_MIT<< " Best Pa=["; for(int k=0;k<best_nPa;k++) {cout<<best_Pa[k]%dim+1<<"("<<floor(double(best_Pa[k])/dim)+1 <<") ";}cout<<"]"<<endl;
      best_score_arr[i]=best_s_MIT;
      for(int k=0;k<best_nPa;k++) {
	net[best_Pa[k]%dim][i]=1;
	net_order[best_Pa[k]%dim][i]= floor(double(best_Pa[k])/dim)+1;
      }
		
      delete[] score;
    }// of i loop


    double best_score=0;
    for(int i=0;i<dim;i++) best_score+=2*Ne*HX_arr[i]-best_score_arr[i];
    cout<<"Final network best score (original MIT) = " << best_score<<endl;

    //Now we return the obtained network (using 0-based indices), without writing it to file.
    int n_links=0;
    for(int i=i_start;i<=i_end;i++){
      for(int j=0;j<dim;j++){
	if (net_order[j][i]>0) n_links++;
      }
    }
    if(out_links) {
      links* Ls = obtain_links(n_links);
      n_links = 0;
      for(int i=i_start;i<=i_end;i++){
	for(int j=0;j<dim;j++){
	  if (net_order[j][i]>0) {
	    output_one_link(Ls, n_links, j,i, net_order[j][i]);
	    n_links++;
	  }
	}
      }
      *out_links = Ls;
    }

    //cleaning up
    cleanUp();
    delete[] best_score_arr;
    delete[] HX_arr;
    delete[] dataLength;
    delete[] g_score;

    return n_links;
  }
};

#endif
