// ================================================================
// Problem: Many-to-many hub location routing problem
//
// Date:
//
// Authors: Bruno Nonato Gomes
//
// Obs: A hybrid heuristic for MMHLP
//
// ================================================================
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <math.h>
#include <vector>
#include <cstring>
#include <sstream>
#include <ctime>
#include <cfloat>
#include <deque>
#include <algorithm>

#define ADJUSTMENT 1000
#define ZERO 0.000000000001
#define EPSILON 0.0001
using namespace std;


// ================================================================
//
// ================================================================

#include <iostream>
#include <math.h>
#include <list>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <qsopt.h>
extern "C" {
	#include <concorde.h>
}

void solving_tsp_concorde(vector< vector<int> > * dist, vector<int> * rota);

int main(int argc, char* argv[]) {
	//Dataset
	int initialize[5][5]=
    {
        {0,3,4,2,7},//A
        {3,0,4,6,3},//B
        {4,4,0,5,8},//C
        {2,6,5,0,6},//D
        {7,3,8,6,0}//E
    };

    //copying initial dataset for the distance structure
    vector< vector<int> > * dist = new vector< vector<int> >( 5, vector<int>(5) );
    for(int i = 0; i < dist->size(); i++){
        for(int j = 0; j < dist->size(); j++){
            dist->at(i)[j] = initialize[i][j];
        }
    }

    //Output tour
	vector<int> * tour = new vector<int>(5, 0);
    cout<<"TOUR = "<<tour->size()<<endl;
    //Concorde TSP proccess
    solving_tsp_concorde(dist, tour);

    //printing final tour generated
    cout<<"TOUR = ";
    for(int i = 0; i < tour->size(); i++){
        cout<<tour->at(i)<<" ";
    }
    cout<<endl;
    return 0;
}


void solving_tsp_concorde(vector< vector<int> > * dist, vector<int> * tour) {

    for(int i = 0; i < tour->size(); i++){
        tour->at(i) = i;
    }
    if(dist->size() > 4 ){//TSP calculado apenas acima de 4 elementos
        int rval = 0;
        int semente = rand();
        double szeit, optval, *mybnd, *mytimebound;
        int ncount, success, foundtour, hit_timebound = 0;
        int *in_tour = (int *) NULL;
        int *out_tour = (int *) NULL;
        CCrandstate rstate;
        char *probname = (char *) NULL;
        static int run_silently = 1;
        CCutil_sprand(semente, &rstate);
        mybnd = (double *) NULL;
        mytimebound = (double *) NULL;
        ncount = dist->size();
        int ecount = (ncount * (ncount - 1)) / 2;
        int *elist = new int[ecount * 2];
        int *elen = new int[ecount];
        int edge = 0;
        int edge_peso = 0;
        for (int kk = 0; kk < ncount; kk++) {
            for (int m = kk + 1; m < ncount; m++) {
                if (kk != m) {
                    elist[edge] = kk;
                    elist[edge + 1] = m;
                    elen[edge_peso]	= dist->at(kk)[m];
                    edge_peso++;
                    edge = edge + 2;
                }
            }
        }
        out_tour = CC_SAFE_MALLOC (ncount, int);
        probname = CCtsp_problabel(" ");
        rval = CCtsp_solve_sparse(ncount, ecount, elist, elen, in_tour,
                out_tour, mybnd, &optval, &success, &foundtour, probname,
                mytimebound, &hit_timebound, run_silently, &rstate);
        for (int kk = 0; kk < ncount; kk++) {
            tour->at(kk) = out_tour[kk];
        }
        szeit = CCutil_zeit();
        CC_IFFREE (elist, int);
        CC_IFFREE (elen, int);
        CC_IFFREE (out_tour, int);
        CC_IFFREE (probname, char);
    }
}
