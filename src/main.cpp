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

//Tipos
typedef struct {
	int *hubs;
	int *individuo;
} solucao;
typedef struct {
	std::list<int> hub_e;
} hub;

struct Aresta { //aresta contém rotas e custos das mesmas
	double custo;
	vector<int> rota;
};

struct SolucaoVRP {
	vector<vector<int> > rotas;
	vector<double> custos;
};

typedef struct {
	char cmd_line[250];
	char exec_name[250];
	char file_name[250];
	double ex;
	double omega;
	double beta;
	double gamma;
	double alpha;
	double T;
	int nhot;
	int n;
	int nv;
	double cq;
	double custo_transp;
	vector<double> a;
	vector<double> o;
	vector<double> d;
	vector<vector<double> > w;
	vector<vector<double> > c;
	vector<double> xcoord;
	vector<double> ycoord;
	clock_t inicio, fim;
	double tempo;

	vector<int> xx;
	vector<int> hubs;
	vector<vector<int> > tour;
	vector<double> tour_length;
	double custo_solucao;

} DATA;

// ================================================================
// Cabeçalhos
// ================================================================

double custo_con_ij(int i, int j, solucao sol, DATA *ptr);
double custo_fixo(solucao sol, DATA *ptr);
double calcula_fo_total(solucao sol, DATA *ptr);
double calcula_fo_vrp_hlp(DATA *ptr, vector<SolucaoVRP> solucaoGeral);
double custo_vrps(SolucaoVRP& vrp);
double custoRota(vector<int>& rota, DATA *ptr);
void Calcula_FO(double *&FO, solucao *s, int *numero_hubs, hub *quem_e_hub,
		double &FOMax, double &FOMin, int &bestpos, DATA *ptr);

solucao *alocar_espaco(DATA *ptr);
solucao cria_vetor_solucao(solucao vetor, DATA *ptr, int p, int q);
void eliminaRota(vector<SolucaoVRP>& solucaoGeral);
void read_data(DATA *ptr);
void help(DATA *ptr);

void associar_spoke_a_hub(solucao *&s, int *numero_hubs, hub *quem_e_hub,
		int i, DATA *ptr);
solucao *inicia_populacao(double alvo, double *&OF, int *&numero_hubs,
		hub *&quem_e_hub, DATA *ptr);
void alocacao(DATA *ptr);
void solving_tsp_concorde(DATA *ptr);
vector<SolucaoVRP> split_tsp2vrp(DATA *ptr);

//buscas atribuicao
solucao busca_local(double *&FO, int pos, solucao &s, int numero_hubs,
		hub quem_e_hub, double &FOMin, double alvo, DATA *ptr);
solucao busca_local_troca_funcao(double *&FO, int pos, solucao &s,
		int numero_hubs, hub &quem_e_hub, double &FOMin, DATA *ptr);
solucao busca_local_remove_hub(double *&FO, int pos, solucao &s,
		int &numero_hubs, hub &quem_e_hub, double &FOMin, DATA *ptr);
solucao busca_local_insere_hub(double *&FO, int pos, solucao &s,
		int &numero_hubs, hub &quem_e_hub, double &FOMin, DATA *ptr);

//meta
solucao VND(double alvo, double *&FO, int pos, solucao &s, int &numero_hubs,
		hub &quem_e_hub, double &FOMin, int *&best, DATA *ptr);
vector<SolucaoVRP> VRP_VND(vector<SolucaoVRP>& solucaoGeral, DATA *&ptr);
vector<SolucaoVRP> VRP_RVND(vector<SolucaoVRP>& solucaoGeral, DATA *&ptr);
vector<SolucaoVRP> VNS_VND(vector<SolucaoVRP>& solucaoGeral, DATA *&ptr);
vector<SolucaoVRP> ILS_VND(vector<SolucaoVRP>& solucaoGeral, DATA *&ptr);

//buscas intra rotas
SolucaoVRP OrOpt(SolucaoVRP& vrp, DATA *ptr);
SolucaoVRP BL2Opt(SolucaoVRP& vrp, DATA *ptr);

//buscas inter rotas
vector<SolucaoVRP> ShiftRota(vector<SolucaoVRP>& vrp, DATA *ptr);
SolucaoVRP crossover(SolucaoVRP& vrp, DATA *ptr);
vector<SolucaoVRP> swap(vector<SolucaoVRP>& solucaoGeral, DATA *ptr);

//outras buscas
vector<SolucaoVRP> tryEmptyRoute(vector<SolucaoVRP>& solucaoGeral, DATA *ptr);
vector<SolucaoVRP> changeDepot(vector<SolucaoVRP>& vrp, DATA *ptr);

//perturbacoes
vector<SolucaoVRP> changeDepotPerturbation(vector<SolucaoVRP>& vrp, DATA *ptr);
vector<SolucaoVRP> insertDepotPerturbation(vector<SolucaoVRP>& vrp, DATA *ptr);

void print_solution(DATA *ptr, vector<SolucaoVRP> solucaoGeral, double custo);
void file_print_solution(int flag, DATA *ptr, vector<SolucaoVRP> solucaoGeral);

struct MyComparator{
	vector<vector<double> > c;
	int source;
	bool operator () (int i, int j){
		return (c[source][i]<c[source][j]);
	}
};

struct MyComparator1{
	vector<double > c;
	bool operator () (int i, int j){
		return (c[i]<c[j]);
	}
};
// ==============================================
//
// ==============================================
int main(int argc, char* argv[]) {
	// input data
	DATA *ptr = NULL;
	ptr = new DATA;
	if (ptr == NULL) {
		cout << "Erro alocando memoria" << endl;
		exit(EXIT_FAILURE);
	}
	ptr->alpha = 0.2;
	ptr->ex = 1.0;
	ptr->beta = 1.0;
	ptr->omega = 1.0;
	ptr->gamma = 1.0;
	ptr->T = 100.0;
	ptr->nhot = 0;
	strcpy(ptr->exec_name, argv[0]);
	for (int i = 0; i < argc; i++) {
		strcat(ptr->cmd_line, argv[i]);
		strcat(ptr->cmd_line, " ");
	}
	if (argc < 2) {
		help(ptr);
		exit(EXIT_FAILURE);
	} else {
		strcpy(ptr->file_name, argv[1]);
		if (argc >= 3)
			ptr->alpha = atof(argv[2]);
		if (argc >= 4)
			ptr->ex = atof(argv[3]);
		if (argc >= 5)
			ptr->beta = atof(argv[4]);
		if (argc >= 6)
			ptr->omega = atof(argv[5]);
		if (argc >= 7)
			ptr->gamma = atof(argv[6]);
		if (argc >= 8)
			ptr->T = atof(argv[7]);
		if (argc >= 9)
			ptr->nhot = atoi(argv[8]);
	}

	// ================================
	// reading the data
	// ================================
	read_data(ptr);
	ptr->nv = 1;



	// ================================
	// testes
	// ================================
//	MyComparator mc;
//	mc.c = ptr->c;
//	vector<int> ord;
//	vector<double> pont;
//	vector<double> aux;
//	pont.resize(ptr->n, 0);
//	vector< vector<int> > mat_ord;
//	for(int i = 0; i < ptr -> n; ++i){
//		ord.push_back(i);
//		aux.push_back(ptr->o[i] + ptr->d[i]);
//	}
//
//	for(int i = 0; i < ptr -> n; ++i){
//		mc.source = i;
//		sort(ord.begin(), ord.end(), mc);
//			for(int j = 0; j < ptr -> n; ++j){
//				pont[ord[j]] += j;
//				cout<<" "<<ord[j];
//			}
//			cout<<"\n";
//		mat_ord.push_back(ord);
//	}
//
//	cout<<"\n \n";
//
//	MyComparator1 mc1;
//	mc1.c = aux;
//	sort(ord.begin(), ord.end(), mc1);
//	for(int j = 0; j < ptr -> n; ++j){
//		pont[ord[j]] += (1-ptr->alpha)*j;
//		cout<<" "<<ord[j];
//	}
//	cout<<"\n";
//
//	mc1.c = ptr->a;
//	sort(ord.begin(), ord.end(), mc1);
//	for(int j = 0; j < ptr -> n; ++j){
//		pont[ord[j]] += 2*j;
//		cout<<" "<<ord[j];
//	}
//	cout<<"\n";
//
//	mc1.c = pont;
//	sort(ord.begin(), ord.end(), mc1);
//	for(int j = 0; j < ptr -> n; ++j){
//		cout<<" "<<ord[j];
//	}
//	cout<<"\n";
//	int aa;
//	cin>>aa;

	// ================================
	// fim testes
	// ================================




	// ================================
	// geting allocation by heuristic
	// ================================
	ptr-> tempo = 0;
	ptr-> inicio = clock();
	alocacao(ptr);

	// ================================
	// geting optimal tsp tour by concorde
	// ================================
	solving_tsp_concorde(ptr);

	// ================================
	// spliting tsp tour to vrp tour's
	// ================================
	vector<SolucaoVRP> solucaoGeral;
	solucaoGeral = split_tsp2vrp(ptr);
	ptr-> custo_solucao = calcula_fo_vrp_hlp(ptr, solucaoGeral);
	print_solution(ptr, solucaoGeral, ptr->custo_solucao);
	file_print_solution(1, ptr, solucaoGeral);

	// ================================
	// Local search VND procedure
	// ================================
//	solucaoGeral = VRP_RVND(solucaoGeral, ptr);
	solucaoGeral = VNS_VND(solucaoGeral, ptr);
	eliminaRota(solucaoGeral);
	ptr-> custo_solucao = calcula_fo_vrp_hlp(ptr, solucaoGeral);
	print_solution(ptr, solucaoGeral, ptr->custo_solucao);
	file_print_solution(2, ptr, solucaoGeral);

	// ================================
	// Local search shiftRota procedures
	// ================================
	// 	for (int i = 0; i < int(solucaoGeral.size()); ++i) {
	// 		bool status = true;
	// 		double custoAtual = custo_vrps(solucaoGeral[i]);
	// 		while (status) {
	// 			SolucaoVRP novaSolucao = ShiftRota(solucaoGeral[i], ptr);
	// 			double custoNova = custo_vrps(novaSolucao);
	// 			if (custoAtual - custoNova >= 0.0001) {
	// 				solucaoGeral[i] = novaSolucao;
	// 				custoAtual = custoNova;
	// 			} else {
	// 				status = false;
	// 			}
	// 		}
	// 	}
	// 	eliminaRota(solucaoGeral);
	// 	ptr-> custo_solucao = calcula_fo_vrp_hlp(ptr, solucaoGeral);
	// 	cout<<" Custo Solucao depois BL Shift: "<<ptr->custo_solucao<<endl;

	// ================================
	// Local search OrOpt procedures
	// ================================
	// 	for (int i = 0; i < int(solucaoGeral.size()); ++i) {
	// 		bool status = true;
	// 		double custoAtual = custo_vrps(solucaoGeral[i]);
	// 		while (status) {
	// 			SolucaoVRP novaSolucao = OrOpt(solucaoGeral[i], ptr);
	// 			double custoNova = custo_vrps(novaSolucao);
	// 			if (custoAtual - custoNova >= 0.0001) {
	// 				solucaoGeral[i] = novaSolucao;
	// 				custoAtual = custoNova;
	// 			} else {
	// 				status = false;
	// 			}
	// 		}
	// 	}
	// 	eliminaRota(solucaoGeral);
	// 	ptr-> custo_solucao = calcula_fo_vrp_hlp(ptr, solucaoGeral);
	// 	cout<<" Custo Solucao depois BL OrOpt: "<<ptr->custo_solucao<<endl;
	//
	//       // ================================
	//       // Local search crossover procedures
	//       // ================================
	// 	for (int i = 0; i < int(solucaoGeral.size()); ++i) {
	// 		bool status = true;
	// 		double custoAtual = custo_vrps(solucaoGeral[i]);
	// 		while (status) {
	// 			SolucaoVRP novaSolucao = crossover(solucaoGeral[i], ptr);
	// 			double custoNova = custo_vrps(novaSolucao);
	// 			if (custoAtual - custoNova >= 0.0001) {
	// 				solucaoGeral[i] = novaSolucao;
	// 				custoAtual = custoNova;
	// 			} else {
	// 				status = false;
	// 			}
	// 		}
	// 	}
	// 	eliminaRota(solucaoGeral);
	// 	ptr-> custo_solucao = calcula_fo_vrp_hlp(ptr, solucaoGeral);
	// 	cout<<" Custo Solucao depois BL Crossover "<<ptr->custo_solucao<<endl;

	// 	ptr -> fim = clock();
	// 	ptr-> tempo = ((double)(ptr->fim-ptr->inicio)/CLOCKS_PER_SEC);
	// 	cout<<"\n Tempo: "<< ptr->tempo<<endl;

	//delete ptr; //as vezes está acontecendo problema ao deletar ptr
	return 0;
}

//===============================
// Print solution procedure
//===============================
void file_print_solution(int flag, DATA *ptr, vector<SolucaoVRP> solucaoGeral) {
	ofstream fp;
	if(flag ==1){
		fp.open("inicial.txt");
	}
	else if(flag ==2){
		fp.open("BL.txt");
	}
	int max_veiculos = 0;
	for (int i = 0; i < int(solucaoGeral.size()); ++i) {
		for (int j = 0; j < int(solucaoGeral[i].rotas.size()); ++j) {
			if(int(solucaoGeral[i].rotas.size())>max_veiculos){
				max_veiculos = int(solucaoGeral[i].rotas.size());
			}
			for (int k = 1; k < int(solucaoGeral[i].rotas[j].size()); ++k) {
				fp<<" "<<solucaoGeral[i].rotas[j][k - 1];
			}
			fp<<" "<<solucaoGeral[i].rotas[j][0]<<endl;
		}
	}

	if (flag==2){
		FILE *fp;
		fp = fopen("res.txt","a");
		fprintf(fp, "%s \t %0.2f \t %0.2lf \t %d \t %d \t %0.4f \t %0.2f \n",
				ptr->file_name, ptr->alpha, ptr->T, int(solucaoGeral.size()),
				max_veiculos, ptr->custo_solucao, ptr->tempo);
		fclose(fp);
	}

}

//===============================
// VNS-VND procedure
//===============================
vector<SolucaoVRP> ILS_VND(vector<SolucaoVRP>& solucaoGeral, DATA *&ptr) {
	int iter = 0;
	double custoAtual;
	double custoNova;
	solucaoGeral = VRP_RVND(solucaoGeral, ptr);
	eliminaRota(solucaoGeral);
	custoAtual = calcula_fo_vrp_hlp(ptr, solucaoGeral);
	ptr->custo_solucao = custoAtual;
	vector<SolucaoVRP> novaSolucao = solucaoGeral;

	while (iter <= 20) {

		//perturbacao change depot
		float pert = rand();
		if(pert<=0.5){
			novaSolucao = changeDepotPerturbation(solucaoGeral, ptr);
			novaSolucao = ShiftRota(novaSolucao, ptr);
			custoNova = calcula_fo_vrp_hlp(ptr, novaSolucao);
		}
		else{
			novaSolucao = insertDepotPerturbation(solucaoGeral, ptr);
			novaSolucao = ShiftRota(novaSolucao, ptr);
			custoNova = calcula_fo_vrp_hlp(ptr, novaSolucao);
		}

		novaSolucao = VRP_RVND(novaSolucao, ptr);
		//eliminaRota(novaSolucao);
		custoNova = calcula_fo_vrp_hlp(ptr, novaSolucao);

		if (custoAtual - custoNova >= 0.0000001) {
			solucaoGeral = novaSolucao;
			custoAtual = custoNova;
			ptr->custo_solucao = custoNova;
		}


		//perturbacao insert depot
//		novaSolucao = insertDepotPerturbation(solucaoGeral, ptr);
//		novaSolucao = ShiftRota(novaSolucao, ptr);
//		custoNova = calcula_fo_vrp_hlp(ptr, novaSolucao);
//		novaSolucao = VRP_VND(novaSolucao, ptr);
		//eliminaRota(novaSolucao);
//		custoNova = calcula_fo_vrp_hlp(ptr, novaSolucao);
//		if (custoAtual - custoNova >= 0.0000001) {
//			solucaoGeral = novaSolucao;
//			custoAtual = custoNova;
//			ptr->custo_solucao = custoNova;
//		}

		iter++;
	}
	return solucaoGeral;
}

//===============================
// VNS-VND procedure
//===============================
vector<SolucaoVRP> VNS_VND(vector<SolucaoVRP>& solucaoGeral, DATA *&ptr) {
	int vizinhanca = 1;
	double custoAtual;
	double custoNova;

	solucaoGeral = VRP_RVND(solucaoGeral, ptr);
	eliminaRota(solucaoGeral);
	custoAtual = calcula_fo_vrp_hlp(ptr, solucaoGeral);
	ptr->custo_solucao = custoAtual;
	vector<SolucaoVRP> novaSolucao = solucaoGeral;
	while (vizinhanca <= 2) {
		switch (vizinhanca) {
		case 3:
			novaSolucao = tryEmptyRoute(solucaoGeral, ptr);
			novaSolucao = ShiftRota(novaSolucao, ptr);
			custoNova = calcula_fo_vrp_hlp(ptr, novaSolucao);
			novaSolucao = VRP_RVND(novaSolucao, ptr);
			//eliminaRota(novaSolucao);
			custoNova = calcula_fo_vrp_hlp(ptr, novaSolucao);

			if (custoAtual - custoNova >= 0.0000001) {
				solucaoGeral = novaSolucao;
				custoAtual = custoNova;
				ptr->custo_solucao = custoNova;
				vizinhanca = 1;
			} else {
				vizinhanca++;
			}
			break;
		case 2:
			novaSolucao = changeDepotPerturbation(solucaoGeral, ptr);
			novaSolucao = ShiftRota(novaSolucao, ptr);
			custoNova = calcula_fo_vrp_hlp(ptr, novaSolucao);

			novaSolucao = VRP_RVND(novaSolucao, ptr);
			//eliminaRota(novaSolucao);
			custoNova = calcula_fo_vrp_hlp(ptr, novaSolucao);

			if (custoAtual - custoNova >= 0.0000001) {
				solucaoGeral = novaSolucao;
				custoAtual = custoNova;
				ptr->custo_solucao = custoNova;
				vizinhanca = 1;
			} else {
				vizinhanca++;
			}
			break;
		case 1:
			novaSolucao = insertDepotPerturbation(solucaoGeral, ptr);
			novaSolucao = ShiftRota(novaSolucao, ptr);
			custoNova = calcula_fo_vrp_hlp(ptr, novaSolucao);

			novaSolucao = VRP_RVND(novaSolucao, ptr);
			//eliminaRota(novaSolucao);
			custoNova = calcula_fo_vrp_hlp(ptr, novaSolucao);

			if (custoAtual - custoNova >= 0.0000001) {
				solucaoGeral = novaSolucao;
				custoAtual = custoNova;
				ptr->custo_solucao = custoNova;
				vizinhanca = 1;
			} else {
				vizinhanca++;
			}
			break;
		}
	}
	return solucaoGeral;
}

//===============================
// RVND procedure
//===============================
vector<SolucaoVRP> VRP_RVND(vector<SolucaoVRP>& solucaoGeral, DATA *&ptr) {

	vector < int > vizinhancas;
	for(int i = 1; i <= 5; ++i){
		vizinhancas.push_back(i);
	}
	random_shuffle (vizinhancas.begin(), vizinhancas.end());

	int vizinhanca = 0;
	double custoAtual;
	double custoNova;
	bool status, melhorou;
	vector<SolucaoVRP> novaSolucao = solucaoGeral;
	custoAtual = calcula_fo_vrp_hlp(ptr, novaSolucao);

	while (vizinhanca < int(vizinhancas.size())) {
		switch (vizinhancas[vizinhanca]) {
		case 5:
			melhorou = false;
			status = true;
			while (status) {
//				cout<<"swap antes"<<endl;
//				print_solution(ptr, solucaoGeral, ptr->custo_solucao);
				novaSolucao = swap(solucaoGeral, ptr);
				custoNova = calcula_fo_vrp_hlp(ptr, novaSolucao);
				if (custoAtual - custoNova >= 0.0000001) {
					solucaoGeral = novaSolucao;
					custoAtual = custoNova;
					melhorou = true;
				} else {
					status = false;
				}
			}
//			cout<<"swap depois"<<endl;
//			print_solution(ptr, solucaoGeral, ptr->custo_solucao);
			if (melhorou) {
				vizinhanca = 0;
				random_shuffle (vizinhancas.begin(), vizinhancas.end());
			} else {
				vizinhanca++;
			}
//			melhorou = false;
//			status = true;
//			while (status) {
//				cout<<"Try antes"<<endl;
//				print_solution(ptr, solucaoGeral, ptr->custo_solucao);
//				novaSolucao = tryEmptyRoute(solucaoGeral, ptr);
//				custoNova = calcula_fo_vrp_hlp(ptr, novaSolucao);
//				if (custoAtual - custoNova >= 0.0000001) {
//					solucaoGeral = novaSolucao;
//					custoAtual = custoNova;
//					melhorou = true;
//				} else {
//					status = false;
//				}
//				cout<<"try depois"<<endl;
//				print_solution(ptr, solucaoGeral, ptr->custo_solucao);
//			}
//			if (melhorou) {
//				vizinhanca = 0;
//				random_shuffle (vizinhancas.begin(), vizinhancas.end());
//			} else {
//				vizinhanca++;
//			}
			break;
		case 2:
			melhorou = false;
//			cout<<"Or antes"<<endl;
//			print_solution(ptr, solucaoGeral, ptr->custo_solucao);
			for (int i = 0; i < int(solucaoGeral.size()); ++i) {
				status = true;
				while (status) {
					novaSolucao[i] = OrOpt(solucaoGeral[i], ptr);
					custoNova = calcula_fo_vrp_hlp(ptr, novaSolucao);
					if (custoAtual - custoNova >= 0.0000001) {
						solucaoGeral = novaSolucao;
						custoAtual = custoNova;
						melhorou = true;
					} else {
						status = false;
					}
				}
			}
//			cout<<"Or depois"<<endl;
//			print_solution(ptr, solucaoGeral, ptr->custo_solucao);
			if (melhorou) {
				vizinhanca = 0;
				random_shuffle (vizinhancas.begin(), vizinhancas.end());
			} else {
				vizinhanca++;
			}
			break;
		case 1:
			melhorou = false;
//			cout<<"2opt antes"<<endl;
//			print_solution(ptr, solucaoGeral, ptr->custo_solucao);
			for (int i = 0; i < int(solucaoGeral.size()); ++i) {
				status = true;
				while (status) {
					novaSolucao[i] = BL2Opt(solucaoGeral[i], ptr);
					custoNova = calcula_fo_vrp_hlp(ptr, novaSolucao);
					if (custoAtual - custoNova >= 0.0000001) {
						solucaoGeral = novaSolucao;
						custoAtual = custoNova;
						melhorou = true;
					} else {
						status = false;
					}
				}
			}
//			cout<<"2opt depois"<<endl;
//			print_solution(ptr, solucaoGeral, ptr->custo_solucao);

			if (melhorou) {
				vizinhanca = 0;
				random_shuffle (vizinhancas.begin(), vizinhancas.end());
			} else {
				vizinhanca++;
			}
			break;
		case 3:
			melhorou = false;
//			cout<<"cross antes"<<endl;
//			print_solution(ptr, solucaoGeral, ptr->custo_solucao);

			for (int i = 0; i < int(solucaoGeral.size()); ++i) {
				status = true;
				while (status) {
					novaSolucao[i] = crossover(solucaoGeral[i], ptr);
					custoNova = calcula_fo_vrp_hlp(ptr, novaSolucao);
					if (custoAtual - custoNova >= 0.0000001) {
						solucaoGeral = novaSolucao;
						custoAtual = custoNova;
						melhorou = true;
					} else {
						status = false;
					}
				}
			}
//			cout<<"cross depois"<<endl;
//			print_solution(ptr, solucaoGeral, ptr->custo_solucao);
			if (melhorou) {
				vizinhanca = 0;
				random_shuffle (vizinhancas.begin(), vizinhancas.end());
			} else {
				vizinhanca++;
			}
			break;
		case 4:
			melhorou = false;
			status = true;
			while (status) {
//				cout<<"shift antes"<<endl;
//				print_solution(ptr, solucaoGeral, ptr->custo_solucao);
				novaSolucao = ShiftRota(solucaoGeral, ptr);
				custoNova = calcula_fo_vrp_hlp(ptr, novaSolucao);
				if (custoAtual - custoNova >= 0.0000001) {
					solucaoGeral = novaSolucao;
					custoAtual = custoNova;
					melhorou = true;
				} else {
					status = false;
				}
			}
//			cout<<"shift depois"<<endl;
//			print_solution(ptr, solucaoGeral, ptr->custo_solucao);
			if (melhorou) {
				vizinhanca = 0;
				random_shuffle (vizinhancas.begin(), vizinhancas.end());
			} else {
				vizinhanca++;
			}
			break;
		}
	}
	return solucaoGeral;
}

//===============================
// VND procedure
//===============================
vector<SolucaoVRP> VRP_VND(vector<SolucaoVRP>& solucaoGeral, DATA *&ptr) {
	int vizinhanca = 1;
	double custoAtual;
	double custoNova;
	bool status, melhorou;
	vector<SolucaoVRP> novaSolucao = solucaoGeral;
	custoAtual = calcula_fo_vrp_hlp(ptr, novaSolucao);

	while (vizinhanca <= 5) {
		switch (vizinhanca) {
		case 5:
			melhorou = false;
			status = true;
			while (status) {
				print_solution(ptr, solucaoGeral, ptr->custo_solucao);
				novaSolucao = swap(solucaoGeral, ptr);
				custoNova = calcula_fo_vrp_hlp(ptr, novaSolucao);
				if (custoAtual - custoNova >= 0.0000001) {
					solucaoGeral = novaSolucao;
					custoAtual = custoNova;
					melhorou = true;
				} else {
					status = false;
				}
			}
			print_solution(ptr, solucaoGeral, ptr->custo_solucao);
			if (melhorou) {
				vizinhanca = 0;
				//random_shuffle (vizinhancas.begin(), vizinhancas.end());
			} else {
				vizinhanca++;
			}
//			melhorou = false;
//			status = true;
//			while (status) {
//				novaSolucao = tryEmptyRoute(solucaoGeral, ptr);
//				custoNova = calcula_fo_vrp_hlp(ptr, novaSolucao);
//				if (custoAtual - custoNova >= 0.0000001) {
//					solucaoGeral = novaSolucao;
//					custoAtual = custoNova;
//					melhorou = true;
//				} else {
//					status = false;
//				}
//			}
//			if (melhorou) {
//				vizinhanca = 1;
//			} else {
//				vizinhanca++;
//			}
//			break;
		case 2:
			melhorou = false;
			for (int i = 0; i < int(solucaoGeral.size()); ++i) {
				status = true;
				while (status) {
					novaSolucao[i] = OrOpt(solucaoGeral[i], ptr);
					custoNova = calcula_fo_vrp_hlp(ptr, novaSolucao);
					if (custoAtual - custoNova >= 0.0000001) {
						solucaoGeral = novaSolucao;
						custoAtual = custoNova;
						melhorou = true;
					} else {
						status = false;
					}
				}
			}
			if (melhorou) {
				vizinhanca = 1;
			} else {
				vizinhanca++;
			}
			break;
		case 1:
			melhorou = false;
			for (int i = 0; i < int(solucaoGeral.size()); ++i) {
				status = true;
				while (status) {
					novaSolucao[i] = BL2Opt(solucaoGeral[i], ptr);
					custoNova = calcula_fo_vrp_hlp(ptr, novaSolucao);
					if (custoAtual - custoNova >= 0.0000001) {
						solucaoGeral = novaSolucao;
						custoAtual = custoNova;
						melhorou = true;
					} else {
						status = false;
					}
				}
			}
			if (melhorou) {
				vizinhanca = 1;
			} else {
				vizinhanca++;
			}
			break;
		case 3:
			melhorou = false;
			for (int i = 0; i < int(solucaoGeral.size()); ++i) {
				status = true;
				while (status) {
					novaSolucao[i] = crossover(solucaoGeral[i], ptr);
					custoNova = calcula_fo_vrp_hlp(ptr, novaSolucao);
					if (custoAtual - custoNova >= 0.0000001) {
						solucaoGeral = novaSolucao;
						custoAtual = custoNova;
						melhorou = true;
					} else {
						status = false;
					}
				}
			}
			if (melhorou) {
				vizinhanca = 1;
			} else {
				vizinhanca++;
			}
			break;
		case 4:
			melhorou = false;
			status = true;
			while (status) {
				novaSolucao = ShiftRota(solucaoGeral, ptr);
				custoNova = calcula_fo_vrp_hlp(ptr, novaSolucao);
				if (custoAtual - custoNova >= 0.0000001) {
					solucaoGeral = novaSolucao;
					custoAtual = custoNova;
					melhorou = true;
				} else {
					status = false;
				}
			}
			if (melhorou) {
				vizinhanca = 1;
			} else {
				vizinhanca++;
			}
			break;
		}
	}
	return solucaoGeral;
}

//===============================
// Elimina rotas não utilizadas
//===============================
void eliminaRota(vector<SolucaoVRP>& solucaoGeral) {
	for (int i = 0; i < int(solucaoGeral.size()); ++i) {
		for (int j = 0; j < int(solucaoGeral[i].rotas.size()); ++j) {
			// Apaga rotas com somente um hub: exemplo 2 - 2
			if (int(solucaoGeral[i].rotas[j].size()) == 2) {
				if(int(solucaoGeral[i].rotas.size()) != 1){
					solucaoGeral[i].custos.erase(solucaoGeral[i].custos.begin() + j);
					solucaoGeral[i].rotas.erase(solucaoGeral[i].rotas.begin() + j);
				}
				else{
					solucaoGeral[i].rotas[j].erase(solucaoGeral[i].rotas[j].end());
					solucaoGeral[i].custos[j] = 0.0;
				}

			}
		}
	}
}

//===============================
// Calcula custo da rotas em cada hub
//===============================
double custo_vrps(SolucaoVRP& vrp) {
	double custo = 0.0;
	for (int r1 = 0; r1 < int(vrp.rotas.size()); ++r1) {
		custo += vrp.custos[r1];
	}
	return custo;
}

//===============================
// Calcular custo de rota
//===============================
double custoRota(vector<int>& rota, DATA *ptr) {
	double custo = 0.0;
	for (int i = 0; i < int(rota.size() - 1); i++) {
		custo += ptr->c[rota[i]][rota[i + 1]];
	}
	return custo;
}



//===============================
// insert depot perturbation
//===============================
vector<SolucaoVRP> insertDepotPerturbation(vector<SolucaoVRP>& solucaoGeral, DATA *ptr) {
	vector<SolucaoVRP> solucaoPerturbada = solucaoGeral;
	double melhorCusto = DBL_MAX;


	//escolhendo aleatório
//	int d = rand() % int(solucaoGeral.size() - 1);
//	int r = rand() % int(solucaoGeral[d].rotas.size() - 1);
//	int v = 1 + rand() % int(solucaoGeral[d].rotas[r].size() - 2);
//
//	vector<SolucaoVRP> teste = solucaoGeral;
//	bool flag = true;
//	// selecionando cliente (vertice) que virará depósito
//	vector<int> vertice;
//	vertice.push_back(teste[d].rotas[r][v]);
//
//	//inserindo cliente (vertice) como depósito
//	SolucaoVRP vrpaux;
//	vrpaux.rotas.push_back(vertice);
//	vrpaux.custos.push_back(0);
//	teste.push_back(vrpaux);
//
//	// apagando cliente (vertice) da rota
//	teste[d].rotas[r].erase(teste[d].rotas[r].begin() + v);
//	teste[d].custos[r] = custoRota(teste[d].rotas[r], ptr);
//	if(teste[d].custos[r] > ptr->T){
//		flag = false;
//	}
//	if(flag){
//		double custoAtual = calcula_fo_vrp_hlp(ptr, teste);
//		melhorCusto	= custoAtual;
//		solucaoPerturbada = teste;
//	}

	//escolhendo melhor inserção
	for (int i = 0; i < int(solucaoGeral.size()); ++i) {
		for (int r1 = 0; r1 < int(solucaoGeral[i].rotas.size()); ++r1) {
			for (int v = 1; v < int(solucaoGeral[i].rotas[r1].size()) - 1; ++v) {
				//voltando solucao original
				vector<SolucaoVRP> teste = solucaoGeral;
				bool flag = true;

				// selecionando cliente (vertice) que virará depósito
				vector<int> vertice;
				vertice.push_back(teste[i].rotas[r1][v]);

				//inserindo cliente (vertice) como depósito
				SolucaoVRP vrpaux;
				vrpaux.rotas.push_back(vertice);
				vrpaux.custos.push_back(0);
				teste.push_back(vrpaux);

				// apagando cliente (vertice) da rota
				teste[i].rotas[r1].erase(teste[i].rotas[r1].begin() + v);
				teste[i].custos[r1] = custoRota(teste[i].rotas[r1], ptr);
				if(teste[i].custos[r1] > ptr->T){
					flag = false;
				}
				if(flag){
					double custoAtual = calcula_fo_vrp_hlp(ptr, teste);
					//Verifica se custo melhor
					if (custoAtual < melhorCusto){
						melhorCusto	= custoAtual;
						solucaoPerturbada = teste;
					}
				}
			}
		}
	}
	// eliminaRota(melhorSolucao);
	// retorna melhor solucao perturbada
	return solucaoPerturbada;
}


//===============================
// change depot perturbation
//===============================
vector<SolucaoVRP> changeDepotPerturbation(vector<SolucaoVRP>& solucaoGeral, DATA *ptr) {
	vector<SolucaoVRP> solucaoPerturbada = solucaoGeral;
	double melhorCusto = DBL_MAX;


	//escolhendo aleatório
//	int d = rand() % int(solucaoGeral.size());
//	int r = rand() % int(solucaoGeral[d].rotas.size());
//	int v = 1 + rand() % int(solucaoGeral[d].rotas[r].size() - 1);
//	vector<SolucaoVRP> teste = solucaoGeral;
//	bool flag = true;
//	// transformando cliente (vertice) em depósito
//	int vertice = teste[d].rotas[r][v];
//	teste[d].rotas[r][v] = teste[d].rotas[r][0];
//	teste[d].rotas[r][0] = vertice;
//	teste[d].rotas[r][int(teste[d].rotas[r].size())-1] = vertice;
//	for (int r2 = 0; r2 < int(teste[d].rotas.size()); ++r2) {
//		teste[d].rotas[r2][0] = vertice;
//		teste[d].rotas[r2][int(teste[d].rotas[r2].size())-1] = vertice;
//		teste[d].custos[r2] = custoRota(teste[d].rotas[r2], ptr);
//		if(teste[d].custos[r2] > ptr->T){
//			flag = false;
//			break;
//		}
//	}
//	if(flag){
//		double custoAtual = calcula_fo_vrp_hlp(ptr, teste);
//		melhorCusto	= custoAtual;
//		solucaoPerturbada = teste;
//	}

	// escolhendo melhor troca
	for (int i = 0; i < int(solucaoGeral.size()); ++i) {
		for (int r1 = 0; r1 < int(solucaoGeral[i].rotas.size()); ++r1) {
			for (int v = 1; v < int(solucaoGeral[i].rotas[r1].size()) - 1; ++v) {
				//voltando solucao original
				vector<SolucaoVRP> teste = solucaoGeral;
				bool flag = true;
				// transformando cliente (vertice) em depósito
				int vertice = teste[i].rotas[r1][v];
				teste[i].rotas[r1][v] = teste[i].rotas[r1][0];
				teste[i].rotas[r1][0] = vertice;
				teste[i].rotas[r1][int(teste[i].rotas[r1].size())-1] = vertice;
				for (int r2 = 0; r2 < int(teste[i].rotas.size()); ++r2) {
					teste[i].rotas[r2][0] = vertice;
					teste[i].rotas[r2][int(teste[i].rotas[r2].size())-1] = vertice;
					teste[i].custos[r2] = custoRota(teste[i].rotas[r2], ptr);
					if(teste[i].custos[r2] > ptr->T){
						flag = false;
						break;
					}
				}
				if(flag){
					double custoAtual = calcula_fo_vrp_hlp(ptr, teste);
					//Verifica se custo melhor
					if (custoAtual < melhorCusto){
						melhorCusto	= custoAtual;
						solucaoPerturbada = teste;
					}
				}
			}
		}
	}
	// eliminaRota(melhorSolucao);
	// retorna melhor solucao perturbada
	return solucaoPerturbada;
}


//===============================
// local search change depot
//===============================
vector<SolucaoVRP> changeDepot(vector<SolucaoVRP>& solucaoGeral, DATA *ptr) {
	vector<SolucaoVRP> melhorSolucao = solucaoGeral;
	double melhorCusto = ptr->custo_solucao;

	for (int i = 0; i < int(solucaoGeral.size()); ++i) {
		for (int r1 = 0; r1 < int(solucaoGeral[i].rotas.size()); ++r1) {
			for (int v = 1; v < int(solucaoGeral[i].rotas[r1].size()) - 1; ++v) {
				//voltando solucao original
				vector<SolucaoVRP> teste = solucaoGeral;
				bool flag = true;
				int vertice = teste[i].rotas[r1][v];
				teste[i].rotas[r1][v] = teste[i].rotas[r1][0];
				teste[i].rotas[r1][0] = vertice;
				teste[i].rotas[r1][int(teste[i].rotas[r1].size())-1] = vertice;

				for (int r2 = 0; r2 < int(teste[i].rotas.size()); ++r2) {
					teste[i].rotas[r2][0] = vertice;
					teste[i].rotas[r2][int(teste[i].rotas[r2].size())-1] = vertice;

					//busca local OrOpt na troca de deposito
					bool status = true;
					double custoAtual = custo_vrps(teste[i]);
					while (status) {
						SolucaoVRP novaSolucao = teste[i];
						novaSolucao = OrOpt(teste[i], ptr);
						double custoNova = custo_vrps(novaSolucao);
						if (custoAtual - custoNova >= 0.0000001) {
							teste[i] = novaSolucao;
							custoAtual = custoNova;
						} else {
							status = false;
						}
					}

					teste[i].custos[r2] = custoRota(teste[i].rotas[r2], ptr);
					if(teste[i].custos[r2] > ptr->T){
						flag = false;
						break;
					}
				}
				if(flag){
					double custoAtualGeral = calcula_fo_vrp_hlp(ptr, teste);
					//Verifica se custo melhor
					if (custoAtualGeral < melhorCusto){
						melhorCusto	= custoAtualGeral;
						melhorSolucao = teste;
					}
				}
			}
		}
	}
	//eliminaRota(melhorSolucao);
	return melhorSolucao;
}

//===============================
// Busca try empty a route
//===============================
vector<SolucaoVRP> tryEmptyRoute(vector<SolucaoVRP>& solucaoGeral, DATA *ptr) {
	double custoMenorRota = DBL_MAX;
	int menor_i, menor_j;
	vector<SolucaoVRP> teste = solucaoGeral;

	//geting menor rota
	for (int i = 0; i < int(teste.size()); ++i) {
		for (int j = 0; j < int(teste[i].rotas.size()); ++j) {
			if (i == 0 && j == 0) {
				custoMenorRota = teste[i].custos[j];
				menor_i = i;
				menor_j = j;
			} else if (teste[i].custos[j] < custoMenorRota) {
				custoMenorRota = teste[i].custos[j];
				menor_i = i;
				menor_j = j;
			}
		}
	}

	//escolhendo aleatorio
//	menor_i = rand() % int(solucaoGeral.size());
//	menor_j = rand() % int(solucaoGeral[menor_i].rotas.size());
//	custoMenorRota = solucaoGeral[menor_i].custos[menor_j];



	deque<int> verticesRemovidos;
	//apagando rota

	if( int(teste[menor_i].rotas.size()) == 1 ){
		//remove depot
		verticesRemovidos.push_back(teste[menor_i].rotas[menor_j][0]);
		for (int i = 1; i < int(teste[menor_i].rotas[menor_j].size())-1; ++i) {
			verticesRemovidos.push_back(teste[menor_i].rotas[menor_j][i]);
		}
		//apaga vector da solucao
		teste[menor_i].rotas.erase(teste[menor_i].rotas.begin() + menor_j);
		teste[menor_i].custos.erase(teste[menor_i].custos.begin() + menor_j);
		teste.erase(teste.begin() + menor_i);
	}
	else{
		for (int i = 1; i < int(teste[menor_i].rotas[menor_j].size())-1; ++i) {
			verticesRemovidos.push_back(teste[menor_i].rotas[menor_j][i]);
		}
		teste[menor_i].rotas.erase(teste[menor_i].rotas.begin() + menor_j);
		teste[menor_i].custos.erase(teste[menor_i].custos.begin() + menor_j);
	}

	while (!verticesRemovidos.empty()) {
		int hubSelecionado = -1;
		int rotaSelecionada = -1;
		int posicaoSelecionada = -1;
		double custoSelecionada = DBL_MAX;
		for (int h = 0; h < int(teste.size()); ++h) {
			for (int r = 0; r < int(teste[h].rotas.size()); ++r) {
				for (int v1 = 1; v1 < int(teste[h].rotas[r].size()); ++v1) {
					teste[h].rotas[r].insert(teste[h].rotas[r].begin() + v1, verticesRemovidos[0]);
					double custoInsercao = ptr->c[teste[h].rotas[r][v1 - 1]][teste[h].rotas[r][v1]]
										   + ptr->c[teste[h].rotas[r][v1]][teste[h].rotas[r][v1 + 1]]
										   - ptr->c[teste[h].rotas[r][v1 - 1]][teste[h].rotas[r][v1 + 1]];
					if ((teste[h].custos[r] + custoInsercao <= ptr->T) && (custoInsercao < custoSelecionada)){
						custoSelecionada = custoInsercao;
						hubSelecionado = h;
						rotaSelecionada = r;
						posicaoSelecionada = v1;
					}
					teste[h].rotas[r].erase(teste[h].rotas[r].begin() + v1);
				}
			}
		}
		if(hubSelecionado != -1){
			teste[hubSelecionado].rotas[rotaSelecionada]
								 .insert(teste[hubSelecionado]
								 .rotas[rotaSelecionada].begin() + posicaoSelecionada, verticesRemovidos[0]);
			teste[hubSelecionado].custos[rotaSelecionada] += custoSelecionada;
			verticesRemovidos.pop_front();
		}
		else{
			return solucaoGeral;
		}
	}
	return teste;
}

//===============================
// local Search Crossover - change two arcs between two routes
//===============================
SolucaoVRP crossover(SolucaoVRP& vrp, DATA *ptr) {
	SolucaoVRP melhorSolucao = vrp;
	double melhorCusto = custo_vrps(vrp);
	for (int r1 = 0; r1 < int(vrp.rotas.size()); ++r1) {
		for (int r2 = 0; r2 < int(vrp.rotas.size()); ++r2) {
			if (r1 != r2) {
				for (int v = 1; v < int(vrp.rotas[r1].size()) - 2; ++v) {
					for (int v1 = 1; v1 < int(vrp.rotas[r2].size()) - 2; ++v1) {
						SolucaoVRP teste = vrp;
						teste.rotas[r1].erase(teste.rotas[r1].begin() + v + 1,
								teste.rotas[r1].end());
						teste.rotas[r2].erase(teste.rotas[r2].begin() + v1 + 1,
								teste.rotas[r2].end());
						teste.rotas[r1].insert(teste.rotas[r1].begin() + v + 1,
								vrp.rotas[r2].begin() + v1 + 1,
								vrp.rotas[r2].end());
						teste.rotas[r2].insert(
								teste.rotas[r2].begin() + v1 + 1,
								vrp.rotas[r1].begin() + v + 1,
								vrp.rotas[r1].end());
						teste.custos[r1] = custoRota(teste.rotas[r1], ptr);
						teste.custos[r2] = custoRota(teste.rotas[r2], ptr);
						if ((teste.custos[r1] <= ptr->T) && (teste.custos[r2]
								<= ptr->T)) {
							double custoNova = custo_vrps(teste);
							//Verifica capacidade da rota e se custo melhor
							if (custoNova < melhorCusto) {
								melhorCusto = custoNova;
								melhorSolucao = teste;
							}
						}
					}
				}
			}
		}
	}

	return melhorSolucao;
}

//===============================
// local search swap - change two vertex between two routes
//===============================
vector<SolucaoVRP> swap(vector<SolucaoVRP>& solucaoGeral, DATA *ptr) {
	vector<SolucaoVRP> melhorSolucao = solucaoGeral;
	vector<SolucaoVRP> teste = solucaoGeral;
	double melhorCusto = ptr->custo_solucao;
	for (int i = 0; i < int(teste.size()); ++i) {
		for (int j = 0; j < int(teste.size()); ++j) {
			for (int r1 = 0; r1 < int(teste[i].rotas.size()); ++r1) {
				for (int r2 = 0; r2 < int(teste[j].rotas.size()); ++r2) {
					// bloqueia rotas iguas de um mesmo deposito
					if ((i!=j) || ( i==j && r1 != r2)) {
						for (int v = 1; v < int(teste[i].rotas[r1].size()) - 1; ++v) {

							int vertice = teste[i].rotas[r1][v];
							double	custoRemocao = ptr->c[solucaoGeral[i].rotas[r1][v - 1]][solucaoGeral[i].rotas[r1][v	+ 1]]
													- ptr->c[solucaoGeral[i].rotas[r1][v - 1]][solucaoGeral[i].rotas[r1][v]]
													- ptr->c[solucaoGeral[i].rotas[r1][v]][solucaoGeral[i].rotas[r1][v + 1]];
							teste[i].rotas[r1].erase(teste[i].rotas[r1].begin() + v);

							for (int v1 = 1; v1 < int(teste[j].rotas[r2].size()) - 1; ++v1) {
								int vertice1 = teste[j].rotas[r2][v1];
								double	custoRemocao1 = ptr->c[solucaoGeral[j].rotas[r2][v1 - 1]][solucaoGeral[j].rotas[r2][v1	+ 1]]
														- ptr->c[solucaoGeral[j].rotas[r2][v1 - 1]][solucaoGeral[j].rotas[r2][v1]]
														- ptr->c[solucaoGeral[j].rotas[r2][v1]][solucaoGeral[j].rotas[r2][v1 + 1]];
								teste[j].rotas[r2].erase(teste[j].rotas[r2].begin() + v1);

								//trocando vertices entre rotas
								teste[j].rotas[r2].insert(teste[j].rotas[r2].begin() + v1, vertice);
								teste[i].rotas[r1].insert(teste[i].rotas[r1].begin() + v, vertice1);
								double novoCusto = calcula_fo_vrp_hlp(ptr, teste);


								double custoInsercao1 = ptr->c[teste[j].rotas[r2][v1 - 1]][teste[j].rotas[r2][v1]]
								                       + ptr->c[teste[j].rotas[r2][v1]][teste[j].rotas[r2][v1 + 1]]
													   - ptr->c[teste[j].rotas[r2][v1 - 1]][teste[j].rotas[r2][v1 + 1]];

								double custoInsercao = ptr->c[teste[i].rotas[r1][v - 1]][teste[i].rotas[r1][v]]
								                       + ptr->c[teste[i].rotas[r1][v]][teste[i].rotas[r1][v + 1]]
													   - ptr->c[teste[i].rotas[r1][v - 1]][teste[i].rotas[r1][v + 1]];
								//Verifica capacidade da rota e se custo melhor
								if (((teste[j].custos[r2] + custoInsercao1 + custoRemocao1 <= ptr->T)
										&& (teste[i].custos[r1] + custoInsercao + custoRemocao <= ptr->T) )
										&& (novoCusto < melhorCusto)){
									melhorCusto	= novoCusto;
									ptr->custo_solucao = melhorCusto;
									melhorSolucao = teste;
									melhorSolucao[i].custos[r1] = solucaoGeral[i].custos[r1] + custoRemocao + custoInsercao;
									melhorSolucao[j].custos[r2] = solucaoGeral[j].custos[r2] + custoRemocao1 + custoInsercao1;
								}

								teste[j].rotas[r2].erase(teste[j].rotas[r2].begin() + v1);
								teste[i].rotas[r1].erase(teste[i].rotas[r1].begin() + v);

								teste[j].rotas[r2].insert(teste[j].rotas[r2].begin() + v1, vertice1);
							}
							teste[i].rotas[r1].insert(teste[i].rotas[r1].begin() + v, vertice);
						}
					}
				}
			}
		}
	}
	return melhorSolucao;
}



//===============================
// local Search ShiftRota - Insert one vertex of one route in other route
//===============================
vector<SolucaoVRP> ShiftRota(vector<SolucaoVRP>& solucaoGeral, DATA *ptr) {
	vector<SolucaoVRP> melhorSolucao = solucaoGeral;
	vector<SolucaoVRP> teste = solucaoGeral;
	double melhorCusto = ptr->custo_solucao;
	for (int i = 0; i < int(teste.size()); ++i) {
		for (int j = 0; j < int(teste.size()); ++j) {
			for (int r1 = 0; r1 < int(teste[i].rotas.size()); ++r1) {
				for (int r2 = 0; r2 < int(teste[j].rotas.size()); ++r2) {
					// bloqueia rotas iguas de um mesmo deposito
					if ((i!=j) || ( i==j && r1 != r2)) {
						for (int v = 1; v < int(teste[i].rotas[r1].size()) - 1; ++v) {
							int vertice = teste[i].rotas[r1][v];
							double	custoRemocao = ptr->c[solucaoGeral[i].rotas[r1][v - 1]][solucaoGeral[i].rotas[r1][v	+ 1]]
													- ptr->c[solucaoGeral[i].rotas[r1][v - 1]][solucaoGeral[i].rotas[r1][v]]
													- ptr->c[solucaoGeral[i].rotas[r1][v]][solucaoGeral[i].rotas[r1][v + 1]];
							teste[i].rotas[r1].erase(teste[i].rotas[r1].begin() + v);
							for (int v1 = 1; v1 < int(teste[j].rotas[r2].size()); ++v1) {
								teste[j].rotas[r2].insert(teste[j].rotas[r2].begin() + v1, vertice);
								double novoCusto = calcula_fo_vrp_hlp(ptr, teste);
								double custoInsercao = ptr->c[teste[j].rotas[r2][v1 - 1]][teste[j].rotas[r2][v1]]
								                       + ptr->c[teste[j].rotas[r2][v1]][teste[j].rotas[r2][v1 + 1]]
													   - ptr->c[teste[j].rotas[r2][v1 - 1]][teste[j].rotas[r2][v1 + 1]];
								//Verifica capacidade da rota e se custo melhor
								if ((teste[j].custos[r2] + custoInsercao <= ptr->T)
										&& (novoCusto < melhorCusto)){
									melhorCusto	= novoCusto;
									ptr->custo_solucao = melhorCusto;
									melhorSolucao = teste;
									melhorSolucao[i].custos[r1] = solucaoGeral[i].custos[r1] + custoRemocao;
									melhorSolucao[j].custos[r2] = solucaoGeral[j].custos[r2] + custoInsercao;
								}
								teste[j].rotas[r2].erase(teste[j].rotas[r2].begin()	+ v1);
							}
							teste[i].rotas[r1].insert(teste[i].rotas[r1].begin() + v, vertice);
						}
					}
				}
			}
		}
	}
	return melhorSolucao;
}

//===============================
// Local Search 2Opt - change two arcs in one route
//===============================
SolucaoVRP BL2Opt(SolucaoVRP& vrp, DATA *ptr) {
	double melhorCusto = 0.0;
	SolucaoVRP melhorSolucao = vrp;
	// seleciona uma rota para realizar o movimento
	for (int r = 0; r < int(vrp.rotas.size()); ++r) {
		for (int c = 1; c < int(vrp.rotas[r].size()) - 1; ++c) {
			for (int c1 = c + 1; c1 < int(vrp.rotas[r].size()) - 1; ++c1) {
				SolucaoVRP teste = vrp;
				double custoRemocao = ptr->c[vrp.rotas[r][c - 1]][vrp.rotas[r][c]]
                                 + ptr->c[vrp.rotas[r][c1]][vrp.rotas[r][c1 + 1]];
				double custoInsercao = ptr->c[vrp.rotas[r][c - 1]][vrp.rotas[r][c1]]
                                 + ptr->c[vrp.rotas[r][c]][vrp.rotas[r][c1 + 1]];

				if((custoInsercao - custoRemocao) < melhorCusto){
					melhorCusto = custoInsercao - custoRemocao;

					for( int i = c; i <= c1; ++i)
					{
						teste.rotas[r][i] = vrp.rotas[r][c1-i+c];
					}

					teste.custos[r] += melhorCusto;
					melhorSolucao = teste;
				}
			}
		}
	}
	return melhorSolucao;
}

//===============================
// Local Search OrOpt - change the position of one vertex in one route
//===============================
SolucaoVRP OrOpt(SolucaoVRP& vrp, DATA *ptr) {
	double melhorCusto = 0;
	SolucaoVRP melhorSolucao = vrp;
	// seleciona uma rota para realizar o movimento
	for (int r = 0; r < int(vrp.rotas.size()); ++r) {
		for (int c = 1; c < int(vrp.rotas[r].size()) - 1; ++c) {
			SolucaoVRP teste = vrp;
			int vertice = teste.rotas[r][c];
			double custoBase = ptr->c[vrp.rotas[r][c - 1]][vrp.rotas[r][c + 1]]
					- ptr->c[vrp.rotas[r][c - 1]][vrp.rotas[r][c]]
					- ptr->c[vrp.rotas[r][c]][vrp.rotas[r][c + 1]];
			teste.rotas[r].erase(teste.rotas[r].begin() + c);
			for (int i = 1; i < int(teste.rotas[r].size()); ++i) {
				if (i != c) {
					teste.rotas[r].insert(teste.rotas[r].begin() + i, vertice);
					double custo = custoBase
							+ ptr->c[teste.rotas[r][i - 1]][teste.rotas[r][i]]
							+ ptr->c[teste.rotas[r][i]][teste.rotas[r][i + 1]]
							- ptr->c[teste.rotas[r][i - 1]][teste.rotas[r][i
									+ 1]];
					if (custo < melhorCusto) {
						melhorCusto = custo;
						melhorSolucao = teste;
						melhorSolucao.custos[r] = vrp.custos[r] + custo;
					}
					teste.rotas[r].erase(teste.rotas[r].begin() + i);
				}
			}
		}
	}
	return melhorSolucao;
}

//===============================
// Print solution procedure
//===============================
void print_solution(DATA *ptr, vector<SolucaoVRP> solucaoGeral, double custo) {
	for (int i = 0; i < int(solucaoGeral.size()); ++i) {
		cout << "\n Hub: " << solucaoGeral[i].rotas[0][0] << endl;
		for (int j = 0; j < int(solucaoGeral[i].rotas.size()); ++j) {
			//Imprimindo veiculos
			cout << " \n      " << j + 1 << ": ";
			for (int k = 1; k < int(solucaoGeral[i].rotas[j].size()); ++k) {
				cout << solucaoGeral[i].rotas[j][k - 1] << " ";
			}
			cout << solucaoGeral[i].rotas[j][0] << " ";
			cout << " --> " << solucaoGeral[i].custos[j] << endl;
		}
		cout << endl;
	}
	double custo1 = calcula_fo_vrp_hlp(ptr, solucaoGeral);
	printf("\n Custo Solucao = %8.4f ", custo1);
	ptr -> fim = clock();
	ptr-> tempo = ((double) (ptr->fim - ptr->inicio) / CLOCKS_PER_SEC);
	printf("\n Tempo: %1.2f \n", ptr->tempo);
}

//===============================
// Compute FO geral procedure
//===============================
double calcula_fo_vrp_hlp(DATA *ptr, vector<SolucaoVRP> solucaoGeral) {
	double custo_rotas = 0.0;
	ptr->custo_transp = 0.0;
	for (int i = 0; i < int(solucaoGeral.size()); ++i) {
		ptr->custo_transp += ptr->ex * ptr->a[solucaoGeral[i].rotas[0][0]];
		ptr->xx[solucaoGeral[i].rotas[0][0]] = solucaoGeral[i].rotas[0][0];
		for (int j = 0; j < int(solucaoGeral[i].rotas.size()); ++j) {
			for (int k = 1; k < int(solucaoGeral[i].rotas[j].size()); ++k) {
				custo_rotas += ptr->gamma * ptr->c[solucaoGeral[i].rotas[j][k
						- 1]][solucaoGeral[i].rotas[j][k]];

				ptr->xx[solucaoGeral[i].rotas[j][k]]
						= solucaoGeral[i].rotas[j][0];
			}
			if (solucaoGeral[i].rotas[j].size() > 2) {
				// vehicle fixed cost
				custo_rotas += ptr->cq * ptr->omega;
			}
		}
	}

	// Allocation  costs
	for (int i = 0; i < ptr->n; i++) {
		for (int j = i + 1; j < ptr->n; j++) {
			ptr->custo_transp += ptr->w[i][j] * (ptr->c[i][ptr->xx[i]] + 10000
					* max(0.0, 2 * ptr->c[i][ptr->xx[i]] - ptr->T) + ptr->alpha
					* ptr->c[ptr->xx[i]][ptr->xx[j]] + ptr->c[ptr->xx[j]][j]
					+ 10000 * max(0.0, 2 * ptr->c[ptr->xx[j]][j] - ptr->T))
					+ ptr->w[j][i] * (ptr->c[j][ptr->xx[j]] + 10000 * max(0.0,
							2 * ptr->c[j][ptr->xx[j]] - ptr->T) + ptr->alpha
							* ptr->c[ptr->xx[j]][ptr->xx[i]]
							+ ptr->c[ptr->xx[i]][i] + 10000 * max(0.0, 2
							* ptr->c[ptr->xx[i]][i] - ptr->T));
		}
	}
	return ptr->custo_transp + custo_rotas;
}

//================================
// Solving TSP by concorde
//================================
void solving_tsp_concorde(DATA *ptr) {
	ptr->tour.resize(ptr->n);
	for (int m = 0; m < ptr->n; m++)
		ptr->tour[m].resize(ptr->n, -1);
	int cont;
	for (int k = 0; k < int(ptr->hubs.size()); k++) {
		vector<int> aux;
		aux.resize(ptr->n, -1);
		ptr->tour[k][0] = ptr->hubs[k];
		aux[0] = ptr->hubs[k];
		cont = 1;
		for (int i = 0; i < ptr->n; i++) {
			if ((i != ptr->xx[i]) && (ptr->xx[i] == ptr->hubs[k])) {
				ptr->tour[k][cont] = i;
				aux[cont] = i;
				cont++;
			}
		}
		if ((cont - 1) >= 4) { //tsp tour only exits with 3 or more nodes
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
			ncount = cont;
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
						elen[edge_peso]
								= ptr->c[ptr->tour[k][kk]][ptr->tour[k][m]];
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
				ptr->tour[k][kk] = aux[out_tour[kk]];
			}
			szeit = CCutil_zeit();
			CC_IFFREE (elist, int);
			CC_IFFREE (elen, int);
			CC_IFFREE (out_tour, int);
			CC_IFFREE (probname, char);
		}
		if ((cont - 1) != 0) { //iserting the return arc deposit
			ptr->tour[k][cont] = ptr->tour[k][0];
		}
		aux.clear();
	}
}

//===============================
// Menor caminho - Dijkstra
//===============================
deque<int> menorCaminho(vector<vector<Aresta> > &grafo, int origem, int destino) {

	deque<int> solucao;
	deque<int> antecessores(grafo.size(), -1);
	deque<double> custos(grafo.size(), DBL_MAX);
	deque<bool> explorado(grafo.size(), false);

	custos[origem] = 0;
	explorado[origem] = true;

	int ultimoExplorado = origem;

	do {
		explorado[ultimoExplorado] = true;
		int vizinho = 0;
		while (vizinho != int(grafo.size())) {

			if (!explorado[vizinho] && grafo[ultimoExplorado][vizinho].custo
					!= -1) {
				if (custos[vizinho] > custos[ultimoExplorado]
						+ grafo[ultimoExplorado][vizinho].custo) {
					custos[vizinho] = custos[ultimoExplorado]
							+ grafo[ultimoExplorado][vizinho].custo;
					antecessores[vizinho] = ultimoExplorado;
				}
			}

			++vizinho;
		}

		// verifica qual o próximo vértice a ser explorado
		int proximo = -1;
		double custoProximo = DBL_MAX;
		for (int i = 0; i < int(explorado.size()); ++i) {
			if (!explorado[i]) {
				if (custos[i] < custoProximo) {
					proximo = i;
					custoProximo = custos[i];
				}
			}
		}

		ultimoExplorado = proximo;

	} while (ultimoExplorado != -1 && ultimoExplorado != destino);

	// cria a sequencia para alcançar o destino a partir da origem
	int atual = destino;
	while (atual != -1) {
		solucao.push_front(atual);
		atual = antecessores[atual];
	}

	return solucao;
}

//===============================
// Split procedure
//===============================
vector<SolucaoVRP> split_tsp2vrp(DATA *ptr) {
	vector<SolucaoVRP> solucaoGeral;

	double custoMaximo = ptr->T;
	for (int i = 0; i < ptr->n; ++i) {
		vector<int> tsp;
		double custoTSP = 0;

		//salvando rotas TSP e os respectivos custos
		if (ptr->tour[i][0] != -1) {
			tsp.push_back(ptr->tour[i][0]);
			int j = 1;
			while (j < ptr-> n && ptr->tour[i][j] != -1) {
				tsp.push_back(ptr->tour[i][j]);
				custoTSP += ptr->c[tsp[tsp.size() - 2]][tsp[tsp.size() - 1]];
				++j;
			}
			if (tsp[0] != tsp[tsp.size() - 1]) {
				tsp.push_back(tsp[0]);
				custoTSP += ptr->c[tsp[tsp.size() - 1]][tsp[0]];
			}
		}

		if (custoTSP > custoMaximo) {
			int inicio = 1;
			int fim = inicio;
			Aresta arestaNula;
			arestaNula.custo = -1;

			vector<Aresta> vetorVazio(tsp.size() - 1, arestaNula);
			vector<vector<Aresta> > grafo(tsp.size() - 1, vetorVazio);

			//montando grafo em que cada aresta representa uma rota e seu custo
			while (inicio < int(tsp.size()) - 1) {
				Aresta aresta;
				aresta.rota.push_back(tsp[0]);
				aresta.custo = 0;
				for (int j = inicio; j <= fim; ++j) {
					aresta.rota.push_back(tsp[j]);
					aresta.custo
							+= ptr->c[aresta.rota[aresta.rota.size() - 2]][aresta.rota[aresta.rota.size()
									- 1]];
				}
				aresta.rota.push_back(tsp[0]);
				aresta.custo
						+= ptr->c[aresta.rota[aresta.rota.size() - 2]][aresta.rota[aresta.rota.size()
								- 1]];

				if (aresta.custo <= custoMaximo) {
					grafo[inicio - 1][fim] = aresta;
				}
				++fim;

				if (fim == int(tsp.size()) - 1) {
					++inicio;
					fim = inicio;
				}
			}

			//encontra menor caminho no grafo de rotas montado
			deque<int> solucaoDijkstra = menorCaminho(grafo, 0, tsp.size() - 2);

			SolucaoVRP vrp;
			for (int j = 0; j < int(solucaoDijkstra.size()) - 1; ++j) {
				Aresta a = grafo[solucaoDijkstra[j]][solucaoDijkstra[j + 1]];
				vrp.rotas.push_back(a.rota);
				vrp.custos.push_back(a.custo);
			}
			solucaoGeral.push_back(vrp);

		} else if (int(tsp.size()) > 0) {
			SolucaoVRP vrp;
			vrp.rotas.push_back(tsp);
			vrp.custos.push_back(custoTSP);
			solucaoGeral.push_back(vrp);
		}
	}

	// calculando custo total da solucao
	return solucaoGeral;
}

//====================================================
// Allocation heuristic - multipartida
//====================================================
void alocacao(DATA *ptr) {
	int *numero_hubs; //Número de hubs a serem definidos
	double *FO; //Função de custo
	double FOMax = 0;
	double FOMin = 1000000000;
	int bestpos = 0; //Salva posição do melhor individuo
	solucao *s; //struct de soluções contendo os individuos e os hubs
	hub *quem_e_hub; //struct contendo lista de quem são os hubs de cada
	ptr->xx = vector<int> (ptr->n);
	float alvo = 0;
	srand((unsigned) time(NULL));
	FO = new double[ptr->n];
	int *best = new int[ptr->n]; //armazena melhor solução

	//Inicializando População - baseado no GRASP
	s = inicia_populacao(alvo, FO, numero_hubs, quem_e_hub, ptr);
	Calcula_FO(FO, s, numero_hubs, quem_e_hub, FOMax, FOMin, bestpos, ptr);

	for (int i = 0; i < ptr->n; i++) {
		best[i] = s[bestpos].individuo[i];
	}

	//Aplicadno busca Local (VND) nas solucoes construidas
	for (int i = 0; i < ptr->n; i++) {
		s[i] = VND(alvo, FO, i, s[i], numero_hubs[i], quem_e_hub[i], FOMin,
				best, ptr);
		if (FO[i] < FO[bestpos]) {
			bestpos = i;
			FOMin = FO[bestpos];
		}
		//parada pelo alvo
		if (FOMin <= alvo + 0.1) {
			break;
		}
		ptr->fim = clock();
		ptr-> tempo = ((double) (ptr->fim - ptr->inicio) / CLOCKS_PER_SEC);
		//parada pelo tempo
		if (ptr->tempo > 3600) {
			break;
		}
	}

	for (int i = 0; i < ptr->n; i++) {
		ptr->xx[i] = s[bestpos].individuo[i];
		if (ptr->xx[i] == i) {
			ptr->hubs.push_back(i);
		}
	}
}

//====================================================
// VND - Variable Neighborhood Descent
//====================================================
solucao VND(double alvo, double *&FO, int pos, solucao &s, int &numero_hubs,
		hub &quem_e_hub, double &FOMin, int *&best, DATA *ptr) {
	int vizinhanca = 1;
	solucao *s_aux = new solucao[1];
	hub *quem_e_hub_aux = new hub[1];
	int *numero_hubs_aux = new int[1];
	s_aux[0] = cria_vetor_solucao(s_aux[0], ptr, 0, 0);
	double *OF_antes = new double[1];
	double *OF_depois = new double[1];
	for (int j = 0; j < ptr->n; j++) {
		s_aux[0].hubs[j] = s.hubs[j];
		s_aux[0].individuo[j] = s.individuo[j];
		if (s.hubs[j] == 1) {
			quem_e_hub_aux[0].hub_e.push_back(j);
		}
	}
	numero_hubs_aux[0] = numero_hubs;
	OF_antes[0] = FO[pos];
	OF_depois[0] = FO[pos];
	while (vizinhanca <= 4) { //explorando todas vizinhancas
		switch (vizinhanca) {
		case 1: //Busca local - vizinhanca 1
			s = busca_local(OF_depois, 0, s, numero_hubs, quem_e_hub, FOMin,
					alvo, ptr);
			if (OF_depois[0] < OF_antes[0]) {
				OF_antes[0] = OF_depois[0];
				if (OF_depois[0] < FOMin) {
					FOMin = OF_depois[0];
					for (int i = 0; i < ptr->n; i++) {
						best[i] = s.individuo[i];
					}
					if (FOMin <= alvo + 0.1) {
						break;
					}
				}
			}
			vizinhanca = 2;

		case 2: //Busca local - vizinhanca 2
			s = busca_local_insere_hub(OF_depois, 0, s, numero_hubs,
					quem_e_hub, FOMin, ptr);//
			if (OF_depois[0] < OF_antes[0]) {
				OF_antes[0] = OF_depois[0];
				if (OF_depois[0] < FOMin) {
					FOMin = OF_depois[0];
					for (int i = 0; i < ptr->n; i++) {
						best[i] = s.individuo[i];
					}
				}
				vizinhanca = 1;
			} else {
				vizinhanca++;
			}
		case 3://Busca local - vizinhanca 3
			s = busca_local_troca_funcao(OF_depois, 0, s, numero_hubs,
					quem_e_hub, FOMin, ptr); //
			if (OF_depois[0] < OF_antes[0]) {
				OF_antes[0] = OF_depois[0];
				if (OF_depois[0] < FOMin) {
					FOMin = OF_depois[0];
					for (int i = 0; i < ptr->n; i++) {
						best[i] = s.individuo[i];
					}
				}
				vizinhanca = 1;
			} else {
				vizinhanca++;
			}
		case 4://Busca local - vizinhanca 4
			if (numero_hubs > 1) {
				s = busca_local_remove_hub(OF_depois, 0, s, numero_hubs,
						quem_e_hub, FOMin, ptr);//
				if (OF_depois[0] < OF_antes[0]) {
					OF_antes[0] = OF_depois[0];
					if (OF_depois[0] < FOMin) {
						FOMin = OF_depois[0];
						for (int i = 0; i < ptr->n; i++) {
							best[i] = s.individuo[i];
						}
					}
					vizinhanca = 1;
				} else {
					vizinhanca++;
				}
			} else {
				vizinhanca++;
			}
		}
		if (FOMin <= alvo + 0.1) {
			break;
		}
	}

	FO[pos] = OF_depois[0];
	delete[] s_aux;
	delete[] quem_e_hub_aux;
	delete[] numero_hubs_aux;
	return s;
}

//====================================================
// Busca Local - insere hub (vizinhança 4)
//====================================================
solucao busca_local_insere_hub(double *&FO, int pos, solucao &s,
		int &numero_hubs, hub &quem_e_hub, double &FOMin, DATA *ptr) {
	solucao *s_aux = new solucao[1];
	hub *quem_e_hub_aux = new hub[1];
	int *numero_hubs_aux = new int[1];
	s_aux[0] = cria_vetor_solucao(s_aux[0], ptr, 0, 0);
	double *OF_antes = new double[1];
	double *OF_depois = new double[1];
	for (int j = 0; j < ptr->n; j++) {
		s_aux[0].hubs[j] = s.hubs[j];
		s_aux[0].individuo[j] = s.individuo[j];
		if (s.hubs[j] == 1) {
			quem_e_hub_aux[0].hub_e.push_back(j);
		}
	}
	numero_hubs_aux[0] = numero_hubs;
	OF_antes[0] = FO[pos];
	OF_depois[0] = FO[pos];

	for (int n = 0; n < ptr->n; n++) {
		if (s_aux[0].hubs[n] == 0) {
			s_aux[0].hubs[n] = 1;
			s_aux[0].individuo[n] = n;
			numero_hubs_aux[0]++;
			std::list<int>::iterator insert_hub;
			insert_hub = quem_e_hub_aux[0].hub_e.begin();
			quem_e_hub_aux[0].hub_e.insert(insert_hub, n);
			for (int j = 0; j < ptr->n; j++) {
				double dist1 = ptr->c[j][s_aux[0].individuo[j]];
				if (dist1 > ptr->c[j][n]) {
					s_aux[0].individuo[j] = n;
				}
			}

			OF_depois[0] = calcula_fo_total(s_aux[0], ptr);

			if (OF_depois[0] < OF_antes[0]) {
				OF_antes[0] = OF_depois[0];
				FO[pos] = OF_depois[0];
				quem_e_hub.hub_e.clear();
				numero_hubs = 0;
				for (int m = 0; m < ptr->n; m++) {
					s.hubs[m] = s_aux[0].hubs[m];
					s.individuo[m] = s_aux[0].individuo[m];
					if (s_aux[0].hubs[m] == 1) {
						quem_e_hub.hub_e.push_back(m);
						numero_hubs++;
					}
				}
				break;
			} else {
				quem_e_hub_aux[0].hub_e.clear();
				numero_hubs_aux[0] = 0;
				for (int m = 0; m < ptr->n; m++) {
					s_aux[0].hubs[m] = s.hubs[m];
					s_aux[0].individuo[m] = s.individuo[m];
					if (s.hubs[m] == 1) {
						quem_e_hub_aux[0].hub_e.push_back(m);
						numero_hubs_aux[0]++;
					}
				}
			}
		}
	}
	delete[] s_aux;
	delete[] quem_e_hub_aux;
	delete[] numero_hubs_aux;
	return s;
}

//====================================================
// Busca Local - remove hub (vizinhança 3)
//====================================================
solucao busca_local_remove_hub(double *&FO, int pos, solucao &s,
		int &numero_hubs, hub &quem_e_hub, double &FOMin, DATA *ptr) {
	solucao *s_aux = new solucao[1];
	hub *quem_e_hub_aux = new hub[1];
	int *numero_hubs_aux = new int[1];
	s_aux[0] = cria_vetor_solucao(s_aux[0], ptr, 0, 0);
	double *OF_antes = new double[1];
	double *OF_depois = new double[1];
	for (int j = 0; j < ptr->n; j++) {
		s_aux[0].hubs[j] = s.hubs[j];
		s_aux[0].individuo[j] = s.individuo[j];
		if (s.hubs[j] == 1) {
			quem_e_hub_aux[0].hub_e.push_back(j);
		}
	}
	numero_hubs_aux[0] = numero_hubs;
	OF_antes[0] = FO[pos];
	OF_depois[0] = FO[pos];

	for (int n = 0; n < ptr->n; n++) {
		if (s_aux[0].hubs[n] == 1) {
			quem_e_hub_aux[0].hub_e.remove(n);
			s_aux[0].hubs[n] = 0;
			numero_hubs_aux[0]--;
			for (int j = 0; j < ptr->n; j++) {
				double dist1 = 1000000;
				if (s_aux[0].individuo[j] == n) {
					std::list<int>::iterator pont_mut;
					for (pont_mut = quem_e_hub_aux[0].hub_e.begin(); pont_mut
							!= quem_e_hub_aux[0].hub_e.end(); pont_mut++) {
						if (dist1 > ptr->c[j][*pont_mut]) {
							dist1 = ptr->c[j][*pont_mut];
							s_aux[0].individuo[j] = *pont_mut;
						}
					}
				}
			}
			OF_depois[0] = calcula_fo_total(s_aux[0], ptr);

			if (OF_depois[0] < OF_antes[0]) {
				OF_antes[0] = OF_depois[0];
				FO[pos] = OF_depois[0];
				quem_e_hub.hub_e.clear();
				numero_hubs = 0;
				for (int m = 0; m < ptr->n; m++) {
					s.hubs[m] = s_aux[0].hubs[m];
					s.individuo[m] = s_aux[0].individuo[m];
					if (s_aux[0].hubs[m] == 1) {
						quem_e_hub.hub_e.push_back(m);
						numero_hubs++;
					}
				}
				break;
			} else {
				quem_e_hub_aux[0].hub_e.clear();
				numero_hubs_aux[0] = 0;
				for (int m = 0; m < ptr->n; m++) {
					s_aux[0].hubs[m] = s.hubs[m];
					s_aux[0].individuo[m] = s.individuo[m];
					if (s.hubs[m] == 1) {
						quem_e_hub_aux[0].hub_e.push_back(m);
						numero_hubs_aux[0]++;
					}
				}
			}
		}
	}
	delete[] s_aux;
	delete[] quem_e_hub_aux;
	delete[] numero_hubs_aux;
	return s;
}

//====================================================
// Busca Local - troca função (vizinhança 2)
//====================================================
solucao busca_local_troca_funcao(double *&FO, int pos, solucao &s,
		int numero_hubs, hub &quem_e_hub, double &FOMin, DATA *ptr) {
	solucao *s_aux = new solucao[1];
	hub *quem_e_hub_aux = new hub[1];
	int *numero_hubs_aux = new int[1];
	s_aux[0] = cria_vetor_solucao(s_aux[0], ptr, 0, 0);
	double *OF_antes = new double[1];
	double *OF_depois = new double[1];
	for (int j = 0; j < ptr->n; j++) {
		s_aux[0].hubs[j] = s.hubs[j];
		s_aux[0].individuo[j] = s.individuo[j];
		if (s.hubs[j] == 1) {
			quem_e_hub_aux[0].hub_e.push_back(j);
		}
	}
	numero_hubs_aux[0] = numero_hubs;
	OF_antes[0] = FO[pos];
	OF_depois[0] = FO[pos];
	for (int j = 0; j < ptr->n; j++) {
		if (s_aux[0].hubs[j] == 1) {
			for (int n = 0; n < ptr->n; n++) {
				if ((s_aux[0].individuo[n] == j) & (s_aux[0].hubs[n] == 0)) {
					quem_e_hub_aux[0].hub_e.remove(j);
					s_aux[0].hubs[j] = 0;
					s_aux[0].individuo[j] = n;
					s_aux[0].hubs[n] = 1;
					s_aux[0].individuo[n] = n;
					std::list<int>::iterator insert_hub;
					insert_hub = quem_e_hub_aux[0].hub_e.begin();
					quem_e_hub_aux[0].hub_e.insert(insert_hub, n);
					associar_spoke_a_hub(s_aux, numero_hubs_aux,
							quem_e_hub_aux, 0, ptr);
					OF_depois[0] = calcula_fo_total(s_aux[0], ptr);
					if (OF_depois[0] < OF_antes[0]) {
						OF_antes[0] = OF_depois[0];
						FO[pos] = OF_depois[0];
						quem_e_hub.hub_e.clear();
						for (int m = 0; m < ptr->n; m++) {
							s.hubs[m] = s_aux[0].hubs[m];
							s.individuo[m] = s_aux[0].individuo[m];
							if (s_aux[0].hubs[m] == 1) {
								quem_e_hub.hub_e.push_back(m);
							}
						}
						break;
					} else {
						quem_e_hub_aux[0].hub_e.clear();
						for (int m = 0; m < ptr->n; m++) {
							s_aux[0].hubs[m] = s.hubs[m];
							s_aux[0].individuo[m] = s.individuo[m];
							if (s.hubs[m] == 1) {
								quem_e_hub_aux[0].hub_e.push_back(m);
							}
						}
					}
				}
			}
		}
	}
	delete[] s_aux;
	delete[] quem_e_hub_aux;
	delete[] numero_hubs_aux;
	return s;
}

//====================================================
// Busca Local - Shift (vizinhança 1)
//====================================================
solucao busca_local(double *&FO, int pos, solucao &s, int numero_hubs,
		hub quem_e_hub, double &FOMin, double alvo, DATA *ptr) {
	solucao *s_aux = new solucao[1];
	s_aux[0] = cria_vetor_solucao(s_aux[0], ptr, 0, 0);
	double OF_antes = 0;
	double OF_depois = 0;
	int hub_solucao;
	for (int j = 0; j < ptr->n; j++) {
		s_aux[0].hubs[j] = s.hubs[j];
		s_aux[0].individuo[j] = s.individuo[j];
	}
	OF_antes = FO[pos];
	for (int j = 0; j < ptr->n; j++) {
		if (s_aux[0].hubs[j] == 0) {
			hub_solucao = s_aux[0].individuo[j];
			std::list<int>::iterator acesso;
			for (acesso = quem_e_hub.hub_e.begin(); acesso
					!= quem_e_hub.hub_e.end(); acesso++) {
				if (*acesso != hub_solucao) {
					s_aux[0].individuo[j] = *acesso;
					OF_depois = calcula_fo_total(s_aux[0], ptr);
					if (OF_depois < OF_antes) {
						OF_antes = OF_depois;
						FO[pos] = OF_depois;
						s.individuo[j] = *acesso;
						hub_solucao = *acesso;
						if (FO[pos] <= alvo + 0.1) {
							break;
						}
					}
					s_aux[0].individuo[j] = hub_solucao;
				}
			}
		}
		if (FO[pos] <= alvo + 0.1) {
			break;
		}
	}
	delete[] s_aux;
	return s;
}

//====================================================
// Geração da Populaçõa Inicial - GRASP
//====================================================
solucao *inicia_populacao(double alvo, double *&OF, int *&numero_hubs,
		hub *&quem_e_hub, DATA *ptr) {
	solucao *s = new solucao[ptr->n]; //struct solucao
	solucao *temporario = new solucao[1]; //struct solucao temporario
	temporario = alocar_espaco(ptr);
	numero_hubs = new int[ptr->n];
	quem_e_hub = new hub[ptr->n];
	int *numero_hubs_temp;
	numero_hubs_temp = new int[1];
	hub *quem_e_hub_temp;
	quem_e_hub_temp = new hub[1];
	double OF_aux = 0;
	double *custo_marginal; //custo marginal (FOatual - FO)
	custo_marginal = new double[ptr->n];
	double custo_marginal_max;
	double custo_marginal_min;
	int *lista_candidatos; //candidatos a se tornarem hubs
	lista_candidatos = new int[ptr->n + 1];
	int *lista_nao_hubs;
	lista_nao_hubs = new int[ptr->n]; // nós que não serão hubs
	double lambda = 0.05; //Parametro de controle da aleatoriedade
	int i = 0;

	//Criando soluções iniciais
	while (i < ptr->n) {
		for (int j = 0; j < ptr->n; j++) {
			s[i] = cria_vetor_solucao(s[i], ptr, 0, j);
			s[i].hubs[j] = 1;
			quem_e_hub[i].hub_e.push_back(j);
			numero_hubs[i] = 1;
			OF[i] = calcula_fo_total(s[i], ptr);
			i++;
			if (i >= ptr->n) {
				break;
			}
		}
	}

	//Improvement Step - Adiciona hubs e cria uma lista restrita de candidatos. Escolhendo um dos candidatos que possuem custo marginal "negativo". Repete até lista vazia
	for (int i = 0; i < ptr->n; i++) {
		lista_candidatos[0] = ptr->n;
		for (int j = 0; j < ptr->n; j++) {
			lista_nao_hubs[j] = 0;
		}
		while (lista_candidatos[0] != 0) {
			custo_marginal_max = -10000000;
			custo_marginal_min = 0;
			for (int j = 0; j < ptr->n; j++) {
				custo_marginal[j] = 0;
				if ((s[i].hubs[j] == 0) & (lista_nao_hubs[j] == 0)) {
					//Copiar solucao para temp
					for (int n = 0; n < ptr->n; n++) {
						temporario[0].hubs[n] = s[i].hubs[n];
						temporario[0].individuo[n] = s[i].individuo[n];
					}
					quem_e_hub_temp[0].hub_e.clear();
					std::list<int>::iterator acesso;
					for (acesso = quem_e_hub[i].hub_e.begin(); acesso
							!= quem_e_hub[i].hub_e.end(); acesso++) {
						quem_e_hub_temp[0].hub_e.push_back(*acesso);
					}
					numero_hubs_temp[0] = numero_hubs[i];

					//Setar nó j como hub
					temporario[0].hubs[j] = 1;
					temporario[0].individuo[j] = j;
					quem_e_hub_temp[0].hub_e.push_back(j);
					numero_hubs_temp[0] = numero_hubs_temp[0] + 1;

					//Reassociar spokes a hub mais próximo, calcular OF e custo marginal
					associar_spoke_a_hub(temporario, numero_hubs_temp,
							quem_e_hub_temp, 0, ptr);
					OF_aux = calcula_fo_total(temporario[0], ptr);
					custo_marginal[j] = OF_aux - OF[i];

					//Pegando custo margunal max e min
					if (custo_marginal[j] < 0) {
						lista_nao_hubs[j] = 0;
						if (custo_marginal[j] < custo_marginal_min) {
							custo_marginal_min = custo_marginal[j];
						}
						if (custo_marginal[j] > custo_marginal_max) {
							custo_marginal_max = custo_marginal[j];
						}
					} else {
						lista_nao_hubs[j] = -1; //caso custo marginal maior que 0, entao vai pra lista de não hubs
					}
				}
			}

			// Preenchendo lista de candidatos
			lista_candidatos[0] = 0;
			if (custo_marginal_min < 0) {
				for (int j = 0; j < ptr->n; j++) {
					if ((custo_marginal[j] < 0)
							& (custo_marginal[j] <= custo_marginal_min + lambda
									* (custo_marginal_max - custo_marginal_min))) {
						lista_candidatos[0] = lista_candidatos[0] + 1;
						lista_candidatos[lista_candidatos[0]] = j;
					}
				}
			}
			//Escolhendo um candidato para se tornar hub
			if (lista_candidatos[0] > 0) {
				int escolhido = 1 + rand() % lista_candidatos[0];

				//setar escolhido como hub, alocar nós ao hub mais proximo e calcular FO
				s[i].hubs[lista_candidatos[escolhido]] = 1;
				s[i].individuo[lista_candidatos[escolhido]]
						= lista_candidatos[escolhido];
				numero_hubs[i] = numero_hubs[i] + 1;
				std::list<int>::iterator insert_hub;
				insert_hub = quem_e_hub[i].hub_e.begin();
				quem_e_hub[i].hub_e.insert(insert_hub,
						lista_candidatos[escolhido]);
				associar_spoke_a_hub(s, numero_hubs, quem_e_hub, i, ptr);
				OF[i] = calcula_fo_total(s[i], ptr);
			}
		}
		if (OF[i] <= alvo) {
			break;
		}
	}
	delete[] temporario;
	delete[] custo_marginal;
	delete[] lista_candidatos;
	delete[] numero_hubs_temp;
	delete[] quem_e_hub_temp;
	return s;
}

//==============================================================================
//Funções para o calculo do custo da FO
//==============================================================================
//Custo do transporte entre i e j, ida e volta
double custo_con_ij(int i, int j, solucao sol, DATA *ptr) {
	return ptr->w[i][j] * (ptr->c[i][sol.individuo[i]] + 10000 * max(0.0, 2
			* ptr->c[i][sol.individuo[i]] - ptr->T) + ptr->alpha
			* ptr->c[sol.individuo[i]][sol.individuo[j]]
			+ ptr->c[sol.individuo[j]][j] + 10000 * max(0.0, 2
			* ptr->c[sol.individuo[j]][j] - ptr->T)) + ptr->w[j][i]
			* (ptr->c[j][sol.individuo[j]] + 10000 * max(0.0, 2
					* ptr->c[j][sol.individuo[j]] - ptr->T) + ptr->alpha
					* ptr->c[sol.individuo[j]][sol.individuo[i]]
					+ ptr->c[sol.individuo[i]][i] + 10000 * max(0.0, 2
					* ptr->c[sol.individuo[i]][i] - ptr->T));
}

//Custo Fixo Total
double custo_fixo(solucao sol, DATA *ptr) {
	double f = 0;
	for (int i = 0; i < ptr->n; i++) {
		f = f + ptr->ex * sol.hubs[i] * ptr->a[i];
	}
	return f;
}

//Calcula a Fo de uma solução
double calcula_fo_total(solucao sol, DATA *ptr) {
	double fo;
	fo = custo_fixo(sol, ptr);
	for (int i = 0; i < ptr->n; i++) {
		for (int j = i + 1; j < ptr->n; j++) {
			fo += custo_con_ij(i, j, sol, ptr);
		}
	}
	return fo;
}

//====================================================
// Função para calculo da FO de toda população
//====================================================
void Calcula_FO(double *&FO, solucao *s, int *numero_hubs, hub *quem_e_hub,
		double &FOMax, double &FOMin, int &bestpos, DATA *ptr) {
	FOMax = 0;
	for (int i = 0; i < ptr->n; i++) {
		FO[i] = 0;
		std::list<int>::iterator ponteiro;
		for (ponteiro = quem_e_hub[i].hub_e.begin(); ponteiro
				!= quem_e_hub[i].hub_e.end(); ponteiro++) {
			FO[i] = FO[i] + ptr->ex * ptr->a[*ponteiro];
		}
		for (int j = 0; j < ptr->n; j++) {
			FO[i] = FO[i] + ptr->c[j][s[i].individuo[j]] * (ptr->o[j]
					+ ptr->d[j]) + 10000 * max(0.0, 2
					* ptr->c[j][s[i].individuo[j]] - ptr->T);
		}
		for (int j = 0; j < ptr->n; j++) {
			for (int w = j + 1; w < ptr->n; w++) {
				if (s[i].individuo[j] != s[i].individuo[w])
					FO[i]
							= FO[i]
									+ (ptr->w[j][w] * ptr->alpha
											* ptr->c[s[i].individuo[j]][s[i].individuo[w]])
									+ (ptr->w[w][j] * ptr->alpha
											* ptr->c[s[i].individuo[w]][s[i].individuo[j]]);
			}
		}
		if (FO[i] < FOMin) {
			FOMin = FO[i];
			bestpos = i;
		}
		if (FO[i] > FOMax) {
			FOMax = FO[i];
		}
	}
}

//====================================================
// Associação dos spokes ao hub mais próximo
//====================================================
void associar_spoke_a_hub(solucao *&s, int *numero_hubs, hub *quem_e_hub,
		int i, DATA *ptr) {
	if (numero_hubs[i] == 1) {
		std::list<int>::iterator ponteiro;
		ponteiro = quem_e_hub[i].hub_e.begin();
		for (int j = 0; j < ptr->n; j++) {
			s[i].individuo[j] = *ponteiro;
		}
	} else { //testando quais hubs são mais próximos dos nós
		for (int w = 0; w < ptr->n; w++) {
			double dist = 1000000;
			if (s[i].hubs[w] != 1) {
				std::list<int>::iterator ponteiro;
				ponteiro = quem_e_hub[i].hub_e.begin();
				for (int j = 0; j < numero_hubs[i]; j++) {
					if (dist > ptr->c[w][*ponteiro]) {
						dist = ptr->c[w][*ponteiro];
						s[i].individuo[w] = *ponteiro;
					}
					ponteiro++;
				}
			}
		}
	}
}

//====================================================
// Crização de várias structs de solução
//====================================================
solucao *alocar_espaco(DATA *ptr) {
	solucao *s = new solucao[ptr->n];
	for (int i = 0; i < ptr->n; i++) {
		s[i] = cria_vetor_solucao(s[i], ptr, 0, 0);
	}
	return s;
}

//====================================================
// Criação para struct solução
//====================================================
solucao cria_vetor_solucao(solucao vetor, DATA *ptr, int p, int q) {
	try {
		vetor.hubs = new int[ptr->n];
		vetor.individuo = new int[ptr->n];
		for (int i = 0; i < ptr->n; i++) {
			vetor.hubs[i] = p;
			vetor.individuo[i] = q;
		}
	} catch (std::bad_alloc) { // ENTRA AQUI SE NÃO ALOCAR MEMÓRIA PARA A MATRIZ.
		cout << "PROBLEMAS DE MEMÓRIA. Bye ...";
		//system("PAUSE");
		exit(-1);
	}
	return vetor;
}

// ================================================================
// reading procedure
// ================================================================
void read_data(DATA *ptr) {
	double daux;
	ifstream in(ptr->file_name);
	if (!in.is_open()) {
		help(ptr);
		exit(EXIT_FAILURE);
	}
	in >> ptr->n;
	in >> ptr->nv;
	in >> ptr->cq;
	for (int i = 0; i < ptr->n; i++) {
		in >> daux;
		ptr->a.push_back(daux / (ADJUSTMENT * ADJUSTMENT));
	}

	for (int i = 0; i < ptr->n; i++) {
		vector<double> vaux;
		for (int j = 0; j < ptr->n; j++) {
			in >> daux;
			vaux.push_back(daux / ADJUSTMENT);
		}
		ptr->w.push_back(vaux);
	}
	for (int i = 0; i < ptr->n; i++) {
		daux = 0.0;
		for (int j = 0; j < ptr->n; j++) {
			daux += ptr->w[i][j];
		}
		ptr->o.push_back(daux);
		daux = 0.0;
		for (int j = 0; j < ptr->n; j++) {
			daux += ptr->w[j][i];
		}
		ptr->d.push_back(daux);
	}

	for (int i = 0; i < ptr->n; i++) {
		vector<double> vaux;
		for (int j = 0; j < ptr->n; j++) {
			in >> daux;
			vaux.push_back(daux / ADJUSTMENT);
		}
		ptr->c.push_back(vaux);
	}
	for (int i = 0; i < ptr->n; i++) {
		in >> daux;
		ptr->xcoord.push_back(daux);
		in >> daux;
		ptr->ycoord.push_back(daux);
	}
	in.close();
}

// ================================================================
//
// ================================================================
void help(DATA *ptr) {
	cout << "\n\n" << ptr->exec_name
			<< " [data file name] [alpha - discount factor] [ex - fixed factor] [beta - allocation factor] [omega - vehicle factor] [gamma - arc factor] [T - tour time limit] [nhot - number of hot start cycles]"
			<< "\n\n" << endl;
}
// ================================================================
//
// ================================================================
void print_data(DATA *ptr) {
	printf("n : %d\n", ptr->n);
	printf("fixed costs  o   d :\n");
	for (int i = 0; i < ptr->n; i++) {
		printf("%3d  %18.2f %18.2f %18.2f\n", i, ptr->ex * ptr->a[i],
				ptr->o[i], ptr->d[i]);
	}
	printf("demands :\n");
	for (int i = 0; i < ptr->n; i++) {
		for (int j = 0; j < ptr->n; j++) {
			printf("%12.2f ", ptr->w[i][j]);
		}
		printf("\n");
	}
	printf("costs :\n");
	for (int i = 0; i < ptr->n; i++) {
		for (int j = 0; j < ptr->n; j++) {
			printf("%12.2f ", ptr->c[i][j]);
		}
		printf("\n");
	}
	printf("(x,y) coordenates :\n");
	for (int i = 0; i < ptr->n; i++) {
		printf("%3d %12.2f %12.2f\n", i, ptr->xcoord[i], ptr->ycoord[i]);
	}
}
