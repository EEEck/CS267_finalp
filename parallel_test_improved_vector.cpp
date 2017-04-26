#include <iostream>
//#define ARMA_DONT_USE_WRAPPER
#include <armadillo>
#include <omp.h>
#include<cmath>
//#include <chrono>
using namespace std;
using namespace arma;
//using namespace std::chrono;
//Blas level 3 2 1
int main(int argc, char** argv)
{
	//Initi vars
	int nbasis =800; //632;
	int nocc= 73;
	int nvirt =nbasis- nocc ;
	int naux= 4*nbasis;
	int nbatch=5;
	vec orben = vec(nvirt+nocc, fill::ones); orben.randu(nvirt+nocc);
	cube B_iaq(nvirt, naux, nocc, fill::randu); B_iaq.randn(nvirt, naux, nocc);
	double emp2=0.0;


//goal batches insted of just a slice
//
//
// care only one nested loop and get in bact and unique pairs
	double t1 = omp_get_wtime();
	#pragma omp parallel reduction(+:emp2)
	{
		int id = omp_get_thread_num();
		if(id == 0)
		{
			int max_thrds = omp_get_max_threads();
			cout <<"Max Threads used:" << max_thrds << endl;
		}
		double mp2part_per_thr = 0.0;
		#pragma omp for schedule(dynamic) 
		for (int i = 0; i < nocc; ++i)
		{

			//cube B_iaq_batch(B_iaq_batch.slice(i).memptr(), nvirt, naux,nbatch, false, false);
			mat B_i(B_iaq.slice(i).memptr(), nvirt, naux, false, false);

			for (int j = i+1; j < nocc; ++j)
			{
				mat B_j(B_iaq.slice(j).memptr(), nvirt, naux, false, false);
				double nume;
				double d_ijab;
				mat W_ab=B_i*B_j.t();
				for (int a = 0; a < nvirt; ++a)
				{
					const double* wa = W_ab.memptr();
					double dija = orben[i] + orben[j] - orben[a];
					#pragma omp simd reduction(+:mp2part_per_thr)
					for (int b = a+1; b < nvirt; ++b)
					{
						//double d_ijab= dija - orben[b];
						double wijab = wa[b] / (dija - orben[b]);
						//double nume = abs(W_ab(a,b) - W_ab(b,a));
						//mp2part_per_thr+=nume*nume/d_ijab;
						mp2part_per_thr+= wijab * (wa[b] - W_ab[a*nvirt+b]);
						//mp2part_per_thr+=(W_ab(a,b) - W_ab(b,a))*(W_ab(a,b) - W_ab(b,a))/(d_ijab=(orben(i)+orben(j)-orben(b+nocc)-orben(a+nocc)));
					}
				}
			}
		}
			emp2+=mp2part_per_thr;
	}

	double t2 = omp_get_wtime();
	cout << "P time:"<< t2-t1<< endl;
	double zahler = 0.0;
	double nenner = 0.0;
	double mp22= 0.0;
	double t3 = omp_get_wtime();
	/*	
	for (int i = 0; i < nocc; ++i)
	{
		for (int j = 0; j < nocc; ++j)
		{
			for (int a = 0; a < nvirt; ++a)
			{
				for (int b = 0; b < nvirt; ++b)
				{
					zahler = 0.0;
					for (int q = 0; q < naux; ++q) {
						zahler += B_iaq(a,q,i)*B_iaq(b,q,j)-B_iaq(b,q,i)*B_iaq(a,q,j);
					}
					nenner = orben(i) + orben(j) - orben(a+nocc) - orben(b+nocc);
					mp22 += zahler*zahler/nenner; 
				}
			}
		}
	}

	mp22 *= 0.25;
	*/
	double t4 = omp_get_wtime();
	cout << "P time naive:"<< t4-t3<< endl;
	cout << "Compare please: "<< emp2 << "\t"<< mp22 << endl;
	//cube B_iaq=randu
	//n_aux=B_iaQ.n_slices
	//high_resolution_clock::time_point t1 = high_resolution_clock::now();
	//eig_sym(eigval, eigvec, A);
	//high_resolution_clock::time_point t2 = high_resolution_clock::now();
	//auto duration = duration_cast<seconds>( t2 - t1 ).count();
	//cout << duration<< endl;
	//cout <<size(A) << endl;

	return 0; }



