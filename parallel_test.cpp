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
double partial_mp2(const mat& B_i, const mat& B_j,vec& Orb_en, const int& ei , const int& ej){ // int% means passing by reference!!!
   		mat W_ab=B_i*B_j.t();
   		double emp2part =0.0;
   		double d_ijab = 0.0;
		size_t norb = Orb_en.n_elem;
   		size_t nvir = W_ab.n_rows;
        	size_t nocc = norb - nvir;
		double nume;
   		for (int a = 0; a < nvir; ++a)
   		{
   			for (int b = a+1; b < nvir; ++b)
   			{
   				d_ijab=(Orb_en(ei)+Orb_en(ej)-Orb_en(b+nocc)-Orb_en(a+nocc));
   				// emp2part+=W_ab(a,b)*((W_ab(a,b)) -W_ab(b,a))/d_ijab;
   				nume = abs(W_ab(a,b) - W_ab(b,a));
   				emp2part+=nume*nume/d_ijab;
   				// std::cout << "W_ab(a,b)*((W_ab(a,b)) -W_ab(b,a))/d_ijab = " << W_ab(a,b)*((W_ab(a,b)) -W_ab(b,a))/d_ijab << std::endl;
   			}
   		}
   		return emp2part;
}
int main(int argc, char** argv)
{
	//Initi vars
	int nbasis = 800;
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
	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		if(id == 0)
		{
			int max_thrds = omp_get_max_threads();
			cout <<"Max Threads used:" << max_thrds << endl;
		}
		double mp2part_per_thr = 0.0;
		#pragma omp for
		for (int i = 0; i < nocc; ++i)
		{

			//cube B_iaq_batch(B_iaq_batch.slice(i).memptr(), nvirt, naux,nbatch, false, false);
			mat B_i(B_iaq.slice(i).memptr(), nvirt, naux, false, false);

			for (int j = i+1; j < nocc; ++j)
			{
				mat B_j(B_iaq.slice(j).memptr(), nvirt, naux, false, false);
				mp2part_per_thr+=partial_mp2(B_i,B_j,orben,i,j);
			}
		}
		#pragma omp critical //change to atomic for better performance!
		{
			emp2+=mp2part_per_thr;
		} 
	}

	double t2 = omp_get_wtime();
	cout << "P time:"<< t2-t1<< endl;
	double zahler = 0.0;
	double nenner = 0.0;
	double mp22= 0.0;
	double t3 = omp_get_wtime();
	/*for (int i = 0; i < nocc; ++i)
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
	*/

	mp22 *= 0.25;

	double t4 = omp_get_wtime();
	//cout << "P time naive:"<< t4-t3<< endl;
	//cout << "Compare please: "<< emp2 << "\t"<< mp22 << endl;
	//cube B_iaq=randu
	//n_aux=B_iaQ.n_slices
	//high_resolution_clock::time_point t1 = high_resolution_clock::now();
	//eig_sym(eigval, eigvec, A);
	//high_resolution_clock::time_point t2 = high_resolution_clock::now();
	//auto duration = duration_cast<seconds>( t2 - t1 ).count();
	//cout << duration<< endl;
	//cout <<size(A) << endl;

	return 0; }



