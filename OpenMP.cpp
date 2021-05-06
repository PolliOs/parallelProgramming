#include <iostream>
#include <omp.h>
#include <cmath>
#include <chrono>
using namespace std::chrono;

int THREAD_NUM = 4;

using namespace std;

double LeibnizFormula(int l, int r) {
    double sum = 0;
    for (double i = l; i < r; i+=1){
        sum += (double) (1 / (2 * i + 1)) * double(((int) i % 2) ? 1 : -1);
    }
    return sum;
}

int main(int argc, char** argv) {
    int n;
    cout << "Enter num of threads:";
    cin >> THREAD_NUM;
    cout << "Enter n:";
    cin >> n;
    omp_set_num_threads(THREAD_NUM);
    int num_per_thread = ceil(n / THREAD_NUM);
    double res = 0;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
//#pragma omp parallel for
//    for ( int k=0; k<THREAD_NUM; k++) {
//        tempRes = LeibnizFormula(k*num_per_thread, min(n, (k+1)*num_per_thread));
//        std::cout << "tempRes = " << tempRes << "k = " << k << "\n";
////        res =  res + tempRes;
//#pragma omp critical
//        std::cout << "tempRes = " << tempRes << "\n";
//        res =  res + tempRes;
   // }

#pragma omp parallel
    {
        double tempRes = 0;
#pragma omp for
        for (int k = 0; k < THREAD_NUM; k++)
            tempRes += LeibnizFormula(k*num_per_thread, min(n, (k+1)*num_per_thread));
#pragma omp atomic
        res += tempRes;
       // std::cout << "tempRes = " << tempRes << "\n";
    }

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

    double actualRes = 0;
    for (double i = 0; i < n; i+=1){
       // double temp = (double) (1 / (2 * i + 1)) * double(((int) i % 2) ? 1 : -1);
       // std::cout << "temp = " << temp << "\n";
        actualRes = actualRes + (double) (1 / (2 * i + 1)) * double(((int) i % 2) ? 1 : -1);
    }

    std::cout << "Time of execution " << time_span.count() << " seconds.\n";
    printf("res = %.20f\n", res);
    printf("actualRes = %.20f\n", actualRes);
    return 0;
}
