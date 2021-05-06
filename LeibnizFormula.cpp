#include <bits/stdc++.h>

#include <time.h>
#include <mpi.h>

using namespace std;
int ProcNum = 0; // Number of available processes
int ProcRank = 0; // Rank of current process

void ProcessInitialization (double* &pVector, int &N, int &blockSize, double &res) {
    if (ProcRank == 0) {
        while(N < ProcNum){
            printf("\nEnter size of the initial objects: ");
            scanf("%d", &N);
            if (N < ProcNum) {
                printf("Size of the objects must be greater than number of processes! \n ");
            }
        }
    }
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    blockSize = N/ProcNum;
    MPI_Bcast(&blockSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&res, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

// Result vector replication
void ReseiveResult(double& res) {
        double tmp;
        MPI_Recv(&tmp, 1, MPI_DOUBLE, ProcRank, 0, MPI_COMM_WORLD,  MPI_STATUS_IGNORE);
        res += tmp;
}


// Function for sequential matrix-vector multiplication
void SerialResultCalculation(double* pMatrix, double* pVector, double* pResult, int Size) {
    int i, j; // Loop variables
    for (i=0; i<Size; i++) {
        pResult[i] = 0;
        for (j=0; j<Size; j++)
            pResult[i] += pMatrix[i*Size+j]*pVector[j];
    }
}
// Function for calculating partial matrix-vector multiplication
void ParallelResultCalculation(double* pVector, int N, int blockSize) {
    double i; // Loop variables
    double currRes = 0;
    double l = ProcRank*blockSize;
    double r = ProcRank == ProcNum - 1 ? (ProcRank + 1) * blockSize : N;
    for (i = l; i < r; i += 1) {
        currRes = currRes + (double) (1 / (2 * i + 1)) * double(((int) i % 2) ? 1 : -1);
    }
    MPI_Send(&currRes, 1, MPI_DOUBLE, ProcRank, 0, MPI_COMM_WORLD);
  //  printf("send rank %d \n", ProcRank);
   // pVector[ProcRank] = currRes;
}

void consistentResultCalculation(int n, double parallelRes) {
    double res = 0;
    for(double i = 0; i < n; i = i + 1){
        res = res + double((double) (1 / (2 * i + 1)) * double(((int) i % 2) ? 1 : -1));
    }
    if(abs(res - parallelRes) > 10e-6){
        cout << "Assertion failed! " << res << " " << parallelRes << "\n";
    }

}

// Function for computational process termination
void ProcessTermination (double* pMatrix, double* pVector, double* pResult,
                         double* pProcRows, double* pProcResult) {
    if (ProcRank == 0)
        delete [] pMatrix;
    delete [] pVector;
    delete [] pResult;
    delete [] pProcRows;
    delete [] pProcResult;
}
int main(int argc, char* argv[]) {
    double* pVector; // Vector with results of each block
    int N;
    int blockSize;
    double Start, Finish, Duration;
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    double res;
    ProcessInitialization(pVector,  N, blockSize, res);
    Start = MPI_Wtime();
    ParallelResultCalculation(pVector, N, blockSize);
    ReseiveResult(res);
    Finish = MPI_Wtime();
    Duration = Finish-Start;
   // double res = 0;
    MPI_Status status;
    if (ProcRank == 0) {
        cout << "res = " << res << "\n";
        printf("Time of consistent execution = %.20f\n", Duration);

        Start = MPI_Wtime();
        consistentResultCalculation(N, res);
        Finish = MPI_Wtime();

        printf("Time of parallel execution = %.20f\n", Finish - Start);
    }


    MPI_Finalize();
}

