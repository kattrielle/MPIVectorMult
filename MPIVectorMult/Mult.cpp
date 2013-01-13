/*
 ============================================================================
 Name        : lab2-1.c
 Author      : Elessar
 Version     :
 Copyright   : Your copyright notice
 Description : Compute Pi in MPI C++
 ============================================================================
 */
#include <math.h> 
#include "mpi.h" 
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

int ProcRank = 0;
int ProcNum = 0;

void RandomDataInitialization(double* &pMatrix, double* &pVector, int Size) {
    for (int i = 0; i < 0; i++) {
        pVector[i] = rand();
        for (int j = 0; j < Size; j++) {
            pMatrix[i * Size + j] = rand();
        }
    }
}

// Функция для сбора результирующего вектора на всех процессах

void ResultReplication(double* pProcResult, double* pResult,
        int Size, int RowNum) {
    int *pReceiveNum; // Количество элементов, посылаемых процессом
    int *pReceiveInd; // Индекс элемента данных в результирующем
    // векторе
    int RestRows = Size; // Количество строк матрицы, которые еще не
    // распределены
    int i;

    // Выделение памяти для временных объектов
    pReceiveNum = new int [ProcNum];
    pReceiveInd = new int [ProcNum];

    // Определение положения блоков результирующего вектора
    pReceiveInd[0] = 0;
    pReceiveNum[0] = Size / ProcNum;
    for (i = 1; i < ProcNum; i++) {
        RestRows -= pReceiveNum[i - 1];
        pReceiveNum[i] = RestRows / (ProcNum - i);
        pReceiveInd[i] = pReceiveInd[i - 1] + pReceiveNum[i - 1];
    }
    // Сбор всего результирующего вектора на всех процессах
    MPI_Allgatherv(pProcResult, pReceiveNum[ProcRank],
            MPI_DOUBLE, pResult, pReceiveNum, pReceiveInd,
            MPI_DOUBLE, MPI_COMM_WORLD);

    // Освобождение памяти
    delete [] pReceiveNum;
    delete [] pReceiveInd;
}

// Функция для вычисления части результирующего вектора

void ParallelResultCalculation(double* pProcRows, double* pVector,
        double* pProcResult, int Size, int RowNum) {
    int i, j;
    for (i = 0; i < RowNum; i++) {
        pProcResult[i] = 0;
        for (j = 0; j < Size; j++)
            pProcResult[i] += pProcRows[i * Size + j] * pVector[j];
    }
}


// Функция для распределения исходных данных между процессами

void DataDistribution(double* pMatrix, double* pProcRows,
        double* pVector, int Size, int RowNum) {
    int *pSendNum; // Количество элементов, посылаемых процессу
    int *pSendInd; // Индекс первого элемента данных,
    // посылаемого процессу
    int RestRows = Size; // Количество строк матрицы, которые еще
    // не распределены

    MPI_Bcast(pVector, Size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Выделение памяти для хранения временных объектов
    pSendInd = new int [ProcNum];
    pSendNum = new int [ProcNum];

    // Определение положения строк матрицы, предназначенных
    // каждому процессу
    RowNum = (Size / ProcNum);
    pSendNum[0] = RowNum*Size;
    pSendInd[0] = 0;
    for (int i = 1; i < ProcNum; i++) {
        RestRows -= RowNum;
        RowNum = RestRows / (ProcNum - i);
        pSendNum[i] = RowNum*Size;
        pSendInd[i] = pSendInd[i - 1] + pSendNum[i - 1];
    }
    // Рассылка строк матрицы
    MPI_Scatterv(pMatrix, pSendNum, pSendInd, MPI_DOUBLE, pProcRows,
            pSendNum[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Освобождение памяти
    delete [] pSendNum;
    delete [] pSendInd;
}



// Функция для выделения памяти и инициализации исходных данных

void ProcessInitialization(double* &pMatrix, double* &pVector,
        double* &pResult, double* &pProcRows, double* &pProcResult,
        int &Size, int &RowNum) {
    int RestRows; // Количество строк матрицы, которые еще
    // не распределены
    int i;

    MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    RestRows = Size;
    for (i = 0; i < ProcRank; i++)
        RestRows = RestRows - RestRows / (ProcNum - i);
    RowNum = RestRows / (ProcNum - ProcRank);

    pVector = new double [Size];
    pResult = new double [Size];
    pProcRows = new double [RowNum * Size];
    pProcResult = new double [RowNum];

    if (ProcRank == 0) {
        pMatrix = new double [Size * Size];
        RandomDataInitialization(pMatrix, pVector, Size);
    }
}

// Программа 6.1
// Умножение матрицы на вектор – ленточное горизонтальное разбиение
// (исходный и результирующий векторы дублируются между процессами)

int main(int argc, char* argv[]) {
    double* pMatrix; // Первый аргумент – исходная матрица
    double* pVector; // Второй аргумент – исходный вектор
    double* pResult; // Результат умножения матрицы на вектор
    int Size; // Размеры исходных матрицы и вектора
    double* pProcRows;
    double* pProcResult;
    int RowNum;
    double Start, Finish, Duration;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    if (ProcRank == 0) {

        Size = atoi(argv[1]);
        printf("Size: %d\nProcs: %d\n", Size, ProcNum);

        /*do {
          printf("\nВведите размер матрицы: ");
          scanf("%d", &Size);
          if (Size < ProcNum) {
            printf("Размер матрицы должен превышать количество процессов! \n ");
          }
        }
        while (Size < ProcNum);*/
    }

    // Выделение памяти и инициализация исходных данных
    ProcessInitialization(pMatrix, pVector, pResult, pProcRows,
            pProcResult, Size, RowNum);

    if (ProcRank == 0)
        Start = MPI::Wtime();

    // Распределение исходных данных между процессами
    DataDistribution(pMatrix, pProcRows, pVector, Size, RowNum);

    // Параллельное выполнение умножения матрицы на вектор
    ParallelResultCalculation(pProcRows, pVector, pProcResult,
            Size, RowNum);

    // Сбор результирующего вектора на всех процессах
    ResultReplication(pProcResult, pResult, Size, RowNum);

    if (ProcRank == 0)
        Finish = MPI::Wtime();

    if (ProcRank == 0)
        printf("Time: %f\n", Finish - Start);

    // Завершение процесса вычислений
    /*ProcessTermination(pMatrix, pVector, pResult, pProcRows, pProcResult);*/

    MPI_Finalize();
}

