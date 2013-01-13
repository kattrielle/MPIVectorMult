/*
 * File:   VectorMult.cpp
 * Author: kate
 *
 * Created on 13 Январь 2013 г., 17:01
 */

#include <cstdlib>
#include <math.h>
#include <mpich2/mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;
int ProcRank = 0;
int ProcNum = 0;

/**
 *
 * @param pVector - инициализируемый вектор
 * @param Size - размер вектора
 */
void RandomDataGeneration(double* pVector, int Size)
{
	for (int i = 0; i < Size; i++) {
		//pVector[i] = rand();
		pVector[i] = 1;
	}
}

double Multiplication(double *vector1, double* vector2, int size)
{
	double result;
	for (int i = 0; i < size; i++) {
		result += vector1[i] * vector2[i];
	}
	return result;
}

/*
 *
 */
int main(int argc, char** argv)
{
	double* pVector1; // Первый аргумент – исходная матрица
	double* pVector2; // Второй аргумент – исходный вектор
	double* procVector1; //Распараллеленный кусок первого вектора
	double* procVector2; //Распараллеленный кусок второго вектора
	double pResult; // Результат скалярного умножения векторов
	double currentResult; //Результат умножения "маленьких" векторов
	int Size; // Размеры исходных матрицы и вектора
	int currentSize;
	double startTime, endTime;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	
	Size = atoi(argv[1]);

	if (ProcRank == 0) {
		
		printf("Size: %d\nProcs: %d\n", Size, ProcNum);

		pVector1 = new double [Size];
		pVector2 = new double [Size];

		RandomDataGeneration(pVector1, Size);
		RandomDataGeneration(pVector2, Size);

		//        for (int i = 0; i < Size; i++)
		//            printf("%f\n", pVector1[i]);
		//        printf("%d\n", currentSize);
	}

	currentSize = Size / ProcNum;
	printf("curr size: %d\n", currentSize);

	procVector1 = new double[ currentSize ];
	procVector2 = new double[ currentSize ];

	startTime = MPI_Wtime();

	MPI_Scatter(pVector1, currentSize, MPI_DOUBLE,
		procVector1, currentSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(pVector2, currentSize, MPI_DOUBLE,
		procVector2, currentSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	currentResult = Multiplication(procVector1, procVector2, currentSize);
	printf("curr res:%f\n", currentResult);

	MPI_Reduce(&currentResult, &pResult, 1,
		MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	endTime = MPI_Wtime();

	if (ProcRank == 0) {
		printf("calculation result: %f\n", pResult);
		printf("time: %f\n", endTime - startTime);
	}

	MPI_Finalize();
	return 0;
}

