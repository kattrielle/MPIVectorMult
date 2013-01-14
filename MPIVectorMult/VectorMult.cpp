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
#include <mpi/mpi.h>

using namespace std;
int ProcRank = 0;
int ProcNum = 0;

/**
 * Генерация данных для вектора
 * @param pVector - инициализируемый вектор
 * @param Size - размер вектора
 */
void RandomDataGeneration(double* pVector, int Size)
{
	for (int i = 0; i < Size; i++) {
		pVector[i] = rand();
		//pVector[i] = 1;
	}
	for (int i = Size; i < Size + Size % ProcNum; i++ )
		pVector[i] = 0;
}

/**
 * Реализация скалярного произведения
 * @param vector1
 * @param vector2 
 * @param size длины векторов
 * @return результат скалярного произведения
 */
double Multiplication(double *vector1, double* vector2, int size)
{
	double result = 0;
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
	double* pVector1; // Первый вектор для умножения
	double* pVector2; // Второй вектор для умножения
	double* procVector1; //Распараллеленный кусок первого вектора
	double* procVector2; //Распараллеленный кусок второго вектора
	double pResult; // Результат скалярного умножения векторов
	double currentResult; //Результат умножения "маленьких" векторов
	int Size; // Размеры исходных векторов
	int currentSize; // Размер распараллеленных векторов
	double startTime, endTime;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	
	Size = atoi(argv[1]);

	if (ProcRank == 0) {
		
//		printf("Size: %d\nProcs: %d\n", Size, ProcNum);

		pVector1 = new double [Size + Size % ProcNum];
		pVector2 = new double [Size + Size % ProcNum];

		RandomDataGeneration(pVector1, Size);
		RandomDataGeneration(pVector2, Size);
	}

	currentSize = Size / ProcNum;
	//printf("curr size: %d\n", currentSize);

	procVector1 = new double[ currentSize ];
	procVector2 = new double[ currentSize ];

	startTime = MPI_Wtime();

	MPI_Scatter(pVector1, currentSize, MPI_DOUBLE,
		procVector1, currentSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(pVector2, currentSize, MPI_DOUBLE,
		procVector2, currentSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	currentResult = Multiplication(procVector1, procVector2, currentSize);
	//printf("curr res:%f\n", currentResult);

	MPI_Reduce(&currentResult, &pResult, 1,
		MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	endTime = MPI_Wtime();

	if (ProcRank == 0) {
//		printf("calculation result: %f\n", pResult);
		printf(" Size: %d Proc: %d time: %f Precision: %f\n", Size, 
			ProcNum, endTime - startTime, MPI_Wtick());
	}

	MPI_Finalize();
	return 0;
}

