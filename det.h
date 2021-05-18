#pragma once
/*
 * matrix_helpers.h
 *
 *  Created on: Jan 3, 2020
 *      Author: Trotsky
 */

#ifndef DET_H
#define DET_H

typedef struct {
	float* m;
	int n;
} Mat;

void fill_matrix(int, float*);
void print_matrix(int, float*);
void get_min(int, int, Mat, Mat*);
float cofactor(int, int, float);
float det(Mat);

#endif