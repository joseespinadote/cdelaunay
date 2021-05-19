/*
 * matrix_helpers.c
 *
 *  Created on: Jan 3, 2020
 *      Author: Trotsky
 */

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "det.h"

void get_min(int x, int y, Mat m, Mat* new_m) {
	int iCont;

	new_m->n = m.n - 1;
	new_m->m = (float*)malloc(new_m->n * new_m->n * sizeof(float));

	iCont = 0;
	for (int i = 0; i < m.n; i++) {
		if (i != y) {
			for (int j = 0; j < m.n; j++) {
				if (j != x) {
					new_m->m[iCont] = (float)m.m[i * m.n + j];
					iCont++;
				}
			}
		}
	}
}

float det(Mat matrix) {
	float fDet = 0;
	float* m = matrix.m;

	if (matrix.n == 1) {
		fDet = m[0];
	}
	else if (matrix.n == 2) {
		fDet = m[0] * m[3] - m[1] * m[2];
	}
	else if (matrix.n > 2) {
		fDet = 0;
		for (int j = 0; j < matrix.n; j++) {
			Mat mAux;
			get_min(j, 0, matrix, &mAux);
			fDet += cofactor(j, 0, matrix.m[j]) * det(mAux);
			free(mAux.m);
		}
	}
	return fDet;
}

float cofactor(int x, int y, float num) {
	return powf(-1, (float)x + (float)y) * num;
}
