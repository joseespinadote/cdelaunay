#pragma once
#define _CRT_SECURE_NO_WARNINGS 1 
#define _WINSOCK_DEPRECATED_NO_WARNINGS 1 
#define _CRT_SECURE_NO_DEPRECATE 1
#define _CRT_NONSTDC_NO_DEPRECATE 1
#pragma warning (disable: 6262)
#pragma warning (disable: 6031)
#define VERTICE_COMPARTIDO_1 0
#define VERTICE_COMPARTIDO_2 1
#define VERTICE_OPUESTO 2

/* por josé espina */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "det.h"

#define LARGO_MALLA 1024
#define TAMANO_MALLA_X 10
#define TAMANO_MALLA_Y 10
#define BUFFER_SIZE 64

typedef struct Vertex {
	float x, y;
} Vertex;

typedef struct Triangle {
	Vertex* vertices[3];
	struct Triangle* next[3];
} Triangle;

void initMesh(Triangle* triangles, Vertex* vertices);
void getDetsByTriangle(Triangle* triangle, float* dets, float x, float y);
int getVertexIdByVertex(Triangle* triangle, Vertex* vertex);
Vertex* getThirdVertex(Triangle* triangle, Vertex* vertex1, Vertex* vertex2);
int getThirdVertexId(Triangle* triangle, Vertex* vertex1, Vertex* vertex2);
float circleTest(Triangle* triangle, float x, float y);
float circleTestByVertex(Triangle* triangle, Vertex* vertice);
void exportData(Triangle* triangles, int numTotalTriangles, char* fileOutput);
void clearTriangle(Triangle* triangle);
void intercambioDeDiagonal(Triangle* triangleA, Triangle* triangleB);

void initMesh(Triangle* triangles, Vertex* vertices) {
	int i, j;
	for (i = 0; i < LARGO_MALLA + 2; i++) {
		vertices[i] = (Vertex){-1, -1};
	}
	for (i = 0; i < LARGO_MALLA; i++) {
		for (j = 0; j < 3; j++) {
			triangles[i].vertices[j] = NULL;
			triangles[i].next[j] = NULL;
		}
	}
	vertices[0].x = 0;
	vertices[0].y = TAMANO_MALLA_Y;
	vertices[1].x = 0;
	vertices[1].y = 0;
	vertices[2].x = TAMANO_MALLA_X;
	vertices[2].y = 0;
	vertices[3].x = TAMANO_MALLA_X;
	vertices[3].y = TAMANO_MALLA_Y;
	triangles[0].vertices[0] = &vertices[0];
	triangles[0].vertices[1] = &vertices[1];
	triangles[0].vertices[2] = &vertices[2];
	triangles[1].vertices[0] = &vertices[0];
	triangles[1].vertices[1] = &vertices[2];
	triangles[1].vertices[2] = &vertices[3];
	triangles[0].next[1] = &triangles[1];
	triangles[1].next[2] = &triangles[0];
}

void getDetsByTriangle(Triangle* triangle, float* dets, float x, float y) {
	Mat matrix;
	matrix.m = malloc(4 * sizeof(float));
	matrix.n = 2;
	matrix.m[0] = triangle->vertices[1]->x - triangle->vertices[0]->x;
	matrix.m[1] = x - triangle->vertices[0]->x;
	matrix.m[2] = triangle->vertices[1]->y - triangle->vertices[0]->y;
	matrix.m[3] = y - triangle->vertices[0]->y;
	dets[0] = det(matrix);
	matrix.m[0] = triangle->vertices[2]->x - triangle->vertices[1]->x;
	matrix.m[1] = x - triangle->vertices[1]->x;
	matrix.m[2] = triangle->vertices[2]->y - triangle->vertices[1]->y;
	matrix.m[3] = y - triangle->vertices[1]->y;
	dets[1] = det(matrix);
	matrix.m[0] = triangle->vertices[0]->x - triangle->vertices[2]->x;
	matrix.m[1] = x - triangle->vertices[2]->x;
	matrix.m[2] = triangle->vertices[0]->y - triangle->vertices[2]->y;
	matrix.m[3] = y - triangle->vertices[2]->y;
	dets[2] = det(matrix);
}

int getVertexIdByVertex(Triangle* triangle, Vertex* vertex) {
	int i;
	for (i = 0; i < 3; i++)
		if (triangle->vertices[i] == vertex)
			return i;
	return -1;
}

Vertex* getThirdVertex(Triangle* triangle, Vertex* vertex1, Vertex* vertex2) {
	int i;
	for (i = 0; i < 3; i++)
		if (triangle->vertices[i] != vertex1 && triangle->vertices[i] != vertex2)
			return triangle->vertices[i];
	return NULL;
}

/* Función que retorna id del vértice opuesto no compartido por el vecino dado un triángulo */
int getThirdVertexId(Triangle* triangle, Vertex* vertex1, Vertex* vertex2) {
	int i;
	for (i = 0; i < 3; i++)
		if (triangle->vertices[i] != vertex1 && triangle->vertices[i] != vertex2)
			return i;
	return -1;
}

float circleTest(Triangle* triangle, float x, float y) {
	Mat matrix;
	matrix.m = malloc(16 * sizeof(float));
	matrix.n = 4;
	matrix.m[0] = triangle->vertices[0]->x;
	matrix.m[1] = triangle->vertices[0]->y;
	matrix.m[2] = powf(triangle->vertices[0]->x, 2) + powf(triangle->vertices[0]->y, 2);
	matrix.m[3] = 1;
	matrix.m[4] = triangle->vertices[1]->x;
	matrix.m[5] = triangle->vertices[1]->y;
	matrix.m[6] = powf(triangle->vertices[1]->x, 2) + powf(triangle->vertices[1]->y, 2);
	matrix.m[7] = 1;
	matrix.m[8] = triangle->vertices[2]->x;
	matrix.m[9] = triangle->vertices[2]->y;
	matrix.m[10] = powf(triangle->vertices[2]->x, 2) + powf(triangle->vertices[2]->y, 2);
	matrix.m[11] = 1;
	matrix.m[12] = x;
	matrix.m[13] = y;
	matrix.m[14] = powf(x, 2) + powf(y, 2);
	matrix.m[15] = 1;
	return det(matrix);
}

float circleTestByVertex(Triangle* triangle, Vertex* vertice) {
	return circleTest(triangle, vertice->x, vertice->y);
}

void intercambioDeDiagonal(Triangle* triangleA, Triangle* triangleB,
	int Avc1, int Avc2, int Aop, int Bop) {
	/*
	 - determinar los vertices compartidos y opuestos entre A y B, sin perder
	dirección contrarreloj
	 - el triangulo A pasará a ser compuesto por Avc2/Bvc1, Aop y Bop
	 - el triangulo B pasará a ser compuesto por Aop, Avc1/Bvc2 y Bop
	 - se debe actualizar el vecindario :
	  - el vecino asociado a Avc1 no cambia
	  - el vecino asociado a Bvc1 no cambia
	  - el vecino asociado a Bop será el ex Avc2
	  - el vecino asociado a Avc2 será B
	  - el vecino asociado a Aop será el ex Bvc2
	  - el vecino asociado a Bvc2 será A
	  - el vecino asociado, al vertice opuesto del vecino asociado a Aop, será A
	  - el vecino asociado, al vertice opuesto del vecino asociado a Bop, será B
	donde:
		Avc1 = es el ínidce del primer vertice de la arista compartida de A con B
		Avc2 = es el ínidce del primer vertice de la arista compartida de A con B
		Aop = es el ínidce del vertice opuesto a la arista compartida de A con B
		Bvc1 = es el ínidce del primer vertice de la arista compartida de B con A
		Bvc2 = es el ínidce del primer vertice de la arista compartida de B con A
		Bop = es el ínidce del vertice opuesto a la arista compartida de B con A
	*/
	int Bvc1, Bvc2;
}

void exportData(Triangle *triangles, int numTotalTriangles, char *fileOutput) {
	FILE* fpOutput;
	int i;
	fpOutput = fopen(fileOutput, "w");
	for (i = 0; i < numTotalTriangles; i++) {
		fprintf(fpOutput, "%f %f\n",
			triangles[i].vertices[0]->x,
			triangles[i].vertices[0]->y
		);
		fprintf(fpOutput, "%f %f\n",
			triangles[i].vertices[1]->x,
			triangles[i].vertices[1]->y
		);
		fprintf(fpOutput, "%f %f\n",
			triangles[i].vertices[2]->x,
			triangles[i].vertices[2]->y
		);
		fprintf(fpOutput, "%f %f\n\n",
			triangles[i].vertices[0]->x,
			triangles[i].vertices[0]->y
		);
	}
	fclose(fpOutput);
}

void clearTriangle(Triangle *triangle) {
	triangle->next[0] = NULL;
	triangle->next[1] = NULL;
	triangle->next[2] = NULL;
	triangle->vertices[0] = NULL;
	triangle->vertices[1] = NULL;
	triangle->vertices[2] = NULL;
}

int main(int argc, char* argv[]) {
	FILE* fpInput;
	char fileInput[BUFFER_SIZE], fileOutput[BUFFER_SIZE];
	int i, j,
		// idsVerticesTriangulo es un vector con los vertices ordenados contrarreloj
		// que se crea en función de la arista donde cae el nuevo punto
		idsVerticesTriangulo[3],
		// idVerticeOpuestoVecino es el id del vértice opuesto del vecino que no comparte
		// con alguno de los nuevos triángulos n1, n2 o n3. Sirve para hacer test del
		// círculo
		idVerticeOpuestoVecino,
		// similar a idsVerticesTriangulo, pero del vecino con el que comparte arista en caso
		// de que el punto cae en arista compartida entre 2 triángulos
		idsVerticesVecino[3] = { -1,-1,-1 },
		// valores iniciales de contadores de trinagulos y vértcies totales en la malla
		numTotalTriangles = 2, numTotalVertices = 4,
		// se usa como buleano para determinar si el punto cae en borde o dentro de un
		// triángulo
		puntoEnBorde,
		// se usa como apoyo para actualizar vecinos muy lejanos (es decir, vecinos
		// de vecinos respecto al triángulo donde cae el nuevo punto)
		idVerticeOpuestoVecinoLejano;
	float x, y, dets[3];
	Triangle triangles[LARGO_MALLA], copyTriangle, copyNext, *farNext;
	Vertex vertices[LARGO_MALLA + 2];

	strcpy(fileInput, argv[1]);
	strcpy(fileOutput, argv[2]);

	initMesh(triangles, vertices);

	fpInput = fopen(fileInput, "r");
	if (!fpInput) {
		printf("No pude abrir el archivo puntos.txt\n");
		return 1;
	}
	while (!feof(fpInput)) {
		fscanf(fpInput, "%f %f", &x, &y);
		// se buscará el triangulo t en cada triangulo[i] de la malla
		for (i = 0; i < numTotalTriangles; i++) {
			// se calcula los determinantes de cada lado del triangulo t su el nuevo punto
			getDetsByTriangle(&triangles[i], dets, x, y);
			// se limpian las variables copyTriangle, que permite llevar una copia
			// del triangulo t, que será eliminado, y copyNext de un posible
			// vecino de t (el que se necesite en el momento)
			clearTriangle(&copyTriangle);
			clearTriangle(&copyNext);
			// si punto cae en arista, se ordenan vértices enumerados como: el primer
			// que comparte esa arista, el segundo y el que no (al que llamaremos opuesto)
			// todo siempre en sentido contrarreloj
			puntoEnBorde = 0;
			if (dets[0] == 0 && dets[1] > 0 && dets[2] > 0) {
				puntoEnBorde = 1;
				idsVerticesTriangulo[VERTICE_COMPARTIDO_1] = 0;
				idsVerticesTriangulo[VERTICE_COMPARTIDO_2] = 1;
				idsVerticesTriangulo[VERTICE_OPUESTO] = 2;
			} else if (dets[0] > 0 && dets[1] == 0 && dets[2] > 0) {
				puntoEnBorde = 1;
				idsVerticesTriangulo[VERTICE_COMPARTIDO_1] = 1;
				idsVerticesTriangulo[VERTICE_COMPARTIDO_2] = 2;
				idsVerticesTriangulo[VERTICE_OPUESTO] = 0;
			} else if (dets[0] > 0 && dets[1] > 0 && dets[2] == 0) {
				puntoEnBorde = 1;
				idsVerticesTriangulo[VERTICE_COMPARTIDO_1] = 2;
				idsVerticesTriangulo[VERTICE_COMPARTIDO_2] = 0;
				idsVerticesTriangulo[VERTICE_OPUESTO] = 1;
			} else if (dets[0] > 0 && dets[1] > 0 && dets[2] > 0) {
				/*
				Caso en que el punto cae dentro del triángulo t
				 - Se crearán los triángulo n1, n2 y n3 dentro de t, con sus aristas y el nuevo vertice
				 - Se deben configurar la nueva vecindad dado n1, n2 y n3
				 - Los ex vecinos del difunto t, deben actualizarse respecto a n1, n2 y n3
				 - Se debe hacer el test del círculo entre:
				    - n1 y V3, donde V3 es el vecino opuesto al vértice 3 de t
					- n2 y V1, donde V1 es el vecino opuesto al vértice 1 de t
					- n3 y V2, donde V2 es el vecino opuesto al vértice 2 de t
				*/
				// guardo un respaldo del triagulo original t donde cae el nuevo punto
				copyTriangle = triangles[i];
				// se limpia el triangulo t que se reutilizará con el nuevo triángulo n1
				clearTriangle(&triangles[i]);
				// se crea el nuevo vértice
				vertices[numTotalVertices-1].x = x;
				vertices[numTotalVertices-1].y = y;
				// n1
				triangles[i].vertices[0] = copyTriangle.vertices[0];
				triangles[i].vertices[1] = copyTriangle.vertices[1];
				triangles[i].vertices[2] = &vertices[numTotalVertices - 1];
				// n2
				triangles[numTotalTriangles].vertices[0] = &vertices[numTotalVertices - 1];
				triangles[numTotalTriangles].vertices[1] = copyTriangle.vertices[1];
				triangles[numTotalTriangles].vertices[2] = copyTriangle.vertices[2];
				numTotalTriangles++;
				//  n3
				triangles[numTotalTriangles].vertices[0] = copyTriangle.vertices[0];
				triangles[numTotalTriangles].vertices[1] = &vertices[numTotalVertices - 1];
				triangles[numTotalTriangles].vertices[2] = copyTriangle.vertices[2];
				numTotalTriangles++;
				// se actualizan los vecinos locales de n1:
				//  - en el vértice 1 de n1 está n2
				//  - en el vértice 2 de n1 está n3
				// el vértice 3 se actualizará más adelante, cuando se revise si 
				// es que hay vecino en el vertice 3 del difunto t
				triangles[i].next[0] = &triangles[numTotalTriangles - 1];
				triangles[i].next[1] = &triangles[numTotalTriangles - 2];
				triangles[i].next[2] = copyTriangle.next[2];
				// se actualizan los vecinos locales de n2
				triangles[numTotalTriangles - 2].next[0] = copyTriangle.next[0];
				triangles[numTotalTriangles - 2].next[1] = &triangles[numTotalTriangles - 1];
				triangles[numTotalTriangles - 2].next[2] = &triangles[i];
				// se actualizan los vecinos locales de n3
				triangles[numTotalTriangles - 1].next[0] = &triangles[numTotalTriangles - 2];
				triangles[numTotalTriangles - 1].next[1] = copyTriangle.next[1];
				triangles[numTotalTriangles - 1].next[2] = &triangles[i];
				// ¿tenia vecino V1 el difunto t frente al vertice 1?, o, en otras palabras,
				// ¿hay vecino frente al vertice 1 de n2 que se requiera actualizar la referencia?
				if (copyTriangle.next[0] != NULL) {
					// se rescata id del vecino opuesto para actualizar referencia de vecino lejano
					idVerticeOpuestoVecino = getThirdVertexId(
						copyTriangle.next[0],
						copyTriangle.next[0]->vertices[1],
						copyTriangle.next[0]->vertices[2]);
					// se actualiza la referencia del vecino respecto a nuevo triángulo n2
					copyTriangle.next[0]->next[idVerticeOpuestoVecino] = &triangles[numTotalTriangles - 2];
					// se realiza procedimiento del intercambio de diagonal, incluyendo el test del círculo
					intercambioDeDiagonal(
						&triangles[numTotalTriangles - 2],
						copyTriangle.next[0], 1, 2, 0,
						idVerticeOpuestoVecino);
				}
				// ¿tenia vecino V2 el difunto t frente al vertice 2? (o frente a n3)
				if (copyTriangle.next[1] != NULL) {
					idVerticeOpuestoVecino = getThirdVertexId(copyTriangle.next[1], copyTriangle.vertices[2], copyTriangle.vertices[0]);
					copyTriangle.next[1]->next[idVerticeOpuestoVecino] = &triangles[numTotalTriangles - 1];
					// se hace test del círculo e intercambio de diagonal
					intercambioDeDiagonal(
						&triangles[numTotalTriangles - 1],
						copyTriangle.next[1], 2, 0, 1,
						idVerticeOpuestoVecino);
				}
				// ¿tenia vecino V3 el difunto t frente al vertice 3? (o frente a n1)
				if (copyTriangle.next[2] != NULL) {
					// idsVerticesVecino[0] es ahora vértice opuesto a la arista compartida entre n1 y V3
					idVerticeOpuestoVecino = getThirdVertexId(copyTriangle.next[2], copyTriangle.next[2]->vertices[0], copyTriangle.next[2]->vertices[1]);
					// se actualiza referencia del vértice opuesto a la arista compartida en el vecino lejano V3
					copyTriangle.next[2]->next[idVerticeOpuestoVecino] = &triangles[i];
					intercambioDeDiagonal(
						&triangles[i],
						copyTriangle.next[2], 0, 1, 2,
						idVerticeOpuestoVecino);
				}
				numTotalVertices++;
				break;
			}
			if (puntoEnBorde == 1) {
				vertices[numTotalVertices-1].x = x;
				vertices[numTotalVertices-1].y = y;
				copyTriangle = triangles[i];
				clearTriangle(&triangles[i]);
				clearTriangle(&copyNext);
				if (triangles[i].next[VERTICE_OPUESTO] != NULL) {
					copyNext = *triangles[i].next[VERTICE_OPUESTO];
					clearTriangle(triangles[i].next[VERTICE_OPUESTO]);
					idsVerticesVecino[VERTICE_COMPARTIDO_1] = getVertexIdByVertex(&copyNext, triangles[i].vertices[VERTICE_COMPARTIDO_2]);
					idsVerticesVecino[VERTICE_COMPARTIDO_2] = getVertexIdByVertex(&copyNext, triangles[i].vertices[VERTICE_COMPARTIDO_1]);
					idsVerticesVecino[VERTICE_OPUESTO] = getThirdVertexId(&copyNext, triangles[i].vertices[VERTICE_COMPARTIDO_1], triangles[i].vertices[VERTICE_COMPARTIDO_2]);
					printf("tr: %d, %d, %d\n", idsVerticesTriangulo[VERTICE_COMPARTIDO_1], idsVerticesTriangulo[VERTICE_COMPARTIDO_2], idsVerticesTriangulo[VERTICE_OPUESTO]);
					printf("v: %d, %d, %d\n\n", idsVerticesVecino[VERTICE_COMPARTIDO_1],idsVerticesVecino[VERTICE_COMPARTIDO_2],idsVerticesVecino[VERTICE_OPUESTO]);
					triangles[i].vertices[0] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]];
					triangles[i].vertices[1] = &vertices[numTotalVertices - 1];
					triangles[i].vertices[2] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_OPUESTO]];
					triangles[numTotalTriangles].vertices[0] = &vertices[numTotalVertices - 1];
					triangles[numTotalTriangles].vertices[1] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]];
					triangles[numTotalTriangles].vertices[2] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_OPUESTO]];
					numTotalTriangles++;
					triangles[i].next[VERTICE_OPUESTO]->vertices[0] = &vertices[numTotalVertices - 1];
					triangles[i].next[VERTICE_OPUESTO]->vertices[1] = copyNext.vertices[idsVerticesVecino[VERTICE_COMPARTIDO_2]];
					triangles[i].next[VERTICE_OPUESTO]->vertices[2] = copyNext.vertices[idsVerticesVecino[VERTICE_OPUESTO]];
					triangles[numTotalTriangles].vertices[0] = &vertices[numTotalVertices - 1];
					triangles[numTotalTriangles].vertices[1] = copyNext.vertices[idsVerticesVecino[VERTICE_OPUESTO]];
					triangles[numTotalTriangles].vertices[2] = copyNext.vertices[idsVerticesVecino[VERTICE_COMPARTIDO_1]];
					numTotalTriangles++;
					/* actualizacion del vecindario local */
					triangles[i].next[0] = &triangles[numTotalTriangles - 2];
					triangles[i].next[2] = triangles[i].next[VERTICE_OPUESTO];
					triangles[numTotalTriangles - 2].next[1] = &triangles[i];
					triangles[numTotalTriangles - 2].next[2] = &triangles[numTotalTriangles - 1];
					triangles[i].next[VERTICE_OPUESTO]->next[1] = &triangles[numTotalTriangles - 1];
					triangles[i].next[VERTICE_OPUESTO]->next[2] = &triangles[i];
					triangles[numTotalTriangles - 1].next[1] = &triangles[numTotalTriangles - 2];
					triangles[numTotalTriangles - 1].next[2] = triangles[i].next[VERTICE_OPUESTO];
					farNext = copyTriangle.next[VERTICE_COMPARTIDO_1];
					if (farNext != NULL) {
						idVerticeOpuestoVecinoLejano = getThirdVertexId(farNext, copyTriangle.vertices[VERTICE_COMPARTIDO_2], copyTriangle.vertices[VERTICE_OPUESTO]);
						farNext->next[idVerticeOpuestoVecinoLejano] = &triangles[numTotalTriangles - 2];
					}
					farNext = copyTriangle.next[VERTICE_COMPARTIDO_2];
					if (farNext != NULL) {
						idVerticeOpuestoVecinoLejano = getThirdVertexId(farNext, copyTriangle.vertices[VERTICE_COMPARTIDO_1], copyTriangle.vertices[VERTICE_OPUESTO]);
						farNext->next[idVerticeOpuestoVecinoLejano] = &triangles[i];
					}
					farNext = copyNext.next[VERTICE_COMPARTIDO_1];
					if (farNext != NULL) {
						idVerticeOpuestoVecinoLejano = getThirdVertexId(farNext, copyNext.vertices[VERTICE_COMPARTIDO_2], copyNext.vertices[VERTICE_OPUESTO]);
						farNext->next[idVerticeOpuestoVecinoLejano] = triangles[i].next[VERTICE_OPUESTO];
					}
					farNext = copyNext.next[VERTICE_COMPARTIDO_2];
					if (farNext != NULL) {
						idVerticeOpuestoVecinoLejano = getThirdVertexId(farNext, copyNext.vertices[VERTICE_COMPARTIDO_1], copyNext.vertices[VERTICE_OPUESTO]);
						farNext->next[idVerticeOpuestoVecinoLejano] = &triangles[numTotalTriangles - 1];
					}
					if (triangles[i].next[1] != NULL) {
						idVerticeOpuestoVecinoLejano = getThirdVertexId(triangles[i].next[1], triangles[i].vertices[0], triangles[i].vertices[2]);
						printf("c: %.2f\n", circleTest(&triangles[i], triangles[i].next[1]->vertices[idVerticeOpuestoVecinoLejano]->x, triangles[i].next[1]->vertices[idVerticeOpuestoVecinoLejano]->y));
					}
					//triangles[numTotalTriangles - 2]
				} else {
					triangles[i].vertices[0] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]];
					triangles[i].vertices[1] = &vertices[numTotalVertices - 1];
					triangles[i].vertices[2] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_OPUESTO]];
					triangles[numTotalTriangles].vertices[0] = &vertices[numTotalVertices - 1];
					triangles[numTotalTriangles].vertices[1] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]];
					triangles[numTotalTriangles].vertices[2] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_OPUESTO]];
					numTotalTriangles++;
					triangles[i].next[0] = &triangles[numTotalTriangles - 1];
					triangles[numTotalTriangles - 1].next[1] = &triangles[i];
					farNext = copyTriangle.next[VERTICE_COMPARTIDO_2];
					if (farNext != NULL) {
						idVerticeOpuestoVecinoLejano = getThirdVertexId(farNext, copyTriangle.vertices[VERTICE_COMPARTIDO_1], copyTriangle.vertices[VERTICE_OPUESTO]);
						farNext->next[idVerticeOpuestoVecinoLejano] = &triangles[numTotalTriangles - 1];
					}
					farNext = copyTriangle.next[VERTICE_COMPARTIDO_1];
					if (farNext != NULL) {
						idVerticeOpuestoVecinoLejano = getThirdVertexId(farNext, copyTriangle.vertices[VERTICE_COMPARTIDO_2], copyTriangle.vertices[VERTICE_OPUESTO]);
						farNext->next[idVerticeOpuestoVecinoLejano] = &triangles[i];
					}
					/* circleTest */
					if (triangles[i].next[1] != NULL) {
						idVerticeOpuestoVecinoLejano = getThirdVertexId(triangles[i].next[1], triangles[i].vertices[0], triangles[i].vertices[2]);
						printf("a: %.2f\n", circleTest(&triangles[i], triangles[i].next[1]->vertices[idVerticeOpuestoVecinoLejano]->x, triangles[i].next[1]->vertices[idVerticeOpuestoVecinoLejano]->y));
					}
					if (triangles[numTotalTriangles-1].next[0] != NULL) {
						idVerticeOpuestoVecinoLejano = getThirdVertexId(triangles[numTotalTriangles - 1].next[0], triangles[numTotalTriangles - 1].vertices[1], triangles[numTotalTriangles - 1].vertices[2]);
						printf("b: %.2f\n", circleTest(triangles[numTotalTriangles - 1].next[0], triangles[numTotalTriangles - 1].next[1]->vertices[idVerticeOpuestoVecinoLejano]->x, triangles[numTotalTriangles - 1].next[1]->vertices[idVerticeOpuestoVecinoLejano]->y));
					}
				}
				numTotalVertices++;
			}
		}
	}
	exportData(triangles, numTotalTriangles, "salida.txt");
	return 0;
}
