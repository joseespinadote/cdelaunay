/*
to do :

tarea 2
encapsular intercambio diagonal
propagar el intercambio de diag
revisar robustez algoritmo
*/

#pragma once
#define _CRT_SECURE_NO_WARNINGS 1 
#define _WINSOCK_DEPRECATED_NO_WARNINGS 1 
#define _CRT_SECURE_NO_DEPRECATE 1
#define _CRT_NONSTDC_NO_DEPRECATE 1
#pragma warning (disable: 6262)
#pragma warning (disable: 6031)

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "det.h"

#define MAX_VERTICES 1500000
#define MAX_TRIANGLES 1500000
#define TAMANO_MALLA_X 100
#define TAMANO_MALLA_Y 100
#define BUFFER_SIZE 64
#define VERTICE_COMPARTIDO_1 0
#define VERTICE_COMPARTIDO_2 1
#define VERTICE_OPUESTO 2

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
void intercambioDeDiagonal(Triangle* triangleA, Triangle* triangleB, int Avc1, int Avc2, int Aop, int Bop);

void initMesh(Triangle* triangles, Vertex* vertices) {
	int i, j;
	for (i = 0; i < MAX_VERTICES; i++) {
		vertices[i] = (Vertex){-1, -1};
	}
	for (i = 0; i < MAX_TRIANGLES; i++) {
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
	free(matrix.m);
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

/* Funci�n que retorna id del v�rtice opuesto no compartido por el vecino dado un tri�ngulo */
int getThirdVertexId(Triangle* triangle, Vertex* vertex1, Vertex* vertex2) {
	int i;
	for (i = 0; i < 3; i++)
		if (triangle->vertices[i] != vertex1 && triangle->vertices[i] != vertex2)
			return i;
	return -1;
}

float circleTest(Triangle* triangle, float x, float y) {
	float d;
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
	d = det(matrix);
	free(matrix.m);
	return d;
}

float circleTestByVertex(Triangle* triangle, Vertex* vertice) {
	return circleTest(triangle, vertice->x, vertice->y);
}
/*
Intercambia diagonal entre tri�ngulos A y B siempre que pase el test del c�rculo
Avc1 = es el �nidce del primer vertice de la arista compartida de A con B
Avc2 = es el �nidce del primer vertice de la arista compartida de A con B
Aop = es el �nidce del vertice opuesto a la arista compartida de A con B
Bvc1 = es el �nidce del primer vertice de la arista compartida de B con A
Bvc2 = es el �nidce del primer vertice de la arista compartida de B con A
Bop = es el �nidce del vertice opuesto a la arista compartida de B con A
*/
void intercambioDeDiagonal(Triangle* triangleA, Triangle* triangleB,
	int Avc1, int Avc2, int Aop, int Bop) {
	if (triangleA == NULL || triangleB == NULL) return;
	/*
	esta funci�n:
	 - el triangulo A pasar� a ser compuesto por Avc2/Bvc1, Aop y Bop
	 - el triangulo B pasar� a ser compuesto por Aop, Avc1/Bvc2 y Bop
	 - se debe actualizar el vecindario :
	  - el vecino asociado a Avc1 no cambia
	  - el vecino asociado a Bvc1 no cambia
	  - el vecino asociado a Bop ser� el ex Avc2
	  - el vecino asociado a Avc2 ser� B
	  - el vecino asociado a Aop ser� el ex Bvc2
	  - el vecino asociado a Bvc2 ser� A
	  - el vecino asociado, al vertice opuesto del vecino asociado a Aop, ser� A
	  - el vecino asociado, al vertice opuesto del vecino asociado a Bop, ser� B
	*/
	int Bvc1, Bvc2, idVerticeVecinoLejano;
	Bvc1 = getVertexIdByVertex(triangleB, triangleA->vertices[Avc2]);
	Bvc2 = getVertexIdByVertex(triangleB, triangleA->vertices[Avc1]);
	// se hace el test del c�rculo
	if (circleTest(triangleA, triangleB->vertices[Bop]->x, triangleB->vertices[Bop]->y) > 0.001) {
		// si el resultado es positivo, se hace el cambio de diagonal, que,
		// practicamente, es cambiar los v�rtices Avc2 y Bvc1 de A y B respectivamente
		triangleA->vertices[Avc1] = triangleB->vertices[Bop];
		triangleB->vertices[Bvc1] = triangleA->vertices[Aop];
		// ahora se debe actualizar el vecindario:
		// triangleA->next[Avc1] y triangleB->next[Bvc1] se mantienen
		triangleB->next[Bop] = triangleA->next[Avc2];
		// el vecino lejano referenciado por Avc2, cambia (si existe)
		if (triangleA->next[Avc2] != NULL) {
			idVerticeVecinoLejano = getThirdVertexId(triangleA->next[Avc2], triangleB->vertices[Bvc1], triangleB->vertices[Bvc2]);
			triangleA->next[Avc2]->next[idVerticeVecinoLejano] = triangleB;
		}
		triangleA->next[Avc2] = triangleB;
		triangleA->next[Aop] = triangleB->next[Bvc2];
		// el vecino lejano referenciado por Bvc2, cambia (si es que existe)
		if (triangleB->next[Bvc2] != NULL) {
			idVerticeVecinoLejano = getThirdVertexId(triangleB->next[Bvc2], triangleA->vertices[Avc2], triangleA->vertices[Avc1]);
			triangleB->next[Bvc2]->next[idVerticeVecinoLejano] = triangleA;
		}
		triangleB->next[Bvc2] = triangleA;

		// se propaga el intercambio de diags
		if (triangleA->next[Aop] != NULL) {
			idVerticeVecinoLejano = getThirdVertexId(triangleA->next[Aop], triangleA->vertices[Avc1], triangleA->vertices[Avc2]);
			intercambioDeDiagonal(triangleA, triangleA->next[Aop], Avc1, Avc2, Aop, idVerticeVecinoLejano);
		}
		if (triangleB->next[Bop] != NULL) {
			idVerticeVecinoLejano = getThirdVertexId(triangleB->next[Bop], triangleB->vertices[Bvc1], triangleB->vertices[Bvc2]);
			intercambioDeDiagonal(triangleB, triangleB->next[Bop], Bvc1, Bvc2, Bop, idVerticeVecinoLejano);
		}
	}
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

void generateDelaunayNet(char* strFileInput, char* strFileOutput) {

}

int main(int argc, char* argv[]) {
	FILE* fpInput;
	char fileInput[BUFFER_SIZE], fileOutput[BUFFER_SIZE];
	int i, j,
		// idsVerticesTriangulo es un vector con los vertices ordenados contrarreloj
		// que se crea en funci�n de la arista donde cae el nuevo punto
		idsVerticesTriangulo[3],
		// idVerticeOpuestoVecino es el id del v�rtice opuesto del vecino que no comparte
		// con alguno de los nuevos tri�ngulos n1, n2 o n3. Sirve para hacer test del
		// c�rculo
		idVerticeOpuestoVecino,
		// similar a idsVerticesTriangulo, pero del vecino con el que comparte arista en caso
		// de que el punto cae en arista compartida entre 2 tri�ngulos
		idsVerticesVecino[3] = { -1,-1,-1 },
		// valores iniciales de contadores de trinagulos y v�rtcies totales en la malla
		numTotalTriangles = 2, numTotalVertices = 4,
		// se usa como buleano para determinar si el punto cae en borde o dentro de un
		// tri�ngulo
		puntoEnBorde,
		// se usa como buleano para determinar si el punto ya existe en la malla
		puntoExiste;
	float x, y, dets[3];
	Triangle* triangles, copyTriangle, copyNext;
	Vertex* vertices;
	strcpy(fileInput, argv[1]);
	strcpy(fileOutput, argv[2]);
	triangles = malloc(sizeof(Triangle) * MAX_TRIANGLES);
	vertices = malloc(sizeof(Vertex) * MAX_VERTICES);
	if (triangles == NULL || vertices == NULL) {
		printf("Error malloc\n");
		exit(0);
	}
	initMesh(triangles, vertices);
	fpInput = fopen(fileInput, "r");
	if (!fpInput) {
		printf("No pude abrir el archivo puntos.txt\n");
		return 1;
	}
	while (!feof(fpInput)) {
		fscanf(fpInput, "%f %f", &x, &y);
		// se chequea si el punto ya existe
		puntoExiste = 0;
		for (j = 0; j < numTotalVertices;j++) {
			if (vertices[j].x == x && vertices[j].y == y) {
				puntoExiste = 1;
			}
		}
		if (puntoExiste == 1) continue;
		// se buscar� el triangulo t en cada triangulo[i] de la malla
		for (i = 0; i < numTotalTriangles; i++) {
			// se calcula los determinantes de cada lado del triangulo t su el nuevo punto
			getDetsByTriangle(&triangles[i], dets, x, y);
			// si punto cae en arista, se ordenan v�rtices enumerados como: el primer
			// que comparte esa arista, el segundo y el que no (al que llamaremos opuesto)
			// todo siempre en sentido contrarreloj
			puntoEnBorde = 0;
			// Caso en que el punto cae en el primer borde de t (es decir, entre arista 1 y 2)
			if (dets[0] == 0 && dets[1] > 0 && dets[2] > 0) {
				puntoEnBorde = 1;
				idsVerticesTriangulo[VERTICE_COMPARTIDO_1] = 0;
				idsVerticesTriangulo[VERTICE_COMPARTIDO_2] = 1;
				idsVerticesTriangulo[VERTICE_OPUESTO] = 2;
			}
			// Caso en que el punto cae en el segundo borde de t (es decir, entre arista 2 y 3)
			else if (dets[0] > 0 && dets[1] == 0 && dets[2] > 0) {
				puntoEnBorde = 1;
				idsVerticesTriangulo[VERTICE_COMPARTIDO_1] = 1;
				idsVerticesTriangulo[VERTICE_COMPARTIDO_2] = 2;
				idsVerticesTriangulo[VERTICE_OPUESTO] = 0;
			}
			// Caso en que el punto cae en el tercer borde de t (es decir, entre arista 3 y 1)
			else if (dets[0] > 0 && dets[1] > 0 && dets[2] == 0) {
				puntoEnBorde = 1;
				idsVerticesTriangulo[VERTICE_COMPARTIDO_1] = 2;
				idsVerticesTriangulo[VERTICE_COMPARTIDO_2] = 0;
				idsVerticesTriangulo[VERTICE_OPUESTO] = 1;
			}
			// Caso en que el punto cae dentro del tri�ngulo t
			else if (dets[0] > 0 && dets[1] > 0 && dets[2] > 0) {
				/*
				 - Se crear�n los tri�ngulo n1, n2 y n3 dentro de t, con sus aristas y el nuevo vertice
				 - Se deben configurar la nueva vecindad dado n1, n2 y n3
				 - Los ex vecinos del difunto t, deben actualizarse respecto a n1, n2 y n3
				 - Se debe hacer el test del c�rculo entre:
				    - n1 y V3, donde V3 es el vecino opuesto al v�rtice 3 de t
					- n2 y V1, donde V1 es el vecino opuesto al v�rtice 1 de t
					- n3 y V2, donde V2 es el vecino opuesto al v�rtice 2 de t
				*/
				// guardo un respaldo del triagulo original t donde cae el nuevo punto
				copyTriangle = triangles[i];
				// se limpia el triangulo t que se reutilizar� con el nuevo tri�ngulo n1
				clearTriangle(&triangles[i]);
				// se crea el nuevo v�rtice
				vertices[numTotalVertices].x = x;
				vertices[numTotalVertices].y = y;
				// n1, referido como triangles[i]
				triangles[i].vertices[0] = copyTriangle.vertices[0];
				triangles[i].vertices[1] = copyTriangle.vertices[1];
				triangles[i].vertices[2] = &vertices[numTotalVertices];
				// n2, referido, m�s adelante, como triangles[numTotalTriangles - 2]
				triangles[numTotalTriangles].vertices[0] = &vertices[numTotalVertices];
				triangles[numTotalTriangles].vertices[1] = copyTriangle.vertices[1];
				triangles[numTotalTriangles].vertices[2] = copyTriangle.vertices[2];
				numTotalTriangles++;
				//  n3, referido, m�s adelante, como triangles[numTotalTriangles - 1]
				triangles[numTotalTriangles].vertices[0] = copyTriangle.vertices[0];
				triangles[numTotalTriangles].vertices[1] = &vertices[numTotalVertices];
				triangles[numTotalTriangles].vertices[2] = copyTriangle.vertices[2];
				numTotalTriangles++;
				// se actualizan los vecinos locales de n1:
				// el vecino de n1, en su vertice 1, es n2
				triangles[i].next[0] = &triangles[numTotalTriangles - 2];
				// el vecino de n1, en su vertice 2, es n3
				triangles[i].next[1] = &triangles[numTotalTriangles - 1];
				// el vecino de n1, en su vertice 3, es el ex vecino de t en su vertice 2
				triangles[i].next[2] = copyTriangle.next[2];
				// se actualizan los vecinos locales de n2
				triangles[numTotalTriangles - 2].next[0] = copyTriangle.next[0];
				triangles[numTotalTriangles - 2].next[1] = &triangles[numTotalTriangles - 1];
				triangles[numTotalTriangles - 2].next[2] = &triangles[i];
				// se actualizan los vecinos locales de n3
				triangles[numTotalTriangles - 1].next[0] = &triangles[numTotalTriangles - 2];
				triangles[numTotalTriangles - 1].next[1] = copyTriangle.next[1];
				triangles[numTotalTriangles - 1].next[2] = &triangles[i];
				// �tenia vecino V1 el difunto t frente al vertice 1?, o, en otras palabras,
				// �hay vecino frente al vertice 1 de n2 que se requiera actualizar la referencia?
				if (copyTriangle.next[0] != NULL) {
					// se rescata id del vecino opuesto para actualizar referencia de vecino lejano
					idVerticeOpuestoVecino = getThirdVertexId(
						copyTriangle.next[0],
						copyTriangle.vertices[1],
						copyTriangle.vertices[2]);
					// se actualiza la referencia del vecino respecto a nuevo tri�ngulo n2
					copyTriangle.next[0]->next[idVerticeOpuestoVecino] = &triangles[numTotalTriangles - 2];
				}
				// �tenia vecino V2 el difunto t frente al vertice 2? (o frente a n3)
				if (copyTriangle.next[1] != NULL) {
					idVerticeOpuestoVecino = getThirdVertexId(
						copyTriangle.next[1],
						copyTriangle.vertices[2],
						copyTriangle.vertices[0]);
					copyTriangle.next[1]->next[idVerticeOpuestoVecino] = &triangles[numTotalTriangles - 1];
				}
				// �tenia vecino V3 el difunto t frente al vertice 3? (o frente a n1)
				if (copyTriangle.next[2] != NULL) {
					// idsVerticesVecino[0] es ahora v�rtice opuesto a la arista compartida entre n1 y V3
					idVerticeOpuestoVecino = getThirdVertexId(
						copyTriangle.next[2],
						copyTriangle.vertices[0],
						copyTriangle.vertices[1]);
					// se actualiza referencia del v�rtice opuesto a la arista compartida en el vecino lejano V3
					copyTriangle.next[2]->next[idVerticeOpuestoVecino] = &triangles[i];
				}
				if (copyTriangle.next[0] != NULL) {
					idVerticeOpuestoVecino = getThirdVertexId(
						copyTriangle.next[0],
						copyTriangle.vertices[1],
						copyTriangle.vertices[2]);
					intercambioDeDiagonal(
						&triangles[numTotalTriangles - 2],
						copyTriangle.next[0], 1, 2, 0,
						idVerticeOpuestoVecino);
				}
				if (copyTriangle.next[1] != NULL) {
					idVerticeOpuestoVecino = getThirdVertexId(
						copyTriangle.next[1],
						copyTriangle.vertices[2],
						copyTriangle.vertices[0]);
					intercambioDeDiagonal(
						&triangles[numTotalTriangles - 1],
						copyTriangle.next[1], 2, 0, 1,
						idVerticeOpuestoVecino);
				}
				if (copyTriangle.next[2] != NULL) {
					// idsVerticesVecino[0] es ahora v�rtice opuesto a la arista compartida entre n1 y V3
					idVerticeOpuestoVecino = getThirdVertexId(
						copyTriangle.next[2],
						copyTriangle.vertices[0],
						copyTriangle.vertices[1]);
					intercambioDeDiagonal(
						&triangles[i],
						copyTriangle.next[2], 0, 1, 2,
						idVerticeOpuestoVecino);
				}
				numTotalVertices++;
				break;
			}
			// Caso en que el nuevo punto le�do cae en un v�rtice del tri�ngulo t
			if (puntoEnBorde == 1) {
				// se carga nuevo vertice a la malla
				vertices[numTotalVertices].x = x;
				vertices[numTotalVertices].y = y;
				// se saca una copia de respaldo del triangulo t
				// porque morir� pronto...
				copyTriangle = triangles[i];
				// muere t, quien se transformar� en el nuevo n1
				clearTriangle(&triangles[i]);
				// �ten�a el difunto triangulo t ten�a un vecino en el v�rtice
				// opuesto a la arista donde cae el nuevo punto?
				// (y por ende, �vecino de t con el que comparte la arista donde cae el punto?)
				if (copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]] != NULL) {
					// se saca una copia de respaldo del vecino con el que comparte artista donde cae el punto
					copyNext = *copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]];
					// muere vecino de la arista compartida donde cae el punto
					clearTriangle(copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]]);
					// se arma set de indices del vecino, asegur�ndose que se conserve el ordenamiento contrarreloj
					idsVerticesVecino[VERTICE_COMPARTIDO_1] = getVertexIdByVertex(
						&copyNext,
						copyTriangle.vertices[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]]);
					idsVerticesVecino[VERTICE_COMPARTIDO_2] = getVertexIdByVertex(
						&copyNext,
						copyTriangle.vertices[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]]);
					idsVerticesVecino[VERTICE_OPUESTO] = getThirdVertexId(
						&copyNext,
						copyTriangle.vertices[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]],
						copyTriangle.vertices[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]]);
					// se crea nuevo tri�ngulo n1
					triangles[i].vertices[0] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]];
					triangles[i].vertices[1] = &vertices[numTotalVertices];
					triangles[i].vertices[2] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_OPUESTO]];
					// se crea nuevo tri�ngulo n2
					triangles[numTotalTriangles].vertices[0] = &vertices[numTotalVertices];
					triangles[numTotalTriangles].vertices[1] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]];
					triangles[numTotalTriangles].vertices[2] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_OPUESTO]];
					numTotalTriangles++;
					// se crea nuevo tri�ngulo n3
					copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]]->vertices[0] =
						&vertices[numTotalVertices];
					copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]]->vertices[1] =
						copyNext.vertices[idsVerticesVecino[VERTICE_COMPARTIDO_2]];
					copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]]->vertices[2] =
						copyNext.vertices[idsVerticesVecino[VERTICE_OPUESTO]];
					// se crea nuevo tri�ngulo n4
					triangles[numTotalTriangles].vertices[0] = &vertices[numTotalVertices];
					triangles[numTotalTriangles].vertices[1] = copyNext.vertices[idsVerticesVecino[VERTICE_OPUESTO]];
					triangles[numTotalTriangles].vertices[2] = copyNext.vertices[idsVerticesVecino[VERTICE_COMPARTIDO_1]];
					numTotalTriangles++;
					// se actualiza el vecindario
					// n2 pasa a ser vecino de n1 en su vertice 1 (n1 v0)
					triangles[i].next[0] = &triangles[numTotalTriangles - 2];
					// el ex vecino de t en el v�rtice 2 pasa a ser vecino de n1 en su vertice 1
					triangles[i].next[1] = copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]];
					// de existir, se actualiza el ex vecino de t en el v�rtice 2
					if (copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]] != NULL) {
						idVerticeOpuestoVecino = getThirdVertexId(
							copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]],
							triangles[i].vertices[0],
							triangles[i].vertices[2]);
						copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]]->next[idVerticeOpuestoVecino] =
							&triangles[i];
					}
					// n3 pasa a ser el vecino de n1 en el vertice 3
					triangles[i].next[2] = copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]];
					// el ex vecino de t en el v�rtice 1 pasa a ser vecino de n2 en su vertice 0
					triangles[numTotalTriangles - 2].next[0] = copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]];
					// de existir, se actualiza el ex vecino de t en el v�rtice 1
					if (copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]] != NULL) {
						idVerticeOpuestoVecino = getThirdVertexId(
							copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]],
							triangles[numTotalTriangles - 2].vertices[1],
							triangles[numTotalTriangles - 2].vertices[2]);
						copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]]->next[idVerticeOpuestoVecino] =
							&triangles[numTotalTriangles - 2];
					}
					// n1 pasa a ser vecino de n2 en su vertice 2
					triangles[numTotalTriangles - 2].next[1] = &triangles[i];
					// n4 pasa a ser vecino de n2 en su vertice 3
					triangles[numTotalTriangles - 2].next[2] = &triangles[numTotalTriangles - 1];
					// falta actualizar vecindad de n3
					// el ex vecino de copyNext en el vertice compartido 1 pasa a ser el vecino de n3 en su vertice 0
					copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]]->next[0] = copyNext.next[idsVerticesVecino[VERTICE_COMPARTIDO_1]];
					// de existir, se actualiza el ex vecino de copyNext
					if (copyNext.next[idsVerticesVecino[VERTICE_COMPARTIDO_1]] != NULL) {
						idVerticeOpuestoVecino = getThirdVertexId(
							copyNext.next[idsVerticesVecino[VERTICE_COMPARTIDO_1]],
							copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]]->vertices[1],
							copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]]->vertices[2]);
						copyNext.next[idsVerticesVecino[VERTICE_COMPARTIDO_1]]->next[idVerticeOpuestoVecino] = copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]];
					}
					copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]]->next[1] = &triangles[numTotalTriangles - 1];
					copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]]->next[2] = &triangles[i];
					// vecindad de n4
					triangles[numTotalTriangles - 1].next[0] = copyNext.next[idsVerticesVecino[VERTICE_COMPARTIDO_2]];
					// de existir, se actualiza el ex vecino de copyNext
					if (copyNext.next[idsVerticesVecino[VERTICE_COMPARTIDO_2]] != NULL) {
						idVerticeOpuestoVecino = getThirdVertexId(
							copyNext.next[idsVerticesVecino[VERTICE_COMPARTIDO_2]],
							triangles[numTotalTriangles - 1].vertices[1],
							triangles[numTotalTriangles - 1].vertices[2]);
						copyNext.next[idsVerticesVecino[VERTICE_COMPARTIDO_2]]->next[idVerticeOpuestoVecino] = &triangles[numTotalTriangles - 1];
					}
					// n2 pasa a ser vecino de n4 en el vertice 1
					triangles[numTotalTriangles - 1].next[1] = &triangles[numTotalTriangles - 2];
					triangles[numTotalTriangles - 1].next[2] = copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]];
					// ahora se realizan los test del circulo para cada n
					// intercambio de diagonal para:
					// n1
					if (copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]] != NULL) {
						idVerticeOpuestoVecino = getThirdVertexId(
							copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]],
							triangles[i].vertices[0],
							triangles[i].vertices[2]);
						intercambioDeDiagonal(
							&triangles[i],
							copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]], 2, 0, 1,
							idVerticeOpuestoVecino);
					}
					// n2
					if (copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]] != NULL) {
						idVerticeOpuestoVecino = getThirdVertexId(
							copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]],
							triangles[numTotalTriangles - 2].vertices[1],
							triangles[numTotalTriangles - 2].vertices[2]);
						intercambioDeDiagonal(
							&triangles[numTotalTriangles - 2],
							copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]], 1, 2, 0,
							idVerticeOpuestoVecino);
					}
					// n3
					if (copyNext.next[idsVerticesVecino[VERTICE_COMPARTIDO_1]] != NULL) {
						idVerticeOpuestoVecino = getThirdVertexId(
							copyNext.next[idsVerticesVecino[VERTICE_COMPARTIDO_1]],
							copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]]->vertices[1],
							copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]]->vertices[2]);
						intercambioDeDiagonal(
							copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]],
							copyNext.next[idsVerticesVecino[VERTICE_COMPARTIDO_1]], 1, 2, 0,
							idVerticeOpuestoVecino);
					}	
					// n4
					if (copyNext.next[idsVerticesVecino[VERTICE_COMPARTIDO_2]] != NULL) {
						idVerticeOpuestoVecino = getThirdVertexId(
							copyNext.next[idsVerticesVecino[VERTICE_COMPARTIDO_2]],
							triangles[numTotalTriangles - 1].vertices[1],
							triangles[numTotalTriangles - 1].vertices[2]);
						intercambioDeDiagonal(
							&triangles[numTotalTriangles - 1],
							copyNext.next[idsVerticesVecino[VERTICE_COMPARTIDO_2]], 1, 2, 0,
							idVerticeOpuestoVecino);
					}
				}
				// en el caso que t no comparta la arista donde cay� el nuevo punto con otro trinagulo...
				// (en otras palabras, no tiene vecino en la arista donde cay� el nuevo punto...)
				else {
					// se crea el nuevo tri�ngulo n1
					triangles[i].vertices[0] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]];
					triangles[i].vertices[1] = &vertices[numTotalVertices];
					triangles[i].vertices[2] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_OPUESTO]];
					// se crea el nuevo triangulo n2
					triangles[numTotalTriangles].vertices[0] = &vertices[numTotalVertices];
					triangles[numTotalTriangles].vertices[1] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]];
					triangles[numTotalTriangles].vertices[2] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_OPUESTO]];
					numTotalTriangles++;
					// se actualiza el vecindario
					// n2 pasa a ser vecino de n1 en el v�rtice 1
					triangles[i].next[0] = &triangles[numTotalTriangles - 1];
					// se actualiza el vecino lejano de n1
					triangles[i].next[1] = copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]];
					if (copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]] != NULL) {
						idVerticeOpuestoVecino = getThirdVertexId(
							copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]],
							copyTriangle.vertices[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]],
							copyTriangle.vertices[idsVerticesTriangulo[VERTICE_OPUESTO]]);
						copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]]->next[idVerticeOpuestoVecino] = &triangles[i];
					}
					triangles[i].next[2] = NULL;
					// se actualiza vertice 1 de n2
					triangles[numTotalTriangles - 1].next[0] = copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]];
					if (copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]] != NULL) {
						idVerticeOpuestoVecino = getThirdVertexId(
							copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]],
							copyTriangle.vertices[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]],
							copyTriangle.vertices[idsVerticesTriangulo[VERTICE_OPUESTO]]);
						copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]]->next[idVerticeOpuestoVecino] = &triangles[numTotalTriangles - 1];
					}
					// n1 pasa a ser vecino de n2 en el v�rtice 2
					triangles[numTotalTriangles - 1].next[1] = &triangles[i];
					triangles[numTotalTriangles - 1].next[2] = NULL;
					// intercambio de diagonal para n1, si es que tiene vecino al frente de su vertice 2
					if (copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]] != NULL) {
						idVerticeOpuestoVecino = getThirdVertexId(
							copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]],
							copyTriangle.vertices[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]],
							copyTriangle.vertices[idsVerticesTriangulo[VERTICE_OPUESTO]]);
						intercambioDeDiagonal(
							&triangles[i],
							copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]], 2, 0, 1,
							idVerticeOpuestoVecino);
					}
					// intercambio de diagonal para n2
					if (copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]] != NULL) {
						idVerticeOpuestoVecino = getThirdVertexId(
							copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]],
							copyTriangle.vertices[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]],
							copyTriangle.vertices[idsVerticesTriangulo[VERTICE_OPUESTO]]);
						intercambioDeDiagonal(
							&triangles[numTotalTriangles - 1],
							copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]], 1, 2, 0,
							idVerticeOpuestoVecino);
					}
				}
				numTotalVertices++;
			}
		}
	}
	fclose(fpInput);
	exportData(triangles, numTotalTriangles, "salida.txt");
	return 0;
}
