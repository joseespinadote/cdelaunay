/*
to do :

propagar el intercambio de diag
encapsular codigos
realizar cambios de diagonal en el nuevo punto
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

#define MAX_VERTICES 3000
#define MAX_TRIANGLES 3000
#define MAX_SEGMENTS 1000
#define TAMANO_MALLA_X 100
#define TAMANO_MALLA_Y 100
#define BUFFER_SIZE 64
#define VERTICE_COMPARTIDO_1 0
#define VERTICE_COMPARTIDO_2 1
#define VERTICE_OPUESTO 2
#define EPSILON 0.001

typedef struct Vertex {
	float x, y;
} Vertex;

typedef struct Triangle {
	Vertex* vertices[3];
	struct Triangle* next[3];
} Triangle;

typedef struct Segment {
	Vertex v1, v2;
} Segment;

void initMesh(Triangle* triangles, Vertex* vertices, Segment* constraints);
void getDetsByTriangle(Triangle* triangle, float* dets, float x, float y);
int getVertexIdByVertex(Triangle* triangle, Vertex* vertex);
Vertex* getThirdVertex(Triangle* triangle, Vertex* vertex1, Vertex* vertex2);
int getThirdVertexId(Triangle* triangle, Vertex* vertex1, Vertex* vertex2);
float circleTest(Triangle* triangle, float x, float y);
float circleTestByVertex(Triangle* triangle, Vertex* vertice);
void exportData(Triangle* triangles, int numTotalTriangles, char* fileOutput);
void exportDataRestricted(Segment* constraints, int numTotalSegments, char* fileOutput);
void clearTriangle(Triangle* triangle);
int intercambioDeDiagonal(Triangle* triangleA, Triangle* triangleB, int Avc1, int Avc2, int Aop, int Bop, int skipCircleTest, int skipPropagation);
int generateDelaunayNet(char* strFileInput, Triangle* triangles, Vertex* vertices);
float getDetBySegments(float x1, float y1, float x2, float y2, float xp, float yp);
int restrictDelaunayNet(char* strFileInput, Triangle* triangles, int numTotalTriangles, Segment* constraints);
Triangle* getIdTriangleContainsPoint(Triangle* triangle, int* puntoEnBorde, float x, float y);

void initMesh(Triangle* triangles, Vertex* vertices, Segment* constraints) {
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
	for (i = 0; i < MAX_SEGMENTS; i++) {
		constraints[i] = (Segment){ (Vertex) { -1, -1 }, (Vertex) { -1, -1 } };
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
int intercambioDeDiagonal(Triangle* triangleA, Triangle* triangleB,
	int Avc1, int Avc2, int Aop, int Bop, int skipCircleTest, int skipPropagation) {
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
	if (skipCircleTest == 1 ||
		(circleTest(triangleA, triangleB->vertices[Bop]->x, triangleB->vertices[Bop]->y) > EPSILON)
		) {
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

		if (skipPropagation == 0) {
			// se propaga el intercambio de diags
			if (triangleA->next[Aop] != NULL) {
				idVerticeVecinoLejano = getThirdVertexId(triangleA->next[Aop], triangleA->vertices[Avc1], triangleA->vertices[Avc2]);
				if (intercambioDeDiagonal(triangleA, triangleA->next[Aop], Avc1, Avc2, Aop, idVerticeVecinoLejano, 0, 0) == 0) {
					if (triangleA->next[Aop]->next[getVertexIdByVertex(triangleA->next[Aop], triangleA->vertices[Avc2])] != NULL) {
						if (intercambioDeDiagonal(
							triangleA->next[Aop],
							triangleA->next[Aop]->next[getVertexIdByVertex(triangleA->next[Aop], triangleA->vertices[Avc2])],
							getVertexIdByVertex(triangleA->next[Aop], triangleA->vertices[Avc1]),
							getThirdVertexId(triangleA->next[Aop], triangleA->vertices[Avc1], triangleA->vertices[Avc2]),
							getVertexIdByVertex(triangleA->next[Aop], triangleA->vertices[Avc2]),
							getThirdVertexId(
								triangleA->next[Aop]->next[getVertexIdByVertex(triangleA->next[Aop], triangleA->vertices[Avc2])],
								triangleA->next[Aop]->vertices[getVertexIdByVertex(triangleA->next[Aop], triangleA->vertices[Avc1])],
								triangleA->next[Aop]->vertices[getThirdVertexId(triangleA->next[Aop], triangleA->vertices[Avc1], triangleA->vertices[Avc2])]
							), 0, 0) == 0) {
							if (triangleA->next[Aop]->next[getVertexIdByVertex(triangleA->next[Aop], triangleA->vertices[Avc1])] != NULL) {
								intercambioDeDiagonal(
									triangleA->next[Aop],
									triangleA->next[Aop]->next[getVertexIdByVertex(triangleA->next[Aop], triangleA->vertices[Avc1])],
									getThirdVertexId(triangleA->next[Aop], triangleA->vertices[Avc1], triangleA->vertices[Avc2]),
									getVertexIdByVertex(triangleA->next[Aop], triangleA->vertices[Avc2]),
									getVertexIdByVertex(triangleA->next[Aop], triangleA->vertices[Avc1]),
									getThirdVertexId(
										triangleA->next[Aop]->next[getVertexIdByVertex(triangleA->next[Aop], triangleA->vertices[Avc1])],
										triangleA->next[Aop]->vertices[getThirdVertexId(triangleA->next[Aop], triangleA->vertices[Avc1], triangleA->vertices[Avc2])],
										triangleA->next[Aop]->vertices[getVertexIdByVertex(triangleA->next[Aop], triangleA->vertices[Avc2])]
									), 0, 0);
							}
						}
					}
				}
			}
			if (triangleB->next[Bop] != NULL) {
				idVerticeVecinoLejano = getThirdVertexId(triangleB->next[Bop], triangleB->vertices[Bvc1], triangleB->vertices[Bvc2]);
				if (intercambioDeDiagonal(triangleB, triangleB->next[Bop], Bvc1, Bvc2, Bop, idVerticeVecinoLejano, 0, 0) == 0) {
					if (triangleB->next[Bop]->next[getVertexIdByVertex(triangleB->next[Bop], triangleB->vertices[Bvc2])] != NULL) {
						if (intercambioDeDiagonal(
							triangleB->next[Bop],
							triangleB->next[Bop]->next[getVertexIdByVertex(triangleB->next[Bop], triangleB->vertices[Bvc2])],
							getVertexIdByVertex(triangleB->next[Bop], triangleB->vertices[Bvc1]),
							getThirdVertexId(triangleB->next[Bop], triangleB->vertices[Bvc1], triangleB->vertices[Bvc2]),
							getVertexIdByVertex(triangleB->next[Bop], triangleB->vertices[Bvc2]),
							getThirdVertexId(
								triangleB->next[Bop]->next[getVertexIdByVertex(triangleB->next[Bop], triangleB->vertices[Bvc2])],
								triangleB->next[Bop]->vertices[getVertexIdByVertex(triangleB->next[Bop], triangleB->vertices[Bvc1])],
								triangleB->next[Bop]->vertices[getThirdVertexId(triangleB->next[Bop], triangleB->vertices[Bvc1], triangleB->vertices[Bvc2])]
							), 0, 0) == 0) {
							if (triangleB->next[Bop]->next[getVertexIdByVertex(triangleB->next[Bop], triangleB->vertices[Bvc1])] != NULL) {
								intercambioDeDiagonal(
									triangleB->next[Bop],
									triangleB->next[Bop]->next[getVertexIdByVertex(triangleB->next[Bop], triangleB->vertices[Bvc1])],
									getThirdVertexId(triangleB->next[Bop], triangleB->vertices[Bvc1], triangleB->vertices[Bvc2]),
									getVertexIdByVertex(triangleB->next[Bop], triangleB->vertices[Bvc2]),
									getVertexIdByVertex(triangleB->next[Bop], triangleB->vertices[Bvc1]),
									getThirdVertexId(
										triangleB->next[Bop]->next[getVertexIdByVertex(triangleB->next[Bop], triangleB->vertices[Bvc1])],
										triangleB->next[Bop]->vertices[getThirdVertexId(triangleB->next[Bop], triangleB->vertices[Bvc1], triangleB->vertices[Bvc2])],
										triangleB->next[Bop]->vertices[getVertexIdByVertex(triangleB->next[Bop], triangleB->vertices[Bvc2])]
									), 0, 0);
							}
						}
					}
				}
			}
		}
		return 1;
	}
	return 0;
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
	fprintf(fpOutput, "\n");
	fclose(fpOutput);
}

void exportDataRestricted(Segment* constraints, int numTotalSegments, char* fileOutput) {
	FILE* fpOutput;
	int i;
	fpOutput = fopen(fileOutput, "a");
	for (i = 0; i < numTotalSegments; i++) {
		fprintf(fpOutput, "%f %f\n",
			constraints[i].v1.x,
			constraints[i].v1.y
		);
		fprintf(fpOutput, "%f %f\n\n",
			constraints[i].v2.x,
			constraints[i].v2.y
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

/*
test orientacion entre segmento ((x1, y1),(x2, y2)) y un punto (xp, yp)
*/
float getDetBySegments(float x1, float y1, float x2, float y2, float xp, float yp) {
	Mat matrix;
	float d;
	matrix.m = malloc(4 * sizeof(float));
	matrix.n = 2;
	matrix.m[0] = x2 - x1;
	matrix.m[1] = xp - x1;
	matrix.m[2] = y2 - y1;
	matrix.m[3] = yp - y1;
	d = det(matrix);
	return d;
}

Triangle* getIdTriangleContainsPoint(Triangle* triangle, int *puntoEnBorde, float x, float y) {
	float dets[3];
	Triangle* currentTriangle = triangle;
	*puntoEnBorde = 0;
	while(1) {
		getDetsByTriangle(currentTriangle, dets, x, y);
		if (dets[0] >= 0 && dets[1] >= 0 && dets[2] >= 0) {
			if (dets[0] == 0) *puntoEnBorde = 1;
			else if (dets[1] == 0) *puntoEnBorde = 2;
			else if (dets[2] == 0) *puntoEnBorde = 3;
			return currentTriangle;
		}
		if (dets[0] < 0) currentTriangle = currentTriangle->next[2];
		else if (dets[2] < 0) currentTriangle = currentTriangle->next[1];
		else if (dets[1] < 0) currentTriangle = currentTriangle->next[0];
	}
}

int generateDelaunayNet(char* strFileInput, Triangle *triangles, Vertex *vertices) {
	FILE* fpInput;
	int i, j,
		// id del tri�ngulo que contiene el nuevo punto
		id,
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
	// x e y son las coordenadas de los nuevos puntos leido una linea a la vez
	float x, y, dets[3];
	// Copia de triangulos de apoyo para no perder referencias originales al momento de modificar la malla
	Triangle *currentTriangle, copyTriangle, copyNext;

	fpInput = fopen(strFileInput, "r");
	if (!fpInput) {
		printf("No pude abrir el archivo puntos\n");
		return 1;
	}
	while (!feof(fpInput)) {
		fscanf(fpInput, "%f %f", &x, &y);
		// se chequea si el punto ya existe. En caso de existir, pasamos al siguiente
		// si no se hace esto �ltimo, se producen errores
		puntoExiste = 0;
		for (j = 0; j < numTotalVertices;j++) {
			if (vertices[j].x == x && vertices[j].y == y) {
				puntoExiste = 1;
			}
		}
		if (puntoExiste == 1) continue;
		puntoEnBorde = -1;
		currentTriangle = getIdTriangleContainsPoint(&triangles[numTotalTriangles/2], &puntoEnBorde, x, y);
		// si punto cae en arista, se ordenan v�rtices enumerados como: el primer
		// que comparte esa arista, el segundo y el que no (al que llamaremos opuesto)
		// todo siempre en sentido contrarreloj
		// Caso en que el punto cae en el primer borde de t (es decir, entre arista 1 y 2)
		if (puntoEnBorde == 1) {
			idsVerticesTriangulo[VERTICE_COMPARTIDO_1] = 0;
			idsVerticesTriangulo[VERTICE_COMPARTIDO_2] = 1;
			idsVerticesTriangulo[VERTICE_OPUESTO] = 2;
		}
		// Caso en que el punto cae en el segundo borde de t (es decir, entre arista 2 y 3)
		else if (puntoEnBorde == 2) {
			idsVerticesTriangulo[VERTICE_COMPARTIDO_1] = 1;
			idsVerticesTriangulo[VERTICE_COMPARTIDO_2] = 2;
			idsVerticesTriangulo[VERTICE_OPUESTO] = 0;
		}
		// Caso en que el punto cae en el tercer borde de t (es decir, entre arista 3 y 1)
		else if (puntoEnBorde == 3) {
			idsVerticesTriangulo[VERTICE_COMPARTIDO_1] = 2;
			idsVerticesTriangulo[VERTICE_COMPARTIDO_2] = 0;
			idsVerticesTriangulo[VERTICE_OPUESTO] = 1;
		}
		// Caso en que el punto cae dentro del tri�ngulo t
		else if (puntoEnBorde == 0) {
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
			copyTriangle = *currentTriangle;
			// se limpia el triangulo t que se reutilizar� con el nuevo tri�ngulo n1
			clearTriangle(currentTriangle);
			// se crea el nuevo v�rtice
			vertices[numTotalVertices].x = x;
			vertices[numTotalVertices].y = y;
			// n1, referido como currentTriangle
			currentTriangle->vertices[0] = copyTriangle.vertices[0];
			currentTriangle->vertices[1] = copyTriangle.vertices[1];
			currentTriangle->vertices[2] = &vertices[numTotalVertices];
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
			currentTriangle->next[0] = &triangles[numTotalTriangles - 2];
			// el vecino de n1, en su vertice 2, es n3
			currentTriangle->next[1] = &triangles[numTotalTriangles - 1];
			// el vecino de n1, en su vertice 3, es el ex vecino de t en su vertice 2
			currentTriangle->next[2] = copyTriangle.next[2];
			// se actualizan los vecinos locales de n2
			triangles[numTotalTriangles - 2].next[0] = copyTriangle.next[0];
			triangles[numTotalTriangles - 2].next[1] = &triangles[numTotalTriangles - 1];
			triangles[numTotalTriangles - 2].next[2] = currentTriangle;
			// se actualizan los vecinos locales de n3
			triangles[numTotalTriangles - 1].next[0] = &triangles[numTotalTriangles - 2];
			triangles[numTotalTriangles - 1].next[1] = copyTriangle.next[1];
			triangles[numTotalTriangles - 1].next[2] = currentTriangle;
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
				copyTriangle.next[2]->next[idVerticeOpuestoVecino] = currentTriangle;
			}
			/*
			*
				Desde ahora en adelante *NO* se debiese usar copyTriangle
			*
			*/
			/* si n1 tiene vecino al frente, debe testearse para el intercambio de diagonal */
			if (currentTriangle->next[2] != NULL) {
				idVerticeOpuestoVecino = getThirdVertexId(
					currentTriangle->next[2],
					currentTriangle->vertices[0],
					currentTriangle->vertices[1]);
				intercambioDeDiagonal(
					currentTriangle,
					currentTriangle->next[2], 0, 1, 2,
					idVerticeOpuestoVecino, 0, 0);
			}
			/* si n2 tiene vecino al frente, debe testearse para intercambio de diagonal */
			if (triangles[numTotalTriangles-2].next[0] != NULL) {
				idVerticeOpuestoVecino = getThirdVertexId(
					triangles[numTotalTriangles - 2].next[0],
					triangles[numTotalTriangles - 2].vertices[1],
					triangles[numTotalTriangles - 2].vertices[2]);
				intercambioDeDiagonal(
					&triangles[numTotalTriangles - 2],
					triangles[numTotalTriangles - 2].next[0], 1, 2, 0,
					idVerticeOpuestoVecino, 0, 0);
			}
			/* si n3 tiene vecino al frente, debe testearse para el intercambio de diagonal */
			if (triangles[numTotalTriangles - 1].next[1] != NULL) {
				idVerticeOpuestoVecino = getThirdVertexId(
					triangles[numTotalTriangles - 1].next[1],
					triangles[numTotalTriangles - 1].vertices[2],
					triangles[numTotalTriangles - 1].vertices[0]);
				intercambioDeDiagonal(
					&triangles[numTotalTriangles - 1],
					triangles[numTotalTriangles - 1].next[1], 2, 0, 1,
					idVerticeOpuestoVecino, 0, 0);
			}
			numTotalVertices++;
			continue;
		}
		// Caso en que el nuevo punto le�do cae en un v�rtice del tri�ngulo t
		if (puntoEnBorde > 0) {
			// se carga nuevo vertice a la malla
			vertices[numTotalVertices].x = x;
			vertices[numTotalVertices].y = y;
			// se saca una copia de respaldo del triangulo t
			// porque morir� pronto...
			copyTriangle = *currentTriangle;
			// muere t, quien se transformar� en el nuevo n1
			clearTriangle(currentTriangle);
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
				currentTriangle->vertices[0] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]];
				currentTriangle->vertices[1] = &vertices[numTotalVertices];
				currentTriangle->vertices[2] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_OPUESTO]];
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
				currentTriangle->next[0] = &triangles[numTotalTriangles - 2];
				// el ex vecino de t en el v�rtice 2 pasa a ser vecino de n1 en su vertice 1
				currentTriangle->next[1] = copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]];
				// de existir, se actualiza el ex vecino de t en el v�rtice 2
				if (copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]] != NULL) {
					idVerticeOpuestoVecino = getThirdVertexId(
						copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]],
						currentTriangle->vertices[0],
						currentTriangle->vertices[2]);
					copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]]->next[idVerticeOpuestoVecino] = currentTriangle;
				}
				// n3 pasa a ser el vecino de n1 en el vertice 3
				currentTriangle->next[2] = copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]];
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
				triangles[numTotalTriangles - 2].next[1] = currentTriangle;
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
				copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]]->next[2] = currentTriangle;
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
				/*
				*
					Desde ahora en adelante *NO* se debiese usar copyTriangle
				*
				*/
				// ahora se realizan los test del circulo para cada n
				// intercambio de diagonal para:
				// n1
				if (currentTriangle->next[1] != NULL) {
					idVerticeOpuestoVecino = getThirdVertexId(
						currentTriangle->next[1],
						currentTriangle->vertices[2],
						currentTriangle->vertices[0]);
					intercambioDeDiagonal(
						currentTriangle,
						currentTriangle->next[1], 2, 0, 1,
						idVerticeOpuestoVecino, 0, 0);
				}
				// n2
				if (triangles[numTotalTriangles-2].next[0] != NULL) {
					idVerticeOpuestoVecino = getThirdVertexId(
						triangles[numTotalTriangles - 2].next[0],
						triangles[numTotalTriangles - 2].vertices[1],
						triangles[numTotalTriangles - 2].vertices[2]);
					intercambioDeDiagonal(
						&triangles[numTotalTriangles - 2],
						triangles[numTotalTriangles - 2].next[0], 1, 2, 0,
						idVerticeOpuestoVecino, 0, 0);
				}
				// n3
				if (currentTriangle->next[2]->next[0] != NULL) {
					idVerticeOpuestoVecino = getThirdVertexId(
						currentTriangle->next[2]->next[0],
						currentTriangle->next[2]->vertices[1],
						currentTriangle->next[2]->vertices[2]);
					intercambioDeDiagonal(
						currentTriangle->next[2],
						currentTriangle->next[2]->next[0], 1, 2, 0,
						idVerticeOpuestoVecino, 0, 0);
				}
				// n4
				if (triangles[numTotalTriangles - 1].next[0] != NULL) {
					idVerticeOpuestoVecino = getThirdVertexId(
						triangles[numTotalTriangles - 1].next[0],
						triangles[numTotalTriangles - 1].vertices[1],
						triangles[numTotalTriangles - 1].vertices[2]);
					intercambioDeDiagonal(
						&triangles[numTotalTriangles - 1],
						triangles[numTotalTriangles - 1].next[0], 1, 2, 0,
						idVerticeOpuestoVecino, 0, 0);
				}
			}
			// en el caso que t no comparta la arista donde cay� el nuevo punto con otro trinagulo...
			// (en otras palabras, no tiene vecino en la arista donde cay� el nuevo punto...)
			else {
				// se crea el nuevo tri�ngulo n1
				currentTriangle->vertices[0] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]];
				currentTriangle->vertices[1] = &vertices[numTotalVertices];
				currentTriangle->vertices[2] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_OPUESTO]];
				// se crea el nuevo triangulo n2
				triangles[numTotalTriangles].vertices[0] = &vertices[numTotalVertices];
				triangles[numTotalTriangles].vertices[1] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]];
				triangles[numTotalTriangles].vertices[2] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_OPUESTO]];
				numTotalTriangles++;
				// se actualiza el vecindario
				// n2 pasa a ser vecino de n1 en el v�rtice 1
				currentTriangle->next[0] = &triangles[numTotalTriangles - 1];
				// se actualiza el vecino lejano de n1
				currentTriangle->next[1] = copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]];
				if (copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]] != NULL) {
					idVerticeOpuestoVecino = getThirdVertexId(
						copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]],
						copyTriangle.vertices[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]],
						copyTriangle.vertices[idsVerticesTriangulo[VERTICE_OPUESTO]]);
					copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]]->next[idVerticeOpuestoVecino] = currentTriangle;
				}
				currentTriangle->next[2] = NULL;
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
				triangles[numTotalTriangles - 1].next[1] = currentTriangle;
				triangles[numTotalTriangles - 1].next[2] = NULL;
				// intercambio de diagonal para n1, si es que tiene vecino al frente de su vertice 2
				if (currentTriangle->next[1] != NULL) {
					idVerticeOpuestoVecino = getThirdVertexId(
						currentTriangle->next[1],
						currentTriangle->vertices[2],
						currentTriangle->vertices[0]);
					intercambioDeDiagonal(
						currentTriangle,
						currentTriangle->next[1], 2, 0, 1,
						idVerticeOpuestoVecino, 0, 0);
				}
				// intercambio de diagonal para n2
				if (triangles[numTotalTriangles-1].next[0] != NULL) {
					idVerticeOpuestoVecino = getThirdVertexId(
						triangles[numTotalTriangles - 1].next[0],
						triangles[numTotalTriangles - 1].vertices[1],
						triangles[numTotalTriangles - 1].vertices[2]);
					intercambioDeDiagonal(
						&triangles[numTotalTriangles - 1],
						triangles[numTotalTriangles - 1].next[0], 1, 2, 0,
						idVerticeOpuestoVecino, 0, 0);
				}
			}
			numTotalVertices++;
		}
	}
	fclose(fpInput);
	return numTotalTriangles;
}

int restrictDelaunayNet(char* strFileInput, Triangle* triangles, int numTotalTriangles, Segment* constraints) {
	FILE* fpInput;
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
		numTotalSegments = 0,
		// se usa como buleano para determinar si el punto cae en borde o dentro de un
		// tri�ngulo
		constraintFound,
		segmentoExiste,
		done;
	// xs e ys es la recta a restringir
	float x1, y1, x2, y2, detsA[3], detsB[3] = { -1, -1, -1 },
		detArista1, detArista2;
	Triangle *currentTriangle, copyTriangle, copyNext;

	fpInput = fopen(strFileInput, "r");
	if (!fpInput) {
		printf("No pude abrir el archivo de puntos restringidos\n");
		return 1;
	}
	while (!feof(fpInput)) {
		fscanf(fpInput, "%f %f %f %f", &x1, &y1, &x2, &y2);
		segmentoExiste = 0;
		for (j = 0; j < numTotalSegments; j++) {
			if (constraints[j].v1.x == x1 &&
				constraints[j].v1.y == y1 &&
				constraints[j].v2.x == x2 &&
				constraints[j].v2.y == y2) {
				segmentoExiste = 1;
			}
		}
		if (segmentoExiste == 1) continue;
		constraints[numTotalSegments].v1.x = x1;
		constraints[numTotalSegments].v1.y = y1;
		constraints[numTotalSegments].v2.x = x2;
		constraints[numTotalSegments].v2.y = y2;
		numTotalSegments++;
		/* se busca el triangulo donde parte del segmento restringido */
		done = 0;
		while (done == 0) {
			for (i = 0; i < numTotalTriangles; i++) {
				constraintFound = 0;
				currentTriangle = &triangles[i];
				getDetsByTriangle(currentTriangle, detsA, x1, y1);
				if (detsA[0] >= 0 && detsA[1] >= 0 && detsA[2] >= 0) {
					getDetsByTriangle(currentTriangle, detsB, x2, y2);
					if (detsB[0] >= 0 && detsB[1] >= 0 && detsB[2] >= 0) {
						done = 1;
						break;
					}
					/* donde idsVerticesTriangulo[0] = 0 ser� el v�rtice opuesto
					al prox vecino del fin del segmento restringido */
					if (detsA[0] > 0 && detsA[1] == 0 && detsA[2] == 0 &&
						detsB[0] < 0 && detsB[1] >= 0 && detsB[2] >= 0) {
						constraintFound = 1;
						idsVerticesTriangulo[0] = 2;
						idsVerticesTriangulo[1] = 0;
						idsVerticesTriangulo[2] = 1;
					}
					else if (
						detsA[0] == 0 && detsA[1] == 0 && detsA[2] > 0 &&
						detsB[0] >= 0 && detsB[1] >= 0 && detsB[2] < 0) {
						constraintFound = 1;
						idsVerticesTriangulo[0] = 1;
						idsVerticesTriangulo[1] = 2;
						idsVerticesTriangulo[2] = 0;
					}
					else if (
						detsA[0] == 0 && detsA[1] > 0 && detsA[2] == 0 &&
						detsB[0] >= 0 && detsB[1] < 0 && detsB[2] >= 0) {
						constraintFound = 1;
						idsVerticesTriangulo[0] = 0;
						idsVerticesTriangulo[1] = 1;
						idsVerticesTriangulo[2] = 2;
					}
					if (constraintFound == 1) {
						/* caso en que arista restringida sea colinear con arista de currentTriangle */
						if (detsB[idsVerticesTriangulo[0]] == 0) {
							x1 = currentTriangle->vertices[idsVerticesTriangulo[1]]->x;
							y1 = currentTriangle->vertices[idsVerticesTriangulo[1]]->y;
							break;
						}
						else if (detsB[idsVerticesTriangulo[2]] == 0) {
							x1 = currentTriangle->vertices[idsVerticesTriangulo[2]]->x;
							y1 = currentTriangle->vertices[idsVerticesTriangulo[2]]->y;
							break;
						}
						idVerticeOpuestoVecino = getThirdVertexId(
							currentTriangle->next[idsVerticesTriangulo[0]],
							currentTriangle->vertices[idsVerticesTriangulo[1]],
							currentTriangle->vertices[idsVerticesTriangulo[2]]);
						/*
						hacer el test orientacion entre las aristas que no es la que comparte con el prox triangulo y
						triangles[i].next[idsVerticesTriangulo[0]]
						*/
						detArista1 = getDetBySegments(
							currentTriangle->vertices[idsVerticesTriangulo[0]]->x,
							currentTriangle->vertices[idsVerticesTriangulo[0]]->y,
							currentTriangle->vertices[idsVerticesTriangulo[1]]->x,
							currentTriangle->vertices[idsVerticesTriangulo[1]]->y,
							currentTriangle->next[idsVerticesTriangulo[0]]->vertices[idVerticeOpuestoVecino]->x,
							currentTriangle->next[idsVerticesTriangulo[0]]->vertices[idVerticeOpuestoVecino]->y);
						detArista2 = getDetBySegments(
							currentTriangle->vertices[idsVerticesTriangulo[2]]->x,
							currentTriangle->vertices[idsVerticesTriangulo[2]]->y,
							currentTriangle->vertices[idsVerticesTriangulo[0]]->x,
							currentTriangle->vertices[idsVerticesTriangulo[0]]->y,
							currentTriangle->next[idsVerticesTriangulo[0]]->vertices[idVerticeOpuestoVecino]->x,
							currentTriangle->next[idsVerticesTriangulo[0]]->vertices[idVerticeOpuestoVecino]->y);
						if (detArista1 > 0 && detArista2 > 0)
							intercambioDeDiagonal(currentTriangle, currentTriangle->next[idsVerticesTriangulo[0]],
								idsVerticesTriangulo[1], idsVerticesTriangulo[2], idsVerticesTriangulo[0], idVerticeOpuestoVecino, 1, 1);
						else if (detArista1 < 0) {
							/*
							intercambioDeDiagonal(
								currentTriangle->next[idsVerticesTriangulo[0]],
								currentTriangle->next[idsVerticesTriangulo[0]]->next[
									getVertexIdByVertex(
										currentTriangle->next[idsVerticesTriangulo[0]],
										currentTriangle->vertices[idsVerticesTriangulo[2]])
								],
								getVertexIdByVertex(currentTriangle->next[idsVerticesTriangulo[0]], currentTriangle->vertices[idsVerticesTriangulo[1]]),
								getThirdVertex(currentTriangle->next[idsVerticesTriangulo[0]],
									currentTriangle->vertices[idsVerticesTriangulo[1]],
									currentTriangle->vertices[idsVerticesTriangulo[2]]),
								getVertexIdByVertex(currentTriangle->next[idsVerticesTriangulo[0]], currentTriangle->vertices[idsVerticesTriangulo[2]]),
								getThirdVertex(
									currentTriangle->next[idsVerticesTriangulo[0]]->next[
										getVertexIdByVertex(
											currentTriangle->next[idsVerticesTriangulo[0]],
											currentTriangle->vertices[idsVerticesTriangulo[2]])],
										getVertexIdByVertex(currentTriangle->next[idsVerticesTriangulo[0]], currentTriangle->vertices[idsVerticesTriangulo[2]]),
										getThirdVertex(
												currentTriangle->next[idsVerticesTriangulo[0]],
												currentTriangle->vertices[idsVerticesTriangulo[1]],
												currentTriangle->vertices[idsVerticesTriangulo[2]])
								),
								1, 1
							);
							*/
						}
						else if (detArista2 < 0) {
							
						}
					}
				}
			}
		}
	}
	fclose(fpInput);
	return numTotalSegments;
}

int main(int argc, char* argv[]) {
	int numTotalTriangles, numTotalSegments;
	char fileInput[BUFFER_SIZE],
		fileOutput[BUFFER_SIZE],
		fileInputRestrict[BUFFER_SIZE];
	Triangle* triangles;
	Vertex* vertices;
	Segment* constraints;
	strcpy(fileInput, argv[1]);
	strcpy(fileOutput, argv[2]);
	strcpy(fileInputRestrict, argv[3]);
	triangles = malloc(sizeof(Triangle) * MAX_TRIANGLES);
	vertices = malloc(sizeof(Vertex) * MAX_VERTICES);
	constraints = malloc(sizeof(Segment) * MAX_SEGMENTS);
	if (triangles == NULL || vertices == NULL) {
		printf("Error malloc\n");
		exit(0);
	}
	initMesh(triangles, vertices, constraints);
	numTotalTriangles = generateDelaunayNet(fileInput, triangles, vertices);
	numTotalSegments = restrictDelaunayNet(fileInputRestrict, triangles, numTotalTriangles, constraints);
	exportData(triangles, numTotalTriangles, fileOutput);
	exportDataRestricted(constraints, numTotalSegments, fileOutput);
	return 0;
}
