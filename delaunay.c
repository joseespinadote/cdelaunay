#pragma once
#define _CRT_SECURE_NO_WARNINGS 1 
#define _WINSOCK_DEPRECATED_NO_WARNINGS 1 
#define _CRT_SECURE_NO_DEPRECATE 1
#define _CRT_NONSTDC_NO_DEPRECATE 1
#pragma warning (disable: 6262)
#pragma warning (disable: 6031)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "delaunay.h"
#include "det.h"

void initMesh(Triangle* triangles, Vertex* vertices, Segment* constraints) {
	int i, j;
	for (i = 0; i < MAX_VERTICES; i++) {
		vertices[i] = (Vertex){ -1, -1 };
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

/* Función que retorna id del vértice opuesto no compartido por el vecino dado un triángulo */
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
Intercambia diagonal entre triángulos A y B siempre que pase el test del círculo
Avc1 = es el ínidce del primer vertice de la arista compartida de A con B
Avc2 = es el ínidce del primer vertice de la arista compartida de A con B
Aop = es el ínidce del vertice opuesto a la arista compartida de A con B
Bvc1 = es el ínidce del primer vertice de la arista compartida de B con A
Bvc2 = es el ínidce del primer vertice de la arista compartida de B con A
Bop = es el ínidce del vertice opuesto a la arista compartida de B con A
*/
int intercambioDeDiagonal(Triangle* triangleA, Triangle* triangleB,
	int Avc1, int Avc2, int Aop, int Bop, int skipCircleTest, int skipPropagation) {
	if (triangleA == NULL || triangleB == NULL) return 0;
	/*
	esta función:
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
	*/
	int Bvc1, Bvc2, idVerticeVecinoLejano;
	Bvc1 = getVertexIdByVertex(triangleB, triangleA->vertices[Avc2]);
	Bvc2 = getVertexIdByVertex(triangleB, triangleA->vertices[Avc1]);
	// se hace el test del círculo
	if (skipCircleTest == 1 ||
		(circleTest(triangleA, triangleB->vertices[Bop]->x, triangleB->vertices[Bop]->y) > EPSILON)
		) {
		// si el resultado es positivo, se hace el cambio de diagonal, que,
		// practicamente, es cambiar los vértices Avc2 y Bvc1 de A y B respectivamente
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
				idVerticeVecinoLejano = getThirdVertexId(
					triangleA->next[Aop],
					triangleA->vertices[Avc1],
					triangleA->vertices[Avc2]);
				intercambioDeDiagonal(
					triangleA,
					triangleA->next[Aop],
					Avc1,
					Avc2,
					Aop,
					idVerticeVecinoLejano, 0, 0);
			}
			if (triangleB->next[Bop] != NULL) {
				idVerticeVecinoLejano = getThirdVertexId(
					triangleB->next[Bop],
					triangleB->vertices[Bvc1],
					triangleB->vertices[Bvc2]);
				intercambioDeDiagonal(
					triangleB,
					triangleB->next[Bop],
					Bvc1,
					Bvc2,
					Bop,
					idVerticeVecinoLejano, 0, 0);
			}
		}
		return 1;
	}
	return 0;
}

void exportData(Triangle* triangles, int numTotalTriangles, char* fileOutput) {
	FILE *fpOutput;
	int i;
	fpOutput = fopen(fileOutput, "w");
	for (i = 0; i < numTotalTriangles; i++) {
		if (triangles[i].vertices[0] == NULL) continue;
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

void clearTriangle(Triangle* triangle) {
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

Triangle* getIdTriangleContainsPoint(Triangle* triangle, int* puntoEnBorde, float x, float y) {
	float dets[3];
	Triangle* currentTriangle = triangle;
	*puntoEnBorde = 0;
	while (1) {
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

int isPointAlreadyOnNet(Vertex* vertices, int numTotalVertices, float x, float y) {
	int i = 0;
	for (i = 0; i < numTotalVertices; i++) {
		if (vertices[i].x == x && vertices[i].y == y) {
			return 1;
		}
	}
	return 0;
}

void addPointToMesh(Triangle* triangles, Vertex* vertices, int* numTotalTriangles, int* numTotalVertices, float x, float y) {
		// idsVerticesTriangulo es un vector con los vertices ordenados contrarreloj
		// que se crea en función de la arista donde cae el nuevo punto
	int idsVerticesTriangulo[3] = { -1, -1, -1 },
		// idVerticeOpuestoVecino es el id del vértice opuesto del vecino que no comparte
		// con alguno de los nuevos triángulos n1, n2 o n3. Sirve para hacer test del
		// círculo
		idVerticeOpuestoVecino,
		// similar a idsVerticesTriangulo, pero del vecino con el que comparte arista en caso
		// de que el punto cae en arista compartida entre 2 triángulos
		idsVerticesVecino[3] = { -1,-1,-1 },
		// se usa como buleano para determinar si el punto cae en borde o dentro de un
		// triángulo
		puntoEnBorde;
	// Copia de triangulos de apoyo para no perder referencias originales al momento de modificar la malla
	Triangle* currentTriangle, copyTriangle, copyNext;
	/*
	se chequea si el punto ya existe. En caso de existir, pasamos al siguiente
	si no se hace esto último, se producen errores
	*/
	if (isPointAlreadyOnNet(vertices, *numTotalVertices, x, y) == 1) return;
	// si punto cae en arista, se ordenan vértices enumerados como: el primer
	// que comparte esa arista, el segundo y el que no (al que llamaremos opuesto)
	// todo siempre en sentido contrarreloj
	puntoEnBorde = -1;
	currentTriangle = getIdTriangleContainsPoint(&triangles[*numTotalTriangles / 2], &puntoEnBorde, x, y);
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
	// Caso en que el punto cae dentro del triángulo t
	else if (puntoEnBorde == 0) {
		/*
		- Se crearán los triángulo n1, n2 y n3 dentro de t, con sus aristas y el nuevo vertice
		- Se deben configurar la nueva vecindad dado n1, n2 y n3
		- Los ex vecinos del difunto t, deben actualizarse respecto a n1, n2 y n3
		- Se debe hacer el test del círculo entre:
		- n1 y V3, donde V3 es el vecino opuesto al vértice 3 de t
		- n2 y V1, donde V1 es el vecino opuesto al vértice 1 de t
		- n3 y V2, donde V2 es el vecino opuesto al vértice 2 de t
		*/
		// se guarda un respaldo del triagulo original t donde cae el nuevo punto
		copyTriangle = *currentTriangle;
		// se limpia el triangulo t que se reutilizará con el nuevo triángulo n1
		clearTriangle(currentTriangle);
		// se crea el nuevo vértice
		vertices[*numTotalVertices].x = x;
		vertices[*numTotalVertices].y = y;
		// n1, referido como currentTriangle
		currentTriangle->vertices[0] = copyTriangle.vertices[0];
		currentTriangle->vertices[1] = copyTriangle.vertices[1];
		currentTriangle->vertices[2] = &vertices[*numTotalVertices];
		// n2, referido, más adelante, como triangles[*numTotalTriangles - 2]
		triangles[*numTotalTriangles].vertices[0] = &vertices[*numTotalVertices];
		triangles[*numTotalTriangles].vertices[1] = copyTriangle.vertices[1];
		triangles[*numTotalTriangles].vertices[2] = copyTriangle.vertices[2];
		(*numTotalTriangles)++;
		//  n3, referido, más adelante, como triangles[*numTotalTriangles - 1]
		triangles[*numTotalTriangles].vertices[0] = copyTriangle.vertices[0];
		triangles[*numTotalTriangles].vertices[1] = &vertices[*numTotalVertices];
		triangles[*numTotalTriangles].vertices[2] = copyTriangle.vertices[2];
		(*numTotalTriangles)++;
		// se actualizan los vecinos locales de n1:
		// el vecino de n1, en su vertice 1, es n2
		currentTriangle->next[0] = &triangles[*numTotalTriangles - 2];
		// el vecino de n1, en su vertice 2, es n3
		currentTriangle->next[1] = &triangles[*numTotalTriangles - 1];
		// el vecino de n1, en su vertice 3, es el ex vecino de t en su vertice 2
		currentTriangle->next[2] = copyTriangle.next[2];
		// se actualizan los vecinos locales de n2
		triangles[*numTotalTriangles - 2].next[0] = copyTriangle.next[0];
		triangles[*numTotalTriangles - 2].next[1] = &triangles[*numTotalTriangles - 1];
		triangles[*numTotalTriangles - 2].next[2] = currentTriangle;
		// se actualizan los vecinos locales de n3
		triangles[*numTotalTriangles - 1].next[0] = &triangles[*numTotalTriangles - 2];
		triangles[*numTotalTriangles - 1].next[1] = copyTriangle.next[1];
		triangles[*numTotalTriangles - 1].next[2] = currentTriangle;
		// ¿tenia vecino V1 el difunto t frente al vertice 1?, o, en otras palabras,
		// ¿hay vecino frente al vertice 1 de n2 que se requiera actualizar la referencia?
		if (copyTriangle.next[0] != NULL) {
			// se rescata id del vecino opuesto para actualizar referencia de vecino lejano
			idVerticeOpuestoVecino = getThirdVertexId(
				copyTriangle.next[0],
				copyTriangle.vertices[1],
				copyTriangle.vertices[2]);
			// se actualiza la referencia del vecino respecto a nuevo triángulo n2
			copyTriangle.next[0]->next[idVerticeOpuestoVecino] = &triangles[*numTotalTriangles - 2];
		}
		// ¿tenia vecino V2 el difunto t frente al vertice 2? (o frente a n3)
		if (copyTriangle.next[1] != NULL) {
			idVerticeOpuestoVecino = getThirdVertexId(
				copyTriangle.next[1],
				copyTriangle.vertices[2],
				copyTriangle.vertices[0]);
			copyTriangle.next[1]->next[idVerticeOpuestoVecino] = &triangles[*numTotalTriangles - 1];
		}
		// ¿tenia vecino V3 el difunto t frente al vertice 3? (o frente a n1)
		if (copyTriangle.next[2] != NULL) {
			// idsVerticesVecino[0] es ahora vértice opuesto a la arista compartida entre n1 y V3
			idVerticeOpuestoVecino = getThirdVertexId(
				copyTriangle.next[2],
				copyTriangle.vertices[0],
				copyTriangle.vertices[1]);
			// se actualiza referencia del vértice opuesto a la arista compartida en el vecino lejano V3
			copyTriangle.next[2]->next[idVerticeOpuestoVecino] = currentTriangle;
		}
		/*
		Desde ahora en adelante *NO* se debiese usar copyTriangle

		si n1 tiene vecino al frente, debe testearse para el intercambio de diagonal */
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
		if (triangles[*numTotalTriangles - 2].next[0] != NULL) {
			idVerticeOpuestoVecino = getThirdVertexId(
				triangles[*numTotalTriangles - 2].next[0],
				triangles[*numTotalTriangles - 2].vertices[1],
				triangles[*numTotalTriangles - 2].vertices[2]);
			intercambioDeDiagonal(
				&triangles[*numTotalTriangles - 2],
				triangles[*numTotalTriangles - 2].next[0], 1, 2, 0,
				idVerticeOpuestoVecino, 0, 0);
		}
		/* si n3 tiene vecino al frente, debe testearse para el intercambio de diagonal */
		if (triangles[*numTotalTriangles - 1].next[1] != NULL) {
			idVerticeOpuestoVecino = getThirdVertexId(
				triangles[*numTotalTriangles - 1].next[1],
				triangles[*numTotalTriangles - 1].vertices[2],
				triangles[*numTotalTriangles - 1].vertices[0]);
			intercambioDeDiagonal(
				&triangles[*numTotalTriangles - 1],
				triangles[*numTotalTriangles - 1].next[1], 2, 0, 1,
				idVerticeOpuestoVecino, 0, 0);
		}
		(*numTotalVertices)++;
		return;
	}
	// Caso en que el nuevo punto leído cae en un vértice del triángulo t
	if (puntoEnBorde > 0) {
		// se carga nuevo vertice a la malla
		vertices[*numTotalVertices].x = x;
		vertices[*numTotalVertices].y = y;
		// se saca una copia de respaldo del triangulo t
		// porque morirá pronto...
		copyTriangle = *currentTriangle;
		// muere t, quien se transformará en el nuevo n1
		clearTriangle(currentTriangle);
		// ¿tenía el difunto triangulo t tenía un vecino en el vértice
		// opuesto a la arista donde cae el nuevo punto?
		// (y por ende, ¿vecino de t con el que comparte la arista donde cae el punto?)
		if (copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]] != NULL) {
			// se saca una copia de respaldo del vecino con el que comparte artista donde cae el punto
			copyNext = *copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]];
			// muere vecino de la arista compartida donde cae el punto
			clearTriangle(copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]]);
			// se arma set de indices del vecino, asegurándose que se conserve el ordenamiento contrarreloj
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
			// se crea nuevo triángulo n1
			currentTriangle->vertices[0] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]];
			currentTriangle->vertices[1] = &vertices[*numTotalVertices];
			currentTriangle->vertices[2] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_OPUESTO]];
			// se crea nuevo triángulo n2
			triangles[*numTotalTriangles].vertices[0] = &vertices[*numTotalVertices];
			triangles[*numTotalTriangles].vertices[1] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]];
			triangles[*numTotalTriangles].vertices[2] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_OPUESTO]];
			(*numTotalTriangles)++;
			// se crea nuevo triángulo n3
			copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]]->vertices[0] =
				&vertices[*numTotalVertices];
			copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]]->vertices[1] =
				copyNext.vertices[idsVerticesVecino[VERTICE_COMPARTIDO_2]];
			copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]]->vertices[2] =
				copyNext.vertices[idsVerticesVecino[VERTICE_OPUESTO]];
			// se crea nuevo triángulo n4
			triangles[*numTotalTriangles].vertices[0] = &vertices[*numTotalVertices];
			triangles[*numTotalTriangles].vertices[1] = copyNext.vertices[idsVerticesVecino[VERTICE_OPUESTO]];
			triangles[*numTotalTriangles].vertices[2] = copyNext.vertices[idsVerticesVecino[VERTICE_COMPARTIDO_1]];
			(*numTotalTriangles)++;
			// se actualiza el vecindario
			// n2 pasa a ser vecino de n1 en su vertice 1 (n1 v0)
			currentTriangle->next[0] = &triangles[*numTotalTriangles - 2];
			// el ex vecino de t en el vértice 2 pasa a ser vecino de n1 en su vertice 1
			currentTriangle->next[1] = copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]];
			// de existir, se actualiza el ex vecino de t en el vértice 2
			if (copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]] != NULL) {
				idVerticeOpuestoVecino = getThirdVertexId(
					copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]],
					currentTriangle->vertices[0],
					currentTriangle->vertices[2]);
				copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]]->next[idVerticeOpuestoVecino] = currentTriangle;
			}
			// n3 pasa a ser el vecino de n1 en el vertice 3
			currentTriangle->next[2] = copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]];
			// el ex vecino de t en el vértice 1 pasa a ser vecino de n2 en su vertice 0
			triangles[*numTotalTriangles - 2].next[0] = copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]];
			// de existir, se actualiza el ex vecino de t en el vértice 1
			if (copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]] != NULL) {
				idVerticeOpuestoVecino = getThirdVertexId(
					copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]],
					triangles[*numTotalTriangles - 2].vertices[1],
					triangles[*numTotalTriangles - 2].vertices[2]);
				copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]]->next[idVerticeOpuestoVecino] =
					&triangles[*numTotalTriangles - 2];
			}
			// n1 pasa a ser vecino de n2 en su vertice 2
			triangles[*numTotalTriangles - 2].next[1] = currentTriangle;
			// n4 pasa a ser vecino de n2 en su vertice 3
			triangles[*numTotalTriangles - 2].next[2] = &triangles[*numTotalTriangles - 1];
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
			copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]]->next[1] = &triangles[*numTotalTriangles - 1];
			copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]]->next[2] = currentTriangle;
			// vecindad de n4
			triangles[*numTotalTriangles - 1].next[0] = copyNext.next[idsVerticesVecino[VERTICE_COMPARTIDO_2]];
			// de existir, se actualiza el ex vecino de copyNext
			if (copyNext.next[idsVerticesVecino[VERTICE_COMPARTIDO_2]] != NULL) {
				idVerticeOpuestoVecino = getThirdVertexId(
					copyNext.next[idsVerticesVecino[VERTICE_COMPARTIDO_2]],
					triangles[*numTotalTriangles - 1].vertices[1],
					triangles[*numTotalTriangles - 1].vertices[2]);
				copyNext.next[idsVerticesVecino[VERTICE_COMPARTIDO_2]]->next[idVerticeOpuestoVecino] = &triangles[*numTotalTriangles - 1];
			}
			// n2 pasa a ser vecino de n4 en el vertice 1
			triangles[*numTotalTriangles - 1].next[1] = &triangles[*numTotalTriangles - 2];
			triangles[*numTotalTriangles - 1].next[2] = copyTriangle.next[idsVerticesTriangulo[VERTICE_OPUESTO]];
			/*
			Desde ahora en adelante *NO* se debiese usar copyTriangle

			ahora se realizan los test del circulo para cada n
			intercambio de diagonal para:
			n1
			*/
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
			if (triangles[*numTotalTriangles - 2].next[0] != NULL) {
				idVerticeOpuestoVecino = getThirdVertexId(
					triangles[*numTotalTriangles - 2].next[0],
					triangles[*numTotalTriangles - 2].vertices[1],
					triangles[*numTotalTriangles - 2].vertices[2]);
				intercambioDeDiagonal(
					&triangles[*numTotalTriangles - 2],
					triangles[*numTotalTriangles - 2].next[0], 1, 2, 0,
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
			if (triangles[*numTotalTriangles - 1].next[0] != NULL) {
				idVerticeOpuestoVecino = getThirdVertexId(
					triangles[*numTotalTriangles - 1].next[0],
					triangles[*numTotalTriangles - 1].vertices[1],
					triangles[*numTotalTriangles - 1].vertices[2]);
				intercambioDeDiagonal(
					&triangles[*numTotalTriangles - 1],
					triangles[*numTotalTriangles - 1].next[0], 1, 2, 0,
					idVerticeOpuestoVecino, 0, 0);
			}
		}
		// en el caso que t no comparta la arista donde cayó el nuevo punto con otro trinagulo...
		// (en otras palabras, no tiene vecino en la arista donde cayó el nuevo punto...)
		else {
			// se crea el nuevo triángulo n1
			currentTriangle->vertices[0] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]];
			currentTriangle->vertices[1] = &vertices[*numTotalVertices];
			currentTriangle->vertices[2] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_OPUESTO]];
			// se crea el nuevo triangulo n2
			triangles[*numTotalTriangles].vertices[0] = &vertices[*numTotalVertices];
			triangles[*numTotalTriangles].vertices[1] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]];
			triangles[*numTotalTriangles].vertices[2] = copyTriangle.vertices[idsVerticesTriangulo[VERTICE_OPUESTO]];
			(*numTotalTriangles)++;
			// se actualiza el vecindario
			// n2 pasa a ser vecino de n1 en el vértice 1
			currentTriangle->next[0] = &triangles[*numTotalTriangles - 1];
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
			triangles[*numTotalTriangles - 1].next[0] = copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]];
			if (copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]] != NULL) {
				idVerticeOpuestoVecino = getThirdVertexId(
					copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]],
					copyTriangle.vertices[idsVerticesTriangulo[VERTICE_COMPARTIDO_2]],
					copyTriangle.vertices[idsVerticesTriangulo[VERTICE_OPUESTO]]);
				copyTriangle.next[idsVerticesTriangulo[VERTICE_COMPARTIDO_1]]->next[idVerticeOpuestoVecino] = &triangles[*numTotalTriangles - 1];
			}
			// n1 pasa a ser vecino de n2 en el vértice 2
			triangles[*numTotalTriangles - 1].next[1] = currentTriangle;
			triangles[*numTotalTriangles - 1].next[2] = NULL;
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
			if (triangles[*numTotalTriangles - 1].next[0] != NULL) {
				idVerticeOpuestoVecino = getThirdVertexId(
					triangles[*numTotalTriangles - 1].next[0],
					triangles[*numTotalTriangles - 1].vertices[1],
					triangles[*numTotalTriangles - 1].vertices[2]);
				intercambioDeDiagonal(
					&triangles[*numTotalTriangles - 1],
					triangles[*numTotalTriangles - 1].next[0], 1, 2, 0,
					idVerticeOpuestoVecino, 0, 0);
			}
		}
		(*numTotalVertices)++;
	}
}

int generateDelaunayNet(char* strFileInput, Triangle* triangles, Vertex* vertices, int* numTotalTriangles, int* numTotalVertices) {
	/* coordenada de punto candidato a ser vertice */
	float x, y;
	FILE* fpInput;
	/* valores iniciales de contadores de trinagulos y vértcies totales en la malla */
	*numTotalTriangles = 2;
	*numTotalVertices = 4,
		/* se abre el archivo que contiene los puntos */
		fpInput = fopen(strFileInput, "r");
	if (!fpInput) {
		printf("No pude abrir el archivo puntos\n");
		return 1;
	}
	while (!feof(fpInput)) {
		fscanf(fpInput, "%f %f", &x, &y);
		/* se intenta añadir el nuevo punto a la malla */
		addPointToMesh(triangles, vertices, numTotalTriangles, numTotalVertices, x, y);
	}
	fclose(fpInput);
	return *numTotalTriangles;
}

void applyConstraint(Triangle* triangles, int *numTotalTriangles, Segment *constraint) {
	Triangle *currentTriangle, *oldTriangle;
	float detsA[3] = { -1, -1, -1 },
		detsB[3] = { -1, -1, -1 },
		detArista1, detArista2,
		x1, y1;
	int i, done = 0,
		idsVerticesTriangulo[3] = { -1, -1, -1 },
		idVerticeOpuestoVecino, constraintFound;
	x1 = constraint->v1.x;
	y1 = constraint->v1.y;
	while (done == 0) {
		for (i = 0; i < *numTotalTriangles; i++) {
			constraintFound = 0;
			currentTriangle = &triangles[i];
			getDetsByTriangle(currentTriangle, detsA, x1, y1);
			if (detsA[0] >= 0 && detsA[1] >= 0 && detsA[2] >= 0) {
				getDetsByTriangle(currentTriangle, detsB, constraint->v2.x, constraint->v2.y);
				if (detsB[0] >= 0 && detsB[1] >= 0 && detsB[2] >= 0) {
					done = 1;
					break;
				}
				/*
				idsVerticesTriangulo[0] será siempre el vértice opuesto
				al prox vecino del fin del segmento restringido
				*/
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
					/*
					caso en que arista restringida sea colinear con arista de currentTriangle
					*/
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
					else if (detArista1 <= 0) {
						oldTriangle = currentTriangle;
						currentTriangle = currentTriangle->next[idsVerticesTriangulo[0]];
						idsVerticesTriangulo[0] = getVertexIdByVertex(currentTriangle, oldTriangle->vertices[idsVerticesTriangulo[1]]);
						idsVerticesTriangulo[1] = idVerticeOpuestoVecino;
						idsVerticesTriangulo[2] = getVertexIdByVertex(currentTriangle, oldTriangle->vertices[idsVerticesTriangulo[2]]);
						idVerticeOpuestoVecino = getThirdVertexId(
							currentTriangle->next[idsVerticesTriangulo[0]],
							currentTriangle->vertices[idsVerticesTriangulo[1]],
							currentTriangle->vertices[idsVerticesTriangulo[2]]);
						intercambioDeDiagonal(
							currentTriangle,
							currentTriangle->next[idsVerticesTriangulo[0]],
							idsVerticesTriangulo[1],
							idsVerticesTriangulo[2],
							idsVerticesTriangulo[0],
							idVerticeOpuestoVecino,
							1, 1
						);
						break;
					}
					else if (detArista2 <= 0) {
						oldTriangle = currentTriangle;
						currentTriangle = currentTriangle->next[idsVerticesTriangulo[0]];
						idsVerticesTriangulo[0] = getVertexIdByVertex(currentTriangle, oldTriangle->vertices[idsVerticesTriangulo[2]]);
						idsVerticesTriangulo[1] = getVertexIdByVertex(currentTriangle, oldTriangle->vertices[idsVerticesTriangulo[1]]);
						idsVerticesTriangulo[2] = idVerticeOpuestoVecino;
						idVerticeOpuestoVecino = getThirdVertexId(
							currentTriangle->next[idsVerticesTriangulo[0]],
							currentTriangle->vertices[idsVerticesTriangulo[1]],
							currentTriangle->vertices[idsVerticesTriangulo[2]]);
						intercambioDeDiagonal(
							currentTriangle,
							currentTriangle->next[idsVerticesTriangulo[0]],
							idsVerticesTriangulo[1],
							idsVerticesTriangulo[2],
							idsVerticesTriangulo[0],
							idVerticeOpuestoVecino,
							1, 1
						);
						break;
					}
				}
			}
		}
	}
}

void loadConstraintFromM2DFile(char* strFileInput, Segment* constraints, int* numTotalSegments) {
	FILE* fpInput;
	char buffer[BUFFER_SIZE];
	float x1, y1, x2, y2;
	int i, numVertices=0, segmentoExiste,
		idsVerticesTrianglo[3] = { -1, -1, -1 };
	Vertex vertices[BUFFER_SIZE];
	
	for (i = 0; i < BUFFER_SIZE; i++) {
		vertices[i] = (Vertex){ -1, -1 };
	}
	fpInput = fopen(strFileInput, "r");
	if (!fpInput) {
		printf("No pude abrir el archivo de puntos restringidos\n");
		return 1;
	}
	*numTotalSegments = 0;
	while (!feof(fpInput)) {
		if (fgets(buffer, BUFFER_SIZE, fpInput) != NULL) {
			switch (buffer[0]) {
				case 'v' :
					sscanf(buffer, "%*c %*d %f %f",&vertices[numVertices].x, &vertices[numVertices].y);
					numVertices++;
					break;
				case 't' :
					sscanf(buffer, "%*c %*d %d %d %d",
						&idsVerticesTrianglo[0],
						&idsVerticesTrianglo[1],
						&idsVerticesTrianglo[2]);
					constraints[*numTotalSegments].v1 = vertices[idsVerticesTrianglo[0] - 1];
					constraints[*numTotalSegments].v2 = vertices[idsVerticesTrianglo[1] - 1];
					(*numTotalSegments)++;
					constraints[*numTotalSegments].v1 = vertices[idsVerticesTrianglo[1] - 1];
					constraints[*numTotalSegments].v2 = vertices[idsVerticesTrianglo[2] - 1];
					(*numTotalSegments)++;
					constraints[*numTotalSegments].v1 = vertices[idsVerticesTrianglo[2] - 1];
					constraints[*numTotalSegments].v2 = vertices[idsVerticesTrianglo[0] - 1];
					(*numTotalSegments)++;
			}
		}
	}
	fclose(fpInput);
	return *numTotalSegments;
}

void loadConstraintFromFile(char* strFileInput, Segment* constraints, int* numTotalSegments) {
	FILE* fpInput;
	float x1, y1, x2, y2;
	int i, segmentoExiste;

	fpInput = fopen(strFileInput, "r");
	if (!fpInput) {
		printf("No pude abrir el archivo de puntos restringidos\n");
		return 1;
	}
	while (!feof(fpInput)) {
		fscanf(fpInput, "%f %f %f %f", &x1, &y1, &x2, &y2);
		segmentoExiste = 0;
		for (i = 0; i < *numTotalSegments; i++) {
			if (constraints[i].v1.x == x1 &&
				constraints[i].v1.y == y1 &&
				constraints[i].v2.x == x2 &&
				constraints[i].v2.y == y2) {
				segmentoExiste = 1;
			}
		}
		if (segmentoExiste == 1) continue;
		/* se añade segmento a vector de segmentos*/
		constraints[*numTotalSegments].v1.x = x1;
		constraints[*numTotalSegments].v1.y = y1;
		constraints[*numTotalSegments].v2.x = x2;
		constraints[*numTotalSegments].v2.y = y2;
		(*numTotalSegments)++;
	}
	fclose(fpInput);
	return *numTotalSegments;
}

void clearPolygon(Triangle* triangles, int* numTotalTriangles, Segment* constraints, int* numTotalSegments) {
	int i, j, constraintFound, idConstraint,
		idsVerticesTriangulo[3] = { -1, -1, -1 },
		detArista, firstRemoval=0;
	float detsA[3] = { -1, -1, -1 },
		detsB[3] = { -1, -1, -1 };
	Triangle* currentTriangle, *island=NULL;

	for (i = 0; i < *numTotalSegments; i++) {
		for (j = 0; j < *numTotalTriangles; j++) {
			if (triangles[j].vertices[0] == NULL) continue;
			constraintFound = 0;
			currentTriangle = &triangles[j];
			getDetsByTriangle(currentTriangle, detsA, constraints[i].v1.x, constraints[i].v1.y);
			getDetsByTriangle(currentTriangle, detsB, constraints[i].v2.x, constraints[i].v2.y);
			if (detsA[0] >= 0 && detsA[1] >= 0 && detsA[2] >= 0) {
				// idsVerticesTriangulo[0] corresponde al vertice que comparte con v1 del
				// segmento restringido
				if (detsA[0] == 0 && detsA[1] == 0) {
					idsVerticesTriangulo[0] = 1;
					idsVerticesTriangulo[1] = 2;
					idsVerticesTriangulo[2] = 0;
					constraintFound = 1;
				}
				else if (detsA[1] == 0 && detsA[2] == 0) {
					idsVerticesTriangulo[0] = 2;
					idsVerticesTriangulo[1] = 0;
					idsVerticesTriangulo[2] = 1;
					constraintFound = 1;
				}
				else if (detsA[2] == 0 && detsA[0] == 0) {
					idsVerticesTriangulo[0] = 0;
					idsVerticesTriangulo[1] = 1;
					idsVerticesTriangulo[2] = 2;
					constraintFound = 1;
				}
				if (constraintFound == 1) {
					// se corrobora que el triangulo no esté involucrado en segmento anterior
					idConstraint = (i + *numTotalSegments - 1) % *numTotalSegments;
					if ((currentTriangle->vertices[idsVerticesTriangulo[1]]->x == constraints[idConstraint].v1.x &&
						currentTriangle->vertices[idsVerticesTriangulo[1]]->y == constraints[idConstraint].v1.y) ||
						(currentTriangle->vertices[idsVerticesTriangulo[2]]->x == constraints[idConstraint].v1.x &&
							currentTriangle->vertices[idsVerticesTriangulo[2]]->y == constraints[idConstraint].v1.y))
						continue;
					detArista = getDetBySegments(
						constraints[i].v1.x,
						constraints[i].v1.y,
						constraints[i].v2.x,
						constraints[i].v2.y,
						currentTriangle->vertices[idsVerticesTriangulo[2]]->x,
						currentTriangle->vertices[idsVerticesTriangulo[2]]->y);
					if (detArista > 0)
						continue;
					if (firstRemoval == 0) {
						firstRemoval = 1;
						island = currentTriangle->next[idsVerticesTriangulo[0]];
					}
					removeTriangle(currentTriangle, triangles, numTotalTriangles);
					j = 0;
				}
			}
		}
	}
	if (firstRemoval == 1 && island != NULL) {
		removeIsland(island);
	}
}

void restrictDelaunayNet(
	Triangle* triangles,
	Vertex* vertices,
	Segment* constraints,
	int* numTotalTriangles,
	int* numTotalVertices,
	int* numTotalSegments) {
	int i, j, k,
		idConstraint,
		/* idsVerticesTriangulo es un vector con los vertices ordenados contrarreloj
		que se crea en función de la arista donde cae el nuevo punto */
		idsVerticesTriangulo[3] = { -1, -1, -1 },
		/* similar a idsVerticesTriangulo, pero del vecino con el que comparte arista en caso
		de que el punto cae en arista compartida entre 2 triángulos */
		idsVerticesVecino[3] = { -1, -1, -1 },
		/* se usa como buleano para determinar si el punto cae en borde o dentro de un
		triángulo */
		constraintFound,
		segmentoExiste,
		firstRemoval,
		polygon;
	/* xs e ys es la recta a restringir */
	float x1, y1, x2, y2, detsA[3], detsB[3] = { -1, -1, -1 },
		detArista;
	Triangle *currentTriangle, *island=NULL;

	for (i = 0; i < *numTotalSegments; i++) {
		polygon = 0;
		firstRemoval = 0;
		if (i > 0 &&
			constraints[0].v1.x == constraints[i].v2.x &&
			constraints[0].v1.y == constraints[i].v2.y) {
			polygon = 1;
		}
		if (i < 2) {
			/* se añaden vértices de inicio y fin de segmento restringido a la malla */
			addPointToMesh(triangles, vertices, numTotalTriangles, numTotalVertices, constraints[i].v1.x, constraints[i].v1.y);
			addPointToMesh(triangles, vertices, numTotalTriangles, numTotalVertices, constraints[i].v2.x, constraints[i].v2.y);
			/* se busca el triangulo donde parte del segmento restringido */
			applyConstraint(triangles, numTotalTriangles, &constraints[i]);
		}
		if (polygon == 1 && 0) {
			clearPolygon(triangles, numTotalTriangles, constraints, numTotalSegments);
		}
	}
}

void removeIsland(Triangle* triangle) {
	if (triangle == NULL) return;
	Triangle copy = *triangle;
	clearTriangle(triangle);
	removeIsland(copy.next[0]);
	removeIsland(copy.next[1]);
	removeIsland(copy.next[2]);
}

void removeTriangle(Triangle* triangle, Triangle* triangles, int *numTotalTriangles) {
	int i;
	for (i = 0; i < *numTotalTriangles; i++) {
		if (triangles[i].next[0] == triangle) triangles[i].next[0] = NULL;
		if (triangles[i].next[1] == triangle) triangles[i].next[1] = NULL;
		if (triangles[i].next[2] == triangle) triangles[i].next[2] = NULL;
	}
	clearTriangle(triangle);
}