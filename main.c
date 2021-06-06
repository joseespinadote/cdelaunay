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
#include "delaunay.h"

int main(int argc, char* argv[]) {
	int *numTotalTriangles,
		*numTotalVertices,
		*numTotalSegments,
		clearMode = 1;
	char fileInput[BUFFER_SIZE],
		fileOutput[BUFFER_SIZE],
		fileInputRestrict[BUFFER_SIZE];
	Triangle* triangles;
	Vertex* vertices;
	Segment* constraints;
	numTotalTriangles = calloc(1, sizeof(int));
	numTotalVertices = calloc(1, sizeof(int));
	numTotalSegments = calloc(1, sizeof(int));
	strcpy(fileInput, argv[1]);
	strcpy(fileOutput, argv[2]);
	strcpy(fileInputRestrict, argv[3]);
	triangles = malloc(sizeof(Triangle) * MAX_TRIANGLES);
	vertices = malloc(sizeof(Vertex) * MAX_VERTICES);
	constraints = malloc(sizeof(Segment) * MAX_SEGMENTS);
	if (triangles == NULL || vertices == NULL || numTotalTriangles == NULL || numTotalVertices == NULL || numTotalSegments == NULL) {
		printf("mem error\n");
		exit(0);
	}
	initMesh(triangles, vertices, constraints);
	generateDelaunayNet(fileInput, triangles, vertices, numTotalTriangles, numTotalVertices);
	if (strstr(fileInputRestrict, "m2d") != NULL) {
		loadConstraintFromM2DFile(fileInputRestrict, constraints, numTotalSegments);
	}
	else {
		loadConstraintFromFile(fileInputRestrict, constraints, numTotalSegments);
	}
	restrictDelaunayNet(triangles, vertices, constraints, numTotalTriangles, numTotalVertices, numTotalSegments, clearMode);
	exportData(triangles, *numTotalTriangles, fileOutput);
	return 0;
}
