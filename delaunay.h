#ifndef DELAUNAY_H
#define DELAUNAY_H

#define MAX_VERTICES 3000
#define MAX_TRIANGLES 6000
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

void addPointToMesh(
	Triangle* triangles,
	Vertex* vertices,
	int* numTotalTriangles,
	int* numTotalVertices,
	float x, float y);
void applyConstraint(
	Triangle* triangles,
	int* numTotalTriangles,
	Segment* constraint);
float circleTest(Triangle* triangle, float x, float y);
float circleTestByVertex(Triangle* triangle, Vertex* vertice);
void clearPolygon(
	Triangle* triangles,
	int* numTotalTriangles,
	Segment* constraints,
	int* numTotalSegments);
void clearTriangle(Triangle* triangle);
void exportData(Triangle* triangles, int numTotalTriangles, char* fileOutput);
void exportDataRestricted(Segment* constraints, int numTotalSegments, char* fileOutput);
Triangle* getIdTriangleContainsPoint(Triangle* triangle, int* puntoEnBorde, float x, float y);
void getDetsByTriangle(Triangle* triangle, float* dets, float x, float y);
Vertex* getThirdVertex(Triangle* triangle, Vertex* vertex1, Vertex* vertex2);
int getThirdVertexId(Triangle* triangle, Vertex* vertex1, Vertex* vertex2);
int getVertexIdByVertex(Triangle* triangle, Vertex* vertex);
int generateDelaunayNet(
	char* strFileInput,
	Triangle* triangles,
	Vertex* vertices,
	int* numTotalTriangles,
	int* numTotalVertices);
float getDetBySegments(float x1, float y1, float x2, float y2, float xp, float yp);
void initMesh(Triangle* triangles, Vertex* vertices, Segment* constraints);
void intercambioDiagonalRestriccion(int* idsVerticesTriangulo, Triangle* currentTriangle);
int intercambioDeDiagonal(Triangle* triangleA, Triangle* triangleB, int Avc1, int Avc2, int Aop, int Bop, int skipCircleTest, int skipPropagation);
int isPointAlreadyOnNet(Vertex* vertices, int numTotalVertices, float x, float y);
void loadConstraintFromFile(char* strFileInput, Segment* constraints, int* numTotalSegments);
void loadConstraintFromM2DFile(char* strFileInput, Segment* constraints, int* numTotalSegments);
void restrictDelaunayNet(
	Triangle* triangles,
	Vertex* vertices,
	Segment* constraints,
	int* numTotalTriangles,
	int* numTotalVertices,
	int* numTotalSegments,
	int clearMode);
void removeTriangle(Triangle* triangle, Triangle* triangles, int *numTotalTriangles);
void removeIsland(Triangle* triangle);
#endif