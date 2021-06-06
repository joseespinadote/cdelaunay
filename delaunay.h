#ifndef DELAUNAY_H
#define DELAUNAY_H

#define MAX_VERTICES 3000
#define MAX_TRIANGLES 6000
#define MAX_SEGMENTS 1000
#define TAMANO_MALLA_X 24
#define TAMANO_MALLA_Y 24
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
int generateDelaunayNet(char* strFileInput, Triangle* triangles, Vertex* vertices, int* numTotalTriangles, int* numTotalVertices);
float getDetBySegments(float x1, float y1, float x2, float y2, float xp, float yp);
void restrictDelaunayNet(Triangle* triangles, Vertex* vertices, Segment* constraints, int* numTotalTriangles, int* numTotalVertices, int* numTotalSegments);
Triangle* getIdTriangleContainsPoint(Triangle* triangle, int* puntoEnBorde, float x, float y);
int isPointAlreadyOnNet(Vertex* vertices, int numTotalVertices, float x, float y);
void removeTriangle(Triangle* triangle, Triangle* triangles, int *numTotalTriangles);
void removeIsland(Triangle* triangle);
void applyConstraint(Triangle* triangles, int* numTotalTriangles, Segment* constraint);
void loadConstraintFromFile(char* strFileInput, Segment* constraints, int* numTotalSegments);
void loadConstraintFromM2DFile(char* strFileInput, Segment* constraints, int* numTotalSegments);
void clearPolygon(Triangle* triangles, int* numTotalTriangles, Segment* constraints, int* numTotalSegments);
#endif