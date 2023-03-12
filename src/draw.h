#pragma once
#include "scene.h"
#include <framework/mesh.h>
#include <framework/ray.h>
#include <utility> // std::forward

// Debug flags
extern bool enableDebugDraw;
extern bool HideDrawings;
extern bool DrawUnvisitedInDifferentColor;
extern bool DrawIntersectedBVHs;
extern bool DrawHitTriangle;
extern bool DrawInterpolatedNormal;
extern bool DrawVertexNormals;
extern bool DrawInterpolatedNormalColor;
extern bool DrawShadingRay;
extern bool DrawShadowRay;
extern bool DrawDepthOfField;

// Add your own custom visual debug draw functions here then implement it in draw.cpp.
// You are free to modify the example one however you like.
//
// Please add following lines on top inside your custom draw function:
//
// if (!enableDebugDraw)
//     return;
//
void drawExampleOfCustomVisualDebug();
void drawDebugSphere(const glm::vec3& center, float radius, const glm::vec3& color = glm::vec3(1.0f));
void drawDebugAABB(const AxisAlignedBox& box, DrawMode drawMode = DrawMode::Wireframe, const glm::vec3& color = glm::vec3(1.0f), float transparency = 1.0f);
void drawDebugTriangle(const Vertex& v0, const Vertex& v1, const Vertex& v2, glm::vec3 color);

void drawRay(const Ray& ray, const glm::vec3& color = glm::vec3(1.0f));

void drawAABB(const AxisAlignedBox& box, DrawMode drawMode = DrawMode::Filled, const glm::vec3& color = glm::vec3(1.0f), float transparency = 1.0f);

void drawTriangle (const Vertex& v0, const Vertex& v1, const Vertex& v2 );
void drawMesh(const Mesh& mesh);
void drawSphere(const Sphere& sphere);
void drawSphere(const glm::vec3& center, float radius, const glm::vec3& color = glm::vec3(1.0f));
void drawScene(const Scene& scene);

