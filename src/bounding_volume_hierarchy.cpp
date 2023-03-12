#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "interpolate.h"
#include "intersect.h"
#include "scene.h"
#include "texture.h"
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <interpolate.h>
#include <iostream>
#include <limits>
#include <queue>
int BVH_MAX_LEVELS = 15;
int BVH_MIN_TRI_PER_NODE = 10;
int BVH_SPLIT_AMOUNT = 4;
const float MIN_RAY_LENGTH = 1e-6;
bool SAH = false;

BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene, const Features& features)
    : m_pScene(pScene)
{
    Node root {
        .m_isLeaf = false
    };
    int maxLevels = BVH_MAX_LEVELS, levelsLeft = maxLevels;
    Axis axis = x;
    for (int i = 0; i < m_pScene->meshes.size(); i++) {
        Mesh& m = m_pScene->meshes[i];
        for (int k = 0; k < m.triangles.size(); k++) {
            root.m_children.push_back({ i, k });
        }
    }
    getBoundingBox(root);
    if (root.m_children.size() <= BVH_MIN_TRI_PER_NODE) {
        root.m_isLeaf = true;
        m_numLevels = 1;
        m_numLeaves = 1;
        m_nodes.push_back(root);
        return;
    }

    m_nodes.push_back(root);
    if (!root.m_isLeaf) {
        if (features.extra.enableBvhSahBinning) {
            SAH = true;
            levelsLeft = splitBoundingBoxSAH(0, m_pScene->meshes, levelsLeft);
        } else {
            SAH = false;
            levelsLeft = splitBoundingBox(axis, 0, m_pScene->meshes, levelsLeft);
        }
    }

    m_numLevels = maxLevels - levelsLeft;

    std::vector<int> indexes;
    collectLeafNodes(0, indexes);

    m_numLeaves = indexes.size();
}

int BoundingVolumeHierarchy::getBoundingBox(Node& node)
{
    Mesh* mesh;
    Vertex* v;
    glm::vec3 max = glm::vec3(-std::numeric_limits<float>::max());
    glm::vec3 min = glm::vec3(std::numeric_limits<float>::max());
    std::vector<glm::ivec2> indexes;
    // collectLeafNodesVertices(node, indexes);

    for (const glm::ivec2 i : node.m_children) {
        mesh = &m_pScene->meshes[i.x];
        v = &mesh->vertices[mesh->triangles[i.y].x];
        if (max.x < v->position.x) {
            max.x = v->position.x;
        }
        if (min.x >= v->position.x) {
            min.x = v->position.x;
        }
        if (max.y < v->position.y) {
            max.y = v->position.y;
        }
        if (min.y >= v->position.y) {
            min.y = v->position.y;
        }
        if (max.z < v->position.z) {
            max.z = v->position.z;
        }
        if (min.z >= v->position.z) {
            min.z = v->position.z;
        }

        v = &mesh->vertices[mesh->triangles[i.y].y];
        if (max.x < v->position.x) {
            max.x = v->position.x;
        }
        if (min.x >= v->position.x) {
            min.x = v->position.x;
        }
        if (max.y < v->position.y) {
            max.y = v->position.y;
        }
        if (min.y >= v->position.y) {
            min.y = v->position.y;
        }
        if (max.z < v->position.z) {
            max.z = v->position.z;
        }
        if (min.z >= v->position.z) {
            min.z = v->position.z;
        }

        v = &mesh->vertices[mesh->triangles[i.y].z];
        if (max.x < v->position.x) {
            max.x = v->position.x;
        }
        if (min.x >= v->position.x) {
            min.x = v->position.x;
        }
        if (max.y < v->position.y) {
            max.y = v->position.y;
        }
        if (min.y >= v->position.y) {
            min.y = v->position.y;
        }
        if (max.z < v->position.z) {
            max.z = v->position.z;
        }
        if (min.z >= v->position.z) {
            min.z = v->position.z;
        }
    }
    node.m_lower = min;
    node.m_upper = max;
    return 0;
}

int BoundingVolumeHierarchy::splitBoundingBox(Axis axis, int nodeIndex, std::vector<Mesh>& meshes, int level)
{
    level--;

    float median;
    std::vector<float> points;
    std::vector<glm::ivec2> l1, l2;

    // only one triangle inside
    if (m_nodes[nodeIndex].m_children.size() <= BVH_MIN_TRI_PER_NODE || level == 0) {
        m_nodes[nodeIndex].m_isLeaf = true;
        return level;
    }

    // cycle through all of the AABB triangles (x = mesh, y = tri)
    for (const glm::ivec2& indexes : m_nodes[nodeIndex].m_children) {
        const glm::vec3 t = meshes[indexes.x].triangles[indexes.y];
        const glm::vec3 triangleCenter = (meshes[indexes.x].vertices[t.x].position + meshes[indexes.x].vertices[t.y].position + meshes[indexes.x].vertices[t.z].position);
        if (axis == 0) {
            points.push_back(triangleCenter.x / 3.f);
        } else if (axis == 1) {
            points.push_back(triangleCenter.y / 3.f);
        } else {
            points.push_back(triangleCenter.z / 3.f);
        }
    }

    std::vector<float> pointsCopy { points };

    // sort points
    std::sort(points.begin(), points.end());

    // pick the median element of the sorted vector
    if (points.size() % 2 == 1) {
        int x = points.size() / 2;
        median = points[x];
    } else {
        int x = (points.size() + 1) / 2;
        median = points[x];
    }

    // use the original positions of the points to link them to their triangles
    // and put them into the two arrays
    for (int i = 0; i < pointsCopy.size(); i++) {
        if (pointsCopy[i] < median) {
            l1.push_back(m_nodes[nodeIndex].m_children[i]);
        } else {
            l2.push_back(m_nodes[nodeIndex].m_children[i]);
        }
    }

    // in case the split is really weird, put at least one node in each
    //  goes one way because of 'pointsCopy[i] < median'
    if (l1.size() == 0) {
        l1.push_back(l2[l2.size() - 1]);
        l2.pop_back();
    }

    Node newNode1 {
        .m_isLeaf = false,
        .m_children = l1
    };
    Node newNode2 {
        .m_isLeaf = false,
        .m_children = l2
    };

    getBoundingBox(newNode1);
    getBoundingBox(newNode2);

    m_nodes[nodeIndex].m_children.clear();

    // get the next axis
    Axis newAxis = static_cast<Axis>((axis + 1) % 3);

    const int node1index = m_nodes.size();
    // the y coordinate is not used for nodes, only for triangles
    m_nodes[nodeIndex].m_children.push_back({ node1index, 0 });
    m_nodes.push_back(newNode1);
    const int node1ReachedLevel = splitBoundingBox(newAxis, node1index, meshes, level);

    const int node2index = m_nodes.size();
    // the y coordinate is not used for nodes, only for triangles
    m_nodes[nodeIndex].m_children.push_back({ node2index, 0 });
    m_nodes.push_back(newNode2);
    const int node2ReachedLevel = splitBoundingBox(newAxis, node2index, meshes, level);

    return node1ReachedLevel < node2ReachedLevel ? node1ReachedLevel : node2ReachedLevel;
}

int BoundingVolumeHierarchy::splitBoundingBoxSAH(int nodeIndex, std::vector<Mesh>& meshes, int level)
{
    level--;

    float median;
    std::vector<float> points;
    std::vector<glm::ivec2> l1, l2;

    // only one triangle inside
    if (m_nodes[nodeIndex].m_children.size() <= BVH_MIN_TRI_PER_NODE || level == 0) {
        m_nodes[nodeIndex].m_isLeaf = true;
        return level;
    }

    getBestSplit(nodeIndex, l1, l2);

    // in case the split is really weird, put at least one node in each
    //  goes one way because of 'pointsCopy[i] < median'
    if (l1.size() == 0 || l2.size() == 0) {
        m_nodes[nodeIndex].m_isLeaf = true;
        return level;
    }

    Node newNode1 {
        .m_isLeaf = false,
        .m_children = l1
    };
    Node newNode2 {
        .m_isLeaf = false,
        .m_children = l2
    };
    getBoundingBox(newNode1);
    getBoundingBox(newNode2);

    m_nodes[nodeIndex].m_children.clear();

    const int node1index = m_nodes.size();
    // the y coordinate is not used for nodes, only for triangles
    m_nodes[nodeIndex].m_children.push_back({ node1index, 0 });
    m_nodes.push_back(newNode1);
    const int node1ReachedLevel = splitBoundingBoxSAH(node1index, meshes, level);

    const int node2index = m_nodes.size();
    // the y coordinate is not used for nodes, only for triangles
    m_nodes[nodeIndex].m_children.push_back({ node2index, 0 });
    m_nodes.push_back(newNode2);
    const int node2ReachedLevel = splitBoundingBoxSAH(node2index, meshes, level);

    return node1ReachedLevel < node2ReachedLevel ? node1ReachedLevel : node2ReachedLevel;
}

void BoundingVolumeHierarchy::getBestSplit(int nodeIndex, std::vector<glm::ivec2>& l1, std::vector<glm::ivec2>& l2)
{
    Axis axis;
    Node& n = m_nodes[nodeIndex];
    std::vector<Mesh>& meshes = m_pScene->meshes;
    float bins = BVH_SPLIT_AMOUNT;
    float xL = n.m_upper.x - n.m_lower.x;
    float yL = n.m_upper.y - n.m_lower.y;
    float zL = n.m_upper.z - n.m_lower.z;
    float max = glm::max(glm::max(xL, yL), zL);
    if (xL == max) {
        axis = x;
    } else if (yL == max) {
        axis = y;
    } else {
        axis = z;
    }

    float best = std::numeric_limits<float>::max(), pA, pB, pC = 2.f * (xL * yL + xL * zL + yL * zL);
    float countAll = n.m_children.size();
    std::vector<glm::ivec2> ll1, ll2;
    for (int i = 1; i < bins; i++) {
        float bin;
        if (axis == x) {
            bin = ((1.f / bins) * (i)*xL);
            pA = (2.f * (bin * yL + bin * zL + yL * zL)) / pC;
            pB = (2.f * ((xL - bin) * yL + (xL - bin) * zL + yL * zL)) / pC;
        } else if (axis == y) {
            bin = ((1.f / bins) * (i)*yL);
            pA = (2.f * (xL * bin + xL * zL + bin * zL)) / pC;
            pB = (2.f * (xL * (yL - bin) + xL * zL + (yL - bin) * zL)) / pC;
        } else {
            bin = ((1.f / bins) * (i)*zL);
            pA = (2.f * (xL * yL + xL * bin + yL * bin)) / pC;
            pB = (2.f * (xL * yL + xL * (zL - bin) + yL * (zL - bin))) / pC;
        }
        ll1.clear();
        ll2.clear();
        float count = 0;
        // cycle through all of the AABB triangles (x = mesh, y = tri)
        for (const glm::ivec2& indexes : n.m_children) {
            const glm::vec3 t = meshes[indexes.x].triangles[indexes.y];
            const glm::vec3 triangleCenter = (meshes[indexes.x].vertices[t.x].position + meshes[indexes.x].vertices[t.y].position + meshes[indexes.x].vertices[t.z].position);
            if (axis == 0) {
                if ((triangleCenter.x / 3.f) <= (n.m_lower.x + bin)) {
                    ll1.push_back(indexes);
                    count++;
                } else {
                    ll2.push_back(indexes);
                }
            } else if (axis == 1) {
                if ((triangleCenter.y / 3.f) <= (n.m_lower.y + bin)) {
                    ll1.push_back(indexes);
                    count++;
                } else {
                    ll2.push_back(indexes);
                }
            } else {
                if ((triangleCenter.z / 3.f) <= (n.m_lower.z + bin)) {
                    ll1.push_back(indexes);
                    count++;
                } else {
                    ll2.push_back(indexes);
                }
            }
        }
        const float newBest = 1.f + (pA * count) + (pB * (countAll - count));
        if (newBest < best) {
            best = newBest;
            l1 = ll1;
            l2 = ll2;
        }
    }
}

// Return the depth of the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 1.
int BoundingVolumeHierarchy::numLevels() const
{
    return m_numLevels;
}

// Return the number of leaf nodes in the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 2.
int BoundingVolumeHierarchy::numLeaves() const
{
    return m_numLeaves;
}

// Use this function to visualize your BVH. This is useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDrawLevel(int level)
{
    drawNodesAtLevel(m_nodes[0], 0, level);
}

// Use this function to visualize your leaf nodes. This is useful for debugging. The function
// receives the leaf node to be draw (think of the ith leaf node). Draw the AABB of the leaf node and all contained triangles.
// You can draw the triangles with different colors. NoteL leafIdx is not the index in the node vector, it is the
// i-th leaf node in the vector.
void BoundingVolumeHierarchy::debugDrawLeaf(int leafIdx)
{
    // Draw the AABB as a transparent green box.
    // AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    // drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);
    leafIdx--;
    std::vector<int> indexes;
    std::vector<glm::vec3> triangles;
    collectLeafNodes(0, indexes);
    if (leafIdx < indexes.size()) {
        AxisAlignedBox aabb { m_nodes[indexes[leafIdx]].m_upper, m_nodes[indexes[leafIdx]].m_lower };
        drawAABB(aabb, DrawMode::Wireframe, glm::vec3(0.05f, 1.0f, 0.05f), 1);
        for (const auto& inde : m_nodes[indexes[leafIdx]].m_children) {
            const auto& tri = m_pScene->meshes[inde.x].triangles[inde.y];
            drawTriangle(m_pScene->meshes[inde.x].vertices[tri.x], m_pScene->meshes[inde.x].vertices[tri.y], m_pScene->meshes[inde.x].vertices[tri.z]);
        }
    }
    // Draw the AABB as a (white) wireframe box.
    // drawAABB(aabb, DrawMode::Wireframe);

    // once you find the leaf node, you can use the function drawTriangle (from draw.h) to draw the contained primitives
}

void BoundingVolumeHierarchy::drawNodesAtLevel(struct Node& root, int currentLevel, int level)
{
    if (currentLevel == level) {
        AxisAlignedBox aabb { root.m_lower, root.m_upper };
        drawAABB(aabb, DrawMode::Wireframe);
        if (SAH) {
            Axis axis;
            std::vector<Mesh>& meshes = m_pScene->meshes;
            float bins = BVH_SPLIT_AMOUNT;
            float xL = root.m_upper.x - root.m_lower.x;
            float yL = root.m_upper.y - root.m_lower.y;
            float zL = root.m_upper.z - root.m_lower.z;
            float max = glm::max(glm::max(xL, yL), zL);
            if (xL == max) {
                axis = x;
            } else if (yL == max) {
                axis = y;
            } else {
                axis = z;
            }
            for (int i = 1; i < bins; i++) {
                float bin;
                if (axis == x) {
                    bin = root.m_lower.x + ((1.f / bins) * (i)*xL);
                    drawAABB({ { bin, root.m_lower.y, root.m_lower.z }, { bin, root.m_upper.y, root.m_upper.z } }, DrawMode::Wireframe, glm::vec3(1.0f, 0.0f, 0.0f), 1);
                } else if (axis == y) {
                    bin = root.m_lower.y + ((1.f / bins) * (i)*yL);
                    drawAABB({ { root.m_lower.x, bin, root.m_lower.z }, { root.m_upper.x, bin, root.m_upper.z } }, DrawMode::Wireframe, glm::vec3(1.0f, 0.0f, 0.0f), 1);
                } else {
                    bin = root.m_lower.z + ((1.f / bins) * (i)*zL);
                    drawAABB({ { root.m_lower.x, root.m_lower.y, bin }, { root.m_upper.x, root.m_upper.y, bin } }, DrawMode::Wireframe, glm::vec3(1.0f, 0.0f, 0.0f), 1);
                }
            }
        }
        return;
    }

    if (!root.m_isLeaf) {
        drawNodesAtLevel(m_nodes[root.m_children[0].x], currentLevel + 1, level);
        drawNodesAtLevel(m_nodes[root.m_children[1].x], currentLevel + 1, level);
    }
}

void BoundingVolumeHierarchy::collectLeafNodesVertices(Node& root, std::vector<glm::ivec2>& indexes)
{
    // if node is leaf node, print its data
    if (root.m_isLeaf) {
        for (const glm::ivec2& v : root.m_children) {
            indexes.push_back(v);
        }
    } else {
        collectLeafNodesVertices(m_nodes[root.m_children[0].x], indexes);
        collectLeafNodesVertices(m_nodes[root.m_children[1].x], indexes);
    }
}

void BoundingVolumeHierarchy::collectLeafNodes(int rootIndex, std::vector<int>& indexes)
{
    // if node is leaf node, print its data
    if (m_nodes[rootIndex].m_isLeaf) {
        indexes.push_back(rootIndex);
    } else {
        collectLeafNodes(m_nodes[rootIndex].m_children[0].x, indexes);
        collectLeafNodes(m_nodes[rootIndex].m_children[1].x, indexes);
    }
}

bool checkTriangle(Vertex& v0, Vertex& v1, Vertex& v2, Ray& ray, HitInfo& hitInfo, std::vector<Vertex>& hitTriangle, Mesh& mesh, Features features)
{
    // Holds the old t in case the new one is larger.
    float oldT = ray.t;

    if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {

        // Checks that the ray has at least a certain length and that it makes the t smaller.
        if (!(ray.t > MIN_RAY_LENGTH && ray.t < oldT)) {
            ray.t = oldT;
            return false;
        }

        // Sets the material.
        hitInfo.material = mesh.material;

        // Calculates barycentric coordinates.
        hitInfo.barycentricCoord = computeBarycentricCoord(v0.position, v1.position, v2.position, ray.origin + ray.direction * ray.t);

        // If we have interpolated normals enabled
        if (features.enableNormalInterp) {
            hitInfo.normal = interpolateNormal(v0.normal, v1.normal, v2.normal, hitInfo.barycentricCoord);
        } else {
            hitInfo.normal = glm::normalize(glm::cross(v1.position - v0.position, v2.position - v0.position));
        }

        // Reverses the normal if the point is on the other side of the normal.
        if (glm::dot(ray.direction, hitInfo.normal) > 0) {
            hitInfo.normal = -hitInfo.normal;
        }

        // Sets the hit triangle so the calling function can get this triangle back only if it hit.
        hitTriangle = { v0, v1, v2 };

        // Sets the texture coordinates to hitInfo
        hitInfo.texCoord = interpolateTexCoord(v0.texCoord, v1.texCoord, v2.texCoord, hitInfo.barycentricCoord);

        return true;
    }

    return false;
}

bool inBox(glm::vec3 lower, glm::vec3 upper, glm::vec3 point)
{
    return point.x >= lower.x && point.x <= upper.x && point.y >= lower.y && point.y <= upper.y && point.z >= lower.z && point.z <= upper.z;
}
// Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h.
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const
{

    std::vector<Vertex> hitTriangle;
    bool hit = false;

    // If BVH is not enabled, use the naive implementation.
    if (!features.enableAccelStructure) {
        // Intersect with all triangles of all meshes.
        for (auto& mesh : m_pScene->meshes) {
            for (const auto& tri : mesh.triangles) {
                auto& v0 = mesh.vertices[tri[0]];
                auto& v1 = mesh.vertices[tri[1]];
                auto& v2 = mesh.vertices[tri[2]];
                hit |= checkTriangle(v0, v1, v2, ray, hitInfo, hitTriangle, mesh, features);
            }
        }

    } else {

        // Comparator for the pairs of nodes and distances.
        struct comparator {
            auto operator()(std::pair<int, float>& const a, std::pair<int, float>& const b) const
                -> bool
            {
                return a.second > b.second;
            }
        };

        std::priority_queue<std::pair<int, float>, std::vector<std::pair<int, float>>, comparator> intersections;
        auto& start = this->m_nodes.front();
        auto copyRay = Ray(ray.origin, ray.direction);

        if (inBox(start.m_lower, start.m_upper, copyRay.origin)) {
            intersections.push({ 0, 0 });
        } else if (intersectRayWithShape({ start.m_lower, start.m_upper }, copyRay)) {
            intersections.push({ 0, copyRay.t });
        }

        while (!intersections.empty()) {
            auto top = intersections.top();
            auto& node = m_nodes[top.first];
            intersections.pop();

            // Checks if we already met a node and that this one can't be closer.
            if (hit && top.second > ray.t) {

                if (DrawIntersectedBVHs) {
                    drawDebugAABB({ node.m_lower, node.m_upper }, DrawMode::Wireframe,
                        DrawUnvisitedInDifferentColor ? glm::vec3(1, 0, 0) : glm::vec3(0, 0, 1), 1.0f);
                }

                if (DrawUnvisitedInDifferentColor && enableDebugDraw)
                    continue;
                break;
            }

            if (DrawIntersectedBVHs) {
                drawDebugAABB({ node.m_lower, node.m_upper }, DrawMode::Wireframe, glm::vec3(0, 0, 1), 1.0f);
            }

            if (node.m_isLeaf) {
                for (auto& i : node.m_children) {
                    auto& mesh = m_pScene->meshes[i.x];
                    auto& tri = mesh.triangles[i.y];
                    auto& v0 = mesh.vertices[tri[0]];
                    auto& v1 = mesh.vertices[tri[1]];
                    auto& v2 = mesh.vertices[tri[2]];
                    float oldt = ray.t;
                    hit |= checkTriangle(v0, v1, v2, ray, hitInfo, hitTriangle, mesh, features);
                }
            } else {
                auto& left = node.m_children[0].x;
                auto& right = node.m_children[1].x;
                auto& leftNode = m_nodes[left];
                auto& rightNode = m_nodes[right];

                // Adds left node if it gets hit by ray
                copyRay = Ray(ray.origin, ray.direction, std::numeric_limits<float>::max());
                if (inBox(leftNode.m_lower, leftNode.m_upper, copyRay.origin)) {
                    intersections.push({ left, 0 });
                } else if (intersectRayWithShape({ leftNode.m_lower, leftNode.m_upper }, copyRay)) {
                    intersections.push({ left, copyRay.t });
                }

                // Adds right node if it gets hit by ray
                copyRay = Ray(ray.origin, ray.direction, std::numeric_limits<float>::max());
                if (inBox(rightNode.m_lower, rightNode.m_upper, copyRay.origin)) {
                    intersections.push({ right, 0 });
                } else if (intersectRayWithShape({ rightNode.m_lower, rightNode.m_upper }, copyRay)) {
                    intersections.push({ right, copyRay.t });
                }
            }
        }
        if (hit && DrawHitTriangle) {
            drawDebugTriangle(hitTriangle[0], hitTriangle[1], hitTriangle[2], glm::vec3(0, 1, 0));
        }
    }

    // To indicate if the new intersection is with a sphere.
    bool hitSphere = false;

    // Intersect with spheres.
    for (const auto& sphere : m_pScene->spheres) {
        auto oldT = ray.t;
        if (intersectRayWithShape(sphere, ray, hitInfo)) {

            // Checks that the ray has at least a certain length
            if (!(ray.t > MIN_RAY_LENGTH && ray.t < oldT)) {
                ray.t = oldT;
                continue;
            }

            hitInfo.normal = glm::normalize(ray.origin + ray.direction * ray.t - sphere.center);
            hitInfo.material = sphere.material;
            // Reverses the normal if the point is inside the sphere.
            if (glm::distance(ray.origin, sphere.center) < sphere.radius) {
                hitInfo.normal = -hitInfo.normal;
            }
            hit = true;
            hitSphere = true;
        }
    }

    // Debug Drawing of the interpolated normal.
    if (features.enableNormalInterp && hit && DrawInterpolatedNormal) {
        auto normalRay = Ray(ray.origin + ray.direction * ray.t, hitInfo.normal, .5);
        drawRay(normalRay, DrawInterpolatedNormalColor ? hitInfo.barycentricCoord : glm::vec3(1));

        // If hit sphere don't draw vertex normals.
        for (auto& v : hitTriangle) {
            if (DrawVertexNormals)
                drawRay(Ray(v.position, v.normal, .5), glm::vec3(1, 0, 0));
        }
    }
    return hit;
}
