#pragma once
#include "common.h"
#include <array>
#include <framework/ray.h>
#include <vector>

// Forward declaration.
struct Scene;

struct Node {
    bool m_isLeaf;
    glm::vec3 m_lower;
    glm::vec3 m_upper;
    std::vector<glm::ivec2> m_children;
};

enum Axis {
    x = 0,
    y = 1,
    z = 2
};

extern int BVH_MAX_LEVELS;
extern int BVH_MIN_TRI_PER_NODE;
extern int BVH_SPLIT_AMOUNT;

class BoundingVolumeHierarchy {
public:
    // Constructor. Receives the scene and builds the bounding volume hierarchy.
    BoundingVolumeHierarchy(Scene* pScene, const Features& features);

    // Return how many levels there are in the tree that you have constructed.
    [[nodiscard]] int numLevels() const;

    // Return how many leaf nodes there are in the tree that you have constructed.
    [[nodiscard]] int numLeaves() const;

    // Visual Debug 1: Draw the bounding boxes of the nodes at the selected level.
    void debugDrawLevel(int level);

    // Visual Debug 2: Draw the triangles of the i-th leaf
    void debugDrawLeaf(int leafIdx);

    // Return true if something is hit, returns false otherwise.
    // Only find hits if they are closer than t stored in the ray and the intersection
    // is on the correct side of the origin (the new t >= 0).
    bool intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const;

private:
    std::vector<Node> m_nodes;
    int m_numLevels;
    int m_numLeaves;
    Scene* m_pScene;

    // recursively splits the node at nodeIndex into AABBs
    int splitBoundingBox(Axis axis, int nodeIndex, std::vector<Mesh>& meshes, int level);

    // recursively splits the node at nodeIndex into AABBs using SAH + binning
    int splitBoundingBoxSAH(int nodeIndex, std::vector<Mesh>& meshes, int level);

    // creates a AABB from all the triangles contained within the node
    int getBoundingBox(Node& node);

    // DEBUG ONLY, draws all the triangles at a certain level
    void drawNodesAtLevel(Node& root, int currentLevel, int level);
    // gets all the leaf nodes' indeces
    void collectLeafNodes(int nodeIndex, std::vector<int>& indexes);
    // gets all the triangles within all the leaf nodes (used for bounding boxes)
    void collectLeafNodesVertices(Node& root, std::vector<glm::ivec2>& indexes);

    void getBestSplit(int nodeIndex, std::vector<glm::ivec2>& l1, std::vector<glm::ivec2>& l2);
};