#ifndef DTM_GENERICS_KD_TREE_H
#define DTM_GENERICS_KD_TREE_H

#include "point.h"
#include <memory>
#include <vector>
#include <math.h>
#include <iostream>

namespace G_KDTree
{

    constexpr uint32_t INVALID_IDX = std::numeric_limits<uint32_t>::max();
    constexpr const unsigned int MIN_DEPTH{0};
    constexpr const unsigned int MAX_DEPTH{8};
    constexpr const unsigned int MAX_ITEMS_PER_NODE{100};

    namespace axis
    {
        constexpr const unsigned int WIDTH{0};
        constexpr const unsigned int DEPTH{1};
        constexpr const unsigned int HEIGHT{2};
        constexpr const std::array<int, 2> ALL_2D = {WIDTH, DEPTH};
        constexpr const std::array<int, 3> ALL_3D = {WIDTH, DEPTH, HEIGHT};
        constexpr const unsigned int R3{ALL_3D.size()};
    }

    struct KD_ItemMeta
    {
        int itemKey_;
        int xLocation_;
        int yLocation_;
        int zLocation_;
        KD_ItemMeta(const int aItemKey, const Point3D &aLocation)
            : itemKey_(aItemKey),
              xLocation_(aLocation.x_),
              yLocation_(aLocation.y_),
              zLocation_(aLocation.z_) {};
    };

    struct Node
    {
        int myDepth_;
        bool isLeaf_ = false;
        uint32_t leftIdx_ = G_KDTree::INVALID_IDX;
        uint32_t rightIdx_ = G_KDTree::INVALID_IDX;
        std::vector<KD_ItemMeta> children_;
        Point3D partitionPoint_;
        Node(const Point3D aPartitionPoint, const int aCurrentDepth) : myDepth_(aCurrentDepth), partitionPoint_(aPartitionPoint)
        {
            children_.reserve(50);
        }

        [[nodiscard]] inline bool hasLeftBranch() const {
        return leftIdx_ != G_KDTree::INVALID_IDX;
    }

    [[nodiscard]] inline bool hasRightBranch() const {
        return rightIdx_ != G_KDTree::INVALID_IDX;
    }
    };

    /**
     * @brief Represents the bin in tree form.
     *
     * Contains all items currently inside the bin.
     * The items are stored by their furthest possible point in space.
     * This way, when inputting the smallest possible point of an item, we can find all items which are relevant to this particular item.
     * We can then check for intersection, gravity, stackingStyle restrictions etc...
     *
     * This allows for efficient 3D search and improves performance of the algorithm by a huge amount.
     *
     */
    class KdTree
    {
    private:
        std::vector<Node> nodes_;
        uint32_t rootIdx_ = 0;
        Point3D maxDimensions_;
        int maxDepth_ = G_KDTree::MIN_DEPTH;

        /**
         * @brief Spawn a Node.
         *
         * This method creates a Node used as node in the kd-tree.
         * Each node represents a point in 3R.
         *
         * @param aPartitionPoint   - the 3R point that the node will be representating and splitting.
         * @param aMins             - 3R point marking the minimum border of the search area point for which this node will be responsable.
         * @param aMaxs             - 3R point marking the maximum border of the search area point for which this node will be responsable.
         * @param aCurrentDepth     - the tree depth on which this node is created.
         */
        inline uint32_t generateNode(const Point3D &aPartitionPoint, const int aCurrentDepth)
        {
            nodes_.emplace_back(aPartitionPoint, aCurrentDepth);
            return static_cast<uint32_t>(nodes_.size() - 1);
        }

        /**
         * @brief Pre generate a fixed depth balanced kd-tree.
         *
         * This method creates a kd-tree of a certain depth.
         * The nodes in this tree each contain cartesian coordinates indicating a 3R point in a space.
         * Used for 3R space partitioning of the bin to be packed.
         *
         * @param aRoot             - root of the tree
         * @param aDepth            - the current depth of the tree
         * @param aPartitionPoint   - the 3R point that the node will be representing, and splitting.
         * @param aMins             - 3R point marking the minimum border of the search area for which this node will be responsable.
         * @param aMaxs             - 3R point marking the maximum border of the search area for which this node will be responsable.
         * @param aRequestedDepth   - the requested maximum depth of the tree to be generated.
         */
        void generateTree(const uint32_t currentIdx,
                          const int aDepth,
                          const Point3D &aPartitionPoint,
                          const Point3D &aMins,
                          const Point3D &aMaxs)
        {

            // 1. Calculate the new partition point for the PREVIOUS axis
            // This centers the split point based on the bounds (Mins/Maxs)
            Point3D newPartitionPoint = aPartitionPoint;
            if (aDepth > 0)
            {
                const int prevAxis = (aDepth - 1) % G_KDTree::axis::R3;
                if (prevAxis == G_KDTree::axis::WIDTH)
                {
                    newPartitionPoint.x_ = (aMins.x_ + aMaxs.x_) / 2;
                }
                else if (prevAxis == G_KDTree::axis::DEPTH)
                {
                    newPartitionPoint.y_ = (aMins.y_ + aMaxs.y_) / 2;
                }
                else
                {
                    newPartitionPoint.z_ = (aMins.z_ + aMaxs.z_) / 2;
                }
            }

            // 2. Base Case: Leaf Node
            if (aDepth >= maxDepth_)
            {
                nodes_[currentIdx].isLeaf_ = true;
                nodes_[currentIdx].partitionPoint_ = newPartitionPoint;
                return;
            }

            const int axis = aDepth % G_KDTree::axis::R3;
            const uint32_t left = generateNode(newPartitionPoint, aDepth);
            const uint32_t right = generateNode(newPartitionPoint, aDepth);
            nodes_[currentIdx].leftIdx_ = left;
            nodes_[currentIdx].rightIdx_ = right;

            // 4. Recurse Left: Update Maxs for the current split axis
            Point3D leftMaxs = aMaxs;
            if (axis == G_KDTree::axis::WIDTH)
            {
                leftMaxs.x_ = newPartitionPoint.x_;
            }
            else if (axis == G_KDTree::axis::DEPTH)
            {
                leftMaxs.y_ = newPartitionPoint.y_;
            }
            else
            {
                leftMaxs.z_ = newPartitionPoint.z_;
            }
            generateTree(left, aDepth + 1, newPartitionPoint, aMins, leftMaxs);

            // 5. Recurse Right: Update Mins for the current split axis
            Point3D rightMins = aMins;
            if (axis == G_KDTree::axis::WIDTH)
            {
                rightMins.x_ = newPartitionPoint.x_;
            }
            else if (axis == G_KDTree::axis::DEPTH)
            {
                rightMins.y_ = newPartitionPoint.y_;
            }
            else
            {
                rightMins.z_ = newPartitionPoint.z_;
            }
            generateTree(right, aDepth + 1, newPartitionPoint, rightMins, aMaxs);
        }

        /**
         * @brief Find correct place to add itemKey in tree and add it.
         *
         * This method adds a itemKey to the corresponding tree leaf.
         * Each bin to be packed uses a kd-tree for space partitioning and organizing the multidimensional data (items) that it contains.
         *
         * @param aRoot             - Node from which to start searching, normally start at the aRoot of the tree.
         * @param aItemKey          - The key of the item that needs to be placed in the tree.
         * @param aDepth            - The current depth of the tree.
         * @param aItemMaxPosition  - The furthest point in space that the item reaches, ie top right corner.
         */
        void addItemKeyToLeaf(const uint32_t aCurrentIdx,
                              const int aItemKey,
                              const int aDepth,
                              const Point3D &aItemMaxPosition)
        {
            // 1. Get raw pointer to avoid ref-count overhead in recursion
            Node &node = nodes_[aCurrentIdx];

            if (node.isLeaf_)
            {
                // 2. Use the new Point3D based KD_ItemMeta constructor
                node.children_.emplace_back(aItemKey, aItemMaxPosition);
                return;
            }

            const int axis = aDepth % 3;

            // 3. Coordinate branch logic for Point3D members
            bool goLeft = false;
            if (axis == G_KDTree::axis::WIDTH)
            {
                goLeft = aItemMaxPosition.x_ < node.partitionPoint_.x_;
            }
            else if (axis == G_KDTree::axis::DEPTH)
            {
                goLeft = aItemMaxPosition.y_ < node.partitionPoint_.y_;
            }
            else
            {
                goLeft = aItemMaxPosition.z_ < node.partitionPoint_.z_;
            }

            if (goLeft)
            {
                addItemKeyToLeaf(node.leftIdx_, aItemKey, aDepth + 1, aItemMaxPosition);
            }
            else
            {
                addItemKeyToLeaf(node.rightIdx_, aItemKey, aDepth + 1, aItemMaxPosition);
            }
        }

        /**
         * @brief Print tree to console..
         *
         * @param aRoot - Node whose children will be printed.
         */
        void printTree(const uint32_t aCurrentIdx) const
        {
            const Node &aRoot = nodes_[aCurrentIdx];

            if (!aRoot.children_.empty())
            {
                std::cout << "LEAF "
                          << aRoot.partitionPoint_.x_ << " "
                          << aRoot.partitionPoint_.y_ << " "
                          << aRoot.partitionPoint_.z_ << "\n";

                std::cout << "  CHILDREN: " << aRoot.children_.size() << ".\n\n";
            }

            if (aRoot.hasLeftBranch())
            {
                printTree(aRoot.leftIdx_);
            }
            if (aRoot.hasRightBranch())
            {
                printTree(aRoot.rightIdx_);
            }
        };

        /**
         * @brief Adds itemKeys of items which might be intersecting with the item on the provided position.
         *
         * This method adds the itemKey of items which might be intersecting with the item on the provided position to a provided vector.
         *
         * @param aRoot             - Node from which to start searching, normally start at the aRoot of the tree.
         * @param aDepth            - The current depth of the tree
         * @param aStartPoint       - The 3R point from which to start searching.
         * @param aMaxSearchPoint   - The furthest point in space that an item can be in order to still be considered a intersection candidate.
         * @param aPassedNodes      - The vector to which itemKeys will be added.
         */
        void getIntersectCandidates(const uint32_t aIdx,
                                    const int aDepth,
                                    const Point3D &aStartPoint,
                                    const Point3D &aMaxSearchPoint,
                                    std::vector<int> &aPassedNodes) const
        {
            const Node &current = nodes_[aIdx];

            if (current.isLeaf_)
            {
                for (const auto &item : current.children_)
                {
                    // Bounding box intersection check
                    if (item.xLocation_ >= aStartPoint.x_ &&
                        item.yLocation_ >= aStartPoint.y_ &&
                        item.zLocation_ >= aStartPoint.z_)
                    {
                        aPassedNodes.push_back(item.itemKey_);
                    }
                }
                return;
            }

            // 2. Axis Selection: Use the variables we extract to avoid repeat lookups
            const int axis = aDepth % G_KDTree::axis::R3;
            int startCoord;
            int partitionCoord;
            int maxSearchOffset;

            switch (axis)
            {
            case G_KDTree::axis::WIDTH:
                startCoord = aStartPoint.x_;
                partitionCoord = current.partitionPoint_.x_;
                maxSearchOffset = aMaxSearchPoint.x_;
                break;
            case G_KDTree::axis::DEPTH:
                startCoord = aStartPoint.y_;
                partitionCoord = current.partitionPoint_.y_;
                maxSearchOffset = aMaxSearchPoint.y_;
                break;
            default: // axis 2
                startCoord = aStartPoint.z_;
                partitionCoord = current.partitionPoint_.z_;
                maxSearchOffset = aMaxSearchPoint.z_;
                break;
            }

            // Check if we need to search the LEFT branch
            if (startCoord < partitionCoord && current.hasLeftBranch())
            {
                getIntersectCandidates(current.leftIdx_, aDepth + 1, aStartPoint, aMaxSearchPoint, aPassedNodes);
            }

            // Check if we need to search the RIGHT branch
            // If the search area overlaps or is entirely to the right of the partition
            if (partitionCoord < (startCoord + maxSearchOffset) && current.hasRightBranch())
            {
                getIntersectCandidates(current.rightIdx_, aDepth + 1, aStartPoint, aMaxSearchPoint, aPassedNodes);
            }
        }

        uint32_t generateRoot()
        {
            return generateNode(
                Point3D(maxDimensions_.x_ / 2, maxDimensions_.y_ / 2, maxDimensions_.z_ / 2),
                G_KDTree::MIN_DEPTH
            );
        }

    public:
        KdTree(const int aEstimatedNrOfItemFits, const Point3D aMaxDimensions) : maxDimensions_(aMaxDimensions)

        {
            nodes_.reserve(aEstimatedNrOfItemFits);
            // Calculates the depth of the tree that will be generated.
            // This attempts to make the algorithm more efficient by scaling the tree depending of the estimated number of items that will fit in the bin.
            maxDepth_ = static_cast<int>(ceil(sqrt(aEstimatedNrOfItemFits / G_KDTree::MAX_ITEMS_PER_NODE) + 1));
            maxDepth_ = std::min(maxDepth_, static_cast<int>(G_KDTree::MAX_DEPTH));
            initialize();
        };

        void initialize()
        {
            nodes_.clear();
            rootIdx_ = generateRoot();
            generateTree(
                rootIdx_,
                G_KDTree::MIN_DEPTH,
                nodes_[rootIdx_].partitionPoint_,
                Point3D(0, 0, 0),
                maxDimensions_);
        }

        /**
         * @brief Add a new itemKey to the tree.
         *
         * This method provides a simple interface to the method that adds a item to the tree.
         *
         * @param aItemKey          - The itemKey which will be added to the tree.
         * @param aItemMaxPosition  - The top right corner of the item, this is the position which will be used to search for the correct place to add the item.
         */
        void addItemKeyToLeaf(const int aItemKey, const Point3D &aItemMaxPosition)
        {
            addItemKeyToLeaf(rootIdx_, aItemKey, G_KDTree::MIN_DEPTH, aItemMaxPosition);
        }

        void getIntersectCandidates(
            const Point3D &aStartPoint,
            const Point3D &aMaxSearchPoint,
            std::vector<int> &aPassedNodes) const
        {
            getIntersectCandidates(
                rootIdx_,
                G_KDTree::MIN_DEPTH,
                aStartPoint,
                aMaxSearchPoint,
                aPassedNodes);
        };

        void printTree() const { printTree(rootIdx_); };
    };
}

#endif