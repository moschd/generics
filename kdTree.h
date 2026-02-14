#ifndef DTM_GENERICS_KD_TREE_H
#define DTM_GENERICS_KD_TREE_H

#include "point.h"
#include <memory>
#include <vector>
#include <math.h>
#include <iostream>

namespace G_KDTree
{

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
        bool isLeaf_ = false;
        std::shared_ptr<Node> left_ = std::shared_ptr<Node>(nullptr);
        std::shared_ptr<Node> right_ = std::shared_ptr<Node>(nullptr);
        int myDepth_;
        std::vector<KD_ItemMeta> myChildren_;
        Point3D minSearchDimensions_;
        Point3D maxSearchDimensions_;
        Point3D partitionPoint_;

        Node(const Point3D aPartitionPoint, const Point3D aMins, const Point3D aMaxs, const int aCurrentDepth)
            : myDepth_(aCurrentDepth),
              minSearchDimensions_(aMins),
              maxSearchDimensions_(aMaxs),
              partitionPoint_(aPartitionPoint)
        {
            Node::myChildren_.reserve(125);
        }
    };

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
    inline std::shared_ptr<Node> generateNode(const Point3D &aPartitionPoint,
                                              const Point3D &aMins,
                                              const Point3D &aMaxs,
                                              const int aCurrentDepth)
    {
        return std::make_shared<Node>(aPartitionPoint, aMins, aMaxs, aCurrentDepth);
    }

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
        int maxDepth_;
        std::shared_ptr<Node> treeRoot_;
        Point3D minDimensions_ = Point3D(0, 0, 0);
        Point3D maxDimensions_;

        /**
         * @brief Calculates the depth of the tree that will be generated.
         *
         * This method calculates what would be the desired depth of the tree.
         * This attempts to make the algorithm more efficient by scaling the tree depending of the estimated number of items that will fit in the bin.
         *
         * @param aEstimatedNumberOfItemFits - integer indicating the estimated number of items that will be in the bin when packing has finished.
         */
        void calculateMaxDepth(int aEstimatedNumberOfItemFits) { maxDepth_ = ceil(sqrt(aEstimatedNumberOfItemFits / 125) + 1); };

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
        void generateTree(std::shared_ptr<Node> &aRoot,
                          const int aDepth,
                          const Point3D &aPartitionPoint,
                          const Point3D &aMins,
                          const Point3D &aMaxs,
                          const int aRequestedDepth)
        {

            const int axis = aDepth % 3; // R3 is 3

            // 1. Calculate the new partition point for the PREVIOUS axis
            // This centers the split point based on the bounds (Mins/Maxs)
            Point3D newPartitionPoint = aPartitionPoint;
            if (aDepth > 0)
            {
                const int prevAxis = (aDepth - 1) % 3;
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
            if (aDepth > aRequestedDepth)
            {
                aRoot->isLeaf_ = true;
                aRoot->partitionPoint_ = newPartitionPoint;
                return;
            }

            // 3. Generate child nodes
            aRoot->left_ = generateNode(newPartitionPoint, aMins, aMaxs, aDepth);
            aRoot->right_ = generateNode(newPartitionPoint, aMins, aMaxs, aDepth);

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
            generateTree(aRoot->left_, aDepth + 1, newPartitionPoint, aMins, leftMaxs, aRequestedDepth);

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
            generateTree(aRoot->right_, aDepth + 1, newPartitionPoint, rightMins, aMaxs, aRequestedDepth);
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
        void addItemKeyToLeaf(std::shared_ptr<Node> &aRoot,
                              const int aItemKey,
                              const int aDepth,
                              const Point3D &aItemMaxPosition)
        {
            // 1. Get raw pointer to avoid ref-count overhead in recursion
            Node *node = aRoot.get();

            if (node->isLeaf_)
            {
                // 2. Use the new Point3D based KD_ItemMeta constructor
                node->myChildren_.emplace_back(aItemKey, aItemMaxPosition);
                return;
            }

            const int axis = aDepth % 3;

            // 3. Coordinate branch logic for Point3D members
            bool goLeft = false;
            if (axis == G_KDTree::axis::WIDTH)
            {
                goLeft = aItemMaxPosition.x_ < node->partitionPoint_.x_;
            }
            else if (axis == G_KDTree::axis::DEPTH)
            {
                goLeft = aItemMaxPosition.y_ < node->partitionPoint_.y_;
            }
            else
            {
                goLeft = aItemMaxPosition.z_ < node->partitionPoint_.z_;
            }

            if (goLeft)
            {
                addItemKeyToLeaf(node->left_, aItemKey, aDepth + 1, aItemMaxPosition);
            }
            else
            {
                addItemKeyToLeaf(node->right_, aItemKey, aDepth + 1, aItemMaxPosition);
            }
        }

        /**
         * @brief Print tree to console..
         *
         * @param aRoot         - Node whose children will be printed.
         */
        void printTreeImp(const std::shared_ptr<Node> aRoot) const
        {
            if (aRoot == nullptr)
            {
                return;
            };

            if (!aRoot->Node::myChildren_.empty())
            {
                std::cout << "LEAF "
                          << aRoot->Node::partitionPoint_.x_ << " "
                          << aRoot->Node::partitionPoint_.y_ << " "
                          << aRoot->Node::partitionPoint_.z_ << "\n";

                std::cout << "  CHILDREN: " << aRoot->myChildren_.size() << ".\n\n";

                // for (const int child : aRoot->myChildren_)
                // {
                //     std::cout << " " << child << " ";
                // }
            }

            printTreeImp(aRoot->Node::left_);
            printTreeImp(aRoot->Node::right_);
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
        void getIntersectCandidates(const std::shared_ptr<Node> &aRoot,
                                    const int aDepth,
                                    const Point3D &aStartPoint,
                                    const Point3D &aMaxSearchPoint,
                                    std::vector<int> &aPassedNodes) const
        {
            Node *current = aRoot.get();
            if (!current)
                return;

            // 1. Leaf Logic: Standard loop is fastest here
            if (current->isLeaf_)
            {
                for (const auto &item : current->myChildren_)
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
                partitionCoord = current->partitionPoint_.x_;
                maxSearchOffset = aMaxSearchPoint.x_;
                break;
            case G_KDTree::axis::DEPTH:
                startCoord = aStartPoint.y_;
                partitionCoord = current->partitionPoint_.y_;
                maxSearchOffset = aMaxSearchPoint.y_;
                break;
            default: // axis 2
                startCoord = aStartPoint.z_;
                partitionCoord = current->partitionPoint_.z_;
                maxSearchOffset = aMaxSearchPoint.z_;
                break;
            }

            // 3. Branching Logic: Use the local ints we just grabbed!
            // This avoids reaching into the Point3D structs again.

            // Check if we need to search the LEFT branch
            if (startCoord < partitionCoord)
            {
                getIntersectCandidates(current->left_, aDepth + 1, aStartPoint, aMaxSearchPoint, aPassedNodes);
            }

            // Check if we need to search the RIGHT branch
            // If the search area overlaps or is entirely to the right of the partition
            if (partitionCoord < (startCoord + maxSearchOffset))
            {
                getIntersectCandidates(current->right_, aDepth + 1, aStartPoint, aMaxSearchPoint, aPassedNodes);
            }
        }

        void setNewRoot()
        {
            treeRoot_ = generateNode(
                Point3D(maxDimensions_.x_ / 2, maxDimensions_.y_ / 2, maxDimensions_.z_ / 2),
                minDimensions_,
                maxDimensions_,
                0);
        }

    public:
        KdTree(int aEstimatedNumberOfItemFits, Point3D aMaxDimensions) : maxDepth_(8),
                                                                         maxDimensions_(aMaxDimensions)

        {
            calculateMaxDepth(aEstimatedNumberOfItemFits);
            generateTreeHelper();
        };
        void printTreeImpHelper() const { printTreeImp(treeRoot_); };
        const std::shared_ptr<Node> &getRoot() const { return treeRoot_; };
        void generateTreeHelper()
        {
            setNewRoot();
            generateTree(treeRoot_,
                         treeRoot_->Node::myDepth_,
                         treeRoot_->Node::partitionPoint_,
                         treeRoot_->Node::minSearchDimensions_,
                         treeRoot_->Node::maxSearchDimensions_,
                         maxDepth_);
        };

        /**
         * @brief Helper function to add a new itemKey to the tree.
         *
         * This method provides a simple interface to the method that adds a item to the tree.
         *
         * @param aItemKey          - The itemKey which will be added to the tree.
         * @param aItemMaxPosition  - The top right corner of the item, this is the position which will be used to search for the correct place to add the item.
         */
        void addItemKeyToLeafHelper(const int aItemKey, const Point3D &aItemMaxPosition)
        {
            addItemKeyToLeaf(treeRoot_, aItemKey, 0, aItemMaxPosition);
        }

        void getIntersectCandidatesHelper(const Point3D &aStartPoint, const Point3D &aMaxSearchPoint, std::vector<int> &aPassedNodes) const
        {
            getIntersectCandidates(getRoot(), getRoot()->myDepth_, aStartPoint, aMaxSearchPoint, aPassedNodes);
        };
    };
}

#endif