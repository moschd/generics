#ifndef DTM_GENERICS_POINT_H
#define DTM_GENERICS_POINT_H

struct Point3D
{
    int x_ = 0;
    int y_ = 0;
    int z_ = 0;
    Point3D() = default;
    explicit Point3D(int x, int y, int z) : x_(x), y_(y), z_(z) {};
    inline bool operator==(const Point3D &other) const { return x_ == other.x_ && y_ == other.y_ && z_ == other.z_; }
    inline bool operator!=(const Point3D &other) const { return !(*this == other); }
};

/**
 * @brief High-performance 3D spatial hasher for Point3D objects.
 * * This hasher is optimized for 3D bin packing scenarios where coordinates
 * are often clustered together. It uses a linear combination of coordinates
 * followed by a SplitMix64 finalizer to ensure a uniform distribution (avalanche effect).
 */
struct Point3DHasher
{
    /**
     * @brief Computes a 64-bit hash for a 3D point.
     * * @param p The Point3D object to hash (expects x_, y_, z_ members).
     * @return std::size_t A high-entropy hash value.
     */
    inline std::size_t operator()(const Point3D &p) const noexcept
    {
        // Convert coordinates to 64-bit and multiply by
        // large primes to spread the input bits across the 64-bit space.
        uint64_t h = (static_cast<uint64_t>(p.x_) * 0x1f23b) ^
                     (static_cast<uint64_t>(p.y_) * 0x2e35d) ^
                     (static_cast<uint64_t>(p.z_) * 0x3d47f);

        // This sequence (XOR-Shift followed by multiplication) ensures that
        // even a 1-bit change in the input coordinates results in a
        // big change in the output hash bits.
        h ^= h >> 33;               // Mix the high bits into the low bits
        h *= 0xff51afd7ed558ccdULL; // MurmurHash3/SplitMix64 constant
        h ^= h >> 33;               // Further bit diffusion
        h *= 0xc4ceb9fe1a85ec53ULL; // Secondary mixing constant
        h ^= h >> 33;               // Final bit smear

        return static_cast<std::size_t>(h);
    }
};

#endif