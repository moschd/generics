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

struct Point3DHasher
{
    inline std::size_t operator()(const Point3D &p) const noexcept
    {
        uint64_t h = (static_cast<uint64_t>(p.x_) * 0x1f23b) ^
                     (static_cast<uint64_t>(p.y_) * 0x2e35d) ^
                     (static_cast<uint64_t>(p.z_) * 0x3d47f);

        // SplitMix64 finalize step
        h ^= h >> 33;
        h *= 0xff51afd7ed558ccdULL;
        h ^= h >> 33;
        h *= 0xc4ceb9fe1a85ec53ULL;
        h ^= h >> 33;
        return static_cast<std::size_t>(h);
    }
};

#endif