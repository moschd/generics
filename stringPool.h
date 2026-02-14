#ifndef DTM_GENERICS_STRING_POOL_H
#define DTM_GENERICS_STRING_POOL_H

#include <unordered_set>

class StringPool
{
public:
    // This creates the one and only instance of the pool
    static StringPool &instance()
    {
        static StringPool instance;
        return instance;
    }

    const std::string *intern(const std::string &text)
    {
        auto result = pool_.insert(text);
        return &(*result.first);
    }

private:
    StringPool() = default; // Constructor is private so no one else can create it
    std::unordered_set<std::string> pool_;

    // Copying the pool is not allowed.
    StringPool(StringPool const &) = delete;
    void operator=(StringPool const &) = delete;
};

#endif