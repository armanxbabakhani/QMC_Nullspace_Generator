#ifndef Null_Generator  // Include guard
#define Null_Generator

#include <vector>

/*template<typename T>
class Indexed_Array {
public:
    // Constructor
    Indexed_Array(std::size_t size) : data(size) {}

    // Add a vector to the ith element
    void add(std::size_t i, const std::vector<T>& vector) {
        data[i].push_back(vector);
    }

    // Get the vectors at the ith element
    const std::vector<std::vector<T>> & get(std::size_t i) const {
        return data[i];
    }

private:
    std::vector<std::vector<T>> data;
};*/

template<typename T>
class Indexed_Array {
public:
    // Constructor
    Indexed_Array(std::size_t size) : data(size) {}

    // Add a vector to the ith element
    void add(std::size_t i, const std::vector<T>& vector) {
        data[i].push_back(vector);
    }

    // Get the vectors at the ith element
    const std::vector<std::vector<T>>& get(std::size_t i) const {
        return data[i];
    }

private:
    std::vector<std::vector<T>> data;
};


#endif  // Null_Generator