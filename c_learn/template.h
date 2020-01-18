#ifndef TEMPLATE_H
#define TEMPLATE_H

#include <iostream>
#include <vector>

template<typename T>
void print_array(const std::vector<T>& b) {
    for (T item : b)
        std::cout << item << ",";
}

#endif