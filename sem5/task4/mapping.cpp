#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    int size = 512;
    vector<int> vec(size);
    for (int i = 0; i < size; ++i) {
        vec[i] = i;
    }
    random_shuffle(vec.begin(), vec.end());
    for (int i = 0; i < size; ++i) {
        cout << vec[i] / 64 << ' ' << vec[i] % 64 / 8 << ' ' << vec[i] % 8 << " 0"<< endl;
    }
}
