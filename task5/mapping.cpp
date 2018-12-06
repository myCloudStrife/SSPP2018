#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    int size = 128;
    vector<int> vec(size);
    for (int i = 0; i < size; ++i) {
        vec[i] = i;
    }
    random_shuffle(vec.begin(), vec.end());
    for (int i = 0; i < 125; ++i) {
        cout << vec[i] / 32 << ' ' << vec[i] % 32 / 8 << ' ' << vec[i] % 8 << " 0"<< endl;
    }
}
