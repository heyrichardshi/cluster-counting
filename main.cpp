#include <iostream>
#include <fstream>
#include <cmath>

#include "random.hpp"

using namespace std;

const int m = 128;
const int n = m;
const int none = m * n;
int clusterCount = 0;
int space[m][n] = {0};
int check[m][n] = {0};
int labels[m][n] = {0};
int bigoh = 0;

void generateSpace(double p) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            space[i][j] = p > RandomUniform() ? 1 : 0;
        }
    }
}

void printSpace() {
    cout << "\nSPACE:\n" << endl;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (labels[i][j] < 10) {
                cout << " ";
            }
            cout << space[i][j] << " ";
        }
        cout << endl;
    }
}

void printLabels() {
    cout << "\nLABELS:\n" << endl;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (labels[i][j] == none) {
                cout << " ~ ";
                continue;
            }
            if (labels[i][j] < 10) {
                cout << " ";
            }
            cout << labels[i][j] << " ";
        }
        cout << endl;
    }
}

void propagate(int a, int b, int c) { // a = row, b = column, c = labels
//    cout << "            propagate (" << a << ", " << b << ", " << c << ")" << endl;
    bigoh++;
    if (labels[a][b] == none) {
        return;
    }

    labels[a][b] = c;
    check[a][b] = 1;

    if ((a != m - 1) && (labels[a + 1][b] != c)) { // if not on bottom
        propagate(a + 1, b, c);
    }
    if ((b != n - 1) && (labels[a][b + 1] != c)) {
        propagate(a, b + 1, c);
    }
}

int label(int a, int b) {
//    cout << "labeling (" << a << ", " << b << ")" << endl;
    bigoh++;
    int right = none;
    int bottom = none;

    if (check[a][b] == 1) {
        return labels[a][b];
    }

    if ((a == m - 1) && (b == n - 1)) { // bottom right corner
        if (space[a][b] == 1) { // if pore exists
            labels[a][b] = ++clusterCount;
        }
        else {
            labels[a][b] = none;
        }
        check[a][b] = 1;
        return labels[a][b];
    }

    if (a != m - 1) {
        // recurse on bottom neighbor
        bottom = label(a + 1, b);
    }
    if (b != n - 1) {
        // recurse on right neighbor
        right = label(a, b + 1);
    }

    if (space[a][b] == 1) {
        if ((right == none) && (bottom == none)) {
            labels[a][b] = ++clusterCount;
        }
        else if (right == none) {
            labels[a][b] = bottom;
        }
        else if (bottom == none) {
            labels[a][b] = right;
        }
        else if (right == bottom) {
            labels[a][b] = right;
        }
        else {
            if (right < bottom) {
                labels[a][b] = right;
                propagate(a + 1, b, right);
            }
            else {
                labels[a][b] = bottom;
                propagate(a, b + 1, bottom);
            }
        }
    }
    else {
        labels[a][b] = none;
    }

    check[a][b] = 1;

    return labels[a][b];
}

int countClusters() {
    int * counts = new int[m * n];
    int count = 0;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            if (labels[i][j] != m * n) {
                counts[labels[i][j]]++;
            }
        }
    }

    for (int i = 0; i < m * n; i++) {
        if (counts[i]) {
            count++;
        }
    }

    return count;
}

int main() {
    float p = 0.7;
    generateSpace(p);
    printSpace();
    label(0,0);
    cout << "p = " << p << " & clusters = " << countClusters() << endl;
    printLabels();
    cout << "Ran " << bigoh << " times" << endl;
    cout << "Size " << m * n << endl;
    return 0;
}