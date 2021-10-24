#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <cmath>

class frac {
public:
    frac() {}

    frac(int a, int b) : a(a), b(b) {
        norm();
    }

    friend frac operator*(frac f1, frac f2) {
        return frac(f1.a * f2.a, f1.b * f2.b);
    }

    frac &operator/=(frac f) {
        a *= f.b;
        b *= f.a;
        norm();
        return *this;
    }

    friend frac operator/(frac f1, frac f2) {
        return frac(f1.a * f2.b, f1.b * f2.a);
    }

    frac &operator*=(frac f) {
        a *= f.a;
        b *= f.b;
        norm();
        return *this;
    }

    friend frac operator+(frac f1, frac f2) {
        return frac(f1.a * f2.b + f1.b * f2.a, f1.b * f2.b);
    }

    bool isNull() {
        return a == 0;
    }

    frac &operator+=(frac f) {
        a = a * f.b + f.a * b;
        b *= f.b;
        norm();
        return *this;
    }

    friend frac operator-(frac f1, frac f2) {
        return frac(f1.a * f2.b - f1.b * f2.a, f1.b * f2.b);
    }

    friend frac operator-(frac f1) {
        return frac(-f1.a, f1.b);
    }

    frac &operator-=(frac f) {
        a = a * f.b - f.a * b;
        b *= f.b;
        norm();
        return *this;
    }

    bool operator==(frac f) {
        return a == f.a && b == f.b;
    }

    bool operator!=(frac f) {
        return !(*this == f);
    }

    void Print() {
        std::cout << a << '/' << b << ' ';
    }

    void PrintLaTeX() {
        if (b == 1) {
            std::cout << a;
        } else {
            std::cout << "\\frac{" << a << "}{" << b << "}";
        }
    }

private:
    int a, b;

    int gcd(int a, int b) {
        while (b) {
            a %= b;
            std::swap(a, b);
        }
        return a;
    }

    void norm() {
        int g = gcd(abs(a), abs(b));
        a /= g;
        b /= g;
        if (a < 0 && b < 0) {
            a *= -1;
            b *= -1;
        }
        if (b < 0) {
            a *= -1;
            b *= -1;
        }
    }
};

class Matrix {
public:
    Matrix(int n, int m) : n(n), m(m) {
        M.assign(m, std::vector<frac>(n, frac(0, 1)));
    }

    Matrix() {}

    void get() {
        std::cin >> m >> n;
        M.assign(m, std::vector<frac>(n, frac(0, 1)));
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                int x;
                std::cin >> x;
                M[i][j] = frac(x, 1);
            }
        }
    }

    void Print(bool LaTeX = true) {
        if (LaTeX) {
            std::cout << "$$\n";
            std::cout << "\\left(\\begin{matirx}\n";
            for (int i = 0; i < m; ++i) {
                for (int j = 0; j < n; ++j) {
                    M[i][j].PrintLaTeX();
                    std::cout << (j + 1 == n ? " \\\\\n" : " & ");
                }
            }
            std::cout << "\\end{matrix}\\right)\n";
            std::cout << "\\to \n";
            std::cout << "$$\n";
        } else {
            for (int i = 0; i < m; ++i) {
                for (int j = 0; j < n; ++j) {
                    M[i][j].Print();
                }
                std::cout << '\n';
            }
            std::cout << "\n\n";
        }
    }

    void Gauss(bool LaTeX = true) {
        std::vector <std::pair<int, int>> ind;
        int level = 0;
        for (int i = 0; i < n; ++i) {
            int x = RowWithout0(i, level);
            if (x != -1) {
                Swap(level, x);
                Print(LaTeX);
                ind.push_back(std::make_pair(x, i));
            } else {
                continue;
            }
            for (int j = level + 1; j < m; ++j) {
                Add(j, level, -(M[j][i] / M[level][i]));
            }
            level++;
            Print(LaTeX);
            if (level == m) break;
        }
        std::reverse(ind.begin(), ind.end());
        for (auto x: ind) {
            if (M[x.first][x.second].isNull()) {
                continue;
            }
            for (int j = x.first - 1; j >= 0; --j) {
                Add(j, x.first, -(M[j][x.second] / M[x.first][x.second]));
            }
            Print(LaTeX);
        }
        for (auto x: ind) {
            if (!M[x.first][x.second].isNull()) {
                Div(x.first, M[x.first][x.second]);
            }
        }
        Print(LaTeX);
    }

private:
    std::vector <std::vector<frac>> M;
    int m, n;

    void Swap(int a, int b) {
        std::swap(M[a], M[b]);
    }

    void Add(int a, int b, frac k) { // M[a] += M[b] * k
        for (int i = 0; i < n; ++i) {
            M[a][i] += (M[b][i] * k);
        }
    }

    void Mul(int a, frac k) {
        for (int i = 0; i < n; ++i) {
            M[a][i] *= k;
        }
    }

    void Div(int a, frac k) {
        for (int i = 0; i < n; ++i) {
            M[a][i] /= k;
        }
    }

    int RowWithout0(int c, int level) {
        for (int i = level; i < m; ++i) {
            if (!M[i][c].isNull()) {
                return i;
            }
        }
        return -1;
    }
};

// Неболшая документация
/*
    Программа считывает матрицу из файла "input.txt" при локальном запуске,
    для этого добавлена локальная переменная при компиляции.

    Метод Гаусса работает за O(mn), где m, n - размеры матрицы.
 */

int main() {
#ifdef lolipop
    freopen("input.txt", "r", stdin);
#endif
    Matrix A;
    A.get();
    A.Gauss(false);
    return 0;
}
