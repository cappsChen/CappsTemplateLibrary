- [1, Math](#1--math)
  * [MathPackage](#mathpackage)
    + [Sieve](#sieve)
    + [Combinatorics](#combinatorics)
  * [Polynomial](#polynomial)
    + [**FastFourierTransform**](#--fastfouriertransform--)
    + [NumberTheoreticTransform](#numbertheoretictransform)
  * [ExGcd and ExCrt](#exgcd-and-excrt)
    + [ExGcd](#exgcd)
    + [ExCrt](#excrt)
  * [Linear Algebra](#linear-algebra)
    + [GaussianElimination](#gaussianelimination)
    + [MatrixMultiply](#matrixmultiply)
    + [LinearBasis](#linearbasis)
  * [Random Number Algorithm](#random-number-algorithm)
    + [RandomNumber](#randomnumber)
    + [MillerRabin](#millerrabin)
    + [PollardRho](#pollardrho)
  * [Clarketech](#clarketech)
    + [FastInputOutput](#fastinputoutput)
    + [next k-bit](#next-k-bit)
- [2, Dynamic Programming](#2--dynamic-programming)
  * [MultipleBackpacks](#multiplebackpacks)
- [3, Sorting](#3--sorting)
  * [RadixSort](#radixsort)
  * [CountingSort](#countingsort)
- [4, String](#4--string)
    + [MinimalString](#minimalstring)
  * [String F4](#string-f4)
    + [StringHash](#stringhash)
    + [Kmp](#kmp)
    + [Z-Function](#z-function)
    + [Manacher](#manacher)
  * [Advanced Strings](#advanced-strings)
    + [AcAutomaton](#acautomaton)
    + [SuffixAutomaton](#suffixautomaton)
    + [ExSuffixAutomaton](#exsuffixautomaton)
    + [PalindromeAutomanton](#palindromeautomanton)
    + [SuffixArray](#suffixarray)
- [5, DataStructures](#5--datastructures)
  * [BaseDataStructures](#basedatastructures)
    + [HashTable](#hashtable)
    + [HashTableByPointer](#hashtablebypointer)
    + [**SparseTable(RMQ)**](#--sparsetable-rmq---)
    + [DeletableHeap](#deletableheap)
    + [**DisjointSetUnion**](#--disjointsetunion--)
    + [Weighted**DisjointSetUnion**](#weighted--disjointsetunion--)
  * [Segment Tree](#segment-tree)
    + [RangeQuerySegmentTree](#rangequerysegmenttree)
    + [LazySegmentTree](#lazysegmenttree)
    + [SegmentTreeMerge](#segmenttreemerge)
    + [PresidentTree](#presidenttree)
    + [LiChaoTree](#lichaotree)
  * [Binary Indexed Tree](#binary-indexed-tree)
    + [BIT](#bit)
    + [RangeBIT](#rangebit)
    + [MatBIT](#matbit)
    + [RangeMatBIT](#rangematbit)
  * [Block](#block)
- [6, Tree Theory](#6--tree-theory)
  * [TreePre](#treepre)
  * [VirtualTreePre](#virtualtreepre)
  * [DsuOnTree](#dsuontree)
  * [CentroidDecomposition](#centroiddecomposition)
  * [TreeDiameter](#treediameter)
  * [DecareTree](#decaretree)
  * [TreeHash](#treehash)
- [7, Graph Theory](#7--graph-theory)
  * [Dijkstra](#dijkstra)
  * [Bi-GraphColor](#bi-graphcolor)
  * [RingTree](#ringtree)
  * [TopSort](#topsort)
  * [Connectivity](#connectivity)
    + [StronglyConnectedComponent](#stronglyconnectedcomponent)
    + [TwoSat](#twosat)
    + [VertexBiconnectedComponent](#vertexbiconnectedcomponent)
    + [EdgeBiconnectedComponent](#edgebiconnectedcomponent)
  * [Flow](#flow)
    + [MaxFlow](#maxflow)
    + [MinCostFlow](#mincostflow)
- [8, Calculate Geometry](#8--calculate-geometry)
  * [Point](#point)
  * [Line](#line)
    + [Point-Line](#point-line)
    + [Line-Line](#line-line)
  * [Polygon](#polygon)
    + [ConvexHull](#convexhull)
    + [Area](#area)
    + [Point-Polygon](#point-polygon)
    + [Line-Polygon](#line-polygon)
    + [Polygon-Polygon](#polygon-polygon)
  * [Examples](#examples)
    + [ClosestPair](#closestpair)
    + [DelaunayTriangulation](#delaunaytriangulation)





<a name="a6i3b"></a>

# 1, Math
<a name="80d6f6b4"></a>
## MathPackage
<a name="ac7c11a9"></a>
### Sieve
```cpp
struct Sieve {
    std::vector<int> P, v;

    Sieve(int n) : v(n) {
        for (int i = 2; i < n; i++) {
            if (v[i] == 0) {
                P.push_back(i);
                v[i] = i;
            }
            for (int j = 0; j < P.size() and i * P[j] < n; j++) {
                v[i * P[j]] = P[j];
                if (P[j] == v[i])
                    break;
            }
        }
    }

    // 求所有约数
    auto getDiv(int x) const {
        std::vector<int> _div(1, 1);
        while (x > 1) {
            int D = v[x];
            int l = 0, r = _div.size();
            while (x % D == 0) {
                for (int k = l; k < r; k++)
                    _div.push_back(_div[k] * D);
                x /= D, l = r, r = _div.size();
            }
        }
        return _div;
    }
};
```
<a name="69c0b657"></a>
### Combinatorics
```cpp
template <class T, T P>
class Comb {
    static constexpr int multip(const int &a, const int &b) {
        return 1ll * a * b % P;
    }
    static constexpr i64 multip(const i64 &a, const i64 &b) {
        i64 res = a * b - i64(1.L * a * b / P) * P;
        res %= P;
        res += (res < 0 ? P : 0);
        return res;
    }

    int n;
    std::vector<T> _jc, _ijc, _inv;

public:
    constexpr Comb() : n{0}, _jc{1}, _ijc{1}, _inv{0} {}
    Comb(int n) : Comb() {
        init(n);
    }

    static constexpr T powp(T a, i64 mi) {
        T ans = 1;
        for (; mi; mi >>= 1, a = multip(a, a))
            if (mi & 1)
                ans = multip(ans, a);
        return ans;
    }

    void init(int m) {
        m = std::min(m, P - 1);
        if (m <= n)
            return;

        _jc.resize(m + 1);
        _ijc.resize(m + 1);
        _inv.resize(m + 1);

        for (int i = n + 1; i <= m; i++) {
            _jc[i] = multip(i, _jc[i - 1]);
        }
        _ijc.back() = powp(_jc.back(), P - 2);
        for (int i = m; i > n; i--) {
            _ijc[i - 1] = multip(i, _ijc[i]);
            _inv[i] = multip(_ijc[i], _jc[i - 1]);
        }

        n = m;
    }

    T jc(int x) {
        if (x > n)
            init(x << 1);
        return _jc[x];
    }
    T ijc(int x) {
        if (x > n)
            init(x << 1);
        return _ijc[x];
    }
    T inv(int x) {
        if (x > n)
            init(x << 1);
        return _inv[x];
    }

    T A(int a, int b) {
        if (a < b or b < 0)
            return 0;
        return multip(jc(a), ijc(a - b));
    }
    T C(int a, int b) {
        if (a < b or b < 0)
            return 0;
        return multip(A(a, b), ijc(b));
    }
};
constexpr int P = 998244353;
Comb<int, P> comb;
```
<a name="8b9b9392"></a>
## Polynomial
<a name="8e8de95a"></a>
### FastFourierTransform
```cpp
#define x first
#define y second

template <class T>
struct Complex : public std::pair<T, T> {
    constexpr Complex(T x = T(), T y = T()) : std::pair<T, T>(x, y) {}

    constexpr Complex operator+(const Complex &a) const {
        return Complex(this->x + a.x, this->y + a.y);
    }
    constexpr Complex operator-(const Complex &a) const {
        return Complex(this->x - a.x, this->y - a.y);
    }
    constexpr Complex operator*(const Complex &a) const {
        return Complex(this->x * a.x - this->y * a.y, this->x * a.y + this->y * a.x);
    }
};
using ld = double;
using Comp = Complex<ld>;

constexpr ld pi = acosl(-1);

std::vector<Comp> w[2];
std::vector<int> r;

constexpr void init(int _log) {
    if (r.size() == (1 << _log))
        return;

    int n = 1 << _log;
    r.assign(n, 0);
    for (int i = 1; i < n; i++)
        r[i] = (r[i >> 1] >> 1) | ((i & 1) << (_log - 1));

    w[0].assign(n, Comp());
    w[1].assign(n, Comp());
    for (int i = 0; i < n; i++) {
        auto th = 1.L * pi * i / n;
        w[0][i] = Comp(cosl(th), sinl(th));
        w[1][i] = Comp(cosl(th), -sinl(th));
    }
}

void fft(std::vector<Comp> &a, int op) {
    int n = a.size();
    init(std::__lg(n));
    for (int i = 0; i < n; i++)
        if (i < r[i])
            std::swap(a[i], a[r[i]]);
    for (int mid = 1; mid < n; mid <<= 1) {
        const int d = n / mid;
        for (int R = mid << 1, j = 0; j < n; j += R) {
            for (int k = 0; k < mid; k++) {
                Comp x = a[j + k], y = w[op][d * k] * a[j + mid + k];
                a[j + k] = x + y;
                a[j + mid + k] = x - y;
            }
        }
    }
}

template <class T>
struct Polynomial : public std::vector<T> {
    explicit constexpr Polynomial(int n = 0, T v = T()) : std::vector<T>(n, v) {}
    explicit constexpr Polynomial(const std::vector<T> &a) : std::vector<T>(a) {}
    constexpr Polynomial(const std::initializer_list<T> &a) : std::vector<T>(a) {}

    template <class IT, class = std::_RequireInputIter<IT>>
    explicit constexpr Polynomial(IT first, IT last) : std::vector<T>(first, last) {}

    constexpr friend Polynomial operator*(const Polynomial &a, const Polynomial &b) {
        if (a.size() == 0 or b.size() == 0)
            return Polynomial();
        int n = a.size() + b.size() - 1;
        int _log = std::__lg(2 * n - 1);
        int s = 1 << _log;
        if (std::min(a.size(), b.size()) < 128) {
            Polynomial res(n);
            for (int i = 0; i < a.size(); i++)
                for (int j = 0; j < b.size(); j++)
                    res[i + j] += a[i] * b[j];
            return res;
        }

        std::vector<Comp> p(s), q(s);
        for (int i = 0; i < a.size(); i++)
            p[i].x = a[i];
        for (int i = 0; i < b.size(); i++)
            q[i].x = b[i];

        fft(p, 0), fft(q, 0);
        for (int i = 0; i < s; i++)
            p[i] = p[i] * q[i];
        fft(p, 1);

        Polynomial res(n);
        for (int i = 0; i < n; i++)
            res[i] = p[i].x / s; // 默认浮点数
        return res;
    }
};

using Poly = Polynomial<ld>;
```
<a name="5dd0c9c0"></a>
### NumberTheoreticTransform
```cpp
template <class T, T P>
class Polynomial : public std::vector<T> {
    static constexpr int multip(const int &a, const int &b) {
        return 1ll * a * b % P;
    }
    static constexpr i64 multip(const i64 &a, const i64 &b) {
        i64 res = a * b - i64(1.L * a * b / P) * P;
        res %= P;
        res += (res < 0 ? P : 0);
        return res;
    }

    static constexpr T add(T a, const T &b) {
        a += b;
        a -= (a >= P ? P : 0);
        return a;
    }

    static constexpr T sub(T a, const T &b) {
        a -= b;
        a += (a < 0 ? P : 0);
        return a;
    }

    static constexpr T powp(T a, i64 mi) {
        T ans = 1;
        for (; mi; mi >>= 1, a = multip(a, a))
            if (mi & 1)
                ans = multip(ans, a);
        return ans;
    }
    static std::vector<T> w;

    static void initW(int _log) {
        const int r = 1 << _log;
        if (w.size() >= r)
            return;

        w.assign(r, 0);
        w[r >> 1] = 1;
        T s = powp(3, (P - 1) >> _log);
        for (int i = r / 2 + 1; i < r; i++)
            w[i] = multip(w[i - 1], s);
        for (int i = r / 2 - 1; i > 0; i--)
            w[i] = w[i * 2];
    }

public:
    using std::vector<T>::vector;

    constexpr friend void dft(Polynomial &a) {
        const int n = a.size();
        for (int k = n >> 1; k; k >>= 1) {
            for (int i = 0; i < n; i += k << 1) {
                for (int j = 0; j < k; j++) {
                    T v = a[i + j + k];
                    a[i + j + k] = multip(sub(a[i + j], v), w[k + j]);
                    a[i + j] = add(a[i + j], v);
                }
            }
        }
    }

    constexpr friend void idft(Polynomial &a) {
        const int n = a.size();
        for (int k = 1; k < n; k <<= 1)
            for (int i = 0; i < n; i += k << 1)
                for (int j = 0; j < k; j++) {
                    T u = a[i + j];
                    T &&v = multip(a[i + j + k], w[j + k]);
                    a[i + j + k] = sub(u, v);
                    a[i + j] = add(u, v);
                }
        T &&val = P - (P - 1) / n;
        for (int i = 0; i < n; i++)
            a[i] = multip(a[i], val);
        std::reverse(a.begin() + 1, a.end());
    }

    constexpr friend Polynomial operator*(Polynomial a, Polynomial b) {
        if (a.size() == 0 or b.size() == 0)
            return Polynomial();
        int n = a.size() + b.size() - 1;
        int _log = std::__lg(2 * n - 1);
        int s = 1 << _log;
        if (((P - 1) & (s - 1)) != 0 or std::min(a.size(), b.size()) < 128) {
            Polynomial res(n);
            for (int i = 0; i < a.size(); i++)
                for (int j = 0; j < b.size(); j++)
                    res[i + j] = add(res[i + j], multip(a[i], b[j]));
            return res;
        }

        initW(_log);
        a.resize(s), b.resize(s);
        dft(a), dft(b);
        for (int i = 0; i < s; i++)
            a[i] = multip(a[i], b[i]);
        idft(a);
        return a.resize(n), a;
    }

    constexpr friend Polynomial inv(const Polynomial &a) {
        int n = a.size();
        if (n == 1)
            return {powp(a[0], P - 2)};

        Polynomial half(a.begin(), a.begin() + (n + 1) / 2);
        Polynomial &&b = inv(half), c = a * b;
        for (T &x : c)
            x = (x == 0 ? 0 : P - x); // ?
        c[0] = add(c[0], 2);
        c = c * b;
        return c.resize(n), c;
    }

    constexpr friend Polynomial ln(const Polynomial &a) {
        int n = a.size();

        Polynomial b(n, 0);
        for (int i = 1; i < n; i++)
            b[i - 1] = multip(i, a[i]);
        b = b * inv(a);
        b.resize(n);

        std::vector<T> _inv(n);
        _inv[1] = 1;
        for (int i = 2; i < n; i++) {
            _inv[i] = multip(P - P / i, _inv[P % i]);
        }
        for (int i = n - 1; i; i--)
            b[i] = multip(b[i - 1], _inv[i]);
        b[0] = 0;
        return b;
    }

    constexpr friend Polynomial exp(const Polynomial &a) {
        int n = a.size();
        if (n == 1)
            return {1};

        Polynomial half(a.begin(), a.begin() + (n + 1) / 2);
        Polynomial b = exp(half);
        b.resize(n);
        Polynomial c = ln(b);
        for (int i = 0; i < n; i++)
            c[i] = sub(a[i], c[i]);
        c[0] = add(c[0], 1);
        c = c * b;
        return c.resize(n), c;
    }
};
template <class T, T P>
std::vector<T> Polynomial<T, P>::w;

constexpr int P = 998244353;

using Poly = Polynomial<int, P>;
```
```cpp
V <= 1e9 : 1004535809
V <= 1e15 : 1337006139375617
V <= 4e18 : 4179340454199820289
```
<a name="b625e5b8"></a>
## ExGcd and ExCrt
<a name="4bd36b76"></a>
### ExGcd 
```cpp
template <class T>
struct ExGcd {
    T operator()(const T &a, const T &b, T &x, T &y) {
        if (b == 0)
            return (x = 1, y = 0, a);
        T g = (*this)(b, a % b, y, x);
        y -= a / b * x;
        return g;
    }
};
ExGcd<int> exgcd;
```

对于方程 $ax+by=c$ , 调用 **exgcd **, 求出 $x_0$ 和 $y_0$ 使得 $ax_0+by_0=\gcd(a, b)$

则在 $\gcd(a,b)\mid c$ 的情况下有通解

$x = x_0 \times \frac{c}{\gcd(a, b)}+k\times \frac{b}{\gcd(a, b)}
\\ \ 
\\  \\
y = y_0 \times \frac{c}{\gcd(a, b)}-k\times \frac{a}{\gcd(a, b)}$
<a name="1da48602"></a>
### ExCrt
```cpp
template <class T, class G>
struct ExCrt : public ExGcd<T> {
std::vector<std::pair<T, T>> q;
void insert(T a, T mod) {
    q.push_back({a, mod});
}

// 方程组 x ≡ a (模 mod) 返回最小正解
// 无解返回 -1
T get() {
    T res = 0, M = 1;
    for (auto [a, mod] : q) {
        T r = (a - res) % mod;
        r += (r < 0 ? mod : 0);

        T x, y;
        T g = (*this)(M, mod, x, y);
        if (r % g){
            q.clear();
            return -1;
        }

        x = (G(x) * r / g % (mod / g));
        x += (x < 0 ? mod / g : 0);

        T Last = M;
        M = M / g * mod;
        res = (G(x) * Last % M + res) % M;
    }
    q.clear();
    return res;
}
};
ExCrt<i64, __int128> crt;
```
<a name="74659e3e"></a>
## Linear Algebra
<a name="25634846"></a>
### GaussianElimination
```cpp
using ld = double;
constexpr ld eps = 1e-9;
int sgn(const ld &a) {
    return (a < -eps ? -1 : a > eps);
}
std::string gauss(std::vector<std::vector<ld>> &a) { // 传入增广矩阵
    int n = a.size();
    int c = 0, r = 0;
    for (c = 0, r = 0; c < n; c++) { // c列r行，遍历列
        int tmp = r;
        for (int i = r; i < n; i++) // 寻找列主元
            if (sgn(a[i][c]))
                tmp = i;
        if (sgn(a[tmp][c]) == 0) // 当前列全为0
            continue;

        std::swap(a[tmp], a[r]); // 交换列主元

        for (int i = n; i >= c; i--) // 倒序处理
            a[r][i] /= a[r][c];

        for (int i = r + 1; i < n; i++)
            if (sgn(a[i][c]))
                for (int j = n; j >= c; j--)
                    a[i][j] -= a[r][j] * a[i][c];
        r++;
    }
    if (r < n) {
        for (int i = r; i < n; i++)
            if (sgn(a[i][n]))
                return "NoSolution";
        return "InfSolution";
    }

    // 解放在 a[i][n]  (0<= i < n)
    for (int i = n - 1; i >= 0; i--)
        for (int j = i + 1; j < n; j++)
            a[i][n] -= a[j][n] * a[i][j];
    return "OK";
}
```
<a name="D5M3R"></a>
### MatrixMultiply
```cpp
template <class T, T P>
struct Multiply {
    static constexpr int multip(const int &a, const int &b) {
        return 1ll * a * b % P;
    }
    static constexpr i64 multip(const i64 &a, const i64 &b) {
        i64 res = a * b - i64(1.L * a * b / P) * P;
        res %= P;
        res += (res < 0 ? P : 0);
        return res;
    }

    constexpr T operator()(const T &x, const T &y) const {
        return multip(x, y);
    }
};

template <class T, T P>
struct Add {
    constexpr T operator()(const T &x, const T &y) const {
        T res = x + y;
        return res >= P ? res - P : res;
    }
};

template <class T, class Cmp = std::less<T>>
struct Min {
    const Cmp cmp{};
    constexpr T operator()(const T &a, const T &b) const {
        return std::min(a, b, cmp);
    }
};

constexpr int P = 998244353;
constexpr int D = 200;

template <class T>
using Matrix = std::array<std::array<T, D>, D>;

template <class T, class MergeIn = Multiply<T, P>, class MergeOut = Add<T, P>>
Matrix<T> operator*(const Matrix<T> &a, const Matrix<T> &b) {
    static const MergeIn mergeIn{};
    static const MergeOut mergeOut{};
    static const int n = a.size();

    Matrix<T> c = {};
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            for (int j = 0; j < n; j++) {
                c[i][j] = mergeOut(c[i][j], mergeIn(a[i][k], b[k][j]));
            }
        }
    }

    return std::move(c);
}

using Mat = Matrix<int>;
```
<a name="ee0d510c"></a>
### LinearBasis
```cpp
template <class T>
struct LinearBasis {
    std::vector<T> b;
    const int logn;
    int flag = 0, r = 0;
    // flag 表示线性基能不能异或出 0
    // r 表示线性基的秩

    LinearBasis(T n) : logn(std::__lg(n)), b(std::__lg(n) + 1) {}

    void insert(T x) { // 插入一个数到线性基里
        for (int i = logn; i >= 0; i--) {
            if (x >> i & 1) {
                if (!b[i]) { // 线性基里没有第i位的项
                    b[i] = x;
                    r++;
                    return;
                }
                x ^= b[i]; // 线性基里有,则异或掉第i位
            }
        }
        flag = true;
    }

    bool check(T x) { // 询问线性基能不能异或出 x
        for (int i = logn; i >= 0; i--)
            if (x >> i & 1) {
                if (!b[i]) // 线性基里没有第i位的项
                    return false;
                x ^= b[i]; // 线性基里有,则异或掉第i位
            }
        return true;
    }
};
```
<a name="5599756b"></a>
## Random Number Algorithm
<a name="4067b740"></a>
### RandomNumber
```cpp
template <class T>
struct Rand {
    std::mt19937 myrand;
    Rand(const i64 seed = time(0)) : myrand(seed) {}
    T operator()(T l, T r) {
        return std::uniform_int_distribution<T>(l, r)(myrand);
    }
};
Rand<int> rd;

//std::mt19937_64 rng(std::chrono::steady_clock::now().time_since_epoch().count());
```
<a name="51780be9"></a>
### MillerRabin
```cpp
/*
维基百科 :
n < 4e9, Prime = [2, 7, 61]
n < 3e14, Prime = [2, 3, 5, 7, 11, 13, 17]
n < 3e18, Prime = [2, 3, 5, 7, 11, 13, 17, 19, 23]
n < 3e23, Prime = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
*/
template <class T>
struct MillerRabin {
    const std::vector<int> Prime;
    MillerRabin() : Prime({2, 3, 5, 7, 11, 13, 17, 19, 23}) {}

    static constexpr int mulp(const int &a, const int &b, const int &P) {
        return 1ll * a * b % P;
    }
    static constexpr i64 mulp(const i64 &a, const i64 &b, const i64 &P) {
        i64 res = a * b - i64(1.L * a * b / P) * P;
        res %= P;
        res += (res < 0 ? P : 0);
        return res;
    }

    static constexpr T powp(T a, T mi, const T &mod) {
        T ans = 1;
        for (; mi; mi >>= 1) {
            if (mi & 1)
                ans = mulp(ans, a, mod);
            a = mulp(a, a, mod);
        }
        return ans;
    }

    constexpr bool operator()(const T &v) { // 判断v是不是质数
        if (v < 2 or v != 2 and v % 2 == 0)
            return false;
        T s = v - 1;
        while (!(s & 1))
            s >>= 1;
        for (int x : Prime) {
            if (v == x)
                return true;
            T t = s, m = powp(x, s, v);
            while (t != v - 1 and m != 1 and m != v - 1)
                m = mulp(m, m, v), t <<= 1;
            if (m != v - 1 and !(t & 1))
                return false;
        }
        return true;
    }
};
MillerRabin<i64> isp;
```
<a name="2a0e6bc2"></a>
### PollardRho
如果 $n$ 是质数 $(MillerRabbin判)$ 返回 $n$<br />否则返回 $n$ 的随机一个 $[2,n-1]$ 的因子<br />复杂度理论 $O(n^{\frac{1}{4}}\log n)$  但实际跑得快, 可以按 $O(n^{\frac{1}{4}})$ 算
```cpp
template <class T>
struct PollardRho : public MillerRabin<T> {
    std::mt19937 myrand;
    PollardRho(const i64 seed = time(0)) : myrand(seed) {}

    constexpr T rd(T l, T r) {
        return std::uniform_int_distribution<T>(l, r)(myrand);
    }

    constexpr T operator()(T n) { // 返回 n 的随机一个[2, n-1]内的因子, 或者判定是质数
        if (n == 4)
            return 2;
        MillerRabin<T> &super = *this;
        if (super(n)) // Miller-Rabbin 判质数
            return n; // 如果 n 是质数直接返回 n
        while (true) {
            T c = rd(1, n - 1);
            auto f = [&](T x) { return (super.mulp(x, x, n) + c) % n; };
            T t = 0, r = 0, p = 1, q;
            do {
                for (int i = 0; i < 128; i++) {
                    t = f(t), r = f(f(r));
                    if (t == r || (q = super.mulp(p, std::abs(t - r), n)) == 0)
                        break;
                    p = q;
                }
                T d = std::gcd<T>(p, n);
                if (d > 1)
                    return d;
            } while (t != r);
        }
    }
};
PollardRho<i64> findDiv;
```
<a name="fabKD"></a>
## Clarketech
<a name="I9LvY"></a>
### FastInputOutput
```cpp
#define gc_() (tt == ss and (tt = (ss = In) + fread(In, 1, 1 << 20, stdin), ss == tt) ? EOF : *ss++)
char In[1 << 20], *ss = In, *tt = In;

template <class T>
inline void _read(T &res) {
    bool fu = 0;
    res = 0;
    char ch = gc_();
    while (!isdigit(ch)) {
        if (ch == '-')
            fu = 1;
        ch = gc_();
    }
    while (isdigit(ch))
        res = (res << 3) + (res << 1) + ch - '0', ch = gc_();
    if (fu)
        res = -res;
}

template <class T>
inline std::istream &operator>(std::istream &in, T &a) {
    _read(a);
    return in;
}

template <class T> // 用了快写就只能用putchar 取代 <<
inline std::ostream &operator<(std::ostream &out, const T &a) {
    if (a < 0) {
        putchar('-');
        if (-a > 9)
            out < (-a / 10);
        putchar(-a % 10 + '0');
        return out;
    }
    if (a > 9)
        out < (a / 10);
    putchar(a % 10 + '0');
    return out;
}
```
<a name="himK8"></a>
### next k-bit
```cpp
// end with return unsigned(-1)
unsigned getNext(unsigned x) {
    unsigned &&b = x & -x;
    unsigned &&t = x + b;
    unsigned &&c = t & -t;
    unsigned &&m = (c / b >> 1) - 1;
    return t | m;
}
```
<a name="c46e8dea"></a>
# 2, Dynamic Programming
<a name="bc338d83"></a>
## MultipleBackpacks
```cpp
auto multiBag(std::vector<std::array<int, 3>>& goods, int S) {
    // S 总背包大小
    std::vector<i64> f(S + 1, 0);
    std::vector<int> q(S + 2, 0);
    for (auto [v, w, m] : goods) {
        // v价值, w体积, m数量
        auto calc = [&](int i, int j) {
            return f[j] + 1ll * (i - j) / w * v;
        };
        for (int up = S; up + w > S; up--) {
            int l = 1, r = 0, k = up;
            for (int x = up; x > 0; x -= w) {
                for (; k >= std::max(0ll, x - 1ll * m * w); k -= w) {
                    while (l <= r and calc(x, k) > calc(x, q[r]))
                        r--;
                    q[++r] = k;
                }
                f[x] = calc(x, q[l]);
                if (q[l] == x)
                    l++;
            }
        }
    }
    return f;
}
```

<a name="jySco"></a>
# 3, Sorting
<a name="Y9Jf8"></a>
## RadixSort
```cpp
template <int B, class T>
void radixSort(std::vector<T> &a) {
    const int mask = (1 << B) - 1, n = a.size();

    std::vector<T> b(n);
    std::vector<int> cnt(1 << B);

    T maxV = *std::max_element(begin(a), end(a));

    for (int i = 0; maxV; i += B, maxV >>= B) {
        std::fill(begin(cnt), end(cnt), 0);
        for (int j = 0; j < n; j++)
            cnt[a[j] >> i & mask] += 1;
        for (int j = 1; j < (1 << B); j++)
            cnt[j] += cnt[j - 1];
        for (int j = n - 1; j >= 0; j--)
            b[--cnt[a[j] >> i & mask]] = a[j];
        std::swap(a, b);
    }
}
```
<a name="cf384e46"></a>
## CountingSort
```cpp
// 按 a 的值返回下标排序
auto countingSort(const std::vector<int> &a) {
    int maxV = *std::max_element(begin(a), end(a));
    std::vector<int> cnt(maxV + 1);

    for (int i = 0; i < a.size(); i++)
        cnt[a[i]] += 1;

    for (int i = 1; i <= maxV; i++)
        cnt[i] += cnt[i - 1];

    std::vector<int> res(a.size());
    for (int i = a.size() - 1; i >= 0; i--)
        res[--cnt[a[i]]] = i;

    return res;
}
```

<a name="cb9e5fb0"></a>
# 4, String
<a name="47057cc1"></a>
### MinimalString
```cpp
std::vector<int> minimalString(std::vector<int> &a) {
    int n = a.size();
    int i = 0, j = 1, k = 0;
    while (k < n and i < n and j < n) {
        if (a[(i + k) % n] == a[(j + k) % n])
            k++;
        else {
            (a[(i + k) % n] > a[(j + k) % n] ? i : j) += k + 1;
            i += (i == j);
            k = 0;
        }
    }
    k = std::min(i, j);
    std::vector<int> ans(n);
    for (int i = 0; i < n; i++)
        ans[i] = a[(i + k) % n];
    return ans;
}
// 直接返回字典序最小循环同构串
```
<a name="FS2nN"></a>
## String F4
<a name="eUNpl"></a>
### StringHash
```cpp
template <int D, std::array<int, D> B, std::array<int, D> P>
struct StringHash {
    std::vector<std::array<int, D>> h;

    template <class T>
    StringHash(const T &s) : h(s.size() + 1) {
        for (int i = 0; i < s.size(); i++) {
            for (int k = 0; k < D; k++) {
                h[i + 1][k] = (1ll * h[i][k] * B[k] + s[i] + 1) % P[k];
            }
        }
    }

    // [l, r)
    std::array<int, D> get(int l, int r) {
        static std::vector<std::array<int, D>> spow(1);
        if (r - l < 0)
            throw -1;

        if (spow.size() < r - l + 1) {
            if (spow[0][0] == 0) {
                spow[0].fill(1);
            }
            int n = spow.size();
            spow.resize(r - l + 1);
            for (int i = n; i < spow.size(); i++) {
                for (int k = 0; k < D; k++) {
                    spow[i][k] = 1ll * spow[i - 1][k] * B[k] % P[k];
                }
            }
        }

        std::array<int, D> res = {};
        for (int k = 0; k < D; k++) {
            res[k] = h[r][k] - 1ll * h[l][k] * spow[r - l][k] % P[k];
            res[k] += (res[k] < 0 ? P[k] : 0);
        }
        return res;
    }
};
using Hash = StringHash<2, {133, 331}, {int(1e9 + 21), int(1e9 + 33)}>;

```
<a name="567dadb5"></a>
### Kmp
```cpp
auto kmp(const std::string &s) {
    int n = s.size();
    std::vector<int> link(n);
    for (int i = 1, j = 0; i < n; i++) {
        while (j and s[i] != s[j])
            j = link[j - 1];
        j += (s[i] == s[j]);
        link[i] = j;
    }
    return link;
}
```
<a name="ed84b7ce"></a>
### Z-Function
```cpp
auto zFunction(const std::string &s) {
    int n = s.size();
    std::vector<int> z(n);
    for (int i = 1, l = 0, r = 0; i < n; i++) {
        if (i < r)
            z[i] = std::min(z[i - l], r - i);
        while (i + z[i] < n and s[i + z[i]] == s[z[i]])
            z[i]++;
        if (i + z[i] > r)
            l = i, r = i + z[i];
    }
    return z;
}
```
<a name="237ffca0"></a>
### Manacher
```cpp
struct Manacher {
    const int n;
    std::vector<int> r, f;

    Manacher(const std::string &t)
        : n(t.size()), r(2 * t.size() + 3), f(2 * t.size() + 3) {
        std::string s = "[-";
        for (int i = 0; i < n; i++) {
            s += t[i];
            s += '-';
        }
        s.push_back(']');

        int mid = 1, far = 1;
        for (int i = 1; i < s.size(); i++) {
            r[i] = std::min(r[2 * mid - i], far - i);
            while (s[i + r[i]] == s[i - r[i]])
                r[i] += 1;
            if (far < i + r[i])
                mid = i, far = i + r[i];
            f[i + r[i] - 1] = std::max(f[i + r[i] - 1], r[i]);
        }
        for (int i = f.size() - 2; i; i--)
            f[i] = std::max(f[i], f[i + 1] - 1);
    }

    // 下标, 是否要以 +0.5 为中心
    int getPalinLenFromCenter(int i, int center) const {
        assert(!center and 0 <= i and i < n or
               center and 0 <= i and i < n - 1);

        return r[2 * (i + 1) + center] - 1;
    }

    int getPalinLenFromTail(int i) const {
        assert(0 <= i and i < n);
        return f[2 * (i + 1)];
    }
};
```

<a name="M9rbx"></a>
## Advanced Strings
<a name="9407a1e0"></a>
### AcAutomaton
```cpp
template <int Z, char Base>
struct AcAutomaton {
    std::vector<std::array<int, Z>> son;
    std::vector<std::vector<int>> ID;
    std::vector<int> link;
    int SIZE = 0, tot = 0;

    AcAutomaton(const std::vector<std::string> &s) {
        for (auto t : s)
            SIZE += t.size();
        son.resize(SIZE + 1);
        ID.resize(SIZE + 1);
        link.resize(SIZE + 1);

        for (int i = 0; i < s.size(); i++)
            insert(i, s[i]);
        build();
    }

    void insert(int id, const std::string &s) {
        int p = 0;
        for (char c : s) {
            c -= Base;
            if (!son[p][c])
                son[p][c] = ++tot;
            p = son[p][c];
        }
        ID[p].push_back(id);
    }

    void build() {
        std::queue<int> q;
        for (int &y : son[0])
            if (y) {
                q.push(y);
            }
        while (!q.empty()) {
            int x = q.front();
            q.pop();

            for (int c = 0; int &y : son[x]) {
                if (y) {
                    link[y] = son[link[x]][c];
                    q.push(y);
                } else
                    y = son[link[x]][c];
                c++;
            }
        }
    }
};
```
<a name="f0ec7728"></a>
### SuffixAutomaton
```cpp
template <int Z, char Base>
struct Sam {
    std::vector<std::array<int, Z>> son;
    std::vector<int> link, len;
    // std::vector<i64> cnt;
    int last, tot;

    Sam(int n) : son(n << 1), link(n << 1), len(n << 1) {
        // cnt.assign(n << 1, 0);
        last = tot = 0;
        link[0] = -1;
    }

    Sam(const std::string &s) : SAM(s.size()) {
        for (int i = 0; i < s.size(); i++) {
            add(s[i]);
        }
    }

    void add(char c) {
        c -= Base;
        int cur = ++tot;
        // 在此对 cur 加信息
        len[cur] = len[last] + 1;
        int v = last;
        while (v != -1 and !son[v][c])
            son[v][c] = cur, v = link[v];

        if (v == -1) link[cur] = 0;
        else {
            int q = son[v][c];
            if (len[v] + 1 == len[q]) link[cur] = q;
            else {
                int cl = ++tot;
                len[cl] = len[v] + 1;
                son[cl] = son[q];
                link[cl] = link[q];

                while (v != -1 and son[v][c] == q)
                    son[v][c] = cl, v = link[v];

                link[q] = link[cur] = cl;
            }
        }
        last = cur;
    }
};
```
<a name="4330dc52"></a>
### ExSuffixAutomaton
```cpp
template <int Z, char Base>
struct ExSam {
    std::vector<std::array<int, Z>> son;
    std::vector<int> len, link;
    // std::vector<std::set<int>> ID;
    int tot = 0, last = 0;

    ExSam(const std::vector<std::string> &S) {
        int LEN = 0;
        for (const auto &str : S) {
            LEN += size(str);
        }

        len.resize(LEN << 1);
        link.resize(LEN << 1);
        link[0] = -1;
        son.resize(LEN << 1);
        // ID.resize(LEN << 1);

        for (int i = 0; i < S.size(); i++) {
            last = 0;
            for (char c : S[i]) {
                exadd(c, i);
            }
        }
    }

private:
    void assign(int cur, int id) {
        // ID[cur].insert(id);
    }
    void exadd(char c, int id) {
        c -= Base;
        if (son[last][c]) {
            int v = last, q = son[v][c];
            if (len[q] != len[v] + 1) {
                int cl = ++tot;
                len[cl] = len[v] + 1;
                link[cl] = link[q];
                son[cl] = son[q];

                while (v != -1 and son[v][c] == q)
                    son[v][c] = cl, v = link[v];
                link[q] = cl;
                q = cl;
            }
            int cur = last = q;
            assign(cur, id);
            return;
        }

        int cur = ++tot;
        len[cur] = len[last] + 1;
        assign(cur, id);

        int v = last;
        while (v != -1 and son[v][c] == 0)
            son[v][c] = cur, v = link[v];
        if (v == -1)
            link[cur] = 0;
        else {
            int q = son[v][c];
            if (len[q] == len[v] + 1)
                link[cur] = q;
            else {
                int cl = ++tot;
                len[cl] = len[v] + 1;
                link[cl] = link[q];
                son[cl] = son[q];

                while (v != -1 and son[v][c] == q)
                    son[v][c] = cl, v = link[v];
                link[q] = link[cur] = cl;
            }
        }
        last = cur;
    }
};
```
<a name="m6UQg"></a>
### PalindromeAutomanton
```cpp
template <int Z, char Base>
struct Pam {
    std::vector<std::array<int, Z>> son;
    std::vector<int> link, len, dep, cnt;
    std::string s;

    int cur = 0, tot = 1;

    Pam(int n)
        : son(n + 2), link(n + 2), len(n + 2),
          dep(n + 2), cnt(n + 2) {
        link[0] = 1;
        len[1] = -1;
    }

    Pam(const std::string &s) : Pam(s.size()) {
        for (int i = 0; i < s.size(); i++) {
            add(i, s[i]);
        }
    }

    void assign(int cur, int id) {
        cnt[cur]++;
    }

    int getLink(int x, int i) {
        while (i - len[x] - 1 < 0 or s[i - len[x] - 1] != s[i])
            x = link[x];
        return x;
    }

    void add(int i, char c) {
        c -= Base;
        s.push_back(c);

        int v = getLink(cur, i);

        if (!son[v][c]) {
            link[++tot] = son[getLink(link[v], i)][c];
            son[v][c] = tot;
            len[tot] = len[v] + 2;
            dep[tot] = dep[link[tot]] + 1;
        }

        cur = son[v][c];
        assign(cur, i);
    }

    // Pam 的 linkTree 是 1 为根的
    auto getLinkTree() const {
        std::vector e(tot + 1, std::vector<int>());
        for (int i = 0; i <= tot; i++) {
            if (i != 1)
                e[link[i]].emplace_back(i);
        }
        return e;
    }
};
```
<a name="SuNAg"></a>
### SuffixArray
```cpp
struct SuffixArray {
    std::vector<int> sa, rk, h;

    template <class T>
    SuffixArray(const T &s)
        : n(s.size()), sa(s.size()), rk(s.size()), id(s.size()), tmp(s.size()) {

        std::iota(begin(id), end(id), 0);
        for (int i = 0; i < n; i++)
            rk[i] = s[i];

        countSort();

        for (int w = 1;; w <<= 1) {
            std::iota(begin(id), begin(id) + w + 1, n - w);
            for (int i = 0, p = w; i < n; i++)
                if (sa[i] >= w)
                    id[p++] = sa[i] - w;

            countSort();
            oldrk = rk;

            rk[sa[0]] = 0;
            for (int i = 1, p = 0; i < n; i++)
                rk[sa[i]] = equal(sa[i], sa[i - 1], w) ? p : ++p;

            if (rk[sa.back()] + 1 == n)
                break;
        }

        calcHeight(s);
    }

private:
    const int n;
    std::vector<int> oldrk, id, tmp, cnt;

    template <class T>
    inline void calcHeight(const T &s) {
        h.assign(n, 0);
        for (int i = 0, k = 0; i < n; i++) {
            if (rk[i] == 0)
                continue;
            k -= bool(k);
            while (s[i + k] == s[sa[rk[i] - 1] + k])
                k += 1;
            h[rk[i]] = k;
        }
    }

    // 计数排序
    inline void countSort() {
        int m = *std::max_element(begin(rk), end(rk));
        cnt.assign(m + 1, 0);
        for (int i = 0; i < n; i++)
            cnt[tmp[i] = rk[id[i]]] += 1;
        for (int i = 1; i < cnt.size(); i++)
            cnt[i] += cnt[i - 1];
        for (int i = n - 1; i >= 0; i--)
            sa[--cnt[tmp[i]]] = id[i];
    }

    inline bool equal(int x, int y, int w) {
        int rkx = (x + w < n ? oldrk[x + w] : -1);
        int rky = (y + w < n ? oldrk[y + w] : -1);
        return oldrk[x] == oldrk[y] and rkx == rky;
    }
};
```
<a name="b6cbf042"></a>
# 5, DataStructures
<a name="44b3ffcd"></a>
## BaseDataStructures
<a name="5d3b92c1"></a>
### HashTable
```cpp
using u64 = unsigned long long;

template <class T, int Mod>
struct HashTable {
    int hd[Mod], nt[Mod << 1], tot = 0;
    u64 to[Mod << 1];
    T val[Mod << 1];

    void clear() {
        for (int i = 1; i <= tot; i++)
            hd[to[i] % Mod] = 0;
        tot = 0;
    }

    T operator()(u64 x) {
        int u = x % Mod;
        for (int i = hd[u]; i; i = nt[i])
            if (to[i] == x)
                return val[i];
        return T();
    }

    T &operator[](u64 x) {
        int u = x % Mod;
        for (int i = hd[u]; i; i = nt[i])
            if (to[i] == x)
                return val[i];
        to[++tot] = x;
        nt[tot] = hd[u];
        hd[u] = tot;
        return val[tot] = T();
    }
};
HashTable<int, int(1e6) + 3> mp;
```
<a name="39f7fb48"></a>
### HashTableByPointer
```cpp
using u64 = unsigned i64;

template <class T>
struct HashNode {
    HashNode *next = nullptr;
    u64 to = 0;
    T val = T();
};

template <class T>
struct HashTable {
    int tot = 0;

    HashTable(int n) : mod(n) {
        hd = new HashNode<T> *[n]();
    }

    T &operator[](const u64 &x) {
        auto &p = find(x);
        if (!p) {
            tot++;
            p = new HashNode<T>();
            p->to = x;
        }
        return p->val;
    }

    bool erase(const u64 &x) {
        auto &p = find(x);
        if (!p)
            return false;
        p = p->next;
        tot--;
        return true;
    }

    bool count(const u64 &x) const {
        return (bool)find(x);
    }

private:
    HashNode<T> **hd;
    const int mod;
    HashNode<T> *&find(const u64 &x) const {
        int u = x % mod;
        if (!hd[u] or hd[u]->to == x)
            return hd[u];
        auto p = hd[u];
        for (; p->next; p = p->next) {
            if (p->next->to == x)
                break;
        }
        return p->next;
    }
};
using HashT = HashTable<int>;
```
<a name="2d04ab39"></a>
### **SparseTable(RMQ)**
```cpp
template <class T, class Cmp = std::less<T>>
struct RMQ {
    const Cmp cmp = Cmp();
    std::vector<std::vector<T>> ST;

    RMQ(const std::vector<T> &a) {
        int n = a.size(), logn = std::__lg(n);
        ST.assign(n, std::vector<T>(logn + 1));
        for (int i = 0; i < n; i++)
            ST[i][0] = a[i];
        for (int j = 0; j < logn; j++) {
            for (int i = 0; i + (1 << (j + 1)) - 1 < n; i++) {
                ST[i][j + 1] = std::min(ST[i][j], ST[i + (1 << j)][j], cmp);
            }
        }
    }

    // [l, r)
    T operator()(int l, int r) {
        int log = std::__lg(r - l);
        return std::min(ST[l][log], ST[r - (1 << log)][log], cmp);
    }
};
```
<a name="f007d20a"></a>
### DeletableHeap
```cpp
template <class T, class Cmp = std::less<T>>
struct Heap {
    std::priority_queue<T, std::vector<T>, Cmp> qPush, qErase; // Heap=qPush-qErase

    void push(T x) { qPush.push(x); }

    void erase(T x) { qErase.push(x); }

    T top() {
        while (!qErase.empty() && qPush.top() == qErase.top())
            qPush.pop(), qErase.pop();
        return qPush.top();
    }

    void pop() {
        while (!qErase.empty() && qPush.top() == qErase.top())
            qPush.pop(), qErase.pop();
        qPush.pop();
    }

    int size() {
        return qPush.size() - qErase.size();
    }
};
```
<a name="d00ce9af"></a>
### **DisjointSetUnion**
```cpp
struct DSU {
    std::vector<int> f;
    std::vector<int> size;

    DSU(int n) : f(n), size(n) {
        std::iota(f.begin(), f.end(), 0);
        std::fill(size.begin(), size.end(), 1);
    }

    int find(int x) {
        while (x != f[x])
            x = f[x] = f[f[x]];
        return x;
    }

    void Union(int x, int y) {
        x = find(x), y = find(y);
        if (x == y)
            return;

        if (size[x] < size[y])
            std::swap(x, y);

        size[x] += size[y];
        f[y] = x;
    }
};
```
<a name="3b2c3169"></a>
### Weighted**DisjointSetUnion**
```cpp
template <class T>
struct DSU {
    std::vector<int> f;
    std::vector<int> size;
    std::vector<T> w;

    DSU(int n) : f(n), size(n), w(n) {
        std::iota(f.begin(), f.end(), 0);
        std::fill(size.begin(), size.end(), 1);
    }

    int find(int x) {
        if (f[x] == x)
            return x;
        int pr = f[x], anc = find(pr);

        w[x] = w[x] + w[pr];

        return f[x] = anc;
    }

    void Union(int x, int y, const T &z) {
        T road = w[x] + z, lastWy = w[y];
        x = find(x), y = find(y);
        if (x == y)
            return;

        w[y] = road - lastWy;

        size[x] += size[y];
        f[y] = x;
    }
};

struct Info {
    int val;

    Info(int x = 0) : val(x) {}

    bool operator==(const Info &a) const {
        return val == a.val;
    }

    Info operator+(const Info &a) const {
    }

    Info operator-(const Info &a) const {
    }
};
```
<a name="5be3e824"></a>
## Segment Tree
<a name="hJ67o"></a>
### RangeQuerySegmentTree
```cpp
template <class T, class Merge = std::plus<T>>
struct SegT {
    const Merge merge;
    const int n;
    std::vector<T> t;

    SegT(int n) : n(n), t(4 << std::__lg(n)), merge(Merge()) {}

    SegT(const std::vector<T> &a) : SegT(a.size()) {
        std::function<void(int, int, int)> build = [&](int i, int l, int r) {
            if (r - l == 1) {
                t[i] = a[l];
                return;
            }
            int mid = l + r >> 1;
            build(i << 1, l, mid);
            build(i << 1 | 1, mid, r);
            up(i);
        };
        build(1, 0, n);
    }

    void up(int i) {
        t[i] = merge(t[i << 1], t[i << 1 | 1]);
    }

    // 默认单点赋值
    void modify(int x, const T &v) {
        modify(1, 0, n, x, v);
    }

    void modify(int i, int l, int r, int x, const T &v) {
        if (r - l == 1) {
            t[i] = v;
            return;
        }
        int mid = l + r >> 1;
        if (x < mid)
            modify(i << 1, l, mid, x, v);
        else
            modify(i << 1 | 1, mid, r, x, v);
        up(i);
    }

    // [l, r)
    T rangeQuery(int l, int r) {
        return rangeQuery(1, 0, n, l, r);
    }

    T rangeQuery(int i, int l, int r, int tl, int tr) {
        if (tl <= l and r <= tr) {
            return t[i];
        }
        int mid = l + r >> 1;
        return merge((tl < mid ? rangeQuery(i << 1, l, mid, tl, tr) : T()),
                     (mid < tr ? rangeQuery(i << 1 | 1, mid, r, tl, tr) : T()));
    }
};
```
<a name="ec15cfc4"></a>
### LazySegmentTree
```cpp
template <class T, class Tag>
struct LazySegT {
    int n;
    std::vector<T> info;
    std::vector<Tag> tag;

    LazySegT(int n, T v = T())
        : LazySegT(std::vector(n, v)) {}

    template <class G>
    LazySegT(const std::vector<G> &a) : n(a.size()) {
        info.assign(4 << std::__lg(n), T());
        tag.assign(4 << std::__lg(n), Tag());
        std::function<void(int, int, int)> build = [&](int i, int l, int r) {
            if (r - l == 1) {
                info[i] = a[l];
                return;
            }
            int mid = l + r >> 1;
            build(i << 1, l, mid);
            build(i << 1 | 1, mid, r);
            up(i);
        };
        build(1, 0, n);
    }

    void up(int i) {
        info[i] = info[i << 1] + info[i << 1 | 1];
    }
    void apply(int i, const Tag &v) {
        info[i].apply(v);
        tag[i].apply(v);
    }
    void down(int i) {
        apply(i << 1, tag[i]);
        apply(i << 1 | 1, tag[i]);
        tag[i] = Tag();
    }

    // 单点修改
    void modify(int i, const T &v) {
        modify(1, 0, n, i, v);
    }
    void modify(int i, int l, int r, int x, const T &v) {
        if (r - l == 1) {
            info[i] = v;
            return;
        }
        int mid = l + r >> 1;
        down(i);
        if (x < mid) {
            modify(i << 1, l, mid, x, v);
        } else {
            modify(i << 1 | 1, mid, r, x, v);
        }
        up(i);
    }

    // 区间查询 [l, r)
    T rangeQuery(int l, int r) {
        return rangeQuery(1, 0, n, l, r);
    }
    T rangeQuery(int i, int l, int r, int tl, int tr) {

        if (tl <= l and r <= tr)
            return info[i];

        down(i);
        int mid = l + r >> 1;

        return (tl < mid ? rangeQuery(i << 1, l, mid, tl, tr) : T()) +
               (mid < tr ? rangeQuery(i << 1 | 1, mid, r, tl, tr) : T());
    }

    // 区间修改 [l, r)
    void rangeModify(int l, int r, const Tag &v) {
        return rangeModify(1, 0, n, l, r, v);
    }
    void rangeModify(int i, int l, int r, int tl, int tr, const Tag &v) {

        if (tl <= l and r <= tr) {
            apply(i, v);
            return;
        }
        down(i);
        int mid = l + r >> 1;

        if (tl < mid)
            rangeModify(i << 1, l, mid, tl, tr, v);
        if (mid < tr)
            rangeModify(i << 1 | 1, mid, r, tl, tr, v);
        up(i);
    }

    // 区间左边第一个满足条件的下标
    template <class F>
    int findFirst(int l, int r, F pred) {
        return findFirst(1, 0, n, l, r, pred);
    }
    template <class F>
    int findFirst(int i, int l, int r, int tl, int tr, F pred) {
        if (l >= tr || r <= tl || !pred(info[i])) {
            return -1;
        }
        if (r - l == 1) {
            return l;
        }
        int mid = l + r >> 1;
        down(i);
        int res = findFirst(i << 1, l, mid, tl, tr, pred);
        if (res == -1) {
            res = findFirst(i << 1 | 1, mid, r, tl, tr, pred);
        }
        return res;
    }

    // 区间右边第一个满足条件的下标
    template <class F>
    int findLast(int l, int r, F pred) {
        return findLast(1, 0, n, l, r, pred);
    }
    template <class F>
    int findLast(int i, int l, int r, int tl, int tr, F pred) {
        if (l >= tr || r <= tl || !pred(info[i])) {
            return -1;
        }
        if (r - l == 1) {
            return l;
        }
        int mid = l + r >> 1;
        down(i);
        int res = findLast(i << 1 | 1, mid, r, tl, tr, pred);
        if (res == -1) {
            res = findLast(i << 1, l, mid, tl, tr, pred);
        }
        return res;
    }
};

struct Tag {
    int add;

    Tag(const int &add = 0)
        : add(add){}

    void apply(const Tag &tag) {
        if (tag.add) {
        }
    }
};

struct Node {
    i64 val;
    int len;

    Node(const i64 &val = 0, const int &len = 1)
        : val(val), len(len) {}

    void apply(const Tag &tag) {
        if (tag.add) {
            val += 1ll * tag.add * len;
        }
    }

    Node operator+(const Node &a) {
        return Node(val + a.val, len + a.len);
    }
};
```
<a name="27a97257"></a>
### SegmentTreeMerge
```cpp
Node *merge(Node *p, Node *q, int l, int r) {
    if (!p or !q) {
        return p ? p : q;
    }

    // 默认开新点
    Node *s = new Node();

    if (r - l == 1) {
        // 默认为加
        s->val = p->val + q->val;
        return s;
    }

    int mid = l + r >> 1;
    s->l = merge(p->l, q->l, l, mid);
    s->r = merge(p->r, q->r, mid, r);

    return up(s), s;
}
```
<a name="440f0c4b"></a>
### PresidentTree
```cpp
template <class T>
class PresidentTree {
    struct Node {
        Node *l = nullptr;
        Node *r = nullptr;
        int val = 0;
        void extend() {
            if (!l)
                l = new Node();
            if (!r)
                r = new Node();
        }
    };

    const T Start, Last;
    std::vector<Node *> root;

    void up(Node *p) {
        if (!p->l or !p->r) {
            p->val = (p->l ? p->l->val : p->r->val);
        } else
            p->val = p->l->val + p->r->val;
    }

    void modify(Node *&p, T l, T r, T x) {
        if (!p)
            p = new Node();
        if (r - l == 1) {
            p->val++;
            return;
        }
        auto mid = l + r >> 1;
        if (x < mid)
            modify(p->l, l, mid, x);
        else
            modify(p->r, mid, r, x);
        up(p);
    }
    Node *merge(Node *x, Node *y, T l, T r) {
        if (!x or !y)
            return (x ? x : y);
        Node *p = new Node();
        if (r - l == 1) {
            p->val = x->val + y->val;
            return p;
        }
        auto mid = l + r >> 1;
        p->l = merge(x->l, y->l, l, mid);
        p->r = merge(x->r, y->r, mid, r);
        return up(p), p;
    }

    int getRange(Node *&x, Node *&y, T l, T r, T tl, T tr) {
        if (!x)
            x = new Node();
        if (!y)
            y = new Node();

        if (tl <= l and r <= tr) {
            return y->val - x->val;
        }
        auto mid = l + r >> 1;
        return (tl < mid ? getRange(x->l, y->l, l, mid, tl, tr) : 0) +
               (mid < tr ? getRange(x->r, y->r, mid, r, tl, tr) : 0);
    }

    T getKth(Node *x, Node *y, T l, T r, int k) {
        if (r - l == 1)
            return l;
        T mid = l + r >> 1;
        x->extend(), y->extend();
        int L = y->l->val - x->l->val;

        return (L >= k ? getKth(x->l, y->l, l, mid, k)
                       : getKth(x->r, y->r, mid, r, k - L));
    }

public:
    PresidentTree(const std::vector<T> &a, T min, T max)
        : root(a.size() + 1), Start(min), Last(max + 1) {

        root[0] = new Node();
        for (int i = 1; i <= a.size(); i++) {
            modify(root[i], Start, Last, a[i - 1]);
            root[i] = merge(root[i], root[i - 1], Start, Last);
        }
    }
    // [l, r), [tl, tr)
    int getRange(int l, int r, T tl, T tr) {
        return getRange(root[l], root[r], Start, Last, tl, tr);
    }
    // [l, r)
    T getKth(int l, int r, int k) {
        return getKth(root[l], root[r], Start, Last, k);
    }
};
```
<a name="5e8d79b1"></a>
### LiChaoTree
```cpp
template <class T, T LB, T RB>
struct LiChaoTree {
    using ld = double;

    struct Line {
        ld k, b;          // 线段斜率, 截距
        T l = LB, r = RB; // 线段左右边界

        Line() {}
        Line(ld k, ld b, T l = LB, T r = RB) : k(k), b(b), l(l), r(r) {}

        ld at(T x) { return k * x + b; }

        friend T cross(const Line &a, const Line &b) {
            return floor((a.b - b.b) / (b.k - a.k));
        }
    };

    struct Node {
        Line seg;
        Node *l = nullptr;
        Node *r = nullptr;
    };

    LiChaoTree() { rt = new Node(); } // 一定要加括号初始化!!

    void insert(ld k, ld b, T lm = LB, T rm = RB) {
        Line Z = (Line){k, b, lm, rm};
        insert(rt, LB, RB, Z);
    }
    ld get(T x) {
        return get(rt, LB, RB, x);
    }

private:
    const ld eps = 1e-9;
    int sgn(const ld &a) { return (a < -eps ? -1 : a > eps); }

    Node *rt;
    void insert(Node *&p, T l, T r, Line Z) {
        if (!p) // 一定要传node*&
            p = new Node();

        T mid = l + r >> 1;
        if (Z.l <= l && r <= Z.r) {
            T L = sgn(Z.at(l) - p->seg.at(l));
            T R = sgn(Z.at(r) - p->seg.at(r));
            if (L + R == 2) { // 新线段完全在老线段之上
                p->seg = Z;
            } else if (L == 1 or R == 1) { // 新线段不完全在老线段之上
                if (sgn(Z.at(mid) - p->seg.at(mid)) > 0)
                    std::swap(p->seg, Z);

                if (sgn(cross(Z, p->seg) - mid) < 0)
                    insert(p->l, l, mid, Z);
                else
                    insert(p->r, mid + 1, r, Z);
            }
        } else {
            if (Z.l <= mid)
                insert(p->l, l, mid, Z);
            if (mid < Z.r)
                insert(p->r, mid + 1, r, Z);
        }
    }
    ld get(Node *&p, T l, T r, T x) {
        if (!p)
            return 0;
        T mid = l + r >> 1;
        ld ans = p->seg.at(x);
        if (l < r)
            ans = std::max(ans, (x <= mid ? get(p->l, l, mid, x) : get(p->r, mid + 1, r, x)));
        return ans;
    }
};
```
<a name="fXGNd"></a>
## Binary Indexed Tree
<a name="SqFm4"></a>
### BIT
```cpp
template <class T, class Cmp = std::greater<T>>
struct Max {
    const Cmp cmp = Cmp();
    constexpr T operator()(const T &a, const T &b) const {
        return std::min(a, b, cmp);
    }
};

template <class T, class Merge = std::plus<T>>
struct BIT {
    const Merge merge;
    std::vector<T> t;

    BIT(int n) : t(n + 1), merge(Merge()) {}

    // O(n) build BIT
    BIT(const std::vector<T> &a) : BIT(a.size()) {
        int n = a.size();
        for (int i = 1; i <= n; i++) {
            t[i] = merge(t[i], a[i - 1]);
            int j = i + (i & -i);
            if (j <= n)
                t[j] = merge(t[j], t[i]);
        }
    }

    void modify(int i, const T &x) {
        for (i += 1; i < t.size(); i += i & -i)
            t[i] = merge(t[i], x);
    }

    T posQuery(int i)  {
        T res = T();
        for (i += 1; i; i -= i & -i)
            res = merge(res, t[i]);
        return res;
    }

    // [l, r)
    T rangeQuery(int l, int r)  {
        return posQuery(r - 1) - posQuery(l - 1);
    }
};
```
<a name="uyMle"></a>
### RangeBIT
```cpp
template <class T>
struct RangeBIT {
    BIT<T, std::plus<T>> d, s;

    RangeBIT(int n) : d(n), s(n) {}

    // O(n) build RangeBIT
    RangeBIT(std::vector<T> a)
        : d(diff(a)), s(multIndex(diff(a))) {}

    static std::vector<T> diff(std::vector<T> a) {
        std::adjacent_difference(begin(a), end(a), begin(a));
        return a;
    }

    static std::vector<T> multIndex(std::vector<T> a) {
        for (int i = 0; i < a.size(); i++) {
            a[i] *= i;
        }
        return a;
    }

    // [l, r)
    void rangeModify(int l, int r, const T &x) {
        d.modify(l, x), d.modify(r, -x);
        s.modify(l, l * x), s.modify(r, -r * x);
    }

    // [l, r)
    T rangeQuery(int l, int r)  {
        T res1 = r * d.posQuery(r - 1) - s.posQuery(r - 1);
        T res2 = l * d.posQuery(l - 1) - s.posQuery(l - 1);
        return res1 - res2;
    }
};
```
<a name="fMpZB"></a>
### MatBIT
```cpp
template <class T, class Merge = std::plus<T>>
struct MatBIT {
    const Merge merge;
    std::vector<BIT<T, Merge>> t;

    MatBIT(int n, int m)
        : t(n + 1, BIT<T>(m)), merge(Merge()) {}

    void modify(int x, int y, const T &v) {
        for (int i = x + 1; i < t.size(); i += i & -i) {
            t[i].modify(y, v);
        }
    }

    T posQuery(int x, int y) {
        T res = T();
        for (int i = x + 1; i; i -= i & -i) {
            res = merge(res, t[i].posQuery(y));
        }
        return res;
    }

    // [u, d), [l, r)
    T rangeQuery(int u, int l, int d, int r) {
        u -= 1, l -= 1, d -= 1, r -= 1;
        T res1 = posQuery(d, r) + posQuery(u, l);
        T res2 = posQuery(d, l) + posQuery(u, r);
        return res1 - res2;
    }
};
```
<a name="JkMrn"></a>
### RangeMatBIT
```cpp
template <class T>
struct RangeMatBIT {
    MatBIT<T> p, px, py, pxy;

    RangeMatBIT(int n, int m)
        : p(n, m), px(n, m), py(n, m), pxy(n, m) {}

    // [u, d), [l, r)
    void rangeModify(int u, int l, int d, int r, const T &v) {
        modify(u, l, v);
        modify(d, r, v);
        modify(u, r, -v);
        modify(d, l, -v);
    }

    // [u, d), [l, r)
    T rangeQuery(int u, int l, int d, int r) {
        u -= 1, l -= 1, d -= 1, r -= 1;
        return query(u, l) + query(d, r) - query(d, l) - query(u, r);
    }

private:
    void modify(int x, int y, const T &v) {
        p.modify(x, y, v);
        px.modify(x, y, v * x);
        py.modify(x, y, v * y);
        pxy.modify(x, y, v * x * y);
    }

    T query(int x, int y) {
        T res = T();
        res += p.posQuery(x, y) * (x + 1) * (y + 1);
        res -= px.posQuery(x, y) * (y + 1);
        res -= py.posQuery(x, y) * (x + 1);
        res += pxy.posQuery(x, y);
        return res;
    }
};
```
<a name="LXNVb"></a>
## Block
```cpp
// O(sqrt(n)) 区间加, O(1) 单点查
template <class T, class Merge = std::plus<T>>
struct Block {
    const int n, B;
    const Merge merge;
    std::vector<T> a, b;

    Block(int n, const T &v = T()) : Block(std::vector<T>(n, v)) {}

    Block(const std::vector<T> &_init)
        : n(_init.size()), B(sqrt(2 * _init.size())), a(_init), merge(Merge()) {
        b.assign(n / B + 1, T());
    }

    // [l, r)
    void add(int l, int r, const T &v) {
        for (; l / B == (l - 1) / B and l < r; l++) {
            a[l] = merge(a[l], v);
        }
        for (; r / B == (r - 1) / B and l < r; r--) {
            a[r - 1] = merge(a[r - 1], v);
        }
        for (int i = l / B; i < r / B; i++) {
            b[i] = merge(b[i], v);
        }
    }

    T get(int x) {
        return merge(a[x], b[x / B]);
    }
};
```
<a name="LJKxt"></a>
# 6, Tree Theory
<a name="ZVO88"></a>
## TreePre
```cpp
template <class T>
class TreePre {
    static constexpr int endPoint(int x) {
        return x;
    }
    template <class G>
    static constexpr int endPoint(const std::pair<int, G> &pr) {
        return pr.first;
    }

    constexpr void dfs1(int x, int f) {
        size[x] = 1;
        fa[x] = f;

        for (auto &&p : e[x]) {
            int &&y = endPoint(p);
            if (y != f) {
                dep[y] = dep[x] + 1;
                dfs1(y, x);

                size[x] += size[y];
                if (big[x] == -1 or size[y] > size[big[x]])
                    big[x] = y;
            }
        }
    }
    constexpr void dfs2(int x, int top) {
        dfn[x] = cur++;
        idfn[dfn[x]] = x;
        tp[x] = top;
        if (big[x] != -1)
            dfs2(big[x], top);

        for (auto &&p : e[x]) {
            int &&y = endPoint(p);
            if (y != big[x] and y != fa[x])
                dfs2(y, y);
        }
    }
    const std::vector<std::vector<T>> &e;

public:
    std::vector<int> size, dep, big, tp, fa, dfn, idfn;
    // dfn begin from 0
    int cur = 0;

    constexpr TreePre(const std::vector<std::vector<T>> &g, int root)
        : e(g), tp(g.size()), big(g.size(), -1),
          size(g.size()), dep(g.size()), fa(g.size()),
          dfn(g.size()), idfn(g.size()) {
        // dep begin from 0
        dep[root] = 0;
        dfs1(root, -1);
        dfs2(root, root);
    }

    constexpr int getLca(int x, int y) const {
        while (tp[x] != tp[y])
            (dep[tp[x]] > dep[tp[y]] ? x = fa[tp[x]] : y = fa[tp[y]]);
        return (dep[x] < dep[y] ? x : y);
    }

    constexpr auto dist(int x, int y) const {
        int &&lca = getLca(x, y);
        return dep[x] + dep[y] - 2 * dep[lca];
    }

    // x→y路径剖分的dfn号区间[l, r]
    constexpr auto getRoad(int x, int y) const {
        int lca = getLca(x, y);
        std::vector<std::pair<int, int>> vec1, vec2;
        while (tp[x] != tp[lca]) {
            vec1.push_back({dfn[x], dfn[tp[x]]});
            x = fa[tp[x]];
        }

        if (x != lca) {
            vec1.push_back({dfn[x], dfn[lca] + 1});
        }

        vec1.push_back({dfn[lca], dfn[lca]});

        while (tp[y] != tp[lca]) {
            vec2.push_back({dfn[tp[y]], dfn[y]});
            y = fa[tp[y]];
        }

        if (y != lca) {
            vec2.push_back({dfn[lca] + 1, dfn[y]});
            y = fa[tp[y]];
        }

        vec1.insert(end(vec1), rbegin(vec2), rend(vec2));
        return vec1;
    }

    constexpr int kthAncester(int x, int k) const {
        if (dep[x] < k) {
            return -1;
        }

        int &&d = dep[x] - k;

        while (dep[tp[x]] > d) {
            x = fa[tp[x]];
        }

        return idfn[dfn[x] - dep[x] + d];
    }

    // x is y's ancester
    constexpr bool isAncester(int x, int y) const {
        return dfn[x] <= dfn[y] and dfn[y] < dfn[x] + size[x];
    }
};
```
<a name="RDr69"></a>
## VirtualTreePre
```cpp
// G:原树边, T:虚树边
template <class G, class T>
class VirtualTreePre {
    const std::function<bool(int, int)> dfnCmp = [&](int x, int y) {
        return pre.dfn[x] < pre.dfn[y];
    };

public:
    const TreePre<G> pre;
    std::vector<std::vector<T>> vt;

    VirtualTreePre(const std::vector<std::vector<G>> &e, int root)
        : pre(e, root), vt(e.size()) {}

    // 虚树存在vt, 返回vt根节点
    int build(std::vector<int> a) {
        std::sort(begin(a), end(a), dfnCmp);
        for (int i = (int)a.size() - 1; i; i--) {
            a.push_back(pre.getLca(a[i - 1], a[i]));
        }
        std::sort(begin(a), end(a), dfnCmp);
        a.erase(std::unique(begin(a), end(a)), end(a));

        for (int x : a)
            vt[x].clear();

        for (int i = 1; i < a.size(); i++) {
            int lca = pre.getLca(a[i - 1], a[i]);
            vt[lca].emplace_back(a[i], pre.dep[a[i]] - pre.dep[lca]);
        }

        return a.front();
    }
};
```
<a name="O2C4k"></a>
## DsuOnTree
```cpp
//	先递归轻子树
    
//	if(当前节点有重子树:big!=-1){
//		递归重子树 
//	}
    
//	加上所有轻子树贡献 
    
//	加上当前节点贡献 
    
//	统计当前节点答案 
    
//	if(当前节点是自己父节点的重儿子)
//    	return   
//	else
//    	删除当前这棵树的贡献
```
<a name="le5Oz"></a>
## CentroidDecomposition
```cpp
template <class T>
class CentroidDecomposition {
    const std::vector<std::vector<T>> &e;
    std::vector<bool> vis;
    std::vector<int> size;

    static constexpr int endPoint(int x) {
        return x;
    }
    template <class G>
    static constexpr int endPoint(const std::pair<int, G> &pr) {
        return pr.first;
    }

    void dfsSize(int x, int fa) {
        size[x] = 1;
        for (auto &&p : e[x]) {
            int &&y = endPoint(p);
            if (y != fa and !vis[y]) {
                dfsSize(y, x);
                size[x] += size[y];
            }
        }
    }

    int getRoot(int x, int fa, int m) {
        for (auto &&p : e[x]) {
            int &&y = endPoint(p);
            if (y != fa and !vis[y] and 2 * size[y] > m) {
                return getRoot(y, x, m);
            }
        }
        return x;
    }

    void build(int x) {
        vis[x] = true;
        dfsOrder.push_back(x);

        for (auto &&p : e[x]) {
            int y = endPoint(p);
            if (!vis[y]) {
                dfsSize(y, -1);
                y = getRoot(y, -1, size[y]);

                cdt[x].push_back(y);
                build(y);
            }
        }
    }

public:
    std::vector<std::vector<int>> cdt;
    std::vector<int> dfsOrder;
    int root;
    CentroidDecomposition(std::vector<std::vector<T>> &g)
        : e(g), vis(g.size()), size(g.size()), cdt(g.size()) {
        dfsSize(0, -1);
        root = getRoot(0, -1, size[0]);
        build(root);
    }
};
```
<a name="fneZ8"></a>
## TreeDiameter
```cpp
template <class T>
class TreeDiameter {
    static constexpr std::pair<int, int> edge(int x) {
        return {x, 1};
    }
    template <class G>
    static constexpr int edge(const std::pair<int, G> &pr) {
        return pr;
    }

    const std::vector<T> &e;
    std::vector<i64> dis;

    void dfs(int x, int fa) {
        for (auto p : e[x]) {
            auto [y, w] = edge(p);
            if (y != fa) {
                dis[y] = dis[x] + w;
                dfs(y, x);
            }
        }
    }

public:
    int v1 = 0, v2 = 0;
    i64 diameter = 0;

    TreeDiameter(const std::vector<T> &e) : e(e), dis(e.size()) {
        dfs(0, -1);
        v1 = std::max_element(begin(dis), end(dis)) - begin(dis);
        dis[v1] = 0;
        dfs(v1, -1);
        v2 = std::max_element(begin(dis), end(dis)) - begin(dis);
        diameter = dis[v2];
    }
};
```
<a name="NPt8W"></a>
## DecareTree
```cpp
// less 是小根树
template <class T, class Cmp = std::less<T>>
struct DecareTree {
    const Cmp cmp = Cmp();
    std::vector<std::pair<int, int>> t;
    int root;

    // 0-indexed
    DecareTree(const std::vector<T> &a) : t(a.size(), {-1, -1}) {
        int n = a.size(), r = 0;
        std::vector<int> q(n + 1);
        for (int i = 0; i < n; i++) {
            while (r and !cmp(a[q[r]], a[i])) {
                t[i].first = q[r];
                r -= 1;
            }
            if (r > 0) {
                t[q[r]].second = i;
            }
            q[++r] = i;
        }

        root = q[1];
    }
};
```
<a name="OGLT4"></a>
## TreeHash
```cpp
using u64 = unsigned long long;

template <class T>
struct Rand {
    std::mt19937 myrand;
    Rand(const i64 seed = time(0)) : myrand(seed) {}
    T operator()(T l, T r) {
        return std::uniform_int_distribution<T>(l, r)(myrand);
    }
};
Rand<u64> rd;

u64 f(u64 x) {
    const static u64 r1 = rd(1 << 20, 1 << 24);
    const static u64 r2 = rd(1 << 25, 1 << 30);
    const static u64 mask = (1ll << 31) - 1;

    auto h = [&](u64 y) {
        return (u64)y * y * y * r1 + r2;
    };
    return h(x & mask) + h(x >> 31);
}

```
<a name="ez7t0"></a>
# 7, Graph Theory
<a name="ur07d"></a>
## Dijkstra
```cpp
template <class T, class G>
class Dijkstra {
    const std::vector<std::vector<std::pair<int, T>>> &e;
    std::vector<std::vector<G>> dis;

    auto get(int s) {
        std::vector<G> dis(e.size(), std::numeric_limits<G>::max() / 2);

        using pii = std::pair<G, int>;
        std::priority_queue<pii, std::vector<pii>, std::greater<>> q;

        dis[s] = G();
        q.push({dis[s], s});

        while (!q.empty()) {
            auto [D, x] = q.top();
            q.pop();

            if (D > dis[x])
                continue;

            for (auto &&[y, w] : e[x]) {
                if (dis[y] > dis[x] + w) {
                    dis[y] = dis[x] + w;
                    q.push({dis[y], y});
                }
            }
        }
        return dis;
    }

public:
    Dijkstra(const std::vector<std::vector<std::pair<int, T>>> &g)
        : e(g), dis(g.size()) {}

    G operator()(int x, int y) {
        if (dis[x].empty())
            dis[x] = get(x);
        return dis[x][y];
    }
};
```
<a name="n2ZVF"></a>
## Bi-GraphColor
```cpp
struct BiGraphColor {
    std::vector<int> col;
    bool isBiGraph;

    BiGraphColor(const std::vector<std::vector<int>> &e)
        : col(e.size(), -1), isBiGraph(true) {

        int n = e.size();
        std::function<void(int)> dfs = [&](int x) {
            if (!isBiGraph)
                return;
            for (int y : e[x]) {
                if (col[y] == -1) {
                    col[y] = col[x] ^ 1;
                    dfs(y);
                } else if (col[y] == col[x]) {
                    isBiGraph = false;
                    return;
                }
            }
        };

        for (int i = 0; i < n; i++) {
            if (col[i] == -1) {
                col[i] = 0;
                dfs(i);
            }
            if (!isBiGraph)
                return;
        }
    }
};
```
<a name="NzDoU"></a>
## RingTree
```cpp
struct RingTree {
    std::vector<std::vector<int>> g;
    std::vector<int> ring;
    DSU dsu;

    RingTree(const std::vector<std::vector<int>> &e)
        : g(e.size()), dsu(e.size()) {
        int n = e.size();

        std::vector<int> vis(n, 0), deg(n, 0);
        for (int x = 0; x < n; x++) {
            for (int y : e[x])
                deg[y]++;
        }

        std::function<void(int)> dfs1 = [&](int x) {
            vis[x] = 1;
            for (int y : e[x]) {
                if (vis[y])
                    continue;

                dsu.Union(y, x);
                g[y].push_back(x);

                deg[y]--;
                if (deg[y] == 1)
                    dfs1(y);
            }
        };
        for (int i = 0; i < n; i++)
            if (!vis[i] and deg[i] == 1)
                dfs1(i);

        std::function<void(int)> dfs2 = [&](int x) {
            ring.push_back(x);
            vis[x] = 1;
            for (int y : e[x])
                if (!vis[y])
                    dfs2(y);
        };

        for (int i = 0; i < n; i++)
            if (vis[i] == 0)
                dfs2(i);
    }
};
```
<a name="JI2UC"></a>
## TopSort
```cpp
class TopSort {
    static constexpr int endPoint(int x) {
        return x;
    }
    template <class G>
    static constexpr int endPoint(const std::pair<int, G> &pr) {
        return pr.first;
    }

public:
    template <class T>
    std::vector<int> operator()(const std::vector<T> &e) const {
        int n = e.size();
        std::vector<int> ind(n);
        for (int x = 0; x < n; x++) {
            for (auto p : e[x]) {
                ind[endPoint(p)] += 1;
            }
        }

        std::vector<int> q;
        for (int x = 0; x < n; x++) {
            if (ind[x] == 0) {
                q.push_back(x);
            }
        }

        std::vector<int> res;
        while (!q.empty()) {
            int x = q.back();
            res.push_back(x);
            q.pop_back();

            for (auto p : e[x]) {
                int y = endPoint(p);
                ind[y] -= 1;
                if (ind[y] == 0) {
                    q.push_back(y);
                }
            }
        }

        return res;
    }
};
const TopSort topSort;
```
<a name="LySs0"></a>
## Connectivity
<a name="sD3EM"></a>
### StronglyConnectedComponent
```cpp
class SCC {
    const std::vector<std::vector<int>> &e;
    std::vector<int> q; // stack
    int r = 0, cur = 0;

    void dfs(int x) {
        dfn[x] = low[x] = cur++;
        q[++r] = x;

        for (int y : e[x]) {
            if (dfn[y] == -1) {
                dfs(y);
                low[x] = std::min(low[x], low[y]);
            } else if (bel[y] == -1) {
                low[x] = std::min(low[x], dfn[y]);
            }
        }

        if (dfn[x] == low[x]) {
            int y;
            do {
                y = q[r--];
                bel[y] = cntBlock;
            } while (y != x);
            cntBlock += 1;
        }
    }

public:
    // original graph
    std::vector<int> dfn, low, bel;

    // shrinking graph
    std::vector<std::vector<int>> g;
    int cntBlock = 0;

    SCC(const std::vector<std::vector<int>> &e)
        : e(e), dfn(e.size(), -1), low(e.size()), bel(e.size(), -1) {
        int n = e.size();
        q.assign(n + 1, 0);

        for (int i = 0; i < n; i++) {
            if (dfn[i] == -1) {
                dfs(i);
            }
        }

        g.resize(cntBlock);
        for (int x = 0; x < n; x++) {
            for (int y : e[x]) {
                if (bel[x] == bel[y])
                    continue;
                g[bel[x]].push_back(bel[y]);
            }
        }

        // for (int x = 0; x < cntBlock; x++) {
        //     std::sort(begin(g[x]), end(g[x]));
        //     g[x].erase(std::unique(begin(g[x]), end(g[x])), end(g[x]));
        // }
    }
};

```
<a name="gRywk"></a>
### TwoSat
```cpp
class TwoSat {
    const int n;
    std::vector<std::vector<int>> e;

public:
    std::vector<bool> ans;

    TwoSat(int n) : n(n), e(2 * n) {}

    void add(int x, bool f, int y, bool g) {
        x *= 2, y *= 2;
        e[x + !f].push_back(y + g);
        e[y + !g].push_back(x + f);
    }

    bool work() {
        SCC scc(e);
        const auto &bel = scc.bel;
        ans.assign(n, false);
        for (int i = 0; i < n; i++) {
            if (bel[2 * i] == bel[2 * i + 1])
                return false;
            ans[i] = bel[2 * i] > bel[2 * i + 1];
        }
        return true;
    }
};
```
<a name="bKKdd"></a>
### VertexBiconnectedComponent
```cpp
class VertexBC {
    const std::vector<std::vector<int>> &e;
    int cur = 0;

    void dfs(int x, int root) {
        dfn[x] = low[x] = cur++;

        int sonNum = 0;
        for (int y : e[x]) {
            if (dfn[y] == -1) {
                sonNum += 1;
                dfs(y, root);
                low[x] = std::min(low[x], low[y]);

                if (low[y] >= dfn[x] and x != root) {
                    cutDeg[x] += 1;
                }
            } else {
                low[x] = std::min(low[x], dfn[y]);
            }
        }

        if (x == root) {
            cutDeg[x] = sonNum - 1;
        }
    }

public:
    // original graph
    std::vector<int> dfn, low, cutDeg;
    int componentNum = 0;

    VertexBC(const std::vector<std::vector<int>> &e)
        : e(e), dfn(e.size(), -1), low(e.size()), cutDeg(e.size()) {
        int n = e.size();

        for (int i = 0; i < n; i++) {
            if (dfn[i] == -1) {
                componentNum += 1;
                dfs(i, i);
            }
        }
    }
};
```
<a name="AXGB5"></a>
### EdgeBiconnectedComponent
```cpp
class EdgeBC {
    const std::vector<std::vector<int>> &e;
    std::vector<int> q; // stack
    int r = 0, cur = 0;

    void dfs(int x, int fa) {
        dfn[x] = low[x] = cur++;
        q[++r] = x;

        for (int y : e[x]) {
            if (y == fa) {
                fa = ~fa;
                continue;
            }
            if (dfn[y] == -1) {
                dfs(y, x);
                low[x] = std::min(low[x], low[y]);
            } else {
                low[x] = std::min(low[x], dfn[y]);
            }
        }

        if (dfn[x] == low[x]) {
            int y;
            do {
                y = q[r--];
                bel[y] = cntBlock;
            } while (y != x);
            cntBlock += 1;
        }
    }

public:
    // original graph
    std::vector<int> dfn, low, bel, cutDeg;

    // shrinking graph
    std::vector<std::vector<int>> g;
    int cntBlock = 0, componentNum = 0;

    EdgeBC(const std::vector<std::vector<int>> &e)
        : e(e), dfn(e.size(), -1), low(e.size()), bel(e.size(), -1), cutDeg(e.size()) {
        int n = e.size();
        q.assign(n + 1, 0);

        for (int i = 0; i < n; i++) {
            if (dfn[i] == -1) {
                componentNum += 1;
                dfs(i, -1);
            }
        }

        g.resize(cntBlock);
        for (int x = 0; x < n; x++) {
            for (int y : e[x]) {
                if (bel[x] == bel[y])
                    continue;
                g[bel[x]].push_back(bel[y]);
            }
        }
    }
};
```
<a name="vc9GU"></a>
## Flow
<a name="PZfa1"></a>
### MaxFlow
```cpp
template <class T>
struct Flow {
    const int n;

    std::vector<std::pair<int, T>> e;
    std::vector<std::vector<int>> g;
    std::vector<int> cur, dep;

    Flow(int n) : n(n), g(n) {}

    bool bfs(int s, int t) {
        dep.assign(n, -1);
        std::queue<int> q;
        dep[s] = 0;

        q.push(s);
        while (!q.empty()) {
            const int u = q.front();
            q.pop();

            for (int i : g[u]) {
                auto [v, c] = e[i];

                if (c > 0 and dep[v] == -1) {
                    dep[v] = dep[u] + 1;
                    if (v == t)
                        return true;
                    q.push(v);
                }
            }
        }

        return false;
    }

    T dfs(int u, int t, T f) {
        if (u == t) {
            return f;
        }
        T res = f;
        for (int &i = cur[u]; i < g[u].size(); i++) {
            const int j = g[u][i];
            auto [v, c] = e[j];

            if (c > 0 and dep[v] == dep[u] + 1) {
                T out = dfs(v, t, std::min(res, c));
                e[j].second -= out;
                e[j ^ 1].second += out;

                res -= out;
                if (res == 0) {
                    return f;
                }
            }
        }
        return f - res;
    }

    void add(int u, int v, T c) {
        g[u].push_back(e.size());
        e.emplace_back(v, c);
        g[v].push_back(e.size());
        e.emplace_back(u, 0);
    }

    T work(int s, int t) {
        T ans = 0;
        while (bfs(s, t)) {
            cur.assign(n, 0);
            ans += dfs(s, t, std::numeric_limits<T>::max());
        }
        return ans;
    }
};
```
<a name="fqoyP"></a>
### MinCostFlow
```cpp
template <class T = i64>
struct MCFGraph {
    struct Edge {
        int y, c, f;
    };
    const int n;
    std::vector<Edge> e;
    std::vector<std::vector<int>> g;
    std::vector<T> h, dis;
    std::vector<int> pre;

    bool dijkstra(int s, int t) {
        dis.assign(n, std::numeric_limits<T>::max());
        pre.assign(n, -1);
        using pii = std::pair<T, int>;
        std::priority_queue<pii, std::vector<pii>, std::greater<>> q;
        dis[s] = 0;
        q.emplace(0, s);

        while (!q.empty()) {
            auto [D, x] = q.top();
            q.pop();

            if (dis[x] < D)
                continue;
            for (int i : g[x]) {
                const auto &[y, c, f] = e[i];
                if (c > 0 and dis[y] > D + h[x] - h[y] + f) {
                    dis[y] = D + h[x] - h[y] + f;
                    pre[y] = i;
                    q.emplace(dis[y], y);
                }
            }
        }
        return dis[t] != std::numeric_limits<T>::max();
    }
    MCFGraph(int n) : n(n), g(n) {}
    void add(int x, int y, int c, int f) {
        if (f < 0) { // ** 删除 <=> 最大流
            g[x].push_back(e.size());
            e.emplace_back(y, 0, f);
            g[y].push_back(e.size());
            e.emplace_back(x, c, -f);
        } else // **
            g[x].push_back(e.size()),
                e.emplace_back(y, c, f),
                g[y].push_back(e.size()),
                e.emplace_back(x, 0, -f);
    }
    std::pair<int, T> work(int s, int t) {
        int flow = 0;
        T cost = 0;
        h.assign(n, 0);
        while (dijkstra(s, t)) {
            for (int i = 0; i < n; ++i)
                h[i] += dis[i];
            int aug = std::numeric_limits<int>::max();
            for (int i = t; i != s; i = e[pre[i] ^ 1].y)
                aug = std::min(aug, e[pre[i]].c);
            for (int i = t; i != s; i = e[pre[i] ^ 1].y) {
                e[pre[i]].c -= aug;
                e[pre[i] ^ 1].c += aug;
            }
            flow += aug;
            cost += T(aug) * h[t];
        }
        return std::pair(flow, cost);
    }
};
```
<a name="Lvjzp"></a>
# 8, Calculate Geometry
<a name="FPRJh"></a>
## Point
```cpp
#define x first
#define y second

using ld = double;
constexpr ld eps = 1e-9;
int sgn(const ld &a) {
    return (a < -eps ? -1 : a > eps);
}
int sgn(const i64 &a) {
    return (a < 0 ? -1 : a > 0);
}

template <class T, class G>
struct Point : public std::pair<T, T> {
    using std::pair<T, T>::pair;

    Point operator+(const Point &a) const {
        return Point(this->x + a.x, this->y + a.y);
    }
    Point operator-(const Point &a) const {
        return Point(this->x - a.x, this->y - a.y);
    }
    Point operator-() const {
        return Point(-this->x, -this->y);
    }
    G operator*(const Point &a) const {
        return G(this->x) * a.x + G(this->y) * a.y;
    }
    G operator%(const Point &a) const {
        return G(this->x) * a.y - G(this->y) * a.x;
    }
    Point rot() const {
        return Point(-this->y, this->x);
    }
    Point rot(const double &th) const {
        Point b(cosl(th), sinl(th));
        return Point(this->x * b.x - this->y * b.y,
                     this->x * b.y + this->y * b.x);
    }
    friend Point operator*(const T &a, const Point &b) {
        return Point(a * b.x, a * b.y);
    }
    friend std::istream &operator>>(std::istream &in, Point &p) {
        return in >> p.x >> p.y;
    }
    friend std::ostream &operator<<(std::ostream &ot, const Point &p) {
        return ot << '(' << p.x << ", " << p.y << ')';
    }
};

template <class T, class G>
G dis2(const Point<T, G> &a, const Point<T, G> &b = Point<T, G>(0, 0)) {
    Point<T, G> p = a - b;
    return p * p;
}
template <class T, class G>
double dis(const Point<T, G> &a, const Point<T, G> &b = Point<T, G>(0, 0)) {
    return sqrtl(dis2(a, b));
}

using pii = Point<int, i64>;
using PS = std::vector<pii>;
```
<a name="W0C63"></a>
## Line
<a name="VgYCK"></a>
### Point-Line
```cpp
// a on Seg b-c
bool onSeg(const pii &a, pii b, pii c) {
    b = b - a, c = c - a;
    return sgn(b % c) == 0 and sgn(b * c) <= 0;
}

// a disTo Line b-c
ld disToLine(const pii &a, const pii &b, const pii &c) {
    pii v1 = b - c, v2 = a - c;
    return 1.L * std::abs(v1 % v2) / sqrtl(v1 * v1);
}

// a disTo Seg b-c
ld disToSeg(const pii &a, const pii &b, const pii &c) {
    if (sgn((a - b) * (c - b)) <= 0 or sgn((a - c) * (b - c)) <= 0)
        return std::min(dis(a, b), dis(a, c));
    return disToLine(a, b, c);
}

// a project to Line b-c (here pii = Point<ld, ld>)
pii foot(const pii &a, const pii &b, const pii &c) {
    pii u = a - b, v = c - b;
    return b + 1.L * (u * v) / (v * v) * v;
}

// a symmetry to Line b-c (here pii = Point<ld, ld>)
pii symmetry(const pii &a, const pii &b, const pii &c) {
    pii ft = foot(a, b, c);
    return 2 * ft - a;
}
```
<a name="SQ6bB"></a>
### Line-Line
```cpp
// Line a-b cross Line c-d (here pii = Point<ld, ld>)
pii cross(const pii &a, const pii &b, const pii &c, const pii &d) {
    pii v = c - d;
    ld sa = v % (a - d), sb = (b - d) % v;
    return 1.L * sa / (sa + sb) * (b - a) + a;
}

const ld pi = acosl(-1);
// a-b 和 a-c 的夹角 signed
ld getAngle(const pii &a, const pii &b, const pii &c) {
    auto v1 = b - a, v2 = c - a;
    return atan2l(v1 % v2, v1 * v2); // ∠bac
}

// 对四个不同的点判断四点共圆
// d在abc外接圆外return 1, 内return -1
int inCircle(const pii &a, const pii &b, const pii &c, const pii &d) {
    auto ag1 = getAngle(a, b, c), ag2 = getAngle(d, c, b);
    if (sgn(ag1) == sgn(ag2)) {
        return sgn(pi - std::abs(ag1 + ag2));
    } else {
        return sgn(std::abs(ag1) - std::abs(ag2));
    }
}
```
<a name="fIG8o"></a>
## Polygon
<a name="VSIDx"></a>
### ConvexHull
```cpp
PS convex(PS ps) {
    std::sort(ps.begin(), ps.end());
    ps.erase(std::unique(ps.begin(), ps.end()), ps.end());

    if (ps.size() <= 1)
        return ps; // 防止空点集RE

    auto cmp = [&](const pii &a, const pii &b) {
        pii A = a - ps[0], B = b - ps[0];
        return (sgn(A % B) == 0 ? dis2(A) < dis2(B) : sgn(A % B) > 0);
    };
    std::sort(ps.begin() + 1, ps.end(), cmp);

    PS q = ps;
    int r = 0;
    for (int i = 1; i < ps.size(); q[++r] = ps[i], i++)
        while (r and sgn((q[r] - q[r - 1]) % (ps[i] - q[r])) <= 0)
            r -= 1;
    return q.resize(r + 1), q;
}
```
<a name="MSSFS"></a>
### Area
```cpp
ld area(const PS &ps) {
    ld res = 0;
    for (int i = 0; i < ps.size(); i++) {
        int &&j = (i + 1) % ps.size();
        res += ps[i] % ps[j];
    }
    return res / 2;
}
```
<a name="kasy8"></a>
### Point-Polygon
```cpp
// #include "onSeg"

class InPoly {
    // Seg c-d is Cross Seg a-b
    bool isCross(const pii &a, const pii &b, const pii &c, const pii &d) const {
        pii ab = b - a, cd = d - c;
        if (sgn(ab % cd) == 0)
            return 0; // 共线,寄
        int r1 = sgn(ab % (c - a)), r2 = sgn(ab % (d - a));
        int g1 = sgn(cd % (a - c)), g2 = sgn(cd % (b - c));
        return !(r1 * r2 > 0 or r1 + r2 == -1 or g1 * g2 > 0);
        // 真相交或者 c-d 贴着 a-b 左边
    }

public:
    // a IN/OUT/ON Polygon
    std::string operator()(const pii &a, const PS &ps) const {
        int res = 0;
        pii b = {1 << 30, a.y};
        for (int i = 0; i < ps.size(); i++) {
            int j = (i + 1) % ps.size();
            if (onSeg(a, ps[i], ps[j]))
                return "ON";
            res += isCross(a, b, ps[i], ps[j]);
        }
        return (res % 2 ? "IN" : "OUT");
    }
};
const InPoly inPoly;

// a IN/OUT/ON Convex
std::string inConvex(const pii &a, const PS &ps) {
    if (a == ps[0])
        return "ON";
    if (ps.size() <= 1)
        return "OUT";
    if (ps.size() == 2)
        return onSeg(a, ps[0], ps[1]) ? "ON" : "OUT";
    auto v = a - ps[0];
    if ((ps[1] - ps[0]) % v < 0 or (ps.back() - ps[0]) % v > 0)
        return "OUT";
    int l = 1, r = ps.size() - 1;
    while (l + 1 < r) {
        auto mid = l + r >> 1;
        auto res = (ps[mid] - ps[0]) % v;
        if (res == 0)
            return (ps[mid] == a ? "ON" : (onSeg(a, ps[0], ps[mid]) ? "IN" : "OUT"));
        (res > 0 ? l : r) = mid;
    }
    auto res = (ps[r] - ps[l]) % (a - ps[l]);
    if (res == 0 or onSeg(a, ps[0], ps[l]) or onSeg(a, ps[0], ps[r]))
        return "ON";
    return (res > 0 ? "IN" : "OUT");
}
```
<a name="y2fdt"></a>
### Line-Polygon
```cpp
PS cutPoly(const PS &ps, pii a, pii b) {
    // 返回多边形 ps 在有向直线 a->b 左边的部分
    PS v;
    auto c = b - a;
    for (int i = 0; i < ps.size(); i++) {
        int j = (i + 1) % ps.size();
        auto cr1 = c % (ps[i] - a), cr2 = c % (ps[j] - a);
        if (sgn(cr1) >= 0)
            v.push_back(ps[i]);
        if (sgn(cr1) * sgn(cr2) == -1)
            v.push_back(cross(a, b, ps[i], ps[j]));
    }
    return v;
}


// find point of tangency
class BinarySearchOnConvex {
    PS vs;

    constexpr static bool sgn(const pii &a) {
        return (a.x > 0 or a.x == 0 and a.y > 0);
    }

    constexpr static bool cmp(const pii &a, const pii &b) {
        auto sa = sgn(a), sb = sgn(b);
        if (sa != sb)
            return sa > sb;
        return a % b > 0;
    }

public:
    constexpr BinarySearchOnConvex(const PS &ps) : vs(ps.size()) {
        int n = size(ps);
        for (int i = 0; i < n; i++) {
            int j = (i + 1) % n;
            vs[i] = ps[j] - ps[i];
        }
    }

    int operator()(const pii &a) const {
        int res = std::lower_bound(begin(vs), end(vs), a, cmp) - begin(vs);
        return res % size(vs);
    }
};
```
<a name="oTAG9"></a>
### Polygon-Polygon
```cpp
PS minkowski(const PS &a, const PS &b) {
    int n = a.size(), m = b.size();
    if (!n or !m)
        return PS();

    PS ps(1, a[0] + b[0]);
    int ap = 0, bp = 0;
    while (ap < n and bp < m) {
        auto va = a[(ap + 1) % n] - a[ap];
        auto vb = b[(bp + 1) % m] - b[bp];
        auto res = va % vb;
        if (res > 0)
            ps.push_back(ps.back() + va), ap++;
        if (res < 0)
            ps.push_back(ps.back() + vb), bp++;
        if (res == 0)
            ps.push_back(ps.back() + va + vb), ap++, bp++;
    }
    while (ap < n) {
        auto va = a[(ap + 1) % n] - a[ap];
        ps.push_back(ps.back() + va), ap++;
    }
    while (bp < m) {
        auto vb = b[(bp + 1) % m] - b[bp];
        ps.push_back(ps.back() + vb), bp++;
    }
    return ps.pop_back(), ps;
}
```
<a name="tejXx"></a>
## Examples
<a name="mck3h"></a>
### ClosestPair
```cpp
// 得到平面最近点对距离的平方
template <class T, class G>
G closestPair(std::vector<Point<T, G>> ps) {
    auto pf = [&](const T &x) {
        return 1ll * x * x;
    };

    G ans = 1ll << 60;
    std::function<void(int, int)> find = [&](int l, int r) {
        if (r - l <= 1)
            return;
        int mid = l + r >> 1;
        find(l, mid), find(mid, r);

        std::vector<Point<T, G>> b; // points in caught
        for (int i = l; i < r; i++)
            if (pf(std::fabs(ps[i].x - ps[mid].x)) < ans)
                b.push_back(ps[i]);

        std::sort(b.begin(), b.end(), [&](pii u, pii v) {
            return u.y < v.y; // sort by y
        });

        for (int i = 0; i < b.size(); i++)
            for (int j = i + 1; j < b.size() and pf(b[j].y - b[i].y) < ans; j++)
                ans = std::min(ans, dis2(b[j] - b[i]));
    };
    std::sort(begin(ps), end(ps)); // sort by x
    find(0, ps.size());

    return ans;
}
// O(nlog2n) but almost O(nlogn)
```
<a name="qUopI"></a>
### DelaunayTriangulation
```cpp
const ld pi = acosl(-1);
// a-b 和 a-c 的夹角 signed
ld getAngle(const pii &a, const pii &b, const pii &c) {
    auto v1 = b - a, v2 = c - a;
    return atan2l(v1 % v2, v1 * v2); // ∠bac
}

// 对四个不同的点判断四点共圆
// d在abc外接圆外return 1, 内return -1
int inCircle(const pii &a, const pii &b, const pii &c, const pii &d) {
    auto ag1 = getAngle(a, b, c), ag2 = getAngle(d, c, b);
    if (sgn(ag1) == sgn(ag2)) {
        return sgn(pi - std::abs(ag1 + ag2));
    } else {
        return sgn(std::abs(ag1) - std::abs(ag2));
    }
}

struct Delaunay {
    struct Edge {
        Edge *l, *r, *inv;
        int face, from;

        Edge() : face{0}, from{-1} {
            l = r = inv = nullptr;
        }
    };

    void addFace(Edge *e0, Edge *e1, Edge *e2, int fid) {
        e0->r = e2->l = e1;
        e1->r = e0->l = e2;
        e2->r = e1->l = e0;

        e0->face = e1->face = e2->face = fid;
        face[fid] = e0;
    }

    void init() {
        bucket.resize(2 * (n + 1));
        vertex.resize(n + 3);

        // 加三个点构造全局三角形
        ps.emplace_back(Inf * 10, Inf);
        ps.emplace_back(-Inf, Inf);
        ps.emplace_back(-Inf, -Inf * 10);

        std::vector<Edge *> e(6);
        for (auto &pt : e) {
            pt = new Edge();
        }

        e[0]->from = e[5]->from = n;
        e[1]->from = e[3]->from = n + 1;
        e[2]->from = e[4]->from = n + 2;

        for (int i = 0; i < 3; i++) {
            e[i]->inv = e[i + 3];
            e[i + 3]->inv = e[i];
        }

        face.resize(2);
        addFace(e[0], e[1], e[2], 0);
        addFace(e[5], e[4], e[3], 1);

        for (int i = 0; i < 3; i++)
            vertex[i + n] = e[i];

        bucket[1].resize(n);
        std::iota(begin(bucket[1]), end(bucket[1]), 0);
        belong.assign(n, 1);
    }

    std::vector<Edge *> vertex, face;
    std::vector<std::vector<int>> bucket;
    std::vector<int> belong;
    PS ps;
    const int n;

    static constexpr int Inf = int(1e8);
    Delaunay(const PS &initPs) : n{(int)initPs.size()}, ps(initPs) {
        init();

        for (int i = 0; i < n; i++) {
            addPoint(i);

            Edge *e = vertex[i];
            for (int j = 0; j < 3; j++) {
                if (flip(i, e->r->inv)) {
                    e = e->l->inv;
                    j -= 2;
                } else {
                    e = e->inv->r;
                }
            }
        }
    }

    bool inFace(pii p, int fid) {
        int cnt = 0;
        Edge *e = face[fid];
        do {
            auto a = ps[e->from];
            e = e->r;
            auto b = ps[e->from];

            cnt += (b - a) % (p - a) > 0;
        } while (e != face[fid]);

        return cnt == 3 or cnt == 0;
    }

    void addPoint(int pid) {
        int fid[3] = {belong[pid], face.size(), face.size() + 1};

        auto pointS = bucket[fid[0]];
        bucket[fid[0]].clear();

        face.resize(face.size() + 2);

        std::vector<Edge *> e(6);
        for (int i = 0; i < 6; i++) {
            e[i] = new Edge();
        }
        for (int i = 0; i < 6; i++) {
            e[i]->inv = e[i ^ 1];
        }

        e[1]->from = e[3]->from = e[5]->from = pid;
        Edge *tmpe = face[fid[0]];
        int r = 0, l = 5;
        for (int j = 0; j < 3; j++) {
            Edge *tmpr = tmpe->r;

            int b = tmpr->from;
            e[r]->from = b;

            addFace(e[l], tmpe, e[r], fid[j]);

            r = (r + 2) % 6, l = (l + 2) % 6;
            tmpe = tmpr;
        }

        vertex[pid] = e[1];

        for (int p : pointS) {
            if (p == pid)
                continue;
            for (int j = 0; j < 3; j++) {
                if (inFace(ps[p], fid[j])) {
                    belong[p] = fid[j];
                    bucket[fid[j]].push_back(p);
                    break;
                }
            }
        }
    }

    bool flip(int pid, Edge *e) {
        if (e->face == 0)
            return false;

        auto fir = e->inv->l;

        int inCircleRes = inCircle(ps[pid], ps[fir->r->from],
                                   ps[fir->l->from], ps[e->l->from]);

        if (inCircleRes >= 0)
            return false;

        int fid[2] = {fir->face, e->face};

        auto pointS = bucket[fid[0]];
        pointS.insert(pointS.end(), bucket[fid[1]].begin(), bucket[fid[1]].end());

        bucket[fid[0]].clear();
        bucket[fid[1]].clear();

        auto ev = e->inv;
        auto sec = e->r, thir = e->l, four = fir->l;

        addFace(fir, sec, ev, fid[0]);
        addFace(e, thir, four, fid[1]);
        e->from = pid;
        ev->from = thir->from;

        for (int p : pointS) {
            for (int j = 0; j < 2; j++) {
                if (inFace(ps[p], fid[j])) {
                    belong[p] = fid[j];
                    bucket[fid[j]].push_back(p);
                    break;
                }
            }
        }
        return true;
    }
};
```
