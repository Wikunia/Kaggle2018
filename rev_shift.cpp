#include <bits/stdc++.h>
using namespace std;

constexpr int maxShift = 5; //5
constexpr int cntNeighbors = 30;
constexpr int cntCities = 197769;

class GlobalIndex {
public:
    array<double, 10> penalty_coefficients{1.1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    vector<vector<int>> neighbors;
    bitset<cntCities> primes;
    array<double, cntCities> x, y;
    
    GlobalIndex()
        : neighbors(cntCities+1, vector<int>(100)) {
        primes.set();
        primes[0] = primes[1] = false;
        for (int p = 2; p * p <= cntCities; ++p)  if (primes[p])
            for (int i = p*2; i < cntCities; i += p)
                primes[i] = false;
        // Read coordinates
        ifstream file("cities.csv");
        assert(file && "cities.csv");
        file.ignore(numeric_limits<streamsize>::max(), '\n');
        char delimiter;
        for (int i; file >> i >> delimiter >> x[i] >> delimiter >> y[i];);
        // Read neighbors for not depends from libboost
        ifstream file2("neighbors100.bin", ios::binary);
        assert(file2 && "neighbors100.bin");
        for (auto& arr : neighbors)
            file2.read(reinterpret_cast<char*>(arr.data()), arr.size()*sizeof(arr[0]));
    }
    inline double dist(int l, int r) const {
        const double dx = x[l] - x[r];
        const double dy = y[l] - y[r];
        return sqrt(dx*dx + dy*dy);
    }
    
    inline double dist(int l, int r, int i) const {
        if (primes[l])  return dist(l, r);
        return penalty_coefficients[(i+1)%10] * dist(l, r);
    }
};
const GlobalIndex gi;

class Index {
public:
    explicit Index(const vector<int>& ids) {
        positions.resize(ids.size());
        for (auto& diff_distance : diff_distances)
            diff_distance.resize(ids.size());
        for (auto& diff_distance : rdiff_distances)
            diff_distance.resize(ids.size());
        pure_dist.resize(ids.size()-1);
        Update(ids);
    }
    void Update(const vector<int>& ids) {
        for (int i = 0; i < ids.size(); ++i)
            positions[ids[i]] = i;
        
        for (int i = 0; i < pure_dist.size(); ++i)
            pure_dist[i] = gi.dist(ids[i], ids[i+1]);
        
        for (int j = 0; j < diff_distances.size(); ++j)
            for (int i = 1; i < diff_distances[0].size(); ++i) {
                double value = pure_dist[i-1];
                if (!gi.primes[ids[i-1]])
                    value *= gi.penalty_coefficients[(i+j)%10];
                diff_distances[j][i] = diff_distances[j][i-1] + value;
            }
        
        for (int j = 0; j < rdiff_distances.size(); ++j)
            for (int i = 1, sz = rdiff_distances[0].size(); i < sz; ++i) {
                double value = pure_dist[i-1];
                if (!gi.primes[ids[i]])
                    value *= gi.penalty_coefficients[(sz-i+j+1)%10];
                rdiff_distances[j][i] = rdiff_distances[j][i-1] + value;
            }
    }
    array<vector<double>, 10> diff_distances;
    array<vector<double>, 10> rdiff_distances;
    vector<int> positions;
    vector<double> pure_dist;
    
    inline double getScore(const vector<int>& ids, int beg, int end) const {
        return diff_distances[0][end] - diff_distances[0][beg];
    }
    inline double getScoreReverse(const vector<int>& ids, int beg, int end, int shift = 0) const {
        int i = (beg + end + shift)%10;
        return rdiff_distances[i][end] - rdiff_distances[i][beg];
    }
    inline double getScore(const vector<int>& ids, int beg, int end, int shift) const {
        if (beg <= end) {
            int i = shift ? (cntCities + shift + 1)%10 : 0;
            return diff_distances[i][end] - diff_distances[i][beg];
        }
        return getScoreReverse(ids, end, beg, shift);
    }
    double getScore() const {
        return diff_distances[0].back();
    }
};
auto loadIds(const string& name) {
    vector<int> result;
    ifstream file(name);
    assert(file);
    file.ignore(numeric_limits<streamsize>::max(), '\n');
    int id;
    while (file >> id)
        result.emplace_back(id);
    return result;
}

void saveIds(const vector<int>& ids, const string& name) {
    ofstream file(name);
    file << "Path";
    for (const auto& id : ids)
        file << endl << id;
}

auto compareReverse = [](const Index& index, const vector<int>& ids, int beg, int end) {
    if (end - beg < 3)
        return numeric_limits<double>::lowest();
    double sum = index.getScore(ids, beg, end);
    sum -= index.getScoreReverse(ids, beg + 1, end - 1);
    sum -= gi.dist(ids[beg], ids[end - 1], beg);
    sum -= gi.dist(ids[beg + 1], ids[end], end - 1);
    return sum;
};

auto compareShiftRightN = [](const int n) {
    return [n](const Index& index, const vector<int>& ids, int beg, int end) {
        if (end >= ids.size() - 1 || end - beg < n + 1)
            return numeric_limits<double>::lowest();
        double sum = index.getScore(ids, beg, end + 1);
        sum -= index.getScore(ids, beg + 1, end - n, n);
        for (auto i = 0; i < n - 1; ++i)
            sum -= gi.dist(ids[end - i], ids[end - i - 1], beg + i + 1);
        sum -= gi.dist(ids[beg], ids[end], beg);
        sum -= gi.dist(ids[end - n + 1], ids[beg + 1], beg+n);
        sum -= gi.dist(ids[end - n], ids[end + 1], end);
        return sum;
    };
};

auto compareShiftLeftN = [](const int n) {
    return [n](const Index& index, const vector<int>& ids, int beg, int end) {
        if (beg == 0 || end - beg < n + 1)
            return numeric_limits<double>::lowest();
        double sum = index.getScore(ids, beg - 1, end);
        sum -= index.getScore(ids, beg+n, end - 1, -n);
        for (auto i = 0; i < n - 1; ++i)
            sum -= gi.dist(ids[beg + 1 + i], ids[beg + i], end - 2 - i);
        sum -= gi.dist(ids[beg - 1], ids[beg + n]  , beg - 1);
        sum -= gi.dist(ids[end - 1], ids[beg + n - 1], end - n - 1);
        sum -= gi.dist(ids[beg]  , ids[end]    , end - 1);
        return sum;
    };
};

struct Result {
    Result(int l, int r) : beg(min(l ,r)), end(max(l, r)) {}
    double dist;
    int beg, end;
    friend bool operator<(const Result& l, const Result& r) {return l.dist < r.dist;}
};

auto bestPairNeighbors = [](const Index& index, const vector<int>& ids, auto comp) {
    Result best(0, cntCities);
    best.dist = comp(index, ids, best.beg, best.end);
    for (int it = 0; it < ids.size(); ++it)
        for (int i = 0; i < cntNeighbors; ++i ) {
            Result tmp(it, index.positions[gi.neighbors[ids[it]][i]]);
            tmp.dist = comp(index, ids, tmp.beg, tmp.end);
            best = max(tmp, best);
        }
    return best;
};

auto shiftLeftN = [](const int n) {
    return [n](vector<int>& ids, int beg, int end) {
        rotate(begin(ids) + beg, begin(ids) + beg + n, begin(ids)+end);
        reverse(begin(ids) + end - n, begin(ids) + end);
    };
};

auto shiftRightN = [](const int n) {
    return [n](vector<int>& ids, int beg, int end) {
        auto it = rbegin(ids) + ids.size() - 1;
        rotate(it - end, it - end + n, it - beg);
        reverse(begin(ids) + beg + 1, begin(ids) + beg + n + 1);
    };
};

void reverseMover(vector<int>& ids, int beg, int end) {
    reverse(begin(ids) + beg + 1, begin(ids) + end);
}

template<class Func>
void fastMover(vector<int>& ids, Func func) {
    cout << fixed << setprecision(4);
    vector<function<double(const Index&, const vector<int>&, int, int)>> comparators{compareReverse};
    vector<function<void(vector<int>&, int, int)>>  movers{reverseMover};
    for (int i = 1; i <= maxShift; ++i) {
        comparators.emplace_back(compareShiftLeftN(i));
        comparators.emplace_back(compareShiftRightN(i));
        movers.emplace_back(shiftLeftN(i));
        movers.emplace_back(shiftRightN(i));
    }
    
    Index index(ids);
    while (true) {
        vector<Result> results;
        vector<future<Result>> futures;
        for (auto& comp : comparators)
            futures.emplace_back(async(func, index, ids, comp));
        for (auto& item : futures)
            results.emplace_back(item.get());
        
        auto best = max_element(begin(results), end(results)) - begin(results);
        if (results[best].dist < 0.0001)
            return;
        movers[best](ids, results[best].beg, results[best].end);
        index.Update(ids);
        cout << "Score: " << index.getScore() << setw(10) << -results[best].dist << endl;
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <input.tour> <output.tour>" << endl;
        return 1;
    }
    auto ids = loadIds(argv[1]);
    fastMover(ids, bestPairNeighbors);
    saveIds(ids, argv[2]);
    return 0;
}