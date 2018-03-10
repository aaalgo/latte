#include <fstream>
#include <string>
#include <vector>
#include <utility>

namespace latte {
    using std::string;
    using std::vector;
    using std::ifstream;
    using std::pair;

    class Ref: public vector<string> {
    public:
        static constexpr size_t GENES = 54675;
        Ref (string const &path) {
            ifstream is(path.c_str());
            string line, gene;
            float value;
            getline(is, line);
            while (is >> gene >> value) {
                push_back(gene);
            }
            if (size() != GENES) throw 0;
        }
    };

    template <typename T>
    float auc (vector<pair<T, int>> const &all, size_t n0, size_t n1) {
        if (n0 + n1 != all.size()) throw 0;
        size_t width = 0, height = 0;
        float area = 0;
        for (auto const &p: all) {
            if (p.second == 0) {
                height += 1;
            }
            else {
                area += height;
                width += 1;
                if (width == n1) break;
            }
        }
        return area / (n0 * n1);
}


}
