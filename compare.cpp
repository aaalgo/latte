#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <boost/progress.hpp>
#include <boost/program_options.hpp>
#include <boost/python/numpy.hpp>
#include "latte.h"

using namespace std;
using namespace latte;
namespace py = boost::python;
namespace np = boost::python::numpy;

#if 0
void test_auc () {
    vector<pair<float, int>> all;
    size_t n0 = 113;
    size_t n1 = 119;
    for (unsigned i = 0; i < n0; ++i) {
        all.emplace_back(0, 0);
    }
    for (unsigned i = 0; i < n1; ++i) {
        all.emplace_back(0, 1);
    }
    sort(all.begin(), all.end(), [](pair<float, int> const &p1, pair<float, int> const &p2) { return p1.second < p2.second;});
    if (auc(all, n0, n1) != 1) throw 0;
    cerr << "Test 1 OK" << endl;
    sort(all.begin(), all.end(), [](pair<float, int> const &p1, pair<float, int> const &p2) { return p1.second > p2.second;});
    if (auc(all, n0, n1) != 0) throw 0;
    cerr << "Test 2 OK" << endl;
    random_shuffle(all.begin(), all.end());
    float v = auc(all, n0, n1);
    if (v < 0.3) throw 1;
    if (v > 0.7) throw 1;
    cerr << "Test 3 OK" << endl;
}
#endif

template <typename T>
void compare (np::ndarray v0, np::ndarray v1, vector<float> *ft) {
    size_t n0 = v0.shape(0);
    size_t n1 = v1.shape(0);
    size_t ns = n0 + n1;
    cerr << n0 << '\t' << n1 << endl;
    if (v0.shape(1) != Ref::GENES) throw 0;
    if (v1.shape(1) != Ref::GENES) throw 0;

    T const *p0 = (T const *)(v0.get_data());
    T const *p1 = (T const *)(v1.get_data());

    ft->resize(Ref::GENES);
    boost::progress_display progress(Ref::GENES, cerr);
#pragma omp parallel
    {
        vector<pair<T, int>> all;
        all.reserve(ns);
        // for each gene
#pragma omp for
        for (size_t gene = 0; gene < Ref::GENES; ++gene) {
            all.clear();
            for (size_t i = 0; i < n0; ++i) {
                all.emplace_back(p0[Ref::GENES * i + gene], 0);
            }
            for (size_t i = 0; i < n1; ++i) {
                all.emplace_back(p1[Ref::GENES * i + gene], 1);
            }
            sort(all.begin(), all.end());
            ft->at(gene) = auc(all, n0, n1);
#pragma omp critical
            ++progress;
        }
    }
}


int main (int argc, char *argv[]) {
    string input1_path;
    string input2_path;
    string output_path("output");
    {
        namespace po = boost::program_options;
        po::options_description desc_visible("Allowed options");
        desc_visible.add_options()
            ("help,h", "produce help message.")
            ("input1", po::value(&input1_path), "")
            ("input2", po::value(&input2_path), "")
            ("output", po::value(&output_path), "")
            ;

        po::options_description desc("Allowed options");
        desc.add(desc_visible);

        po::positional_options_description p;
        p.add("input1", 1);
        p.add("input2", 1);
        p.add("output", 1);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).
                         options(desc).positional(p).run(), vm);
        po::notify(vm);
        if (vm.count("help") || input1_path.empty() || input2_path.empty()) {
            cout << "Usage:" << endl;
            cout << desc_visible;
            cout << endl;
            return 0;
        }
	}
    Py_Initialize();
    np::initialize();
    Ref genes("data/ref");
    py::object np_load = py::import("numpy").attr("load");
    np::ndarray rank0 = np::array(np_load(input1_path));
    np::ndarray rank = np::array(np_load(input2_path));

    vector<float> ft;
    if (np::dtype::get_builtin<float>() == rank0.get_dtype()) {
        compare<float>(rank0, rank, &ft);
    }
    else if (np::dtype::get_builtin<uint16_t>() == rank0.get_dtype()) {
        compare<uint16_t>(rank0, rank, &ft);
    }
    else {
        cerr << "dtype not supported." << endl;
        throw 0;
    }


    if (ft.size() != genes.size()) throw 0;
    ofstream os(output_path);
    for (unsigned i = 0; i < Ref::GENES; ++i) {
        os << genes[i] << '\t' << ft[i] << endl;
    }

    return 0;
}
