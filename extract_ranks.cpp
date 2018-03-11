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

template <typename T>
void extract_ranks (np::ndarray v0, np::ndarray v1, unsigned gene, vector<float> *ranks) {
    if (np::dtype::get_builtin<T>() != v1.get_dtype()) throw 0;
    size_t n0 = v0.shape(0);
    size_t n1 = v1.shape(0);
    size_t ns = n0 + n1;
    if (v0.shape(1) != Ref::GENES) throw 0;
    if (v1.shape(1) != Ref::GENES) throw 0;

    T const *p0 = (T const *)(v0.get_data());
    T const *p1 = (T const *)(v1.get_data());

    vector<pair<T, int>> all;
    for (size_t i = 0; i < n0; ++i) {
        all.emplace_back(p0[Ref::GENES * i + gene], 0);
    }
    for (size_t i = 0; i < n1; ++i) {
        all.emplace_back(p1[Ref::GENES * i + gene], 1);
    }
    sort(all.begin(), all.end());
    ranks->clear();
    for (unsigned i = 0; i < all.size(); ++i) {
        if (all[i].second == 1) {
            ranks->push_back(1.0 * i / all.size());
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

    ifstream is("interest");
    map<string, unsigned> lookup;
    for (unsigned i = 0; i < genes.size(); ++i) {
        lookup[genes[i]] = i;
    }
    string probe, gene;
    while (is >> gene >> probe) {
        auto it = lookup.find(probe);
        if (it == lookup.end()) throw 0;
        unsigned ps_id = it->second;
        vector<float> ranks;
        if (np::dtype::get_builtin<float>() == rank0.get_dtype()) {
            extract_ranks<float>(rank0, rank, ps_id, &ranks);
        }
        else if (np::dtype::get_builtin<uint16_t>() == rank0.get_dtype()) {
            extract_ranks<uint16_t>(rank0, rank, ps_id, &ranks);
        }
        else {
            cerr << "dtype not supported." << endl;
            throw 0;
        }
        ofstream os("ranks/" + gene + "_" + probe);
        for (float r: ranks) {
            os << r << endl;
        }
    }


    return 0;
}
