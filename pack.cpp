#include <cstdint>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <boost/program_options.hpp>
#include <boost/python/numpy.hpp>
#include "latte.h"

namespace py = boost::python;
namespace np = boost::python::numpy;

using namespace std;
using namespace latte;

int main (int argc, char *argv[]) {
    string ref_path;
	string output_prefix;
    {
        namespace po = boost::program_options;
        po::options_description desc_visible("Allowed options");
        desc_visible.add_options()
            ("help,h", "produce help message.")
            ("version,v", "print version.")
            ("ref", po::value(&ref_path)->default_value("data/ref"), "")
            ("output,o", po::value(&output_prefix)->default_value("out"), "output prefix")
            ;

        po::options_description desc("Allowed options");
        desc.add(desc_visible);

        po::positional_options_description p;
        p.add("output", 1);


        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).
                         options(desc).positional(p).run(), vm);
        po::notify(vm);
        if (vm.count("help")) {
            cout << "Usage:" << endl;
            cout << desc_visible;
            cout << endl;
            return 0;
        }
	}

    Ref genes(ref_path);
    vector<string> paths;
    {
        string path;
        while (getline(cin, path)) {
            paths.push_back(path);
        }
        cerr << paths.size() << " paths loaded." << endl;
    }

    Py_Initialize();
    np::initialize();

    np::ndarray exprs = np::zeros(py::make_tuple(paths.size(), genes.size()), np::dtype::get_builtin<float>());
    np::ndarray ranks = np::zeros(py::make_tuple(paths.size(), genes.size()), np::dtype::get_builtin<uint16_t>());

    ofstream index(output_prefix + ".idx");
    vector<pair<float, uint16_t>> pv(genes.size());
    for (unsigned i = 0; i < paths.size(); ++i) {
        float *expr = (float *)(exprs.get_data()) + i * Ref::GENES;
        uint16_t *rank = (uint16_t *)(ranks.get_data()) + i * Ref::GENES;
        index << paths[i] << endl;

        ifstream is(paths[i]);
        string line, gene;
        float value;
        getline(is, line);
        uint16_t j = 0;
        while (is >> gene >> value) {
            if (gene != genes[j]) throw 0;
            expr[j] = value;
            pv[j] = make_pair(value, j);
            ++j;
        }
        if (j != genes.size()) throw 0;
        sort(pv.begin(), pv.end());
        for (unsigned j = 0; j < genes.size(); ++j) {
            rank[pv[j].second] = j;
        }
    }
    py::object np_save = py::import("numpy").attr("save");
    np_save(output_prefix + ".expr", exprs);
    np_save(output_prefix + ".rank", ranks);
    return 0;
}

