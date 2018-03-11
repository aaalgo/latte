#include <string>
#include <boost/python/numpy.hpp>
#include <file/CELFileData.h>

namespace py = boost::python;
namespace np = boost::python::numpy;

namespace {
    np::ndarray read_cel (std::string const &path) {
        affxcel::CCELFileData cel;
        cel.SetFileName(path.c_str());
        cel.Read();
        int cols = cel.GetCols();
        int rows = cel.GetRows();
        np::ndarray array = np::zeros(py::make_tuple(rows, cols), np::dtype::get_builtin<float>());
        float *ptr = (float *)array.get_data();
        affxcel::CELFileEntryType e;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                cel.GetEntry(j, i, e);
                *ptr = e.Intensity;
                ++ptr;
            }
        }
        return array;
    }
}

BOOST_PYTHON_MODULE(cLatte) {
    np::initialize();
    def("read_cel", ::read_cel);
}
