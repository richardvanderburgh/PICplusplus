// ExampleWrapper.cpp

#include <boost/python.hpp>
#include "init.hpp"

BOOST_PYTHON_MODULE(init)
{
    using namespace boost::python;

    def("initialize", init::initialize);
}