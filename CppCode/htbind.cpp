#include <string>
#include <map>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <exprtk.hpp>
using namespace std;

namespace py = pybind11;
#include "htHeader.h"

PYBIND11_MAKE_OPAQUE(std::vector<double>);
PYBIND11_MAKE_OPAQUE(std::map<string, MolInfo>);
PYBIND11_MAKE_OPAQUE(std::map<string, ReacInfo>);
PYBIND11_MAKE_OPAQUE(std::map<string, EqnInfo>);

PYBIND11_MODULE(ht, m) {
	py::bind_vector<std::vector<double>>(m, "VectorDouble");

	/////////////////////////////////////////////////////////////////////
    py::class_<MolInfo>(m, "MolInfo")
        .def( 
			py::init<const std::string &, const std::string &, int, double, Model* >(), py::arg("name"), py::arg("grp"), py::arg("isSub") = 0, py::arg( "concInit" ), py::arg( "model" ) )
		.def_readwrite("name", &MolInfo::name)
		.def_readwrite("grp", &MolInfo::grp)
		.def_readwrite("order", &MolInfo::order)
		.def_property("concInit", &MolInfo::getConcInit, &MolInfo::setConcInit)
		.def_readonly("isSub", &MolInfo::isSub)
		.def_readwrite("index", &MolInfo::index);

	/////////////////////////////////////////////////////////////////////
    py::class_<ReacInfo>(m, "ReacInfo")
        .def( 
			py::init<const std::string &, const std::string &, const vector< string >&, const map< string, double>&, Model* >() )
		.def_readwrite("name", &ReacInfo::name)
		.def_readwrite("grp", &ReacInfo::grp)
		.def_property("KA", &ReacInfo::getKA, &ReacInfo::setKA)
		.def_readwrite("tau", &ReacInfo::tau)
		.def_readwrite("tau2", &ReacInfo::tau2)
		.def_readwrite("Kmod", &ReacInfo::Kmod)
		.def_readwrite("Amod", &ReacInfo::Amod)
		.def_readwrite("Nmod", &ReacInfo::Nmod)
		.def_readwrite("gain", &ReacInfo::gain)
		.def_readwrite("baseline", &ReacInfo::baseline)
		.def_readwrite("inhibit", &ReacInfo::inhibit)
		.def_readwrite("prdIndex", &ReacInfo::prdIndex)
		.def_readwrite("kh", &ReacInfo::kh)
		.def_readonly("HillCoeff", &ReacInfo::HillCoeff)
		.def_property("concInit", &ReacInfo::getConcInit, &ReacInfo::setConcInit)
		.def_readwrite("isBuffered", &ReacInfo::isBuffered)
		.def_readonly("subs", &ReacInfo::subs)
		.def( "eval", &ReacInfo::eval, "Evaluator for Reacs" )
		.def( "getReacOrder", &ReacInfo::getReacOrder, "Returns 1+largest of substrate mol orders, and updates ReacMol accordingly." )
		.def( "concInf", &ReacInfo::concInf, "Computes steady-state value of reaction output" );

	/////////////////////////////////////////////////////////////////////
    py::class_<EqnInfo>(m, "EqnInfo")
        .def( 
			py::init<const std::string &, const std::string &, const std::string&, const vector< string >&, Model* >() )
		.def_readwrite("name", &EqnInfo::name)
		.def_readwrite("grp", &EqnInfo::grp)
		.def_readwrite("eqnStr", &EqnInfo::eqnStr)
		.def_property("concInit", &EqnInfo::getConcInit, &EqnInfo::setConcInit)
		.def_readwrite("isBuffered", &EqnInfo::isBuffered)
		.def_readonly("subs", &EqnInfo::subs)
		.def( "eval", &EqnInfo::eval, "Evaluator for Eqns" );

	/////////////////////////////////////////////////////////////////////
    py::class_<Model>(m, "Model")
        .def(py::init())
		.def_readwrite("molInfo", &Model::molInfo)
		.def_readwrite("reacInfo", &Model::reacInfo)
		.def_readwrite("eqnInfo", &Model::eqnInfo)
		.def_readwrite("grpInfo", &Model::grpInfo)
		.def_readonly("plotNames", &Model::plotNames)
		.def_readonly("plotTimePoints", &Model::plotTimePoints)
		.def_readwrite("namedConsts", &Model::namedConsts)
		.def_readonly("currentTime", &Model::currentTime)
		.def_readwrite("dt", &Model::dt)
		.def_readwrite("internalDt", &Model::internalDt)
		.def_readonly("minTau", &Model::minTau)
		.def_readwrite("conc", &Model::conc)
		.def_readwrite("concInit", &Model::concInit)
		.def_readonly("plotvec", &Model::plotvec)
		.def( "makeMol", &Model::makeMol, "Create MolInfo object.", py::arg("name"), py::arg("grp"), py::arg("concInit") = -1.0, py::arg("isSub") = 0 )
		.def( "makeReac", &Model::makeReac, "Create ReacInfo object.", py::arg("name"), py::arg("grp"), py::arg("subs"), py::arg("reacParms"))
		.def( "makeEqn", &Model::makeEqn, "Create EqnInfo object.", py::arg("name"), py::arg("grp"), py::arg("expr"), py::arg( "eqnSubs" ) )
		.def( "addGrp", &Model::addGrp, "Append grpname string to grpInfo vector.", py::arg("grpname") )
		.def( "addPlot", &Model::addPlot, "Append plot name string to plotNames vector.", py::arg("plotname") )
		.def( "setReacSeqDepth", &Model::setReacSeqDepth, "Defines how deep is the sequence of reactions, that is, the size of sortedReacInfo.")
		.def( "assignReacSeq", &Model::assignReacSeq, "Builds up sortedReacOrder vectors.")
		.def( "advance", &Model::advance, "Advances the simulation", py::arg( "runtime" ), py::arg( "settle" ) = 0 )
		.def( "reinit", &Model::reinit, "Reinits all conc values" )
		.def( "allocConc", &Model::allocConc, "Allocates and initializes conc vectors" )
		.def( "getConcVec", &Model::getConcVec, "Returns vector of doubles of conc as a function of time for specified mol index." )
		.def( "getMolOrder", &Model::getMolOrder, "Returns order of named molecule.", py::arg( "molName" ) )
		.def( "updateMolOrder", &Model::updateMolOrder, "Checks if order of named molecule is <0, if so updates it and returns True.", py::arg( "maxOrder"), py::arg( "molName" ) )
		.def( "modifySched", &Model::modifySched, "Modifies scheduling to retain/eliminate subsets of reactions and groups.", py::arg("saveList"), py::arg("deleteList") )
		;
}

