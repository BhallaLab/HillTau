/************************************************************************
 * This program is part of HILLTAU, a framework for fast compact
 * abstractions of signaling events.
 * Copyright	(C) 2021	Upinder S. Bhalla and NCBS
 * It is made available under the terms of the
 * GNU Public License version 3 or later.
 * See the file COPYING.LIB for the full notice.
************************************************************************/

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <cmath>
#include <exprtk.hpp>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl_bind.h>
using namespace std;
#include "htHeader.h"
PYBIND11_MAKE_OPAQUE(std::vector<double>);
PYBIND11_MAKE_OPAQUE(std::map<string, MolInfo>);
PYBIND11_MAKE_OPAQUE(std::map<string, ReacInfo>);
PYBIND11_MAKE_OPAQUE(std::map<string, EqnInfo>);

const double INTERNAL_DT_SCALE = 0.02;

MolInfo::MolInfo( const std::string& name_, const std::string& grp_, int isSub_ = 0, double concInit_ = -1, Model* model_ = 0 ):
	name(name_),
	grp( grp_ ),
	order( 0 ),
	index( 0 ),
	isSub( isSub_ ),
	explicitConcInit( false ),
	model( model_ )
{
	// The local value of concInit is only used in reconcileConcInit.
	// All access to concInit should be through the set/get functions
	// which refer to the system array of concInits.
	if ( concInit_ < 0.0 ) {
		concInit = 0.0;
	} else {
		concInit = concInit_;
		explicitConcInit = true;
	}
};

void MolInfo::setConcInit( double val ) {
	model->concInit[ index ] = val;
}

double MolInfo::getConcInit() const {
	return model->concInit[ index ];
}

void MolInfo::reconcileConcInit() {
	model->concInit[ index ] = concInit;
}

void MolInfo::privateSetConcInit( double val ) {
	concInit = val;
}

//////////////////////////////////////////////////////////////////////

ReacInfo::ReacInfo( const string& name_, const string& grp_, 
	const vector< string >& subs_, 
	const map< string, double>& reacObj, 
	Model* model_ ):
	name(name_),
	grp( grp_ ),
	KA( reacObj.at("KA") ),
	tau( reacObj.at("tau") ),
	tau2( 1.0 ),
	Kmod( 1.0 ),	/// Default halfmax of modifier
	Amod( 4.0 ),	/** Default alpha factor of modifier
					** Values of Amod < 1 make it an inhibitory modifier
					** Values of Amod > 1 make excitatory. 
					** Amod == 1 does not modify
					*/
	Nmod( 1.0 ),	/// Order of modifier
	gain( 1.0 ),
	baseline( 0.0 ),
	inhibit( 0 ),
	prdIndex( 0 ),
	kh( 1.0 ),
	HillCoeff( 1.0 ),
	isBuffered( 0 ),
	overrideConcInit( false ),
	subs( subs_ ),
	hillIndex( 0 ),
	reagIndex( 0 ),
	modIndex( ~0U ),
	oneSub( false ),
	model( model_)
{
	tau2 = tau;
	prdIndex = model->molInfo.at(name)->index;
	if ( subs.size() == 0 ) {
		throw std::range_error( "Error: Reaction " + name + " has zero reagents\n" );
	}
	reagIndex = model->molInfo.at( subs[0] )->index;
	hillIndex = model->molInfo.at( subs.back() )->index;
	overrideConcInit = !model->molInfo.at( name )->explicitConcInit;
	int numUnique = 1;
	if ( reagIndex != hillIndex ) { // At least two subs
		if (subs.size() == 2) {
			numUnique = 2;
		} else if ( subs.size() > 2 ) {
			if ( subs.back() != subs[1] ) {	// A modifier too
				numUnique = 3;
				modIndex = model->molInfo.at(subs[1])->index;
			} else { // We have a reagent and multiple Hill ligands, no mod
				numUnique = 2;
			}
		}
	}
	oneSub = ( numUnique == 1);
	HillCoeff = subs.size() + 1 - numUnique;
	kh = pow( KA, HillCoeff);

	auto t = reacObj.find( "tau2" );
	if ( t != reacObj.end() ) {
		tau2 = t->second;
	}
	auto inh = reacObj.find( "inhibit" );
	if ( inh != reacObj.end() ) {
		inhibit = inh->second;
	}
	auto i = reacObj.find( "Kmod" );
	if ( i != reacObj.end() ) {
		if ( numUnique != 3) {
			cout << "Warning: Reaction " << name << " has <3 reagents but Kmod has been specified. Ignoring.\n";
			Kmod = 1.0;
		} else {
			Kmod = i->second;
		}
	} else if ( numUnique == 3)  {
		throw std::range_error( "Error: Reaction " + name + " has 3 reagents but no Kmod specified for modifier.\n");
	}
	auto a = reacObj.find( "Amod" );
	if ( a != reacObj.end() ) {
		Amod = a->second;
	}
	auto n = reacObj.find( "Nmod" );
	if ( n != reacObj.end() ) {
		Nmod = n->second;
	}
	auto b = reacObj.find( "baseline" );
	if ( b != reacObj.end() ) {
		baseline = b->second;
	}
	auto g = reacObj.find( "gain" );
	if ( g != reacObj.end() ) {
		gain = g->second;
	}
	auto ib = reacObj.find( "isBuffered" );
	if ( ib != reacObj.end() ) {
		isBuffered = ib->second;
	}
	auto ci = reacObj.find( "concInit" );
	if ( ci != reacObj.end() ) {
		setConcInit( ci->second );
		// cout << "Initialize CONCINIT1 = " << getConcInit()  << ", buff = " << isBuffered << endl;
	}
}

void ReacInfo::setKA( double val ) {
	KA = val;
	kh = pow( KA, HillCoeff);
}

double ReacInfo::getKA() const {
	return KA;
}

void ReacInfo::setConcInit( double val ) {
	model->molInfo[name]->privateSetConcInit( val );
	model->concInit[ prdIndex ] = val;
	if ( isBuffered )
		model->conc[ prdIndex ] = val;
}

double ReacInfo::getConcInit() const {
	return model->concInit[ prdIndex ];
}

double ReacInfo::eval( double dt ) const
{
	if ( isBuffered ) {
		// cout << "CONCINIT2 = " << getConcInit() << endl;
		return model->conc[ prdIndex ] = getConcInit();
		// return getConcInit();
	}
	double orig = model->conc[ prdIndex ] - baseline;
	double delta = concInf( model->conc ) - orig;
	if ( delta >= 0.0 ) {
		delta *= 1.0 - exp( -dt/tau );
	} else {
		delta *= 1.0 - exp( -dt/tau2 );
	}
	double ret = baseline + orig + delta;
	if (ret < 0.0 ) {
		throw std::range_error( "Error: negative value on: " + name );
	}
	model->conc[ prdIndex ] = ret;
	return ret;
}

double ReacInfo::concInf( const vector< double >& conc ) const
{
	double h = pow( conc[ hillIndex ], HillCoeff );
	double mod = 1.0;
	// cout << name << ":	vec = " << conc[reagIndex] << 	"	" << conc[hillIndex] << "	" << conc[modIndex] << endl;
	if ( modIndex != ~0U ) {
		// mod = conc[ modIndex ] / Kmod;
		double x = pow( conc[ modIndex ] / Kmod, Nmod );
		mod = ( 1.0 + x ) / ( 1.0 + Amod * x );
	}
	if ( oneSub ) {
		return h / KA;
	}

	double s = conc[ reagIndex ] * gain;
	if ( inhibit ) {
		return s * (1.0 - h / (h + kh * mod ) );
	} else {
		return s * h / (h + kh * mod);
	}
}

int ReacInfo::getReacOrder() const
{
	int ret = 0;
	for (auto s = subs.begin(); s != subs.end(); s++ ) {
		int mo = model->molInfo.at( *s )->order;
		if (mo < 0)
			return -1;
		if ( ret < mo )
			ret = mo;
	}
	ret += 1;
	model->molInfo.at(name)->order = ret;
	return ret;
}

////////////////////////////////////////////////////////////////////

/*
vector<unsigned int> EqnInfo::findMolTokens( const string& eqn )
{
	int isInMol = 0;
	vector< unsigned int > ret;

	for ( unsigned int i = 0; i < eqn.length(); i++ ) {
		if (!isInMol ) {
			if (eqn[i] == 'e' || eqn[i] == 'E' ) {
				if ( i < eqn.length() - 1 ) {
					if ( eqn[i] == '+' || eqn[i] == '-' || isalnum(eqn[i]) )
						continue;
				}
			}
			if ( isalpha( eqn[i] ) ) {
				isInMol = 1;
				ret.push_back( i );
			} else {
				continue;
			}
		} else {
			if ( isalnum( eqn[i] ) || eqn[i] == '_' ) {
				continue;
			} else {
				ret.push_back( i );
				isInMol = 0;
			}
		}
	}

	if (isInMol)
		ret.push_back( eqn.length() );
	if ( ret.size() % 2 != 0 )
		throw "Error: equation token not ended? " + eqn;

	return ret;
}
*/

EqnInfo::EqnInfo( const string& name_, const string& grp_, 
			const string& eqnStr_, const vector< string >& eqnSubs, 
			Model* model_ ):
	name(name_),
	grp( grp_ ),
	eqnStr( eqnStr_ ),
	isBuffered( 0 ),
	subs( eqnSubs ),
	molIndex( 0 ),
	model( model_ )

{
	for ( const auto& s: subs ) {
		auto mi = model->molInfo.find( s );
		if ( mi != model->molInfo.end() ) {
			symbol_table.add_variable( s, model->conc[ mi->second->index ]);
		} else {
			throw std::invalid_argument( "Error: Unable to find variable '" + s + "' in equation " + eqnStr );
		}
	}
	/**
	auto tokens = findMolTokens( eqnStr_ );
	for ( auto i = tokens.begin(); i != tokens.end(); i += 2 ) {
		string sstr = eqnStr_.substr( *i, *(i+1) - *i );
		auto mi = molInfo.find( sstr );
		if ( mi != molInfo.end() ) {
			symbol_table.add_variable( sstr, conc[ mi->second->index ] );
			subs.push_back( sstr );
		}
	}
	*/
	if ( eqnStr.substr(0,9) == "concInit=" ) {
		isBuffered = 1;
		if ( eqnStr.size() > 9 ) {
			setConcInit( stod( eqnStr.substr(9,string::npos) ) );
		} else {
			setConcInit( 0.0 );
			throw std::invalid_argument( "Error: missing float value in '" + eqnStr + "'" );
		}

	} else {
		symbol_table.add_constants();

		expression.register_symbol_table( symbol_table );
		exprtk::parser< double > parser;
		parser.compile( eqnStr, expression );
	}
	molIndex = model->molInfo.at( name )->index;
};

double EqnInfo::eval() const
{
	if ( isBuffered ) {
		model->conc[molIndex] = model->concInit[ molIndex ];
	} else {
		model->conc[molIndex] = expression.value();
	}
	return model->conc[molIndex];
}

void EqnInfo::setConcInit( double val ) {
	model->molInfo[name]->privateSetConcInit( val );
	model->conc[ molIndex ] = model->concInit[ molIndex ] = val;
}

double EqnInfo::getConcInit() const {
	return model->concInit[ molIndex ];
}

////////////////////////////////////////////////////////////////////
Model::Model()
	: 
			currentTime( 0.0 ),
			step( 0 ),
			dt( 1.0 )
{;}

void Model::setReacSeqDepth( int maxDepth )
{
	if ( maxDepth < 1 )
		throw std::range_error( "Error: maxDepth must be >= 1" );
	sortedReacInfo.clear();
	sortedReacInfo.resize( maxDepth );
	sortedEqnInfo.clear();
	for ( auto eri = eqnInfo.begin(); eri != eqnInfo.end(); eri++ ) {
		sortedEqnInfo.push_back( eri->second );
	}
}

void Model::assignReacSeq( const string& name, int seq )
{
	auto ri = reacInfo.at( name ); // Assume it is good.
	sortedReacInfo[seq].push_back( ri );
}

bool onnit( vector< const ReacInfo* >::const_iterator ri, const vector< string >& saveList )
{
	for ( auto f = saveList.begin(); f != saveList.end(); f++ ) {
		if ( *f == (*ri)->name || *f == (*ri)->grp )
			return true;
	}
	return false;
}

bool eonnit( const EqnInfo* eri, const vector< string >& saveList )
{
	for ( auto f = saveList.begin(); f != saveList.end(); f++ ) {
		if ( *f == eri->name || *f == eri->grp ) {
			return true;
		}
	}
	return false;
}


void Model::modifySched( const vector< string >& saveList, const vector< string >& deleteList )
{
	unsigned int numSeq = sortedReacInfo.size();
	vector< vector< const ReacInfo* > > newsri(numSeq);
	if ( saveList.size() > 0 ) {
		for (unsigned int seq = 0; seq < numSeq; seq++ ) {
			auto sri = &(sortedReacInfo[seq] );
			for( auto ri = sri->begin(); ri != sri->end(); ri++ ) {
				if ( deleteList.size() > 0 ) {
					if ( onnit( ri, saveList ) && !onnit(ri, deleteList) ){
						newsri[seq].push_back( *ri );
					} else {
						reacInfo.erase( (*ri)->name );
					}
				} else {
					if ( onnit( ri, saveList ) ) {
						newsri[seq].push_back( *ri );
					} else {
						reacInfo.erase( (*ri)->name );
					}
				}
			}
		}
        sortedReacInfo.clear();
		for (auto sri = newsri.begin(); sri != newsri.end(); sri++ ) {
            if ( sri->size() > 0 ) {
                sortedReacInfo.push_back( *sri );
			}
        }
		sortedEqnInfo.clear();
		for ( auto eri = eqnInfo.begin(); eri != eqnInfo.end(); eri++ ) {
			if ( deleteList.size() > 0 ) {
				if ( eonnit( eri->second, saveList ) && !eonnit( eri->second, deleteList ) ) {
					sortedEqnInfo.push_back( eri->second );
				} else {
					eqnInfo.erase( eri->second->name );
				}
			} else {
				if ( eonnit( eri->second, saveList ) ) {
					sortedEqnInfo.push_back( eri->second );
				} else {
					eqnInfo.erase( eri->second->name );
				}
			}
		}
	} else if ( deleteList.size() > 0 ) {
		for (unsigned int seq = 0; seq < numSeq; seq++ ) {
			auto sri = &(sortedReacInfo[seq] );
			for( auto ri = sri->begin(); ri != sri->end(); ri++ ) {
				if ( !onnit( ri, deleteList ) ) {
					newsri[seq].push_back( *ri );
				}
			}
		}
        sortedReacInfo.clear();
		for (auto sri = newsri.begin(); sri != newsri.end(); sri++ ) {
            if ( sri->size() > 0 ) {
                sortedReacInfo.push_back( *sri );
			}
        }
		sortedEqnInfo.clear();
		for ( auto eri = eqnInfo.begin(); eri != eqnInfo.end(); eri++ ) {
			if ( !eonnit( eri->second, deleteList ) ) {
					sortedEqnInfo.push_back( eri->second );
			}
		}
	}
	// If both lists are empty, just retain original sortedReacInfo.
}

void Model::advance( double runtime, int settle )
{
	if (runtime < 10e-6) return;
	if (settle) {
		/// double newdt = runtime / 10.0;
		double newdt = runtime * INTERNAL_DT_SCALE;
		innerAdvance( runtime, newdt );
	} else {
		double newdt = min( dt, internalDt);
		double adv = max( minTau * 10.0, dt );
		if ( newdt >= runtime / 2.0 ) {
			newdt = pow( 10.0, floor( log10( runtime / 2.0 ) ) );
			innerAdvance( runtime, newdt );
		} else if ( 2.0 * adv < runtime  )  {
			// Advance a few small timesteps, then switch to longer ones
			innerAdvance( adv, newdt );
			innerAdvance( runtime - adv, dt );
		} else { // All small dt
			innerAdvance( runtime, newdt );
		}
	}
}

void Model::advancePlot( double t )
{
	if ( plotNames.size() == 0 ) {
		plotvec.push_back( conc ); // This will need to be transposed.
	} else {
		for (unsigned int ii = 0; ii < plotIdx.size(); ii++ ) {
			plotvec[ii].push_back( conc[plotIdx[ii]] );
		}
	}
	plotTimePoints.push_back( t );
}

void Model::innerAdvance( double runtime, double newdt )
{
	for (double t = 0.0; t < runtime; t += newdt ) {
		if ( newdt > (runtime - t) )
			newdt = runtime - t;
		for (auto r = sortedReacInfo.begin(); r != sortedReacInfo.end(); 
						r++) {
			for (auto ri = r->begin(); ri != r->end(); ri++ ) {
				(*ri)->eval( newdt );
			}
		}
		for (auto e = sortedEqnInfo.begin(); e != sortedEqnInfo.end(); ++e ) {
			(*e)->eval();
		}

		if ( floor( (currentTime + t + newdt ) / dt ) > step ) {
			advancePlot( currentTime + t + newdt );
			step += 1;
		}
	}
	currentTime += runtime;
}

void Model::allocConc()
{
	concInit.resize( molInfo.size(), 0.0 );
	for ( auto m = molInfo.begin(); m != molInfo.end(); m++ ) {
		m->second->reconcileConcInit();
	}
	conc = concInit;
}

double neatRound( double x )
{
	if (x <= 0.0)
		return 0.0;
	double y = pow( 10.0, floor( log10( x ) ) );
	if ( (x/y) >= 5.0 )
		return y * 5.0;
	else if ( (x/y) >= 2.0 )
		return y * 2.0;
	return y;
}

void Model::reinit()
{
	// Logic: Any explicitly defined initialization value is to be used
	// as is. This is happens if ReacInfo::overrideConcInit is false.
	// Any others need to be estimated from the steady-state 
	// value of the reacns.
	currentTime = 0.0;
	step = 0;
	internalDt = dt;
	minTau = 1e20; // dt should be < 0.25x smallest tau at input.
	for (auto r = sortedReacInfo.begin(); r != sortedReacInfo.end(); r++) {
		for (auto ri = r->begin(); ri != r->end(); ri++) {
			if ( (*ri)->isBuffered ) { // ConcInit is used unchanged.
				// cout << "concinit for " << (*ri)->name << " in reinit = " << (*ri)->getConcInit() << endl;
				continue;
			}
			minTau = min( min( minTau, (*ri)->tau ), (*ri)->tau2 );
			if ( (*ri)->overrideConcInit ) {
				unsigned int j = (*ri)->prdIndex;
				if ((*ri)->inhibit ) {
					concInit[j] = (*ri)->concInf( concInit ) + (*ri)->baseline;
					if ( concInit[j] < 0.0 )
						concInit[j] = 0.0;
				} else {
					concInit[j] = (*ri)->baseline;
				}
			}
		}
	}
	if ( dt > INTERNAL_DT_SCALE * minTau ) {
		internalDt = neatRound( INTERNAL_DT_SCALE * minTau );
	}
	auto ci = concInit.begin();
	for (auto c = conc.begin(); c < conc.end(); c++, ci++ ) {
		*c = *ci;
	}

	plotvec.clear();
	if ( plotIdx.size() > 0) {
		assert( plotIdx.size() <= molInfo.size() );
		plotvec.resize( plotIdx.size() );
	}
	plotTimePoints.clear();
	advancePlot( 0.0 );
}

void Model::makeReac( const string & name, const string & grp, 
				const vector< string >& subs, 
				const map< string, double >& reacObj )
{
	auto r = new ReacInfo( name, grp, subs, reacObj, this );
	reacInfo[ name ] = r;
	// If it is a reac, then by definition we don't yet know its order
	molInfo[name]->order = -1;
	// Override group of product mol it is == grp of reac.
	molInfo[name]->grp = grp;
}

void Model::makeMol( const string & name, const string & grp, double concInit = -1.0, int isSub = 0 )
{
	auto mi = molInfo.find( name );
	if ( mi == molInfo.end() ) { // Make new one.
		auto m = new MolInfo( name, grp, isSub, concInit, this );
		m->index = molInfo.size();
		molInfo[ name ] = m;
	} else {
		if ( !isSub ) {
			if ( mi->second->isSub ) { // OK, override grp assignment of sub
				mi->second->grp = grp;
				mi->second->isSub = 0;
			} else if ( mi->second->grp != grp ){ //BAD: same obj in new grp
				throw std::invalid_argument( "Error: Molecule '" + name + "' redefined in new grp '" + grp + "', already in '" + mi->second->grp + "'.\n" );
			}
		}
		if ( concInit >= 0.0 ) {
			// Could be second pass assignment of species with concInits.
			// Unfortunately the Model::concInit array hasn't yet been built
			mi->second->privateSetConcInit( concInit );
			mi->second->explicitConcInit = true;
		}
	}
}

void Model::makeEqn( const string & name, const string & grp, const string& expr, const vector< string >& eqnSubs )

{
	auto e = new EqnInfo( name, grp, expr, eqnSubs, this );
	eqnInfo[ name ] = e;
	molInfo[ name ]->order = 0; // We assume that eqns do not cascade. 
	// We evaluate all eqns after all the reacs are done, so 0 is good
	// Override group of output mol it is == grp of reac.
	molInfo[name]->grp = grp;
}

void Model::addGrp( const string& grpname )
{
	grpInfo.push_back( grpname );
}

void Model::addPlot( const string& plotname )
{
	try {
		auto mi = molInfo.at( plotname );
		plotNames.push_back( plotname );
		plotIdx.push_back( mi->index );
    } catch ( const std::out_of_range& e ) {
		cerr << "Error: HillTau::addPlot: Molecule '" + plotname + "' not found.\n" ;
		std::exit( EXIT_FAILURE );
    }
}

vector< double > Model::getConcVec( int index ) const
{
	vector< double > ret( plotvec.size() );
	if ( plotvec.size() == 0 || plotvec[0].size() <= unsigned( index ) )
		return ret;
	for( unsigned int i = 0; i < ret.size(); ++i )
		ret[i] = plotvec[i][index];
	return ret;
}

int Model::getMolOrder( const string& molName ) const
{
	return molInfo.at(molName)->order;
}

bool Model::updateMolOrder( int maxOrder, const string& molName ) const
{
	auto mi = molInfo.at( molName );
	if (mi->order < 0) {
		mi->order = maxOrder;
		return true;
	}
	return false;
}
