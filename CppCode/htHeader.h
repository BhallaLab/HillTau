/************************************************************************
 * This program is part of HILLTAU, a framework for fast compact
 * abstractions of signaling events.
 * Copyright	(C) 2021	Upinder S. Bhalla and NCBS
 * It is made available under the terms of the
 * GNU Public License version 3 or later.
 * See the file COPYING.LIB for the full notice.
************************************************************************/

class Model;

class MolInfo
{
	public:
			MolInfo( const string& name, const string& grp, double concInit );
			string name;
			string grp;
			int order;
			double concInit;
			unsigned int index;
			bool explicitConcInit;
};	

class ReacInfo
{
	public:
			ReacInfo( const string& name, const string& grp, 
			const vector< string >& subs, 
			const map< string, double>& reacObj, 
			const map< string, MolInfo* >& molInfo );
			string name;
			string grp;
			double KA;
			double tau;
			double tau2;
			double Kmod;
			double gain;
			double baseline;
			int inhibit;
			int prdIndex;
			double kh;
			double HillCoeff;
			double concInf( const vector< double >& conc ) const;
			double eval( Model* model, double dt ) const;
			bool overrideConcInit;
			vector< string > subs;

	private:
			unsigned int hillIndex;
			unsigned int reagIndex;
			unsigned int modIndex;
			bool oneSub;

};

class EqnInfo
{
	public:
			EqnInfo( const string& name, const string& grp, const string& eqnStr, const vector< string >& eqnSubs, const map< string, MolInfo* >& molInfo, vector< double >& conc );
			string name;
			string grp;
			string eqnStr;
			double eval( vector<double>& conc ) const;
			static vector< unsigned int > findMolTokens(const string& eqn);
			vector< string > subs;
	private:
			// Stuff for parser
			unsigned int molIndex;
			exprtk::symbol_table<double> symbol_table;
			exprtk::expression<double> expression;
};

class Model
{
	public:
			Model();
			map< string, MolInfo* > molInfo;
			map< string, ReacInfo* > reacInfo;
			map< string, EqnInfo* > eqnInfo;
			double currentTime;
			double dt;
			vector< double > conc;
			vector< double > concInit;
			vector< vector< double > > plotvec;
			
			void makeMol( const string & name, const string & grp, double concInit );
			void makeReac( const string & name, const string & grp, const vector< string >& subs, const map< string, double >& reacObj );
			void makeEqn( const string & name, const string & grp, const string& expr, const vector< string >& eqnSubs );
			void setReacSeqDepth( int order );
			void assignReacSeq( const string& name, int seq );
			void advance( double runtime, int settle );
			void allocConc();
			void parseEqns();
			void reinit();
	private:
			vector< vector< const ReacInfo* > > sortedReacInfo;
};
