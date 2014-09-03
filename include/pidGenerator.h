

#ifndef PID_GENERATOR
#define PID_GENERATOR 


#include "tofGenerator.h"
#include "dedxGenerator.h"

#include "TF1.h"

#include "iostream"
#include "string"
using namespace std;

class pidGenerator
{
protected:
	TRandom3 * rGen;
	tofGenerator * tof;
	dedxGenerator * dedx;

	double ptCutLow, ptCutHigh;

	/**
	 * pT distributions for each species
	 */
	TF1 * ptPi;		// 0
	TF1 * ptK;		// 1
	TF1 * ptP;		// 2

	static const int nSpecies = 3;
	static const double masses[  ];
	double ptTotal[ nSpecies ];
	double ptSum;

	/**
	 * Output data file
	 */
	string output;
	TFile * fOutput;

	/**
	 * Histograms
	 */

public:
	pidGenerator( double ptCutLow, double ptCutHigh, string output = "pidSimEvents.root" ){
		cout << "[pidGenerator.pidGenerator] " << endl;
		cout << "\tpT Range [ " << ptCutLow << " -> " << ptCutHigh << " ] " << endl;

		this->ptCutLow = ptCutLow;
		this->ptCutHigh = ptCutHigh;
		this->output = output;
		cout << "\tWriting to : " << output << endl;

		// Create the pt distribution functions
		// They each use the same functional form, with different parameters
		string formula = "abs([0])*exp(-(x/abs([1])))*(1.0+[2]*pow(x,[3]))";
		ptPi = new TF1( "ptPi", formula.c_str(), ptCutLow, ptCutHigh );
		ptK = new TF1( "ptK", formula.c_str(), ptCutLow, ptCutHigh );
		ptP = new TF1( "ptP", formula.c_str(), ptCutLow, ptCutHigh );

		/**
		 * From Evan Sangaline
		 * I do not know how these parameters where determined
		 */
		ptPi->SetParameters(515, 0.2329, 0.001179, 7.171);
		ptK->SetParameters(4.414, 0.3156, 10.25, 0.2371);
		ptP->SetParameters(0.2682, 0.2892, 262.5, 1.289);

		ptTotal[ 0 ] = ptPi->Integral( ptCutLow, ptCutHigh );
		ptTotal[ 1 ] = ptK->Integral( ptCutLow, ptCutHigh );
		ptTotal[ 2 ] = ptP->Integral( ptCutLow, ptCutHigh );
		ptSum = ptTotal[ 0 ] + ptTotal[ 1 ] + ptTotal[ 2 ];

		cout << endl << "\t pT distribution composition" << endl;
		for ( int i = 0; i < nSpecies; i++ ){
			cout << "\t" << species( i ) << ": " << (100.0 * (ptTotal[ i ] / ptSum)) << " %" << endl;
		}

		// setup the random number generator
		// Use the unique seed provided by ROOT
		rGen = new TRandom3( 0 );

		// open the output file
		fOutput = new TFile( output.c_str(), "RECREATE" );

	}
	
	~pidGenerator(){
		cout << "pidGenerator.~pidGenerator]" << endl;
	}

	
	

	void prepareHistograms( double pBinWidth, double tofBinWidth, double dedxBinWidth, double padding = 5.0 ) {

	}


/**
 * Protected Implementations
 */
protected:

	/**
	 * Maps a string particle species name to the integer index
	 * @param  name Name of particle species
	 * @return      integer index for the given particle species or -1 if no match was found.
	 */
	int species( string name ){
		if ( "P" == name || "p" == name )
			return 2;
		if ( "K" == name || "k" == name )
			return 1;
		if ( "Pi" == name || "pi" == name || "PI" == name )
			return 0;
		return -1;
	}
	/**
	 * Maps an integer particle species index to the human readable name
	 * @param  index Index of particle species
	 * @return      Human readable string for the given particle index
	 */
	string species( int index ){
		
		if ( 0 == index )
			return "Pi";
		else if ( 1 == index )
			return "K";
		else if ( 2 == index )
			return "P";
		return "";
	}

	/**
	 * Generates a random particle species id. Uses the pt distributions to generate the 
	 * correct relative # of each species
	 * @return return the integer index for the chosen species
	 */
	int randomID() {
		const double u = rGen->Uniform( 0, ptSum );
		if ( u < ptTotal[ 0 ] )
			return 0;
		else if ( u < ( ptTotal[ 0 ] + ptTotal[ 1] ) )
			return 1;
		else 
			return 2;
	}

	/**
	 * Samples the pt distribution for the given species to produce a random pt value in the pre-specified range
	 * @param  id Particle species id
	 * @return    returns a random momentum value for the species.
	 */
	double randomPt( int id ){
		if ( 0 == id ) 
			return ptPi->GetRandom(ptCutLow, ptCutHigh);
    	else if ( 1 == id ) 
    		return ptK->GetRandom(ptCutLow, ptCutHigh);
    	else if( 2 == id ) 
    		return ptP->GetRandom(ptCutLow, ptCutHigh);

    	return 0;
	}

	/**
	 * Produces the effects of momentum smearing
	 * @param  pT    Momentum value before smearing
	 * @param  sigma sigma value for smearing
	 * @return       Randomly smeared momentum value
	 */
	double blurPt( double pT, double sigma ){
		return 1.0 / rGen->Gaus( 1.0 / pT, sigma/pT );
	}
};

/**
 * Particle Species Masses in order of index
 */
const double pidGenerator::masses[] = { 0.139, 0.493, 0.938 };


#endif