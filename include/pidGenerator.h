

#ifndef PID_GENERATOR
#define PID_GENERATOR 


#include "tofGenerator.h"
#include "dedxGenerator.h"

#include "pidGeneratorConfig.h"

#include "TF1.h"

#include "iostream"

using namespace std;

#include "jdbUtils.h"
using namespace jdbUtils;

class pidGenerator
{

/**
 * Protected member properties
 */
protected:
	TRandom3 * rGen;
	tofGenerator * tof;
	dedxGenerator * dedx;

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
	// one for each species and the total, for each of those one for each p bin
	vector< TH2D * > hDedxTof[ nSpecies + 1 ];
	int nBinsP;

	// Tof vs. P Hists
	TH2D* hTof[ nSpecies + 1 ];

	// Dedx vs. P Hists
	TH2D* hDedx[ nSpecies + 1 ];
	
	pidGeneratorConfig config;

/**
 * Publicly visible actions
 */
public:
	pidGenerator( pidGeneratorConfig con ){
		this->config = con;

		cout << "[pidGenerator.pidGenerator] " << endl;
		cout << "\tpT Range [ " << config.pCutLow << " -> " << config.pCutHigh << " ] " << endl;

		// Create the pt distribution functions
		// They each use the same functional form, with different parameters
		string formula = "abs([0])*exp(-(x/abs([1])))*(1.0+[2]*pow(x,[3]))";
		ptPi = new TF1( "ptPi", formula.c_str(), config.pCutLow, config.pCutHigh );
		ptK = new TF1( "ptK", formula.c_str(), config.pCutLow, config.pCutHigh );
		ptP = new TF1( "ptP", formula.c_str(), config.pCutLow, config.pCutHigh );

		/**
		 * From Evan Sangaline
		 * I do not know how these parameters where determined
		 */
		ptPi->SetParameters(515, 0.2329, 0.001179, 7.171);
		ptK->SetParameters(4.414, 0.3156, 10.25, 0.2371);
		ptP->SetParameters(0.2682, 0.2892, 262.5, 1.289);

		ptTotal[ 0 ] = ptPi->Integral( config.pCutLow, config.pCutHigh );
		ptTotal[ 1 ] = ptK->Integral( config.pCutLow, config.pCutHigh );
		ptTotal[ 2 ] = ptP->Integral( config.pCutLow, config.pCutHigh );
		ptSum = ptTotal[ 0 ] + ptTotal[ 1 ] + ptTotal[ 2 ];

		cout << endl << "\tpT distribution composition" << endl;
		for ( int i = 0; i < nSpecies; i++ ){
			cout << "\t\t" << species( i ) << ": " << (100.0 * (ptTotal[ i ] / ptSum)) << " %" << endl;

		}

		// setup the random number generator
		// Use the unique seed provided by ROOT
		rGen = new TRandom3( 0 );

		// open the output file
		fOutput = new TFile( config.outFilename.c_str(), "RECREATE" );


		// Initialize the dedx and tof generators
		dedx = new dedxGenerator( );
		tof = new tofGenerator( config.oneBetaSigma, config.nTofSamples ); 


		// prapare the histograms
		prepareHistograms();
	}
	
	~pidGenerator(){
		cout << "[pidGenerator.~pidGenerator]" << endl;

		fOutput->Write();
	}

	
	
	void generate( ){

		cout << "[pidGenerator." << __FUNCTION__ << "]" << endl;
		cout << "Generating " << config.nEvents << " Events" << endl;

		for ( unsigned long long iEvent = 0; iEvent < config.nEvents; iEvent++ ){
			progressBar( iEvent, config.nEvents );


			const int id = randomID();
			const double realP = randomPt( id );

			const double measuredP = blurPt( realP, config.smearP );
			const int pBin = int( TMath::Floor( (measuredP - config.pCutLow ) / config.pBinWidth ) );
			//const double centerP = config.pCutLow + config.pBinWidth * pBin + config.pBinWidth * 0.5;
			if ( pBin < 0 || pBin > nBinsP ){
				iEvent--;
				continue;
			} // momentum outside of range

			double measuredDedx = dedx->random( measuredP, masses[ id ] );
			double measuredTof = tof->random( measuredP, masses[ id ] );

			// mismatch some events
			if ( rGen->Uniform() < config.tofMismatch ){
				int newID = randomID();
				// TODO: add momentum smearing here also?
				measuredTof = tof->random( randomPt( newID ), masses[ newID ] );
			}

			// dedx merging
			if ( rGen->Uniform() * 3.1415926 < (3.1415926 * config.phiMismatch ) ){
				int newID = randomID();
				const double tempP = randomPt( newID );
				if ( TMath::Abs( ( tempP / measuredP ) - 1.0 ) < config.pMismatch ){
					measuredDedx += dedx->random( tempP, masses[ newID ] );
				}
			}

			
			

			// Fill the Histograms
			// the tof vs. P
			hTof[ id ]->Fill( measuredP, measuredTof );
			hTof[ nSpecies ]->Fill( measuredP, measuredTof );
			//dedx vs. P
			hDedx[ id ]->Fill( measuredP, measuredDedx );
			hDedx[ nSpecies ]->Fill( measuredP, measuredDedx );

			// in centered coords
			measuredDedx = TMath::Log10( measuredDedx );
			
			double dedxMeanK = TMath::Log10( dedx->mean( measuredP, masses[ config.centerSpecies ] )  );
			double tofMeanK = tof->mean( measuredP, masses[ config.centerSpecies ] );

			measuredDedx = ( measuredDedx - dedxMeanK ) / config.centerSigmaDedx ;
			measuredTof = ( measuredTof - tofMeanK ) / config.centerSigmaTof ;

			hDedxTof[ id ][ pBin ]->Fill( measuredDedx, measuredTof );
			hDedxTof[ nSpecies ][ pBin ]->Fill( measuredDedx, measuredTof );



		} // end loop events


	}


/**
 * Protected Implementations
 */
protected:
	void prepareHistograms(  ) {

		cout << "[pidGenerator." << __FUNCTION__ << "]" << endl;

		// determine the number of p bins to use
		double pRange = config.pCutHigh - config.pCutLow;
		nBinsP = (int)(pRange / config.pBinWidth) + 1;

		double lowP = config.pCutLow;
		double highP = lowP + config.pBinWidth;

		for ( int i = 0; i < nBinsP; i++ ){

			double p = (lowP + highP) / 2.0; // take the center of the p bin

			// compute the histogram limits
			double kDedxMean = TMath::Log10( dedx->mean( p, masses[ 1 ] ) );
			double kTofMean = tof->mean( p, masses[ 1 ] );

			double 	lowDedxMean = 0,
					highDedxMean = 0,
					lowTofMean = 0,
					highTofMean = 0;

			// loop over the species and determine the bounds
			for ( int ip = 0; ip < nSpecies; ip ++ ){
				double nDedxMean = ( TMath::Log10( dedx->mean( p, masses[ ip ] )) - kDedxMean ) / config.centerSigmaDedx;
				double nTofMean = ( tof->mean( p, masses[ ip ]) - kTofMean) / config.centerSigmaTof ;

				if ( nDedxMean < lowDedxMean ) lowDedxMean = nDedxMean;
				if ( nDedxMean > highDedxMean ) highDedxMean = nDedxMean;

				if ( nTofMean < lowTofMean ) lowTofMean = nTofMean;
				if ( nTofMean > highTofMean ) highTofMean = nTofMean;
			}

		
			double dedxPad = ( highDedxMean - lowDedxMean ) * config.paddingScaleDedx;
			double tofPad = ( highTofMean - lowTofMean ) * config.paddingScaleTof;
			lowDedxMean -= (config.paddingDedx + dedxPad);
			highDedxMean += (config.paddingDedx + dedxPad);

			lowTofMean -= (config.paddingTof + tofPad);
			highTofMean += (config.paddingTof + tofPad);

			double nBinsDedx = ( highDedxMean - lowDedxMean ) / config.dedxBinWidth;
			double nBinsTof = ( highTofMean - lowTofMean ) / config.tofBinWidth;

			//cout << "dedx: " << lowDedxMean << " -> " << highDedxMean << endl;

			for ( int ip = 0; ip < nSpecies + 1; ip++ ){
				
				string name = "h_dedx_tof_p" + ts(ip) + "_b"+ts(i);
				string title = species( ip ) + " : " + ts( lowP, 4 ) + " <= P < " + ts( highP, 4 ) + "; N_{#sigma} dedx ; N_{#sigma} 1/#beta";

				TH2D* temp = new TH2D( name.c_str(), title.c_str(), 
								nBinsDedx, lowDedxMean, highDedxMean, 
								nBinsTof, lowTofMean, highTofMean );

				hDedxTof[ ip ].push_back( temp );
			}

			lowP += config.pBinWidth;
			highP += config.pBinWidth;
		} // loop on p Bins


		// for the vs P plots
		double nBinsTof = ( config.tofCutHigh - config.tofCutLow ) / config.tofVsPBinWidth;
		double nBinsDedx = ( config.dedxCutHigh - config.dedxCutLow ) / config.dedxVsPBinWidth;

		// create the tof and dedx histos
		for ( int ip = 0; ip < nSpecies + 1; ip ++ ){

			string name = "h_tof_p" + ts(ip);
			string title = species( ip ) + " : tof ; P [GeV]; 1/#beta";

			hTof[ ip ] = new TH2D( name.c_str(), title.c_str(),
									nBinsP, config.pCutLow, config.pCutHigh,
									nBinsTof, config.tofCutLow, config.tofCutHigh );

			name = "h_dedx_p" + ts(ip);
			title = species( ip ) + " : dedx ; P [GeV] ; dedx [KeV/cm]";

			hDedx[ ip ] = new TH2D( name.c_str(), title.c_str(),
									nBinsP, config.pCutLow, config.pCutHigh,
									nBinsDedx, config.dedxCutLow, config.dedxCutHigh );

		}

	}


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
		if ( nSpecies == index)
			return "All";
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
			return ptPi->GetRandom(config.pCutLow, config.pCutHigh);
    	else if ( 1 == id ) 
    		return ptK->GetRandom(config.pCutLow, config.pCutHigh);
    	else if( 2 == id ) 
    		return ptP->GetRandom(config.pCutLow, config.pCutHigh);

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