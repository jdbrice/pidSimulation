
#include <iostream>
#include "allroot.h"
#include "constants.h"
#include "histoBook.h"
#include "xmlConfig.h"

#include "pidGenerator.h"
#include "pidGeneratorConfig.h"

pidGeneratorConfig fromConfig( xmlConfig * config );

int main( int argc, char* argv[] ) {


	if ( argc >= 2 ){
		xmlConfig config( argv[ 1 ] );
		config.report();

		// load the simulation options from the config file
		pidGeneratorConfig pidGenConfig = fromConfig( &config );

		// create a new pidSimulation generator from the config
		pidGenerator *pidGen = new pidGenerator( pidGenConfig );
		
		// generate the events
		pidGen->generate(  );

		// generate a quick report pdf
		pidGen->report();

		delete pidGen;
	}


	return 0;
}

/**
 * Loads the pidGenerator configuration into the class from an xmlConfig object
 * @param  config xmlConfig object with the config options specified
 * @return        the pidGeneratorConfig with the specified values loaded in
 */
pidGeneratorConfig fromConfig( xmlConfig * config ){

	pidGeneratorConfig pgc;

	pgc.outFilename 		= config->getString( "output", pgc.outFilename );
	pgc.nEvents 			= config->getInt( "nEvents", pgc.nEvents );
	
	pgc.oneBetaSigma 		= config->getDouble( "tof.oneBetaSigma", pgc.oneBetaSigma );
	pgc.nTofSamples 		= config->getDouble( "tof.nSamples", pgc.nTofSamples );
	
	pgc.paddingScaleDedx 	= config->getDouble( "paddingScale:dedx", pgc.paddingScaleDedx );
	pgc.paddingScaleTof 	= config->getDouble( "paddingScale:tof", pgc.paddingScaleTof );
	pgc.paddingDedx 		= config->getDouble( "padding:dedx", pgc.paddingDedx );
	pgc.paddingTof 			= config->getDouble( "padding:tof", pgc.paddingTof );

	pgc.pCutLow 			= config->getDouble( "binning.pRange:low", pgc.pCutLow );
	pgc.pCutHigh 			= config->getDouble( "binning.pRange:high", pgc.pCutHigh );
	pgc.pBinWidth 			= config->getDouble( "binning.dedxTofBinWidth:p", pgc.pBinWidth );
	pgc.tofBinWidth 		= config->getDouble( "binning.dedxTofBinWidth:tof", pgc.tofBinWidth );
	pgc.dedxBinWidth 		= config->getDouble( "binning.dedxTofBinWidth:dedx", pgc.dedxBinWidth );

	pgc.tofCutLow 			= config->getDouble( "binning.tof:low", pgc.tofCutLow );
	pgc.tofCutHigh 			= config->getDouble( "binning.tof:high", pgc.tofCutHigh );
	pgc.tofVsPBinWidth  	= config->getDouble( "binning.tof:binWidth", pgc.tofVsPBinWidth );

	pgc.dedxCutLow 			= config->getDouble( "binning.dedx:low", pgc.dedxCutLow );
	pgc.dedxCutHigh 		= config->getDouble( "binning.dedx:high", pgc.dedxCutHigh );
	pgc.dedxVsPBinWidth  	= config->getDouble( "binning.dedx:binWidth", pgc.dedxVsPBinWidth );

	pgc.centerSpecies 		= config->getDouble( "centering:species", pgc.centerSpecies );
	pgc.centerSigmaDedx 	= config->getDouble( "centering:sigmaDedx", pgc.centerSigmaDedx );
	pgc.centerSigmaTof 		= config->getDouble( "centering:sigmaTof", pgc.centerSigmaTof );

	pgc.smearP 				= config->getDouble( "detectorEffects.smearP", pgc.smearP );
	pgc.tofMismatch 		= config->getDouble( "detectorEffects.tofMismatch", pgc.tofMismatch );
	pgc.phiMismatch 		= config->getDouble( "detectorEffects.phiMismatch", pgc.phiMismatch );
	pgc.pMismatch 			= config->getDouble( "detectorEffects.pMismatch", pgc.pMismatch );


	return pgc;

}
















