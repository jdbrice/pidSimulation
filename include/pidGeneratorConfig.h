#ifndef PID_GENERATOR_CONFIG
#define PID_GENERATOR_CONFIG 

#include "string"

class pidGeneratorConfig
{
public:
	pidGeneratorConfig(){}
	~pidGeneratorConfig(){}


	unsigned long long nEvents = 1e6;
	/**
	 * For the dedx vs. tof centered histos
	 */
	double pCutLow = 0.2;
	double pCutHigh = 4.0;
	double pBinWidth = 0.05;

	double paddingScaleTof = 0.1;
	double paddingScaleDedx = 0.1;
	double paddingTof = 5.0;
	double paddingDedx = 5.0;

	double dedxBinWidth = 0.25;
	double tofBinWidth = 0.25;


	/**
	 * For the tof v. P and the dedx v. P.
	 */
	double tofCutLow = 0.5; 
	double tofCutHigh = 4.0;
	double tofVsPBinWidth = 0.005;

	double dedxCutLow = 0.0;
	double dedxCutHigh = 0.014;
	double dedxVsPBinWidth = 0.001;


	string outFilename = "pidSimEvents.root";

	/**
	 * Tof Generator
	 */
	double oneBetaSigma = 0.012;
	int nTofSamples = 10;

	/**
	 * For centering coords
	 */
	int centerSpecies = 1;			// k
	double centerSigmaDedx = 0.06;	// k
	double centerSigmaTof = 0.012;	// k

	/**
	 * Detector Effects
	 */
	double smearP = 0.02;
	double tofMismatch = 0.08;
	double phiMismatch = 0.03;
	double pMismatch = 0.05;

	
};


#endif