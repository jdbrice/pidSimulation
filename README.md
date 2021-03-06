pidSimulation
=============

Project for simulating events for use in combined tof + dedx pid studies.
Beased on a methods provided by Evan Sangaline. 

1. Checkout project into working directoy
2. cd pidSimulation/bin
3. ./fullbuild
4. This will produce the executable "simulate"
5. Run the simulation with the default values with "./simulate" or run with a config file "./simulate config.xml"



#Configuration

A default configuration is hard coded, but xml configuration files can be used to override any value. 



```xml

<config>

	<output>simukatedPidEvents.root</output>
	<nEvents>50000000</nEvents>

	<tof>
		<oneBetaSigma>0.012</oneBetaSigma>
		<nSamples>10</nSamples>
	</tof>

	
	<paddingScale tof="0.1" dedx="0.1" />
	<padding tof="5.0" dedx="5.0" />
	
	<binning>		
		<!-- DeDx vs. Tof in slices of P -->
		<pRange low="0.2" high="4.0" />
		<dedxTofBinWidth p="0.05" tof="0.10" dedx="0.10" />

		<!-- tof vs. P. Uses p binning from dedxTof-->
		<tof low="0.5" high="4.0" binWidth="0.001"/>
		<!-- dedx vs. P. Uses p binning from dedxTof-->
		<dedx low="0.0" high="0.014" binWidth="0.0001" />

	</binning>

	<centering species="1" sigmaDedx="0.06" sigmaTof="0.012" />


	<detectorEffects>
		
		<smearP>0.02</smearP>
		<tofMismatch>0.08</tofMismatch>
		<phiMismatch>0.03</phiMismatch>
		<pMismatch>0.05</pMismatch>

	</detectorEffects>


</config>

```
