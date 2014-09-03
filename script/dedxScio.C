

#include "../include/dedxGenerator.h"
#include "iostream"
using namespace std;

void dedxScio() {

	dedxGenerator *dedx = new dedxGenerator(  );

	cout << " dedx = " << dedx->random( 1, 0.134 ) << endl;

	dedxGenerator *dedx1 = new dedxGenerator(  );

	cout << " dedx = " << dedx1->random( 1, 0.134 ) << endl;

	dedxGenerator *dedx2 = new dedxGenerator(  );

	cout << " dedx = " << dedx2->random( 1, 0.134 ) << endl;

	
}