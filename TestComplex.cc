#include <iostream>
#include <TComplex.h>


void TestComplex(){

	TComplex z(0.,TMath::Pi()/2.);
	
	double x = z.Rho();

	cout << TComplex::Exp(z) << endl;
	cout << x << endl;
	
}