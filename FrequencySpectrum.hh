void FrequencySpectrum(unsigned nSamples,double tMin, double tMax,Float_t * time, Float_t * voltage
		, TH1 *hist){
		
	const int range = 0; // extension of samples to be used beyond [tMin, tMax]
	double deltaT = (time[nSamples - 1] - time[0])/(double)(nSamples - 1); // sampling time interval
	double fCut = 0.5/deltaT; // cut frequency = 0.5 * sampling frequency
	int iSmin = floor(tMin/deltaT) - range; // first sample to use
	int iSmax = ceil(tMax/deltaT) + range; // last sample to use
	iSmin = max(iSmin,0); // check low limit
	iSmax = min(iSmax, (int)nSamples - 1); // check high limit
	int iS0 = (iSmin + iSmax)/2;
	
	int nBins = hist->GetNbinsX();
	double f1 = hist->GetXaxis()->GetXmin();
	double f2 = hist->GetXaxis()->GetXmax();
		
	for(int iBin = 1; iBin <= nBins; ++iBin){
		double freq = (hist->GetXaxis()->GetBinLowEdge(iBin) + hist->GetXaxis()->GetBinLowEdge(iBin))/2.;
		
		TComplex s(0.,0.); // Fourier transform at freq
		TComplex I(0.,1.); // i
		
		for(int iS = iSmin; iS <= iSmax; ++iS){
			//s += (double)voltage[iS]*TComplex::Exp(-I*(TComplex)((iS-iS0)*TMath::Pi()*freq/fCut));
			s += (double)voltage[iS]*TComplex::Exp(-I*((iS-iS0)*TMath::Pi()*freq/fCut));
		}
		
		double content = hist->GetBinContent(iBin) + s.Rho();
		hist->SetBinContent(iBin,content);
		//cout << "SetBinContent" << endl;
	}
}