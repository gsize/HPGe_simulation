
Double_t fun_FWHM( Double_t *energy,Double_t *par)
{
	Double_t x0=energy[0];
	return (par[0]+par[1]*TMath::Sqrt(x0 +par[2]* x0 * x0));
}

void plotE(TString file_name("data.root")){
	gROOT->Reset();

	TFile *f =new  TFile(Form("build/%s",file_name.Data()));
	if(!f->IsOpen())
	{
		cout<<file_name<<"isn't opend"<<endl;
		return 0;
	}

	TCanvas* c1 = new TCanvas("c1", "  ");
	c1->Divide(1,2);
	TDirectory* dir = f->Get("histo");
	c1->cd(1);
	TH1D* hist1 = (TH1D*)dir->Get("1");      
	hist1->Draw("HIST");
	c1->cd(2);
	TH1D* hist2 = (TH1D*)dir->Get("2"); 
	hist2->Draw("HIST");
	/*
	   int nbinsx=hist1->GetNbinsX();
	   int nbinsy=nbinsx;
	   Float_t * source = new float[nbinsx];
	   Float_t ** response = new float *[nbinsy];
	   TF1 *fun_FWHM_calb= new TF1("fun_FWHM",fun_FWHM,0.,2.,3);
	   fun_FWHM_calb->SetParameters(0.54027,0.0219425,0.00158056);

	   TH2F *decon_unf_resp = new TH2F("decon_unf_resp","Root File",nbinsy,ymin,ymax,nbinsx,xmin,xmax);
	   for (Int_t i=0;i<nbinsy;i++) response[i]=new float[nbinsx];
	   for (Int_t i = 0; i < nbinsx; i++) source[i] = 1000. * (hist1->GetBinContent(i + 1));
	   for (Int_t i = 0; i < nbinsy; i++){
	   double energy = hist1->GetBinCenter(i) ;
	   double sig = fun_FWHM_calb->Eval(energy)/2.355;
	   for (Int_t j = 0; j< nbinsx; j++){
	   double energy1 = hist1->GetBinCenter(j) ;
	   response[i][j] = TMath::Gaus(energy1,energy,sig,kTRUE);
	   }
	   }
	   TSpectrum *s = new TSpectrum();
	   s->Unfolding(source,response,nbinsx,nbinsy,1000,1,1);
	   TH1D *d = new TH1D("Spectrum"," Gamma Spectrum",nbinsx,0.,2);
	   for (i = 0; i < nbinsy; i++) d->SetBinContent(i + 1,source[i]);
	   TCanvas* c2 = new TCanvas("c2", "  ");
	   d->Draw();

*/
}
