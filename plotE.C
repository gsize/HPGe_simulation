
TH1D *h_init ;
TH1D *h_edep ;

TH1D *h0;
TH1D *h1; 
Double_t fun_FWHM( Double_t *energy,Double_t *par)
{
	Double_t x0=energy[0];
	return (par[0]+par[1]*TMath::Sqrt(x0 +par[2]* x0 * x0));
}

double eff_fun(double *x,double *par)
{
	double eff=0;
	for(int i=1;i<7;i++)
	{
		eff += par[i-1]*TMath::Power(x[0],2-i);
	}
	return (TMath::Exp(eff));
}

void plotE(TString file_name="HPGe_data")
{
	read_data(file_name);

	int num = 23;
	Double_t energy_0[]={
	0.053, 0.06202  ,  0.08602 , 0.09802 , 
 0.121781 ,0.2446981, 0.295941 , 0.344281 ,
 0.3677891, 0.4111161, 0.4439651, 0.563991 ,
 0.688671 , 0.78891  , 0.867371 , 0.9640791,
 1.085871 , 1.08971  , 1.1120691, 1.2129481,
 1.299141 , 1.4081   , 1.52811   
	};
	Double_t data_0[]={ 4. , 4.,4.  ,4.,            
28.668 ,7.61   ,0.448  ,26.558 ,  
0.8618 ,  2.2370 , 3.15760 , 0.491 ,   
0.859  ,  12.96   , 4.258   ,14.65   , 
10.238  , 1.729   ,13.69   ,1.426   , 
1.625   , 21.069  , 0.28 };
int num1=h_init->GetNbinsX();
cout<<num1<<endl;
double *data_init= new double[num];
double *data_edep= new double[num];
int j=0;
	for(Int_t i=0; i<num1; i++)
	{
		if(h_init->GetBinContent(i+1)>5)
		{
			double ntmp=h_init->GetBinContent(i+1);
			data_init[j]=ntmp;
			j++;
	}
	}
	for(int i=0;i<num;i++){
Int_t btmp=h_edep->FindBin(energy_0[i]);
data_edep[i]=get_area(btmp,h_edep);
printf("%6.3lf\t%8.2lf\t%8.2lf\t%8.5lf\n",energy_0[i],data_init[i],data_edep[i],data_edep[i]/data_init[i]);
	}
	TF1 *fun_eff=  new TF1("fun_eff",eff_fun,0.050,1.6,6);
	fun_eff->SetParameters(0.109574,-7.255860, 1.839348,  -0.462271, 0.060579, -0.003032);
	
TGraph *g_eff = new TGraph(num);
g_eff->SetTitle("data_eff");
for(int i=0; i<num;i++)
{
	g_eff->SetPoint(i,energy_0[i],data_edep[i]/data_init[i]);
}
	TCanvas* c4 = new TCanvas("ce", "  ");
	g_eff->Fit("fun_eff","R+");
	g_eff->Draw("A*");
	gPad->SetLogy(1);
/*
TGraph *g_0 = new TGraph(num);
g_0->SetTitle("data_0");
TGraph *g_init = new TGraph(num);
g_init->SetTitle("data_init");
double d0_max =TMath::MaxElement(num,data_0);
double dinit_max =TMath::MaxElement(num,data_init);
for(int i=0; i<num;i++)
{
	g_0->SetPoint(i,energy_0[i],100.0*data_0[i]/d0_max);
	g_init->SetPoint(i,energy_0[i],100.0*data_init[i]/dinit_max);
}
	TCanvas* c3 = new TCanvas("c3", "  ");
	c3->Divide(1,2);
	c3->cd(1);
	g_0->Draw();
	c3->cd(2);
	g_init->Draw();
*/	
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
	//f->Close();
}

double get_area(int bin, TH1D *h)
{
   double sum=0;
   Int_t dn=8;
   double bg=0;
   Int_t k=0;
	for(int i=-dn;i<dn;i++)
	{
		sum += h->GetBinContent(bin+i);
		if(i>3 ||i<-3)
		{
		bg += h->GetBinContent(bin+i);
		k++;
		}	
	}
	return (sum-(1.*dn*bg/k));
}

Int_t read_data(TString file_name)
{
	//gROOT->Reset();
TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
   dir.ReplaceAll("plotE.C","");
   dir.ReplaceAll("/./","/");
/*
	TFile *f1= TFile::Open(Form("%sbuild/%s_t0.root",dir.Data(),file_name.Data()));
	if(!f1->IsOpen())
	{
		cout<<file_name<<"isn't opend"<<endl;
		return ;
	}
	TTree *tr; f1->GetObject("HPGe_detector;1",tr);
	Int_t nentries = (Int_t)tr->GetEntries();
	Double_t edep_init;
	Double_t edep;
	//cout<<"tree:"<<nentries<<endl;
	tr->SetBranchAddress("Edep_init",&edep_init);
	tr->SetBranchAddress("Edep",&edep);
	h0 = new TH1D("edep_i","edep_0",8192,0,2.);
	h1 = new TH1D("edep_1","edep_1",8192,0,2.);
	for(Int_t i=0; i<nentries; i++)
	{
		tr->GetEntry(i);
		h0->Fill(edep_init);
		h1->Fill(edep);
	}
	TCanvas* c2 = new TCanvas("c2", "  ");
	c2->Divide(1,2);
	c2->cd(1);
	h0->Draw("HIST");
	c2->cd(2);
	h1->Draw();
*/

	TFile *f2= TFile::Open(Form("%sbuild/%s.root",dir.Data(),file_name.Data()));
  f2->GetObject("1;1",h_init);
  f2->GetObject("2;1",h_edep);
	TCanvas* c1 = new TCanvas("c1", "  ");
	c1->Divide(1,2);
	c1->cd(1);
	h_init->SetTitle("Source of gamma spectrum");
	h_init->GetXaxis()->SetTitle("Energy/MeV");
	h_init->GetXaxis()->CenterTitle();
	h_init->GetYaxis()->SetTitle("Count");
	h_init->GetYaxis()->CenterTitle();
	h_init->Draw("HIST");
	gPad->SetLogy(1);
	c1->cd(2);
	h_edep->SetTitle("Respond");
	h_edep->GetXaxis()->SetTitle("Energy/MeV");
	h_edep->GetXaxis()->CenterTitle();
	h_edep->GetYaxis()->SetTitle("Count");
	h_edep->GetYaxis()->CenterTitle();
	h_edep->Draw("HIST");
	gPad->SetLogy(1);
}
    
        
        
