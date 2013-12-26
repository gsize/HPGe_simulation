
Double_t fun_unfold_gause(Double_t energy)
{
        Double_t par[3];
        Double_t ene_tmp=energy*1000.0 ;
        
        par[0]=0.54;par[1]=0.022;par[2]=0.00158;
        return 0.001*(ene_tmp + fun_FWHM(&ene_tmp,par) * gRandom->Gaus()/2.355);
}
Double_t fun_FWHM( Double_t *energy,Double_t *par)
{
	Double_t x0=energy[0];
    return (par[0]+par[1]*TMath::Sqrt(x0 +par[2]* x0 * x0));
}

void plot_FWHM()
{
	Int_t n=18;
	Double_t *en= new Double_t[n];
	Double_t *chn= new Double_t[n];
	Double_t *fwhm= new Double_t[n];

	string dir("calb_energy.Ent");
	ifstream in;
	in.open(dir.c_str());
	if (!in || in.bad())
	{
		printf("bad file!\n");
		return 0;
	}
	for(int i=0;i<n;i++)
	{
		in>>en[i];
		in>>chn[i];
		in>>fwhm[i];
	}
	in.close();
	
    TCanvas* c2 = new TCanvas("Plot FWHM", "FWHM");
    c2->Divide(1,2);
	TGraph *energy_cal=new TGraph(n,en,chn);
	energy_cal->SetTitle("Energy calibration");
	energy_cal->Fit("pol2");
	energy_cal->GetXaxis()->CenterTitle(true);
    energy_cal->GetYaxis()->CenterTitle(true);
	energy_cal->GetHistogram()->GetXaxis()->SetTitle("Energy/keV");
	energy_cal->GetHistogram()->GetYaxis()->SetTitle("Channal");
	TF1 *fun_energy_calb= energy_cal->GetFunction("pol2");
	c2->cd(1);
	energy_cal->Draw("AC*");
	
	for(int i=0;i<n;i++)
	{
		fwhm[i] /= fun_energy_calb->Derivative(en[i]);
		//cout<<fwhm[i]<<endl;
	}
	TGraph *FWHM_cal=new TGraph(n,en,fwhm);
	FWHM_cal->SetTitle("FWHM calibration");
	FWHM_cal->GetXaxis()->CenterTitle(true);
    FWHM_cal->GetYaxis()->CenterTitle(true);
	FWHM_cal->GetHistogram()->GetXaxis()->SetTitle("Energy/keV");
	FWHM_cal->GetHistogram()->GetYaxis()->SetTitle("FWHM/keV");
	TF1 *fun_FWHM_calb= new TF1("fun_FWHM",fun_FWHM,0.,2.,3);
	fun_FWHM_calb->SetParameters(0.54027,0.0219425,0.00158056);
	FWHM_cal->Fit("fun_FWHM");
	c2->cd(2);
	FWHM_cal->Draw("A*");
}

double GetRateOfPeakComputom(TH1D *hist)
{
Int_t binmin=0,binmax=0;

binmin=hist->FindFixBin(1.040);
binmax=hist->FindFixBin(1.096);
printf("d_bin=%d - %d = %d\n",binmax,binmin,binmax-binmin);
Double_t sum=hist->Integral(binmin,binmax);
return sum/(binmax-binmin);
}

void plotE_L(TString file_name="test_gamma.root"){
    gROOT->Reset();

    // Draw histos filled by Geant4 simulation
    TCanvas* c1 = new TCanvas("c1", "  ");
    TString tmp_file(file_name);
    Int_t bins = 4096;
    TH1D *h1 = new TH1D("1.33MeV Gamma Spectrum","1.33MeV Gamma Spectrum",bins,0.,1.5);
    TH1D *h2 = new TH1D("e1.33MeV Gamma length","e1.33MeV Gamma length",bins,0.,250.);
    TH2D *h3 = new TH2D("e1.33MeV Gamma","e1.33MeV Gamma",bins,0.,1.5,bins,0.,250.);
    TFile f(tmp_file.Data());

    if(f.IsOpen())
    {
        TTree *t1;
        f.GetObject("ntuple/test",t1);
        Int_t entrise = t1->GetEntries();

        Double_t e_e;
        Double_t l_e;
        Double_t energy_tmp=0.0;
        Int_t number=0;
        t1->SetBranchAddress("e_g",&e_e);
        t1->SetBranchAddress("l_g",&l_e);

        for(int j=0; j<entrise; j++)
        {
            t1->GetEntry(j);
            if(e_e>1.E-5)
            {
                if(energy_tmp<e_e)
                    energy_tmp=e_e;
                number++;
            }
            h1->Fill(fun_unfold_gause(e_e));
            h2->Fill(l_e);
            h3->Fill(e_e,l_e);
        }
    }

    tmp_file.Clear();
    f.Close();

  TH1D* hist1;
  f.GetObject("histo/1;1",hist1);
  
    c1->Divide(1,3);
    c1->cd(1);
    h1->SetTitle("1.332MeV gamma spectrum");
    h1->GetXaxis()->CenterTitle(true);
    h1->GetYaxis()->CenterTitle(true);
    h1->GetXaxis()->SetTitle("energy /MeV");
    h1->GetYaxis()->SetTitle("Count");
    gPad->SetLogy(1);
    h1->Draw();
    double xxx=GetRateOfPeakComputom(h1);
    double kkk=0;
	kkk=h1->GetBinContent(h1->FindFixBin(1.332));
    printf("coputom:%lf\tPeak:%lf\tPeakComputom:%lf\n",xxx,kkk,kkk/xxx);

    c1->cd(2);
    h2->SetTitle("1.332MeV gamma length");
    h2->GetXaxis()->CenterTitle(true);
    h2->GetYaxis()->CenterTitle(true);
    h2->GetXaxis()->SetTitle("length /mm");
    h2->GetYaxis()->SetTitle("Count");
    h2->Draw();

    c1->cd(3);
    h3->SetTitle("1.332MeV gamma energy_length histo");
    h3->GetXaxis()->CenterTitle(true);
    h3->GetYaxis()->CenterTitle(true);
    h3->GetXaxis()->SetTitle("energy /MeV");
    h3->GetYaxis()->SetTitle("length/mm");
    h3->Draw("FB");

    c1->cd();

plot_FWHM();
}
