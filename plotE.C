#include <vector>
#include <iostream>
#include <stdio.h>

#include "TPaveText.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TString.h"
#include "TF1.h"

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
	std::vector<double> en;
	std::vector<double> chn;
	std::vector<double> fwhm;

	TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
	dir.ReplaceAll("plotE.C","");
	dir.ReplaceAll("/./","/");

	ifstream in;
	in.open(Form("%scalb_energy.Ent",dir.Data()));
	while (1) {
		double v_en,v_chn,v_fwhm;
		in >>  v_en >> v_chn >> v_fwhm;
		if (!in.good()) break;
		en.push_back(v_en);
		chn.push_back(v_chn);
		fwhm.push_back(v_fwhm);
	}//while
	in.close();

	TCanvas* c2 = new TCanvas("Plot FWHM", "FWHM");
	c2->Divide(1,2);
	TGraph *energy_cal=new TGraph(en.size(),&en[0],&chn[0]);
	energy_cal->SetTitle("Energy calibration");
	energy_cal->Fit("pol2");
	energy_cal->GetXaxis()->CenterTitle(true);
	energy_cal->GetYaxis()->CenterTitle(true);
	energy_cal->GetHistogram()->GetXaxis()->SetTitle("Energy/keV");
	energy_cal->GetHistogram()->GetYaxis()->SetTitle("Channal");
	TF1 *fun_energy_calb= energy_cal->GetFunction("pol2");
	c2->cd(1);
	energy_cal->Draw("AC*");

	for(int i=0;i<en.size();i++)
	{
		fwhm[i] /= fun_energy_calb->Derivative(en[i]);
		//cout<<fwhm[i]<<endl;
	}
	TGraph *FWHM_cal=new TGraph(en.size(),&en[0],&fwhm[0]);
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

double eff_fun(double *x,double *par)
{
	double eff=0;
	for(int i=1;i<7;i++)
	{
		eff += par[i-1]*TMath::Power(x[0],2-i);
	}
	return (TMath::Exp(eff));
}

void plot_eff_std()
{
	TF1 *fun_eff_0=  new TF1("fun_eff_0",eff_fun,0.039,1.6,6);
	TF1 *fun_eff_1=  new TF1("fun_eff_1",eff_fun,0.039,1.6,6);
	TF1 *fun_eff_2=  new TF1("fun_eff_2",eff_fun,0.039,1.6,6);
	fun_eff_0->SetParameters(-0.552,-5.687, 0.434, -0.0404, 0.0013, -0.00003);
	fun_eff_1->SetParameters(-0.452590,-5.901407, 0.539222 ,-0.059246, 0.002599 ,-0.000048);
	fun_eff_2->SetParameters(-0.439793, -5.822942, 0.488969, -0.048183, 0.001682 ,-0.000024);

	//TCanvas* c_eff = new TCanvas("Canvas_eff", "Canvas_eff");
	fun_eff_0->SetLineColor(kBlack);
	fun_eff_0->Draw("SAME");
	fun_eff_1->SetLineColor(kRed);
	fun_eff_1->Draw("SAME");
	fun_eff_2->SetLineColor(kBlue);
	fun_eff_2->Draw("SAME");
}

void reb(TF1 *fun1,TF1 *fun2,int bins, double xMin, double xMax)
{
	TCanvas *c1 = new TCanvas("reb","c1",800,1000);
	c1->Divide(1,2);
	c1->cd(1);

	fun1->Draw();
	fun2->Draw("same");

	TH1F *hfun1 = new TH1F("hfun1","fun1",bins,xMin,xMax);
	TH1F *hfun2 = new TH1F("hfun2","fun2",bins,xMin,xMax);
	hfun1->Eval(fun1);
	hfun2->Eval(fun2);

	c1->cd(2);
	hfun1->DrawCopy();
	hfun2->Draw("same");
	//to bypass the warning
	//hfun2->GetXaxis()->SetLimits(0,5);
	hfun1->Add(hfun2, -1);
	//hfun2->GetXaxis()->SetLimits(0,3);
	hfun1->SetLineColor(kBlue);
	hfun1->Draw("same");
}

TGraphErrors * plotEFF_exp(const TString &fname,TH1D *hSource,TH1D *hHPGe)
{
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
	int num1 = hSource->GetNbinsX();
	//cout<<num1<<endl;
	std::vector<double> xList;
	double *data_init= new double[num];
	double *data_edep= new double[num];
	int j=0;
	for(Int_t i=0; i<num1; i++)
	{
		if(hSource->GetBinContent(i+1)>5)
		{
			double ntmp=hSource->GetBinContent(i+1);
			data_init[j]=ntmp;
			j++;
		}
	}
	for(int i=0;i<num;i++){
		Int_t btmp=hHPGe->FindBin(energy_0[i]);
		data_edep[i]=get_area(btmp,hHPGe);
	}
	TF1 *fun_eff=  new TF1("fun_eff",eff_fun,0.039,1.6,6);
	fun_eff->SetParameters(-0.552,-5.687, 0.434, -0.0404, 0.0013, -0.00003);
	//TF1 *fun_eff1=  new TF1("fun_eff_0",eff_fun,0.051,1.6,6);
	//fun_eff1->SetParameters(-0.552,-5.687, 0.434, -0.0404, 0.0013, -0.00003);

	TGraphErrors *g_eff = new TGraphErrors(num);
	g_eff->SetTitle(fname.Data());
	for(int i=0; i<num;i++)
	{
		g_eff->SetPoint(i,energy_0[i],data_edep[i]/data_init[i]);
		g_eff->SetPointError(i,0,TMath::Sqrt(1./data_edep[i]+1./data_init[i])*(data_edep[i]/data_init[i]));
		printf("%d\t%6.3lf\t%8.lf\t%8.2lf\t%6.5lf\n",i,energy_0[i],data_init[i],data_edep[i],data_edep[i]/data_init[i]);
	}
	//TCanvas* c4 = new TCanvas("ce", "  ");
	g_eff->Fit("fun_eff","R+");
	//g_eff->SetMarkerStyle(20);
	g_eff->SetTitle(" ");
	g_eff->GetXaxis()->SetTitle("Energy/MeV");
	g_eff->GetXaxis()->CenterTitle();
	g_eff->GetYaxis()->SetTitle("Efficiency");
	g_eff->GetYaxis()->CenterTitle();
	/* 
	   fun_eff1->SetLineColor(kBlack);
	   fun_eff1->Draw("SAME");
	   */
	//	reb(g_eff->GetFunction("fun_eff"),fun_eff1, 819, 0.051,1.6 );
	return g_eff;
}

double get_area(int bin, TH1D *h)
{
	double sum=0;
	Int_t dn=10;
	Int_t dbg=5;
	double bg=0;
	Int_t k=0;
	for(int i=-dn;i<dn;i++)
	{
		sum += h->GetBinContent(bin+i);
		if(i>dbg ||i<-dbg)
		{
			bg += h->GetBinContent(bin+i);
			k++;
		}
	}
	return (sum-(1.*dn*bg/k));
}

void PlotSpectra(const TString &fileName,TH1D *hSource,TH1D *hHPGe)
{
	//gROOT->Reset();
	TCanvas* c1 = new TCanvas(fileName.Data(), fileName.Data());
	c1->Divide(1,2);
	c1->cd(1);
	hSource->SetTitle("Source of gamma spectrum");
	hSource->GetXaxis()->SetTitle("Energy/MeV");
	hSource->GetXaxis()->CenterTitle();
	hSource->GetYaxis()->SetTitle("Count");
	hSource->GetYaxis()->CenterTitle();
	hSource->Draw("HIST");
	gPad->SetLogy(1);
	gPad->SetGridy(1);
	gPad->SetGridx(1);
	c1->cd(2);
	hHPGe->SetTitle("Respond");
	hHPGe->GetXaxis()->SetTitle("Energy/MeV");
	hHPGe->GetXaxis()->CenterTitle();
	hHPGe->GetYaxis()->SetTitle("Count");
	hHPGe->GetYaxis()->CenterTitle();
	hHPGe->Draw("HIST");
	gPad->SetLogy(1);
	gPad->SetGridy(1);
	gPad->SetGridx(1);

	return ;
}

void PlotEfficiency(const std::vector<TString> &fileList,TObjArray *sourceList, TObjArray *HPGeList)
{
	TMultiGraph *mg = new TMultiGraph();
	TLegend *leg = new TLegend(0.1,0.8,0.3,0.9);
	for(int i=0;i<sourceList->GetEntries();i++)
	{
		TH1D *s1 =(TH1D* )sourceList->At(i);
		TH1D *s2 =(TH1D* )HPGeList->At(i);
		PlotSpectra(fileList[i], s1, s2 );
		TGraphErrors *gr = 0;
		gr = plotEFF_exp(fileList[i],s1,s2);
		gr->SetMarkerStyle(20+i);
		gr->SetLineColor(2+i);
		//gr->Set(20+i);
		mg->Add(gr);
		leg->AddEntry(gr,fileList[i].Data(),"lep");
	}
	TPaveText *pt = new TPaveText(0.6,0.7,0.98,0.98,"brNDC");
	pt->SetFillColor(18);
	pt->SetTextAlign(12);
	pt->SetTextFont(12);
	pt->AddText("Polynomial Fit");
	pt->AddText(" #varepsilon = exp(#sum^{6}_{i=1} a_{i} E^{2-i})");

	TCanvas* c4 = new TCanvas("effGr", "  ");
	mg->Draw("ALP");
	//pt->Draw();
	if(1)
		plot_eff_std();
	
	leg->Draw();
	gPad->SetLogy(1);
	gPad->SetGridy(1);
	gPad->SetGridx(1);

}

int ReadFile(std::vector<TString> &fileList,TObjArray *sourceList, TObjArray *HPGeList)
{
	TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
	dir.ReplaceAll("plotE.C","");
	dir.ReplaceAll("/./","/");
	for(int i=0;i<fileList.size();i++) 
	{
		TString fname;
		fname = fileList[i];
		//TFile *f2= TFile::Open(Form("%sdata/entries10MH1/%s.root",dir.Data(),fname.Data()));
		TFile *f2= TFile::Open(Form("%sdata/livemore/%s.root",dir.Data(),fname.Data()));
		if(!(f2->IsOpen())){
			cout<<"file: "<<file_name<<" isn't opened!"<<endl;
			return 0;
		}
		//f2->GetObject("source;1",h_init);
		//f2->GetObject("HPGe;1",h_edep);
		TDirectory* dire = (TDirectory*)f2->Get("histo");
		TH1D *hSource = (TH1D*)dire->Get("source"); 
		hSource->SetDirectory(0);
		TH1D * hHPGe = (TH1D*)dire->Get("HPGe"); 
		hHPGe->SetDirectory(0);
		sourceList->Add(hSource);
		HPGeList->Add(hHPGe);
		cout<<"Read file "<<fname<<" success!"<<endl;
		f2->Close();
	}
}

void plotE()
{
	std::vector<TString> fileList;
	TObjArray *sourceList = new TObjArray();
	TObjArray *HPGeList = new TObjArray();

	TString fileName;
	fileName = "output_point_80mm";
	fileList.push_back(fileName);
	/*fileName = "output_plane_circle_5mm";
	fileList.push_back(fileName);
	fileName = "output_plane_circle_10mm";
	fileList.push_back(fileName);
	fileName = "output_plane_circle_15mm";
	fileList.push_back(fileName);
	fileName = "output_plane_circle_20mm";
	fileList.push_back(fileName);
	fileName = "output_plane_circle_30mm";
	fileList.push_back(fileName);
	fileName = "output_plane_circle_50mm";
	fileList.push_back(fileName);
*/
	ReadFile(fileList,sourceList,HPGeList);
	PlotEfficiency(fileList,sourceList,HPGeList);
	
	

	if(0)
		plot_FWHM();
}

