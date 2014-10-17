#include <vector>
#include <iostream>
#include <stdio.h>

#include "TPaveText.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TString.h"
#include "TF1.h"

Double_t fun_unfold_gause(Double_t energy);
Double_t fun_FWHM( Double_t *energy,Double_t *par);
double eff_fun(double *x,double *par);

class DataAnalysis {
	public:
		DataAnalysis();
		virtual ~DataAnalysis();

		int ReadFile(const std::vector<TString> fileList);
		void ReadMacfile(TString &macFile);
		void PlotFWHM();
		TGraphErrors* PlotEffExperiment();
		TGraphErrors *PlotEffMC(int i);
		void PlotAllEfficiency();
		void PlotSpectra(int i);
		void PlotAllSpectra();

	private:
		double GetRateOfPeakComputom(TH1D *hist);
		double eff_fun(double *x,double *par);
		double GetArea(TH1D *h, int i);
		double GetArea( TH1D *h, double e);
		void AnalyzeSpectra(TH1D *h,std::vector<double> &peakAddr,std::vector<double> &peakArea);

	private:
		std::vector<TString> fileList;
		TObjArray *sourceList;
		TObjArray *HPGeList;
};

DataAnalysis::DataAnalysis()
{
	sourceList = new TObjArray();
	HPGeList = new TObjArray();
}
DataAnalysis::~DataAnalysis()
{
	//Cleanup
	delete sourceList;
	delete HPGeList;
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

void DataAnalysis::PlotFWHM()
{
	std::vector<double> en;
	std::vector<double> chn;
	std::vector<double> fwhm;

	TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
	dir.ReplaceAll("dataAnalysis.C","");
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

double DataAnalysis::GetRateOfPeakComputom(TH1D *hist)
{
	Int_t binmin=0,binmax=0;

	binmin=hist->FindFixBin(1.040);
	binmax=hist->FindFixBin(1.096);
	printf("d_bin=%d - %d = %d\n",binmax,binmin,binmax-binmin);
	Double_t sum=hist->Integral(binmin,binmax);
	return sum/(binmax-binmin);
}

void DataAnalysis::ReadMacfile(TString &macFile)
{

}


double DataAnalysis::GetArea( TH1D *h,int bin)
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
double DataAnalysis::GetArea( TH1D *h, double energy)
{
	int bin = h->FindBin(energy);
	return (GetArea(h, bin));
}

void DataAnalysis::AnalyzeSpectra(TH1D *h,std::vector<double> &peakAddr,std::vector<double> &peakArea)
{

	TSpectrum *sp = new TSpectrum();
	Int_t nfound = sp->Search(h,2,"nodraw");
	printf("Found %d peaks to fit\n",nfound);

	for(int i=0; i < nfound;i++)
	{
		peakAddr.push_back((double)( (sp->GetPositionX())[i]));
		double area = 0;
		area = GetArea(h,(double)( (sp->GetPositionX())[i]));
		peakArea.push_back(area);
	}

	delete sp;
}

TGraphErrors* DataAnalysis::PlotEffMC(int index )
{
	TH1D *hSource =(TH1D* )sourceList->At(index);
	TH1D *hHPGe =(TH1D* )HPGeList->At(index);

	std::vector<double> speakAddr;
	std::vector<double> speakArea;
	AnalyzeSpectra(hSource,speakAddr,speakArea);

	std::vector<double> peakAddr;
	std::vector<double> peakArea;
	AnalyzeSpectra(hHPGe,peakAddr,peakArea);

	TF1 *fun_eff=  new TF1("fun_eff",eff_fun,0.039,1.6,6);
	fun_eff->SetParameters(-0.552,-5.687, 0.434, -0.0404, 0.0013, -0.00003);

	TGraphErrors *g_eff = new TGraphErrors(peakAddr.size());

	for(int i=0; i<speakAddr.size();i++)
	{
		int k=0;
		double shield = 0.0015;
		while(k<peakAddr.size())
		{
			if(TMath::Abs(peakAddr[k] - speakAddr[i]) < shield )
			{
				g_eff->SetPoint(i,peakAddr[k],peakArea[k]/speakArea[i]);
				g_eff->SetPointError(i,0,TMath::Sqrt(1./peakArea[k]+1./speakArea[i])*(peakArea[k]/speakArea[i]));
				double dif=( (g_eff->GetY())[i] - fun_eff->Eval(peakAddr[k]) )/fun_eff->Eval(peakAddr[k]) ;
				printf("%d\t%6.3lf\t%8.2lf\t%8.2lf\t%8.5lf\t%8.2lf\n",i,peakAddr[k],speakArea[i],peakArea[k],(g_eff->GetY())[i], 100*dif);
				break;
			}
			k++;
		}

	}
	g_eff->Fit("fun_eff","R+");
	g_eff->SetMarkerStyle(20 + index);
	g_eff->SetLineColor(2+index);
	g_eff->SetTitle(fileList[index].Data());
	g_eff->GetXaxis()->SetTitle("Energy/MeV");
	g_eff->GetXaxis()->CenterTitle();
	g_eff->GetYaxis()->SetTitle("Efficiency");
	g_eff->GetYaxis()->CenterTitle();

	return g_eff;
}

void  DataAnalysis::PlotSpectra(int i)
{
	int entries = sourceList->GetEntries();
	if(entries < 1 || i >= entries) return;

	TH1D *hSource =(TH1D* )sourceList->At(i);
	TH1D *hHPGe =(TH1D* )HPGeList->At(i);

	TCanvas* c1 = new TCanvas(fileList[i].Data(), fileList[i].Data());
	c1->Divide(1,2);
	c1->cd(1);
	hSource->SetTitle("gamma Source");
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

void  DataAnalysis::PlotAllSpectra()
{
	for(int i=0;i<sourceList->GetEntries();i++)
	{
		PlotSpectra(i);
	}
}

TGraphErrors* DataAnalysis::PlotEffExperiment()
{
	TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
	dir.ReplaceAll("dataAnalysis.C","");
	dir.ReplaceAll("/./","/");
ifstream inf;
	inf.open(Form("%seff_ba133_08_20130502.Txt",dir.Data()));
	if(!(inf.is_open()))
	{
		printf("open eff file failed!\n");
		return 0;
	}

	int i=0; 
	std::vector<double> xList;
	std::vector<double> yList;
	std::vector<double> yErrorsList;
	string str_tmp;
	while(getline(inf,str_tmp)) {
		i++;
		double x,y,yFit,yError,yDelta; 
		if(i>5){
			sscanf(str_tmp.c_str(),"%lf%lf%lf%lf",&x, &y, &yFit, &yError);
			xList.push_back(x*0.001);
			yList.push_back(y);
			yDelta =TMath::Abs( y * yError *0.01);
			//cout<<yError<<"\t"<<yDelta<<endl;
			yErrorsList.push_back(yDelta);
		}
	}
	inf.close();

	TGraphErrors *gr = new TGraphErrors((Int_t)xList.size(),&xList[0],&yList[0],0,&yErrorsList[0]);
	TF1 *fun_fit = new TF1("eff_fun",eff_fun,0.039,1.5,6);
	fun_fit->SetParameters(-0.552,-5.687, 0.434, -0.0404, 0.0013, -0.00003);
	gr->Fit(fun_fit,"R+");

	return gr;
}
void DataAnalysis::PlotAllEfficiency()
{
	TMultiGraph *mg = new TMultiGraph();
	TLegend *leg = new TLegend(0.1,0.8,0.3,0.9);

	for(int i=0;i<sourceList->GetEntries();i++)
	{
		TGraphErrors *gr = 0;
		gr = PlotEffMC( i);
		mg->Add(gr);
		leg->AddEntry(gr,fileList[i].Data(),"lep");
	}

	TCanvas* c4 = new TCanvas("effGr", "  ");

	if(1)
	{
		TGraphErrors *gr = 0;
		gr= PlotEffExperiment();
		mg->Add(gr);
	}
	mg->Draw("AP");
	leg->Draw();
	gPad->SetLogy(1);
	gPad->SetGridy(1);
	gPad->SetGridx(1);
}

int  DataAnalysis::ReadFile(const std::vector<TString> fList)
{
	TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
	dir.ReplaceAll("dataAnalysis.C","");
	dir.ReplaceAll("/./","/");

	for(int i=0;i<fList.size();i++) 
	{
		TString fname;
		fileList.push_back(fList[i]);
		fname = fileList[i];

		TFile *f2= TFile::Open(Form("%sdata/20141016/%s.root",dir.Data(),fname.Data()));
		if(!(f2->IsOpen())){
			cout<<"file: "<<fname<<" isn't opened!"<<endl;
			return 0;
		}
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
	return 1;
}

void dataAnalysis()
{
	DataAnalysis *da = new DataAnalysis();
	std::vector<TString> fileList;
	TString fileName;
/*	fileName = "output_point_80mm";
	fileList.push_back(fileName);
	fileName = "output_point_85mm";
fileList.push_back(fileName);
fileName = "output_point_90mm";
fileList.push_back(fileName);
fileName = "output_point_95mm";
fileList.push_back(fileName);
fileName = "output_point_100mm";
fileList.push_back(fileName);
*/
fileName = "dl07";
	fileList.push_back(fileName);
	fileName = "dl09";
fileList.push_back(fileName);
fileName = "dl12";
fileList.push_back(fileName);
fileName = "dl16";
fileList.push_back(fileName);
fileName = "dl19";
fileList.push_back(fileName);
	/*
	  fileName = "output_plane_circle_5mm";
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
	da->ReadFile(fileList);
	da->PlotAllEfficiency();
	if(0)
		da->PlotAllSpectra();
	if(0)
		da->PlotEffExperiment();

	if(0)
		da->PlotFWHM();
}

