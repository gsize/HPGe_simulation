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
double eff_fun1(double *x,double *par);
Double_t correct(Double_t *x, Double_t *par);

class DataAnalysis {
	public:
		DataAnalysis();
		virtual ~DataAnalysis();

		int ReadFile();
		void ReadMacfile(TString &macFile);
		void PlotFWHM();
		TGraphErrors* PlotEffExperiment(int index);
		TGraphErrors *PlotEffMC(int i);
		TGraphErrors* PlotEffMCNP();
		void PlotAllEfficiency();
		void PlotSpectra(int i);
		void PlotAllSpectra();
		void SetDirName(TString d){dataDirName = d;};
		void SetDir(TString d){dataDir = d;};
		void ReadEffFile(TString filename,std::vector<double> &xList,std::vector<double> &yList,std::vector<double> &yErrorsList);
		void GetFileNameList(bool  flag=0);

	private:
		double GetRateOfPeakComputom(TH1D *hist);
		double GetBackground( TH1D *h, int bin,int winbin);
		double GetBackground( TH1D *h, double energy,int winbin);
		double GetArea( TH1D *h,  int bin,int winbin);
		double GetNetArea( TH1D *h,double energy,int winbin=20);
		void AnalyzeSpectra(TH1D *h,std::vector<double> &peakAddr,std::vector<double> &peakArea,std::vector<double> &bgArea);

		TF1* getFitFunction(TString funname);
		void PrintChi2(TGraph *gr,TF1 *f );

	private:
		std::vector<TString> fileList;
		std::vector<TString> effFileList;
		TObjArray *sourceList;
		TObjArray *HPGeList;
		TString dataDir;
		TString dataDirName;
		Bool_t flagEffFit;
};

DataAnalysis::DataAnalysis()
{
	sourceList = new TObjArray();
	HPGeList = new TObjArray();
	flagEffFit =0;

	TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
	dir.ReplaceAll("dataAnalysis.C","");
	dir.ReplaceAll("/./","/");
	dataDir = dir + "data/";

}
DataAnalysis::~DataAnalysis()
{
	//Cleanup
	delete sourceList;
	delete HPGeList;
}
TF1* DataAnalysis::getFitFunction(TString funname)
{
	TF1* fun_eff = 0;
	if(flagEffFit)
	{
		fun_eff= new TF1(funname.Data(),eff_fun,0.045,1.500,6);
		//fun_eff->SetParameters(-0.552,-5.687, 0.434, -0.0404, 0.0013, -0.00003);
		fun_eff->SetParameters(-0.349521, -6.056041, 0.605647, -0.068800, 0.003151, -0.000059);
	}
	else
	{
		fun_eff= new TF1(funname.Data(),eff_fun1,0.045,1.500,5);
		fun_eff->SetParameters(-5.869, -0.9255, -0.22389, -0.27728,-0.09987);
	}
	return fun_eff;
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
double eff_fun1(double *x, double *par)
{
	double eff=0.;
	for(int i=0;i<5;i++)
	{
		eff += par[i] * TMath::Power(TMath::Log(x[0]),i);
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

Double_t correct(Double_t *x, Double_t *par) {
	Double_t x0 = TMath::Tan(TMath::Pi()*x[0]/180);
	Double_t result = 0.5*(1. - 1/TMath::Sqrt(1+x0*x0));
	return result;
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

double DataAnalysis::GetArea( TH1D *h,int bin,int winbin)
{
	double area=0.;
	int halfWin = TMath::CeilNint(0.5*winbin);
	area = h->Integral(bin-halfWin ,bin+halfWin );
	return (area );
}
double DataAnalysis::GetNetArea( TH1D *h, double energy,int winbin)
{
	int bin = h->FindBin(energy);
	double bgArea=0.;
	double grossArea=0.;
	bgArea =this->GetBackground(h,bin,winbin);
	grossArea =this->GetArea(h,bin,winbin);
	return (grossArea -bgArea);
}

double DataAnalysis::GetBackground( TH1D *h, int bin,int winbin)
{	
	double bgArea=0.;
	TH1* bg = h->ShowBackground(20,"compton");
	bgArea = this->GetArea((TH1D*)bg,bin,winbin);
	return (bgArea );
}

double DataAnalysis::GetBackground( TH1D *h, double energy,int winbin)
{
	int bin = h->FindBin(energy);
	double bgArea =0.;
	bgArea = this->GetBackground(h,bin,winbin);
	return (bgArea);
}

void DataAnalysis::AnalyzeSpectra(TH1D *h,std::vector<double> &peakAddr,std::vector<double> &peakArea,std::vector<double> &bgArea)
{
	TSpectrum *sp = new TSpectrum();
	Int_t nfound = sp->Search(h,2,"nodraw");
	printf("%s found %d peaks to fit\n",h->GetTitle(),nfound);
	TH1* bg = h->ShowBackground(20,"compton");

	for(int i=0; i < nfound;i++)
	{
		peakAddr.push_back((double)( (sp->GetPositionX())[i]));
	}
	std::sort(peakAddr.begin(),peakAddr.end());
	for(int i=0; i < nfound;i++)
	{
		double area = 0;
		double bgr = 0;
		int winbin =10;
		area = GetNetArea(h,peakAddr[i],winbin);
		peakArea.push_back(area);
		bgr = GetBackground(h,peakAddr[i],winbin);
		bgArea.push_back(bgr);
	}

	delete sp;
}

void DataAnalysis::PrintChi2(TGraph *gr,TF1 *f )
{
	double chi2=0.;
	int n=gr->GetN();
	double* y= gr->GetY();
	double* x= gr->GetX();
	double chi2_1 =gr->Chisquare(f);
	double rms1= TMath::RMS(n,y);
	double rms2= gr->GetRMS(2);
	double R2 =1-TMath::Power( chi2/rms1,2);
	double R2_1 =1-TMath::Power( chi2/rms2,2);
	for(int i=0;i<n;i++)
	{
		chi2 +=TMath::Power(y[i] - f->Eval(x[i]),2);
	}
	printf("chi2:%lg\tChisquare:%lg\trms1:%lg\tR^2:%lg\trms2:%lg\tR^2:%lg\n",chi2,chi2_1,rms1,R2,rms2,R2_1);
	printf("TG->chi2:%lg\n",gr->Chisquare(f,"R"));
}
TGraphErrors* DataAnalysis::PlotEffMC(int index )
{
	TH1D *hSource =(TH1D* )sourceList->At(index);
	TH1D *hHPGe =(TH1D* )HPGeList->At(index);

	std::vector<double> speakAddr;
	std::vector<double> speakArea;
	std::vector<double> sbgArea;
	AnalyzeSpectra(hSource,speakAddr,speakArea,sbgArea);

	std::vector<double> peakAddr;
	std::vector<double> peakArea;
	std::vector<double> bgArea;
	AnalyzeSpectra(hHPGe,peakAddr,peakArea,bgArea);

	int lineColor=1;
	lineColor = 1+index;
	if(index > 8) lineColor++;
	TString titleName(fileList[index]);
	titleName.ReplaceAll(".root","");
	TString effName;
	effName = "eff_" +titleName;
	printf("%s\n",fileList[index].Data());
	TF1 *fun_eff=0;
	fun_eff = getFitFunction(effName);
	fun_eff->SetLineColor(lineColor);
	TGraphErrors *g_eff = new TGraphErrors(peakAddr.size());
	TF1 fun_correct("fun_correct",correct,0,90,0);
	double ang=30.;
	if(titleName.Contains("d11mm") ||titleName.Contains("d10mm")||titleName.Contains("d50mm")||titleName.Contains("90deg"))
		ang=90.;
	double gemCorrect=fun_correct.Eval(ang);
	if(titleName.Contains("NG") || titleName.Contains("ISO"))
		gemCorrect=1.0;
	for(int i=0,j=0; i<speakAddr.size();i++)
	{
		int k=0;
		double shield = 0.0015;
		while(k<peakAddr.size())
		{
			if(TMath::Abs(peakAddr[k] - speakAddr[i]) < shield )
			{
				g_eff->SetPoint(j,peakAddr[k],gemCorrect *peakArea[k]/speakArea[i]);
				double err = gemCorrect*TMath::Sqrt(1./peakArea[k]+1./speakArea[i])*(peakArea[k]/speakArea[i]);
				g_eff->SetPointError(j,0,err);
				double dif=( (g_eff->GetY())[i] - fun_eff->Eval(peakAddr[k]) )/fun_eff->Eval(peakAddr[k]) ;
				printf("%d\t%8.4lf\t%8.3lf\t%8.3lf\t%8.3lf\t%8.3lf\t%8.5lf\t%8.2lf\n",j,peakAddr[k],speakArea[i],sbgArea[i],peakArea[k],bgArea[i],(g_eff->GetY())[i], 100*dif);
				j++;
				break;
			}
			k++;
		}

	}
	g_eff->Fit(fun_eff,"R");
	PrintChi2(g_eff,fun_eff );
	g_eff->SetMarkerStyle(20 + index);
	g_eff->SetLineColor(lineColor);
	g_eff->SetTitle(titleName);
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
	hSource->ShowBackground(20,"compton same");
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
	hHPGe->ShowBackground(20,"compton same");
	gPad->SetLogy(1);
	gPad->SetGridy(1);
	gPad->SetGridx(1);
	c1->cd();
	return ;
}

void  DataAnalysis::PlotAllSpectra()
{
	for(int i=0;i<sourceList->GetEntries();i++)
	{
		PlotSpectra(i);
	}
}

TGraphErrors* DataAnalysis::PlotEffMCNP()
{
	ifstream inf;
	inf.open(Form("%seff_80mm_point_MCNP.txt",dataDir.Data()));
	if(!(inf.is_open()))
	{
		printf("open eff file failed!\n");
		return 0;
	}

	int i=0; 
	std::vector<double> xList;
	std::vector<double> yList;
	string str_tmp;
	while(getline(inf,str_tmp)) {
		i++;
		double x,yExp,yMC; 
		if(i>0){
			sscanf(str_tmp.c_str(),"%lf%lf%lf",&x, &yExp, &yMC);
			xList.push_back(x*0.001);
			yList.push_back(yMC);
		}
	}
	inf.close();

	TGraphErrors *gr = new TGraphErrors((Int_t)xList.size(),&xList[0],&yList[0]);
	TF1 *fun_fit = 0;
	TString effName("eff_funMC");
	printf("%s\n",effName.Data());
	fun_fit = getFitFunction(effName);
	gr->Fit(fun_fit,"R+");
	gr->SetTitle("eff MCNP");
	gr->GetXaxis()->SetTitle("Energy/MeV");
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->SetTitle("Eff");
	gr->GetYaxis()->CenterTitle();
	return gr;
}
TGraphErrors* DataAnalysis::PlotEffExperiment(int index)
{
	TString effFile;
	effFile = dataDir + effFileList[index];
	TString titleName(effFileList[index]);
	titleName.ReplaceAll(".Txt","");

	std::vector<double> xList;
	std::vector<double> yList;
	std::vector<double> yErrorsList;
	ReadEffFile(effFile ,xList,yList,yErrorsList);

	TGraphErrors *gr = new TGraphErrors((Int_t)xList.size(),&xList[0],&yList[0],0,&yErrorsList[0]);
	TF1 *fun_fit = 0;
	TString effName("eff_funExp");
	printf("%s\n",effName.Data());
	fun_fit = getFitFunction(effName);
	fun_fit->SetLineColor(1);
	gr->Fit(fun_fit,"R+");
	printf("Chisquare:%lg\n",fun_fit->GetChisquare());
	gr->SetMarkerStyle(5);
	gr->SetTitle(titleName.Data());

	return gr;
}
void DataAnalysis::PlotAllEfficiency()
{
	TMultiGraph *mg = new TMultiGraph();
	TLegend *leg = new TLegend(0.1,0.8,0.3,0.9);

	for(int i=0;i<sourceList->GetEntries();i++)
	{
	//	if(i==6|| i==8|| i==11|| i==10)
		//if(i==0|| i==3|| i==5|| i==6|| i==9)
	//	if( i==3|| i==5|| i==6)
//		if( i==8)
		//if(i==0) continue;
		{
			TGraphErrors *gr = 0;
			gr = PlotEffMC( i);
			mg->Add(gr);
			leg->AddEntry(gr,gr->GetTitle(),"lep");
		}
	}
	if(1)
	{
		for(int i=0; i< effFileList.size();i++)
		{
			//if(i == 1)
{
				TGraphErrors *gr = 0;
				gr= PlotEffExperiment(i);
				mg->Add(gr);
				leg->AddEntry(gr,gr->GetTitle(),"lep");
			}
		}
	}
	TCanvas* c4 = new TCanvas("effGr", "  ");
	mg->Draw("AP");
	mg->GetXaxis()->SetTitle("Energy/MeV");
	mg->GetXaxis()->CenterTitle();
	mg->GetYaxis()->SetTitle("Efficiency");
	mg->GetYaxis()->CenterTitle();
	leg->Draw();
	//gPad->SetLogy(1);
	gPad->SetGridy(1);
	gPad->SetGridx(1);
	gPad->Update();
}

int  DataAnalysis::ReadFile()
{

	for(int i=0;i<fileList.size();i++) 
	{
		TString fname;
		fname = fileList[i];

		TFile *f2= TFile::Open(Form("%s%s/%s",dataDir.Data(),dataDirName.Data(),fname.Data()));
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

void DataAnalysis::ReadEffFile(TString filename,std::vector<double> &xList,std::vector<double> &yList,std::vector<double> &yErrorsList)
{
	ifstream inf;
	inf.open(filename.Data());
	if(!(inf.is_open()))
	{
		printf("open file failed! %s\n",filename.Data());
		return;
	}

	cout<<"Read efficiency file "<<filename<<" success!"<<endl;
	int i=0; 
	string str_tmp;
	while(getline(inf,str_tmp)) {
		i++;
		double x,y,yFit,yError,yDelta; 
		if(i>5){
			sscanf(str_tmp.c_str(),"%lf%lf%lf%lf",&x, &y, &yFit, &yError);
			xList.push_back(x*0.001);
			yList.push_back(y);
			yDelta =TMath::Abs( y * yError *0.01);
			yErrorsList.push_back(yDelta);
		}
	}
	inf.close();
}

void DataAnalysis::GetFileNameList(bool  flag)
{
	TSystemDirectory dire;
	TString ext;
	if(flag)
	{
		TString dataPath;
		dataPath = dataDir + dataDirName;
		dire.SetDirectory(dataPath.Data());
		ext = ".root";
	}
	else
	{
		dire.SetDirectory(dataDir.Data());
		ext = ".Txt";
	}
	TList *files = dire.GetListOfFiles();
	if(files){
		TSystemFile *file;
		TString fname;
		TIter next(files);
		while((file=(TSystemFile*)next())){
			fname = file->GetName();
			if(!file->IsDirectory() && fname.EndsWith(ext.Data()) )
			{
				if(fname.BeginsWith("eff"))
				{
					effFileList.push_back(fname.Data());
					cout<<"Load spectra file "<<fname<<" success!"<<endl;
				}
				else
				{
					fileList.push_back(fname.Data());
					cout<<"Load efficiency file "<<fname<<" success!"<<endl;
				}
			}
		}
	}
}

void dataAnalysis(TString dirName="point")
{
	DataAnalysis *da = new DataAnalysis();
	da->SetDirName(dirName);
	da->GetFileNameList(1);
	da->GetFileNameList(0);
	da->ReadFile();
	da->PlotAllEfficiency();
	if(0)
		da->PlotAllSpectra();

	if(0)
		da->PlotFWHM();
}

