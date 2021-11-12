#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <RQ_OBJECT.h>
#include <stdlib.h>
#include <string>
#include <TButton.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TDecompSVD.h>
#include <TF1.h>
#include <TFile.h>
#include <TGButton.h>
#include <TGClient.h>
#include <TGFrame.h>
#include "TGLBoxPainter.h"
#include "TGLHistPainter.h"
#include <TGraph.h>
#include "TH1F.h"
#include "TH1.h"
#include <TH2F.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TNtuple.h>
#include <TObjArray.h>
#include <TPaveLabel.h>
#include <TPolyLine.h>
#include <TProfile.h>
#include <TRandom3.h>
#include "TRandom.h"
#include <TRootEmbeddedCanvas.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include "TTree.h"
#include <TVectorD.h>
#include <vector>

//#include "bbgrinch.C"
//#include "bbcal.C"
//#include "bbhodo.C"
//#include "bbgem.C"
#include "Bigbite.C"

std::string user_input;

Int_t gCurrentEntry = -1;
Int_t runnum;
char skipAns;

const Int_t kCanvSize = 100;
const Int_t kNrows = 24;
const Int_t kNcols = 12;
TCanvas *subCanv[2];

TString bbFile;
TChain *fbbChain;
Bigbite *Tbb;

void Setup();
void clicked_displayEntryButton();
void clicked_displayNextButton();
void clicked_displayPrevButton();
void displayEvent();


//Gem
const Int_t nmodmax = 4;   
const Int_t maxch = 60000;
const Int_t nstripsmax = 4000;
  
int Ndata_Ustrips;
double ustriplo[maxch];
int Ndata_Vstrips;
double vstriplo[maxch];
double module[maxch];
double layer[maxch];
double ADCU[maxch];
double ADCV[maxch];
UInt_t evtID;
TH2F *h[5];

//Grinch
const Int_t N_PMT=510;
const Int_t N_ROW=60;
const Int_t nChanADC = 64;
const Int_t nChanVETROC = 510;
const Int_t TDC_array = 1;
Int_t i,j;
Double_t cut = 2000;
Double_t MaxHits = 100;
Int_t next_evt=0;
Int_t usr_evt;
Int_t tdc_event[N_PMT][10];
Int_t tdcTimeArray[N_PMT];
//Array to count hits on each tube for entire run.
double TubeHits[N_PMT];
double multiTubeHits[N_PMT];
//Array to count hits on each tube for individual events.
double EvtTubeHits[N_PMT] = {};
TPaveLabel** Signal_Label;
TEllipse** PMT;
TPaveLabel** PMT_Label;
Color_t hit_color[50]={kRed, kPink+1, kMagenta, kViolet+7, kBlue, kAzure-4, kCyan, kTeal+3, kGreen, kYellow};
// SBS-Offline variables
Int_t ntdc;
Double_t grinch_tdc_pmt[512];
Double_t grinch_tdc_le[200];
Double_t grinch_tdc_te[200];
Double_t grinch_tdc_tot[50];
Double_t grinch_tdc_mult[50];
Double_t grinch_pmt_row[512];
Double_t grinch_pmt_col[512];

//BBcal
//TBox *pspmt[2][26];
//TBox *shpmt[7][27];
TH2F *hPS;
TH2F *hSH;
//Double_t psfactor;
//Double_t psfactorMax;
//Double_t shfactor;
//Double_t ShfactorMax;

//Hodo
TBox* bar[90];
TBox* pmtL[90];
TBox* pmtR[90];
TMarker* hit[90];
double ymin       = 0.00;
double barxmin    = 0.1;
double barxmax    = 0.9;
double height     = 0.0095;
double separation = (1. / 90.);
double pmtLxmin   = 0.0;
double pmtLxmax   = 0.09;
double pmtRxmin   = 0.89;
double pmtRxmax   = 1.00;

//track only events
bool trackBool;

//--------------------------------
// Set up global TCanvas and TChains
//--------------------------------
void Setup(Int_t runNum, char skipAns) {
  //need to get this path from the replay scripts output dir
  TString RootDir = gSystem->Getenv("OUT_DIR");
  bbFile = RootDir+"/gmn_replayed_"+runNum+"_stream0_seg0_76_5.root";
  //bbFile = Form("whodo/combined_%i_-1.root",runNum);
  
  fbbChain = new TChain("T");
  fbbChain->Add(bbFile);
  Tbb = new Bigbite(fbbChain);
  
  //branch address for gem display - it needs these
  fbbChain->SetBranchAddress("fEvtHdr.fEvtNum",&evtID);
  fbbChain->SetBranchAddress("Ndata.bb.gem.hit.ustriplo",&Ndata_Ustrips);
  fbbChain->SetBranchAddress("bb.gem.hit.ustriplo",ustriplo);
  fbbChain->SetBranchAddress("Ndata.bb.gem.hit.vstriplo",&Ndata_Vstrips);
  fbbChain->SetBranchAddress("bb.gem.hit.vstriplo",vstriplo);
  fbbChain->SetBranchAddress("bb.gem.hit.module",module);
  fbbChain->SetBranchAddress("bb.gem.hit.layer",layer);
  fbbChain->SetBranchAddress("bb.gem.hit.ADCmaxsampU",ADCU);
  fbbChain->SetBranchAddress("bb.gem.hit.ADCmaxsampV",ADCV);


  switch( skipAns ){
  case 'Y': case 'y':
    trackBool = true;
    cout << "Skipping events with no gem track " << endl;
    break;
    
  case 'N': case 'n':
    trackBool = false;
    cout << "Displaying all events " << endl; 
    break;
  }
}
 

namespace gui {
  TGMainFrame *main = 0;
  TGHorizontalFrame *frame1 = 0;
  TGTab *fTab;
  TGLayoutHints *fL3;
  TGCompositeFrame *tf;
  TGTextButton *exitButton;
  TGTextButton *displayEntryButton;
  TGTextButton *displayNextButton;
  TGTextButton *displayPrevButton;
  TGNumberEntry *entryInput;
  //TGLabel *ledLabel;
  
  TRootEmbeddedCanvas *canv[1];
  
  TGCompositeFrame* AddTabSub(Int_t sub) {
    tf = fTab->AddTab(Form("Run: %d",runnum));

    TGCompositeFrame *fF5 = new TGCompositeFrame(tf, (12+1)*kCanvSize,(6+1)*kCanvSize , kHorizontalFrame);
    TGLayoutHints *fL4 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX |
        kLHintsExpandY, 5, 5, 5, 5);
    TRootEmbeddedCanvas *fEc1 = new TRootEmbeddedCanvas(Form("hcalSubCanv%d",sub), fF5, 14*kCanvSize,8*kCanvSize);
    //Int_t wid = fEc1->GetCanvasWindowId();
    //subCanv[sub] = new TCanvas(Form("subCanv%d",sub),10,10,wid);
    //subCanv[sub]->Divide(12,6);
    //fEc1->AdoptCanvas(subCanv[sub]);
    canv[sub] = fEc1;
    fF5->AddFrame(fEc1,fL4);
    tf->AddFrame(fF5,fL4);
    return tf;
  }
  
  void SetupGUI() {
    if(!main) {
      main = new TGMainFrame(gClient->GetRoot()); // main window 1000 900
      frame1 = new TGHorizontalFrame(main,1,1,kFixedWidth); //"bottom bar with buttons" 150 20
      
      //frame 1 buttons and text labels
      //ledLabel = new TGLabel(frame1,"LED Bit:    , Count:      ");
      displayEntryButton = new TGTextButton(frame1,"&Display Entry","clicked_displayEntryButton()");
      entryInput = new TGNumberEntry(frame1,0,5,-1,TGNumberFormat::kNESInteger);
      displayPrevButton = new TGTextButton(frame1,"&Prev Entry","clicked_displayPrevButton()");
      displayNextButton = new TGTextButton(frame1,"&Next Entry","clicked_displayNextButton()");
      exitButton = new TGTextButton(frame1, "&Exit", "gApplication->Terminate(0)");
      
      TGLayoutHints *frame1LH = new TGLayoutHints(kLHintsTop|kLHintsLeft|
          kLHintsExpandX,2,2,2,2);
      //frame1->AddFrame(ledLabel,frame1LH);
      frame1->AddFrame(displayEntryButton,frame1LH);
      frame1->AddFrame(entryInput,frame1LH);
      frame1->AddFrame(displayPrevButton,frame1LH);
      frame1->AddFrame(displayNextButton,frame1LH);
      frame1->AddFrame(exitButton,frame1LH);
      frame1->Resize(800, displayNextButton->GetDefaultHeight());
      main->AddFrame(frame1, new TGLayoutHints(kLHintsBottom | kLHintsRight, 2, 2, 5, 1));
      
      // Create the tab widget
      fTab = new TGTab(main, 300, 300);
      fL3 = new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5);
      
      // Create Tabs
      for(Int_t i = 0; i < 1; i++) {
        tf = AddTabSub(i);
      }
      
      main->AddFrame(fTab, new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
					     kLHintsExpandY, 2, 2, 5, 1));
      main->MapSubwindows();
      main->Resize();   // resize to default size
      main->MapWindow();
      
      //gStyle->SetPalette(1);
      gStyle->SetOptStat(0);
      gStyle->SetOptTitle(0);
      for(Int_t i = 0; i < 1; i++) {
        subCanv[i] = canv[i]->GetCanvas();
	subCanv[i]->Divide(2,1,0.001,0.001);
      }
    }
  }
}


struct GEMLayer
{
  int layerID, nmodules, U_strips, V_strips;
  TString GEMtype;


};


// This function creates the histogram frame for each layer
//void DrawModules(TH2F *h, GEMLayer gem_layer){
void DrawModules(GEMLayer gem_layer){
  //Read in number of strips
  int U_tot = gem_layer.U_strips*gem_layer.nmodules;
  int V_tot = gem_layer.V_strips;
  int layer = gem_layer.layerID;
  delete h[layer];
  //Make axis ranges from 0 to 1 for UV layers because they are complicated
  if(gem_layer.GEMtype == "UV"){
    h[layer] = new TH2F(Form("Layer_%i",layer),Form("Layer %i",layer),V_tot,0,1,U_tot,0,1);
    //For XY layers simply make axis ranges the same as the number of strips
  }else{
    h[layer] = new TH2F(Form("Layer_%i",layer),Form("Layer %i",layer),V_tot,0,V_tot,U_tot,0,U_tot);
  }
  h[layer]->GetXaxis()->SetLabelSize(0.0);
  h[layer]->GetXaxis()->SetLabelColor(10); 
  h[layer]->GetYaxis()->SetLabelSize(0.0);
  h[layer]->GetYaxis()->SetLabelColor(10); 
  h[layer]->Draw();
  //Some layers have multiple modules. Draw lines to show divisions between them
  for(int imod = 0; imod < gem_layer.nmodules; imod++){
    if(imod == 0) continue;
    TLine *line = new TLine(0,imod*gem_layer.U_strips,V_tot,imod*gem_layer.U_strips);
    line->SetLineWidth(2);
    line->Draw("same");
  }
}

//Most of the work is done here. We take arrays and define all the strip 
//positions which we will use for later plotting
void StripConfig(GEMLayer gem_layer, double U_strip_line[][4], double V_strip_line[][4]){

  //The UV layers are complicated because they are not perpendicular
  if(gem_layer.GEMtype == "UV"){
    //UV layer has 30% with respect to the horizontal axis
    double angle = 30*3.14159/180;
    // Strips follow y = -tan(30)x + b. The only difference between strips is a different "b" offset
    // The strips end at (x,y) = (1,1) so we calcualte db to be enenly spaced between (0,0) and (1,1)
    double db = (1 + tan(angle))/gem_layer.U_strips;
    
    double int_L = 0;
    double int_R = 1;
    int itop = 0;
    
    //There are some plotting issues when the angled line reaches the top of the plot
    for(int istrip=0; istrip < gem_layer.U_strips; istrip++){
      if(db*istrip > 1) {
	itop = istrip;
	break;
      }
    }
    
    //Loop over strips and get positions
    for(int istrip=0; istrip < itop; istrip++){
      //Increment by db for all strips
      double b = db*(istrip+1);
      //The right side of the strip is this formula
      double right_bound = b/tan(angle);
      //If we are hitting the edge of the canvas switch to 1 instead.
      if(right_bound > 1) right_bound = 1;
      
      //All these values follow from setting left boundary to y = 0 and 
      //bottom boundary to x = 0 and using y = -tan(30)x + b to calculate
      //the intercept points
      U_strip_line[istrip][0] = int_L;
      U_strip_line[istrip][1] = right_bound;
      U_strip_line[istrip][2] = b;
      U_strip_line[istrip][3] = 0;
      
      //The V strips are the mirror image of the U strips
      V_strip_line[istrip][0] = 1 - int_L;
      V_strip_line[istrip][1] = 1 - right_bound;
      V_strip_line[istrip][2] = b;
      V_strip_line[istrip][3] = 0;
    }
    //This loop is for strips terminating on the top instead of the sides.
    for(int istrip=itop; istrip < gem_layer.U_strips; istrip++){
      double b = db*(istrip);

      //Again we use the formula for the line and simply plug in
      U_strip_line[istrip][0] = (b - 1)/tan(angle);
      U_strip_line[istrip][1] = int_R;
      U_strip_line[istrip][2] = 1;
      U_strip_line[istrip][3] = -tan(angle) + b;
      
      V_strip_line[istrip][0] = int_R - (b - 1)/tan(angle);
      V_strip_line[istrip][1] = 0;
      V_strip_line[istrip][2] = 1;
      V_strip_line[istrip][3] = -tan(angle) + b;
    }
    
  }
  //For XY layers we simply use straight vertical and horizontal lines
  //No geometry issues to consider here.
  else{
    for(int istrip=0; istrip < gem_layer.U_strips; istrip++){
      U_strip_line[istrip][0] = 0;
      U_strip_line[istrip][1] = gem_layer.V_strips;
      U_strip_line[istrip][2] = istrip;
      U_strip_line[istrip][3] = istrip;
    }
    for(int istrip=0; istrip < gem_layer.V_strips; istrip++){
      V_strip_line[istrip][0] = istrip;
      V_strip_line[istrip][1] = istrip;
      V_strip_line[istrip][2] = 0;
      V_strip_line[istrip][3] = gem_layer.U_strips;
    }
  }
}


//For a strip we take the positions and draw it on the histogram
//The strips are also colored to denote higher ADC values.
void DrawStrip(GEMLayer gem_layer,double strip,double moduleID,double strip_line[4], double adc){

 
  double x[2] = {0};
  double y[2] = {0};

 
  if(gem_layer.GEMtype == "UV"){
    strip += 40; //First 40 strips do not exist on the electronics

      x[0] = strip_line[0];
      x[1] = strip_line[1];
      y[0] = strip_line[2];
      y[1] = strip_line[3];
  }
  else{

      x[0] = strip_line[0];
      x[1] = strip_line[1];
      y[0] = gem_layer.U_strips*moduleID + strip_line[2];
      y[1] = gem_layer.U_strips*moduleID + strip_line[3];

  }
  
  //Set max ADC color to 600 and set line color depending on ADC value
  int ncolors = gStyle->GetNumberOfColors();
  int max_adc = 600;

  int ADCbin = int(adc*1.0/max_adc*double(ncolors));
  
  TLine *line = new TLine(x[0],y[0],x[1],y[1]);
  line->SetLineColor(gStyle->GetColorPalette(TMath::Max(0,TMath::Min(ncolors-1,ADCbin))));
  line->Draw("same");
}

void displayEvent(Int_t entry = -1){ 
  if(entry == -1) {
    gCurrentEntry++;
  } else {
    gCurrentEntry = entry;
  }

  if(gCurrentEntry<0) {
    gCurrentEntry = 0;
  }
  
  fbbChain->GetEntry(gCurrentEntry);

  if (trackBool == true){
    while ( Tbb->bb_tr_n == 0 ){ // only single track events)
      gCurrentEntry++;
      fbbChain->GetEntry(gCurrentEntry);
    }
  }
  gui::entryInput->SetIntNumber(gCurrentEntry);
  //fgemChain->GetEntry(gCurrentEntry);
  //fgrinchChain->GetEntry(gCurrentEntry);
  //fbbcalChain->GetEntry(gCurrentEntry);
  //fhodoChain->GetEntry(gCurrentEntry);
  fbbChain->GetEntry(gCurrentEntry);
  std::cout << endl;
  std::cout << "Displaying event " << gCurrentEntry << std::endl;
  //gui::ledLabel->SetText(TString::Format("LED Bit: %02d, Count: %5d",2,3));
  
  // --------------------------------------------------------------------------
  // Get tree
  // --------------------------------------------------------------------------
  //Long64_t nentries = fChain->GetEntriesFast();
  //for (Long64_t ev=0; ev<nentries;ev++) {

 
  //Draw Gem (histos)
  subCanv[0]->cd(1);
  gPad->Clear();   
  gPad->Divide(nlayers,1);

  ConfigParser("gem_config.cfg");
  //TH2F *h2[nlayers];
  GEMLayer gem_layer[nlayers];
  double U_strip_line[nlayers][nstripsmax][4]; // [0],[1] is x2, x2 
  double V_strip_line[nlayers][nstripsmax][4]; // [2],[3] is y1, y2
  //Set GEM layer variables and strip positions
  for(int ilayer=0; ilayer < nlayers; ilayer++){

    gem_layer[ilayer].layerID = ilayer;
    gem_layer[ilayer].nmodules = nmodules[ilayer];
    gem_layer[ilayer].U_strips = U_strips[ilayer];
    gem_layer[ilayer].V_strips = V_strips[ilayer];
    gem_layer[ilayer].GEMtype = GEMtype[ilayer];

    StripConfig(gem_layer[ilayer],U_strip_line[ilayer],V_strip_line[ilayer]);    
  }
  
  for(int ilayer=0; ilayer < nlayers; ilayer++){
    subCanv[0]->cd(1)->cd(ilayer+1);
    gPad->SetTopMargin(0.001);
    gPad->SetBottomMargin(0.001); 
    //DrawModules(h2[ilayer], gem_layer[ilayer]);
    DrawModules(gem_layer[ilayer]);
  }

  //// Loop over strips in event
  for(int istrip = 0; istrip < Ndata_Ustrips; istrip++){
         
    subCanv[0]->cd(1)->cd(layer[istrip]+1);
     
    int mod_rel = module[istrip];

    for(int ilayer = layer[istrip]; ilayer > 0; ilayer--)
      mod_rel -= gem_layer[ilayer - 1].nmodules;

    //Draw the strips on the histogram
    DrawStrip(gem_layer[(int)layer[istrip]],ustriplo[istrip],mod_rel,U_strip_line[(int)layer[istrip]][(int)ustriplo[istrip]], ADCU[istrip]);
     
  }
   
  for(int istrip = 0; istrip < Ndata_Vstrips; istrip++){
     
    subCanv[0]->cd(1)->cd(layer[istrip]+1); 
 
    int mod_rel = module[istrip];

    for(int ilayer = layer[istrip]; ilayer > 0; ilayer--)
      mod_rel -= gem_layer[ilayer - 1].nmodules;

    //Draw the strips on the histogram
    DrawStrip(gem_layer[(int)layer[istrip]],ustriplo[istrip],mod_rel,V_strip_line[(int)layer[istrip]][(int)ustriplo[istrip]], ADCV[istrip]);
    
  }
 
  
  //Divide 2nd half of canvas for pmt displays
  subCanv[0]->cd(2);
  gPad->Clear();
  gPad->Divide(4,1);
  
  //Draw Grinch
  subCanv[0]->cd(2)->cd(1);
  gPad->SetTopMargin(0.001);
  gPad->SetBottomMargin(0.001); 
  gPad->SetLeftMargin(0.001);
  
  for (int i = 0; i<510; i++){
    PMT[i]->SetFillColor(18);
    PMT[i]->Draw("same");
    PMT_Label[i]->SetFillColor(0);
    PMT_Label[i]->SetTextSize(1);
    PMT_Label[i]->Draw("same");
  }
    
  for(int r=0; r<N_PMT; r++)  //Reset individual events tube hits counter for each event.
    {
      EvtTubeHits[r] = 0;
      tdcTimeArray[r] = 0;
    }
  
  for(int j=0; j<nChanVETROC; j++)
    {
      for(int co = 0; co < TDC_array; co++)
        {
	  //grinch_tdc_mult[j] = Tbb->bb_grinch_tdc_tdc_mult[j];
	  //if (Tbb->grinch_tdc_mult[j] > 0) 
	  // cout << grinch_tdc_mult[j] << endl;
	  tdcTimeArray[j] = grinch_tdc_mult[j];
	  
      	  if(tdcTimeArray[j] > 0)
	    {
	      Int_t d = j;
	      Int_t color_shift = tdcTimeArray[d]-1;
	      
	      //cout << "PMT: " << d << " value: " << tdcTimeArray[d] << endl;
	      EvtTubeHits[d] = EvtTubeHits[d] + 1;

	      Signal_Label[d]->SetFillColor(hit_color[color_shift]);

	      Signal_Label[d]->SetTextColor(kBlack);
	      PMT[d]->SetFillColor(hit_color[color_shift]);
	      PMT_Label[d]->SetFillColor(hit_color[color_shift]);

	      Signal_Label[d]->SetTextSize(2);
	      //cout<<"tdcTimeArray[j]: "<<tdcTimeArray[d]<<" d: "<<d<<endl;
	      Signal_Label[d]->SetLabel(Form("%i",tdcTimeArray[d]));
	  
	      Signal_Label[d]->Draw("same");
	    
	    }
	}
      /*
    if ( Tbb->bb_grinch_tdc_tdc_mult[j] == 0 ) {	
      if ( Tbb->bb_grinch_tdc_tdc[(Int_t)Tbb->bb_grinch_tdc_tdcelemID[j]] <= 10 ) {
	PMT[(Int_t)Tbb->bb_grinch_tdc_tdcelemID[j]]->SetFillColor(kBlue-10);	
	PMT[(Int_t)Tbb->bb_grinch_tdc_tdcelemID[j]]->Draw("same");
      } else if ( Tbb->bb_grinch_tdc_tdc[(Int_t)Tbb->bb_grinch_tdc_tdcelemID[j]] <= 25 ){
	PMT[(Int_t)Tbb->bb_grinch_tdc_tdcelemID[j]]->SetFillColor(kBlue-7);	
	PMT[(Int_t)Tbb->bb_grinch_tdc_tdcelemID[j]]->Draw("same");
      } else{//if( Tbb->bb_grinch_tdc_tdc[(Int_t)Tbb->bb_grinch_tdc_tdcelemID[j]] <= 50 ){	
	PMT[(Int_t)Tbb->bb_grinch_tdc_tdcelemID[j]]->SetFillColor(kBlue);	
	PMT[(Int_t)Tbb->bb_grinch_tdc_tdcelemID[j]]->Draw("same");
      }
    }
      */
    }
  
    
  //Draw Preshower
  subCanv[0]->cd(2)->cd(2);
  gPad->SetTopMargin(0.001);
  gPad->SetBottomMargin(0.001); 
  
  delete hPS;
  hPS = new TH2F("hPS",
		 "Preshower",
		 2,0,2,
		 26,0,26);
  
  //cout << "Ndata for clus:" << endl;
  //cout << "atime: " << Tbb->Ndata_bb_ps_clus_atime << "  id: " << Tbb->Ndata_bb_ps_clus_id << "  e: " << Tbb->Ndata_bb_ps_clus_e << "  ec: " << Tbb->Ndata_bb_ps_clus_e_c << "  row: " << Tbb->Ndata_bb_ps_clus_row << "  col: " << Tbb->Ndata_bb_ps_clus_col << "  x: " << Tbb->Ndata_bb_ps_clus_x << "  y:" << Tbb->Ndata_bb_ps_clus_y << endl;
  //cout << "Ndata for blk:" << endl;
  //cout << "atime: " << Tbb->Ndata_bb_ps_clus_blk_atime << "  id: " << Tbb->Ndata_bb_ps_clus_blk_id << "  e: " << Tbb->Ndata_bb_ps_clus_blk_e << "  ec: " << Tbb->Ndata_bb_ps_clus_blk_e_c << "  row: " << Tbb->Ndata_bb_ps_clus_blk_row << "  col: " << Tbb->Ndata_bb_ps_clus_blk_col << "  x: " << Tbb->Ndata_bb_ps_clus_blk_x << "  y:" << Tbb->Ndata_bb_ps_clus_blk_y << endl;
  
  //cout << "data for clus:" << endl;

  //if (Tbb->Ndata_bb_ps_clus_id == 0){
    
  //cout << "atime: " << Tbb->bb_ps_clus_atime[0] << "  id: " << Tbb->bb_ps_clus_id[0] << "  e: " << Tbb->bb_ps_clus_e[0] << "  ec: " << Tbb->bb_ps_clus_e_c[0] << "  row: " << Tbb->bb_ps_clus_row[0] << "  col: " << Tbb->bb_ps_clus_col[0] << "  x: " << Tbb->bb_ps_clus_x[0] << "  y:" << Tbb->bb_ps_clus_y[0] << endl;
    
    //hPS->Fill(Tbb->bb_ps_clus_col[0],Tbb->bb_ps_clus_row[0],Tbb->bb_ps_clus_e_c[0]);
    
  //}else{
  for ( Int_t k = 0; k < Tbb->Ndata_bb_ps_clus_id; k++ ){
      
    //cout << "atime: " << Tbb->bb_ps_clus_atime[k] << "  id: " << Tbb->bb_ps_clus_id[k] << "  e: " << Tbb->bb_ps_clus_e[k] << "  ec: " << Tbb->bb_ps_clus_e_c[k] << "  row: " << Tbb->bb_ps_clus_row[k] << "  col: " << Tbb->bb_ps_clus_col[k] << "  x: " << Tbb->bb_ps_clus_x[k] << "  y:" << Tbb->bb_ps_clus_y[k] << endl;
    
    if ( Tbb->bb_ps_clus_e_c[k] > 0) {
      hPS->Fill(Tbb->bb_ps_clus_col[k],Tbb->bb_ps_clus_row[k],Tbb->bb_ps_clus_e_c[k]);
    } 
  }

  
  hPS->GetXaxis()->SetLabelSize(0.0);
  hPS->GetXaxis()->SetLabelColor(10); 
  hPS->GetYaxis()->SetLabelSize(0.0);
  hPS->GetYaxis()->SetLabelColor(10);
  hPS->SetFillColor(kBlue);
  hPS->Draw("box");
  

  //Draw Hodo
  subCanv[0]->cd(2)->cd(3);
  gPad->SetTopMargin(0.001);
  gPad->SetBottomMargin(0.001); 
    
  for(int j=0; j<90; j++){
    // reset bar and pmt colours for each event
    if( j == 0 ) {
      bar[j]->SetFillColor(kGray+1);
      bar[j]->Draw();
      pmtL[j]->SetFillColor(kGray+2);
      pmtL[j]->Draw("same");
      pmtR[j]->SetFillColor(kGray+2);
      pmtR[j]->Draw("same");
    }
    else {
      bar[j]->SetFillColor(kGray+1);
      bar[j]->Draw("same");
      pmtL[j]->SetFillColor(kGray+2);
      pmtL[j]->Draw("same");
      pmtR[j]->SetFillColor(kGray+2);
      pmtR[j]->Draw("same");
    }

    hit[j]->SetMarkerStyle(28);
    hit[j]->SetMarkerColor(1);
    hit[j]->SetMarkerSize(1.0);

    /*
    // use GenericDetector variables for PMT display (blue if a "good" tdc hit)
    if( fabs(Tbb->bb_hodotdc_clus_size[j]) < 10000. ) {
      pmtL[j]->SetFillColor(4);
      pmtL[j]->Draw("same");
    }
    if(fabs(Tbb->bb_hodotdc_clus_size[j+90]) < 10000. ) {
      pmtR[j]->SetFillColor(4);
      pmtR[j]->Draw("same");
    }
    */
  }
    
    
  // use TimingHodoscope variables for bar display (red if a "good" hit)
    for(int j=0; j<(Int_t)Tbb->Ndata_bb_hodotdc_clus_bar_tdc_id; j++) {
    bar[(Int_t)Tbb->bb_hodotdc_clus_bar_tdc_id[j]]->SetFillColor(2);
    bar[(Int_t)Tbb->bb_hodotdc_clus_bar_tdc_id[j]]->Draw("same") ;

    //cout << Tbb->bb_hodotdc_allclus_ymean[j] << endl;

    //if ( fabs(Tbb->bb_hodotdc_allclus_ymean[(Int_t)Tbb->Ndata_bb_hodotdc]) < 0.3 ) {
  
    hit[j]->SetX(0.5+Tbb->bb_hodotdc_clus_bar_tdc_timehitpos[j]);
    hit[j]->SetY(ymin+height/2. + (Int_t)Tbb->bb_hodotdc_clus_bar_tdc_id[j]*separation);
    hit[j]->Draw("same");
    pmtL[(Int_t)Tbb->bb_hodotdc_clus_bar_tdc_id[j]]->SetFillColor(4);
    pmtL[(Int_t)Tbb->bb_hodotdc_clus_bar_tdc_id[j]]->Draw("same");
    pmtR[(Int_t)Tbb->bb_hodotdc_clus_bar_tdc_id[j]]->SetFillColor(4);
    pmtR[(Int_t)Tbb->bb_hodotdc_clus_bar_tdc_id[j]]->Draw("same");
    //}
    }
  

  //Draw Shower
  subCanv[0]->cd(2)->cd(4);
  gPad->SetTopMargin(0.001);
  gPad->SetBottomMargin(0.001); 
    
  delete hSH;
  hSH = new TH2F("bSH",
		 "Shower",
		 7,0,7,
		 27,0,27);
  
  for ( int k = 0; k < Tbb->Ndata_bb_sh_clus_id; k++ ){
    if ( Tbb->bb_sh_clus_e_c[k] > 0) {
      hSH->Fill(Tbb->bb_sh_clus_col[k],Tbb->bb_sh_clus_row[k],Tbb->bb_sh_clus_e_c[k]);
    }
  } 
  
  hSH->GetXaxis()->SetLabelSize(0.0);
  hSH->GetXaxis()->SetLabelColor(10); 
  hSH->GetYaxis()->SetLabelSize(0.0);
  hSH->GetYaxis()->SetLabelColor(10);
  hSH->SetFillColor(kBlue);
  hSH->Draw("box");
  gPad->Update();
  
  //END OF DISPLAYEVENT()
}




void clicked_displayNextButton()
{
  //if(gCurrentEntry>gMaxEntries);
  gui::entryInput->SetIntNumber(++gCurrentEntry);
  displayEvent(gCurrentEntry);
  gSystem->ProcessEvents();
}

void clicked_displayPrevButton()
{
  //if(gCurrentEntry>gMaxEntries);
  gui::entryInput->SetIntNumber(--gCurrentEntry);
  displayEvent(gCurrentEntry);
  gSystem->ProcessEvents();
}


void clicked_displayEntryButton()
{
  gCurrentEntry = gui::entryInput->GetIntNumber();
  displayEvent(gCurrentEntry);
}


bool is_number(const std::string& mystring)
{
  std::string::const_iterator it = mystring.begin();
  while (it != mystring.end() && std::isdigit(*it)) ++it;
  return !mystring.empty() && it == mystring.end();
}


void DisplayBB(){
  std::cout << "Enter Run Number" << endl;
  std::cin >> runnum;
  
  while (skipAns != 'y' && skipAns != 'Y' && skipAns != 'n' && skipAns != 'N'){
  std::cout << "Skip events with no tracks? (y/n)" << endl;
  std::cin >> skipAns;
  }
  
  Setup(runnum,skipAns);
  
  Long64_t Nevents = fbbChain->GetEntries();
  cout << Nevents << " events in run." << endl;

  //--------------------------------
  // Gems - later
  //--------------------------------
  ConfigParser("gem_config.cfg");
  //TH2F* h[nlayers];
  //--------------------------------
  // Grinch
  //--------------------------------
  double PMT_radius = 1./N_ROW; //the whole canvas has a max size of 1. We divide by the number of rows. 
  double x,y;//center of PMT. Draws PMT array.
  double labelsize=0.5;
  Signal_Label = new TPaveLabel*[N_PMT];
  PMT = new TEllipse*[N_PMT];
  PMT_Label = new TPaveLabel*[N_PMT];
  Int_t k=0;
  double y_off = -0.0;
  double xrad = 0.5/9.;
  double yrad = 0.5/60.;
  for (int i=0; i<N_ROW; i++ )
    {
      int N_COL = 8 + i%2;
      for (int j=0; j<N_COL; j++ )
        {
	  if ( i%2==0 )
	    {
              x = ((float)j+1.5) * 2. * xrad - xrad;
              //y = ((float)N_ROW - ((float)i * sin(60*TMath::DegToRad())) + 0.5) * yrad  +  y_off;
	      y = (float)(N_ROW - i) * 2. * yrad - yrad;
	    }
          else
            {
              x = (float)(j+1) * 2. * xrad - xrad;
	      // y = ((float)N_ROW - ((float)i * sin(60*TMath::DegToRad())) + 0.5) * yrad  +  y_off;
	      y = (float)(N_ROW - i) * 2. * yrad - yrad;
	    }
	  
	  Signal_Label[k]=new TPaveLabel(x-1.5*xrad*labelsize,y+yrad*labelsize,x+3*xrad*labelsize,y+2*yrad*labelsize,"");
          Signal_Label[k]->SetTextColor(2);
          Signal_Label[k]->SetFillColor(0);
          Signal_Label[k]->SetLabel("");
	  PMT_Label[k] = new TPaveLabel(x-xrad*labelsize, y-yrad*labelsize, x+xrad*labelsize, y+yrad*labelsize, Form("%d",k));
          PMT[k]=new TEllipse(x, y, xrad, yrad);
	  k++;
	  //cout << k << endl;
	    }
    }
  
  /*
  //--------------------------------
  // Preshower + Shower calorimeters
  //--------------------------------
  double PSymin       = 0.05;
  double PSheight     = 0.03;
  double PSwidth      = 0.41;
  double PSradius = 0.03;
  double PSseperation = 0.0346;
  double PSpmtxmin   = 0.10;
  double PSpmtxmax   = 0.49;
  
  
  
  for(int i = 0; i<2; i++){
    for(int j = 0; j<26; j++){
      pspmt[i][j] = new TBox((PSpmtxmin + i*PSwidth), (PSymin + j*PSseperation), (PSpmtxmax+ i*PSwidth),(PSymin + PSheight + j*PSseperation));
      //x = (i+1) * PMT_radius * 5  + 0.05;
      //y = (j+1) * PMT_radius + 0.05;
      //pspmt[i][j] = new TEllipse(x,y,5*PSradius,PSradius/2);
    }
  }
  

  double SHymin       = 0.05;
  double SHheight     = 0.03;
  double SHwidth      = 0.128;
  double SHseperation = 0.0333;
  double SHpmtxmin   = 0.10;
  double SHpmtxmax   = 0.22;

  TH2F *hSH = new TH2F("pSH",
		       "Shower",
		       35,0,7,
		       27,0,27);
  
  for(int i = 0; i<7; i++){
    for(int j = 0; j<27; j++){
      shpmt[i][j] = new TBox((SHpmtxmin + i*SHwidth), (SHymin + j*SHseperation), (SHpmtxmax+ i*SHwidth),(SHymin + SHheight + j*SHseperation));
      //x = (j+1) * PMT_radius * 3  + 0.05;
      //y = (i+1) * PMT_radius * 0.5 + 0.05;
      //shpmt[i][j] = new TEllipse(x,y,sh,yrad);
    }
  }
  */
  
  //--------------------------------
  // Hodoscope
  //--------------------------------
  

  for(int i = 0; i<90; i++ ) {
    bar[i]  = new TBox(barxmin, (ymin + (float)i*separation), barxmax, (ymin + height + (float)i*separation));
    pmtL[i] = new TBox(pmtLxmin,(ymin + (float)i*separation), pmtLxmax,(ymin + height + (float)i*separation));
    pmtR[i] = new TBox(pmtRxmin,(ymin + (float)i*separation), pmtRxmax,(ymin + height + (float)i*separation));
    hit[i]  = new TMarker(0,0,2);
  }
  
  gui::SetupGUI();
  displayEvent(gCurrentEntry);
  
  /*
  while(!kButtonEngaged) {
    while( user_input != "q" ) {
      if(is_number(user_input)) {
	gCurrentEntry = std::stoi(user_input);
      } else{
	gCurrentEntry++;
      }
      displayEvent(gCurrentEntry);
      gSystem->ProcessEvents();
      std::cout << "Display options: Press <enter> for next event. Enter event number to display specific event,  or q to stop." << std::endl;
      getline(std::cin,user_input);
    }
  }
  */
}
 
