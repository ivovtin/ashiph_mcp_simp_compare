void deltae_snd()
{
TCanvas* mcan=new TCanvas("mcan","",0,0,800,600);
TVectorF vx1(12);
TVectorF vx2(16);
TVectorF npe1(12);
TVectorF npe2(16);

Float_t mmu=106.0; 
Float_t mpi=139.0;
Float_t mka=498.0;

vx1[0]=464./mpi;
vx1[1]=504./mpi;
vx1[2]=544./mpi;
vx1[3]=584./mpi;
vx1[4]=624./mpi;
vx1[5]=664./mpi;
vx2[0]=644./mpi;
vx2[1]=684./mpi;
vx2[2]=724./mpi;
vx2[3]=764./mpi;
vx2[4]=804./mpi;
vx2[5]=844./mpi;
vx2[6]=884./mpi;
vx2[7]=924./mpi;

vx1[6]=464./mmu;
vx1[7]=504./mmu;
vx1[8]=544./mmu;
vx1[9]=584./mmu;
vx1[10]=624./mmu;
vx1[11]=664./mmu;
vx2[8]=644./mmu;
vx2[9]=684./mmu;
vx2[10]=724./mmu;
vx2[11]=764./mmu;
vx2[12]=804./mmu;
vx2[13]=844./mmu;
vx2[14]=884./mmu;
vx2[15]=924./mmu;

npe1[0]=0.0; //mpi;
npe1[1]=0.0;//mpi;
npe1[2]=0.0;//mpi;
npe1[3]=1.2;//mpi;
npe1[4]=3.6;//mpi;
npe1[5]=5.7;//mpi;
npe2[0]=0.0;//mpi;
npe2[1]=0.0;//mpi;
npe2[2]=0.0;//mpi;
npe2[3]=0.0;//mpi;
npe2[4]=0.2;//mpi;
npe2[5]=1.3;//mpi;
npe2[6]=3.6;//mpi;
npe2[7]=5.7;//mpi;

npe1[6]=3.1;//mmu;
npe1[7]=5.8;//mmu;
npe1[8]=8.0;//mmu;
npe1[9]=9.6;//mmu;
npe1[10]=11.0;//mmu;
npe1[11]=12.2;//mmu;
npe2[8]=2.3;//mmu;
npe2[9]=4.4;//mmu;
npe2[10]=6.2;//mmu;
npe2[11]=7.4;//mmu;
npe2[12]=9.1;//mmu;
npe2[13]=10.1;//mmu;
npe2[14]=11.0;//mmu;
npe2[15]=11.9;//mmu;

 gStyle->SetOptStat(0);

   TF1 *nphe1 = new TF1("nphe1","([1]*((x^2-(1/([0]^2-1)))/(x^2))*(1+TMath::Erf((x-abs(sqrt(1/([0]^2-1))))/0.0001))/2)",0,20);
   nphe1->SetLineColor(4);
   nphe1->SetLineWidth(2);
   nphe1->SetParameter(0,1.13);
   nphe1->SetParName(0,"n=1.13");
//   nphe1->SetParameter(0,1.05);
//   nphe1->SetParName(0,"n=1.05");
   nphe1->SetParName(1,"N_{pe}(#beta=1)");
   nphe1->SetParameter(1,8);
 
 TF1* nd1=new TF1("nd1","[1]*[0]*(1+1/(x^2))*(1/[2]-1/(x^2))*(1+TMath::Erf((x-abs(sqrt([2])))/0.0001))/2",0,20);   
 nd1->SetParameter(0,0.1);
 nd1->SetParName(0,"#ksi");
 nd1->SetParameter(1,nphe1->GetParameter(1)/2);
 nd1->SetParName(1,"N_{pe}^{#delta}($beta=1)/2");
 nd1->SetLineStyle(9);
 nd1->SetParameter(2,0.97);
 nd1->SetParName(2,"T_{Ch.thr}^{#delta}(n=1.13)");
// nd1->SetParameter(2,1.6);
// nd1->SetParName(2,"T_{Ch.thr}^{#delta}(n=1.05)");
 nd1->SetLineColor(1);

   TF1 *nphe11 = new TF1("nphe11","([2]+[1]*((x^2-(1/([0]^2-1)))/(x^2))*(1+TMath::Erf((x-abs(sqrt(1/([0]^2-1))))/0.0001))/2)+[1]/2*[3]*(1+1/(x^2))*(1/[4]-1/(x^2))*(1+TMath::Erf((x-abs(sqrt([4])))/0.0001))/2",0,20);
   nphe11->SetLineColor(4);
   nphe11->SetLineWidth(2);
   nphe11->SetLineStyle(10);
   nphe11->SetParameter(0,1.05);
   nphe11->SetParName(0,"n=1.05");
//   nphe11->SetParameter(0,1.03);
//   nphe11->SetParName(0,"n=1.03");
//   nphe11->SetParameter(0,1.13);
//   nphe11->SetParameter(0,1.008);
//   nphe11->SetParName(0,"n=1.13");
//   nphe11->SetParName(0,"n=1.008");
   nphe11->SetParName(1,"N_{pe}(#beta=1)");
//   nphe11->SetParameter(1,8); // for n=1.05
//   nphe11->SetParameter(1,20);  // for n=1.13
   nphe11->SetParameter(1,17.6);  // for n=1.05 6cm Integral SiPM OverVoltage=2V -> PDE(500)=25%
//   nphe11->SetParameter(1,21.1);  // for n=1.03 10+10cm Integral SiPM OverVoltage=2V -> PDE(500)=25%
//   nphe11->SetParameter(1,10.);  // for n=1.03 6cm Integral SiPM OverVoltage=2V -> PDE(500)=25%
//   nphe11->SetParameter(1,21.9);  // for n=1.13 normalized amplitude for SiPM
//   nphe11->SetParameter(1,7.2);  // for n=1.008 normalized amplitude for SiPM
//   nphe11->SetParameter(1,11.9);  // for n=1.13 normalized amplitude for MCP PMT
//   nphe11->SetParameter(1,28.2);  // for n=1.05 10cm Integral SiPM OverVoltage=4V -> PDE(500)=40%
   nphe11->SetParName(2,"Const.lev from SiPM noises");
   nphe11->SetParameter(2,0.351); // for 4.5MHz and window 78ns
//   nphe11->SetParameter(2,0.0078); // for 100KHz and window 78ns
 nphe11->SetParameter(3,0.14);
 nphe11->SetParName(3,"#ksi");
// nphe11->SetParameter([1]/2,nphe2->GetParameter(1)/2);
// nd2->SetParName(1,"N_{pe}^{#delta}($beta=1)/2");
// nd2->SetLineStyle(9);
// nphe11->SetParameter(4,1.0);
// nphe11->SetParName(4,"T_{Ch.thr}^{#delta}(n=1.13)");
// nphe11->SetParameter(4,2.13); // for n=1.03
// nphe11->SetParName(4,"T_{Ch.thr}^{#delta}(n=1.03)");
 nphe11->SetParameter(4,1.6); // for n=1.05
 nphe11->SetParName(4,"T_{Ch.thr}^{#delta}(n=1.05)");
// nphe11->SetParameter(4,0.97); // for n=1.13
// nphe11->SetParameter(4,4.0); // for n=1.008
// nphe11->SetParName(4,"T_{Ch.thr}^{#delta}(n=1.13)");
// nphe11->SetParName(4,"T_{Ch.thr}^{#delta}(n=1.008)");

TH2F* null=new TH2F("null","N_{pe}^{Cher.} and #theta=0",200,0,20,251,-0.05,25.05); // for n=1.05
//TH2F* null=new TH2F("null","N_{pe}^{Cher.} and #theta=0",200,0,20,251,-0.05,15.05); // for n=1.03 and t=6cm
//TH2F* null=new TH2F("null","N_{pe}^{Cher.} and #theta=0",200,0,20,251,-0.05,7.5); // for n=1.008
null->SetXTitle("#beta#gamma");
null->SetYTitle("N_{pe}");
TGraph* gnpe1=new TGraph(vx1,npe1);
TGraph* gnpe2=new TGraph(vx2,npe2);

null->Draw();
nphe11->Draw("same");
nd1->Draw("same");

TLegend* leg0=new TLegend(10,13,19,6,"","");
//leg0->AddEntry(nd1,"n=1.008 N_{pe}^{#delta}","lp");
//leg0->AddEntry(nphe11,"n=1.008 N_{pe}^{#mu,#pi,K}+N_{pe}^{#delta}+N_{pe}^{noi}","lp");
//leg0->AddEntry(nd1,"n=1.13 N_{pe}^{#delta}","lp");
//leg0->AddEntry(nphe11,"n=1.13 N_{pe}^{#mu,#pi,K}+N_{pe}^{#delta}+N_{pe}^{noi}","lp");
leg0->AddEntry(nd1,"n=1.05 N_{pe}^{#delta}","lp");
leg0->AddEntry(nphe11,"n=1.05 N_{pe}^{#mu,#pi,K}+N_{pe}^{#delta}+N_{pe}^{noi}","lp");
//leg0->AddEntry(nd1,"n=1.03 N_{pe}^{#delta}","lp");
//leg0->AddEntry(nphe11,"n=1.03 N_{pe}^{#mu,#pi,K}+N_{pe}^{#delta}+N_{pe}^{noi}","lp");
leg0->Draw("same");

mcan->SaveAs("mcan1.png");

TCanvas* mcan2=new TCanvas("mcan2","",200,100,1200,600);
/*
TProfile* npek=new TProfile("npek","",2000,-0.5,1999.5); // npe vs momentum for Kaons
TProfile* npepi=new TProfile("npepi","",2000,-0.5,1999.5); // npe vs momentum for pions
TProfile* npemu=new TProfile("npemu","",2000,-0.5,1999.5); // npe vs momentum for muons
*/
TProfile* npek=new TProfile("npek","",20000,-0.5,19990.5); // npe vs momentum for Kaons
TProfile* npepi=new TProfile("npepi","",20000,-0.5,19990.5); // npe vs momentum for pions
TProfile* npemu=new TProfile("npemu","",20000,-0.5,19990.5); // npe vs momentum for muons

TH2F* null2=new TH2F("null2","",200,0,20000,251,-0.05,25.505); // for n=1.05 and t=6cm
//TH2F* null2=new TH2F("null2","",200,0,20000,251,-0.05,15.505); // for n=1.03 and t=6cm
//TH2F* null2=new TH2F("null2","",200,0,20000,251,-0.05,7.505); // for n=1.008
null2->SetXTitle("P, #frac{MeV}{c}");
null2->SetYTitle("N_{pe}");
//null2->SetYTitle("N_{pe}/N_{pe}(#beta=1)");

for (int i=0; i<20000; i++)
	{
	for (int l=0; l<100; l++) { npek->Fill((float)i,gRandom->Poisson(nphe11->Eval((float)i/mka)));}
//npek->Fill((float)i,nphe11->Eval((float)i/mka));}
	for (int l=0; l<100; l++) { npepi->Fill((float)i,gRandom->Poisson(nphe11->Eval((float)i/mpi)));}
//npepi->Fill((float)i,nphe11->Eval((float)i/mpi));}
	for (int l=0; l<100; l++) { npemu->Fill((float)i,gRandom->Poisson(nphe11->Eval((float)i/mmu)));}
//npemu->Fill((float)i,nphe11->Eval((float)i/mmu));}
	}
null2->Draw();
npek->SetLineColor(kBlack);
npek->SetMarkerColor(kBlack);
npek->Draw("same");
npepi->SetLineColor(kRed);
npepi->SetMarkerColor(kRed);
npepi->Draw("same");
npemu->SetLineColor(kBlue);
npemu->SetMarkerColor(kBlue);
npemu->Draw("same");
TLegend* leg1=new TLegend(200,2.3,500,5.,"","");
TLegendEntry *le1=leg1->AddEntry(npek,"K","elp");
le1->SetTextColor(kBlack);
//leg1->GetEntry()->SetTextColor(kBlack);
//gPad->Update();
//TLegend* leg2=new TLegend(200,1.2,500,1.1,"","");
TLegendEntry *le2=leg1->AddEntry(npepi,"#pi","elp");
le2->SetTextColor(kRed);
//TLegend* leg3=new TLegend(200,1.1,500,1.0,"","");
TLegendEntry *le3=leg1->AddEntry(npemu,"#mu","elp");
le3->SetTextColor(kBlue);
//leg1->GetEntry()->SetTextColor(kBlack);
leg1->Draw("same");
//leg2->Draw("same");
//leg3->Draw("same");
gPad->Update();

mcan2->SaveAs("mcan2.png");

}

//--------END-----//
