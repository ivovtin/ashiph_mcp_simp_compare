void n113_npe_sigma()
{
   gStyle->SetOptStat(0);

   TString dir_out = "out";

   TCanvas *c1 = new TCanvas();
   c1->cd();	

   TF1 *nphe11 = new TF1("nphe11","([2]+[1]*((x^2-(1/([0]^2-1)))/(x^2))*(1+TMath::Erf((x-abs(sqrt(1/([0]^2-1))))/0.0001))/2)+[1]/2*[3]*(1+1/(x^2))*(1/[4]-1/(x^2))*(1+TMath::Erf((x-abs(sqrt([4])))/0.0001))/2",0,20);
   nphe11->SetLineColor(4);
   nphe11->SetLineWidth(2);
   nphe11->SetLineStyle(10);

   nphe11->SetParameter(0,1.13);
   nphe11->SetParName(0,"n=1.13");
   nphe11->SetParName(1,"N_{pe}(#beta=1)");
   //
   //float mu0 = 0.11;      // t15 u=53.5v 
   //float max_mu = 12.71;  // t15 u=53.5v p2
   //float max_mu = 8.17;  // t15 u=53.5v p9
   //
   float mu0 = 0.80;       // t45 u=57v 
   //float max_mu = 11.97;   // t45 u=57v p2
   float max_mu = 7.65;  // t45 u=57v p9
   //
   nphe11->SetParameter(1,max_mu);    // for n=1.05 6cm Integral SiPM OverVoltage=2V -> PDE(500)=25%
   nphe11->SetParName(2,"Const.lev from SiPM noises");
   nphe11->SetParameter(2,mu0);  // mcp - snd
   //nphe11->SetParameter(3,0.14);  
   nphe11->SetParameter(3,mu0);   // mcp - snd  
   nphe11->SetParName(3,"#ksi"); 
   //nphe11->SetParameter(4,1.6); // for n=1.05
   //nphe11->SetParameter(4,1.0); // for n=1.05
   nphe11->SetParameter(4,2.4); // for n=1.05
   nphe11->SetParName(4,"T_{Ch.thr}^{#delta}(n=1.13)");

   Float_t me=0.5;
   Float_t mmu=106.0; 
   Float_t mpi=139.0;
   Float_t mka=498.0;

   const int ncounts = 6;
   float npethr_step = 1.0;
   float P1 = 400;
   float P2 = 900;

   TProfile* npee=new TProfile("npee","",1000,-0.5,990.5); // npe vs momentum for electrons
   TProfile* npemu=new TProfile("npemu","",1000,-0.5,990.5); // npe vs momentum for muons
   TProfile* npepi=new TProfile("npepi","",1000,-0.5,990.5); // npe vs momentum for pions
   TProfile* npek=new TProfile("npek","",1000,-0.5,990.5); // npe vs momentum for pions
   TProfile* pr1 = new TProfile("sigma-sipm",TString::Format("%.1f < P, MeV/c < %.1f",P1,P2).Data(),100,0,6.5,0,7);  // sigma vs npe_trh 

   TProfile* Ksigma[ncounts];
   for (int k=0; k<ncounts; k++)
   {
      Ksigma[k] = new TProfile(TString::Format("Ksigma%d",k).Data(),TString::Format("Nphe_{thr}=%.1f",k*npethr_step).Data(),250,-0.5,990.5);
   }
   
   //TH2F* null2=new TH2F("null2","",200,0,1000,55,-0.05,max_mu+0.5); 
   TH2F* null2=new TH2F("null2","",200,0,1000,55,-0.05,max_mu+3.5); 
   null2->SetXTitle("P, #frac{MeV}{c}");
   null2->SetYTitle("N_{pe}");

   for (int i=0; i<1000; i++)
   {
      int N0ee[ncounts]={0}, N0pi[ncounts]={0}, N0k[ncounts]={0}, Ntotal=0;   
      float Npethr[ncounts] = {0.};

      for (int l=0; l<50000; l++) 
      { 
         float mu_e = nphe11->Eval((float)i/me);
         float mu_pi = nphe11->Eval((float)i/mpi);
         float mu_mu = nphe11->Eval((float)i/mmu);
         float mu_k = nphe11->Eval((float)i/mka);

         float Nphe_e = gRandom->Poisson(mu_e);
         float Nphe_pi = gRandom->Poisson(mu_pi);
         float Nphe_mu = gRandom->Poisson(mu_mu);
         float Nphe_k = gRandom->Poisson(mu_k);
      
         npee->Fill((float)i,Nphe_e);
         npemu->Fill((float)i,Nphe_mu);
         npepi->Fill((float)i,Nphe_pi);
         npek->Fill((float)i,Nphe_k);
         
         for( int k=0; k<ncounts; k++ )
         {
            Npethr[k] = k*npethr_step;
            if( Nphe_e<Npethr[k] ) N0ee[k]++;
            if( Nphe_pi<Npethr[k] ) N0pi[k]++;
            if( Nphe_k<Npethr[k] ) N0k[k]++;
         }
         Ntotal++;   
      }
      float underthr_eff[ncounts];
      float eff[ncounts];   
   
      for( int k=0; k<ncounts; k++ )
      {
         underthr_eff[k] = (float)N0k[k]/(float)Ntotal;
         eff[k] = (float)N0pi[k]/(float)Ntotal;  
         //if( k==5 && i<300 ) cout<<i<<"\t"<<Ntotal<<"\t"<<Npethr[k]<<"\t"<<N0ee[k]<<"\t"<<N0pi[k]<<endl;
         Ksigma[k]->Fill((float)i,abs(sqrt(2.)*(TMath::ErfInverse(1-2*underthr_eff[k])+TMath::ErfInverse(1-2*(1-eff[k])))));

         if( i>P1 && i<P2 ) pr1->Fill((float)Npethr[k],abs(sqrt(2.)*(TMath::ErfInverse(1-2*underthr_eff[k])+TMath::ErfInverse(1-2*(1-eff[k])))));
      }
   }

   null2->Draw();

   npee->SetLineColor(kRed);
   npee->SetMarkerColor(kRed);
   npee->Draw("same");

   npemu->SetLineColor(kGreen);
   npemu->SetMarkerColor(kGreen);
   npemu->Draw("same");

   npepi->SetLineColor(kBlue);
   npepi->SetMarkerColor(kBlue);
   npepi->Draw("same");

   npek->SetLineColor(kBlack);
   npek->SetMarkerColor(kBlack);
   npek->Draw("same");

   TLegend* leg1=new TLegend(100,max_mu/3,200,max_mu/3+1.5,"","");
   TLegendEntry *le1=leg1->AddEntry(npee,"e","elp");
   le1->SetTextColor(kRed);
   TLegendEntry *le2=leg1->AddEntry(npemu,"#mu","elp");
   le2->SetTextColor(kGreen);
   TLegendEntry *le3=leg1->AddEntry(npepi,"#pi","elp");
   le3->SetTextColor(kBlue);
   TLegendEntry *le4=leg1->AddEntry(npek,"K","elp");
   le4->SetTextColor(kBlack);
   leg1->Draw("same");

   c1->SaveAs(dir_out + "/" + "npe_momentum.png");
   
   for( int i=0; i<ncounts; i++ )
   {
      TCanvas *c2 = new TCanvas();
      c2->cd();   
      Ksigma[i]->SetLineColor(kRed);
      Ksigma[i]->SetMarkerColor(kRed);
      //Ksigma[i]->SetTitle("; P, MeV/c; #sigma");
      Ksigma[i]->GetXaxis()->SetTitle("P, MeV/c");
      Ksigma[i]->GetYaxis()->SetTitle("#sigma");
      Ksigma[i]->Draw("same");
      c2->SaveAs(dir_out + "/" + TString::Format("Ksigma_%d.png",i).Data());
   }

   TCanvas *c3 = new TCanvas();
   c3->cd();   
   pr1->SetLineColor(kBlue);
   pr1->SetMarkerColor(kBlue);
   pr1->SetMarkerStyle(4);
   //pr1->SetTitle("; Threshold, npe; #sigma");
   pr1->GetXaxis()->SetTitle("Threshold, npe");
   pr1->GetYaxis()->SetTitle("#sigma");
   pr1->Draw("prof");
   c3->SaveAs(dir_out + "/" + "sigma.png");

   TCanvas *c4 = new TCanvas();
   c4->cd();   
   TString dir_in = "out/snd_1.13";
   TString fnamein = dir_in + "/" + "results.root";
   TFile *infile = TFile::Open(fnamein);

   TProfile* pr2 = (TProfile*)infile->Get("sigma");
   pr1->Draw("prof");
   pr2->Draw("same");
       
   //TLegend* leg2=new TLegend(5, 3.5, 6, 4.0,"","");
   TLegend* leg2=new TLegend(5, 1.5, 6, 2.0,"","");
   TLegendEntry *le5=leg2->AddEntry(pr1,"ashiph sipm","elp");
   le5->SetTextColor(kBlue);
   TLegendEntry *le6=leg2->AddEntry(pr2,"ashiph mcp","elp");
   le6->SetTextColor(kRed);
   leg2->Draw("same");
    
   c4->SaveAs(dir_out + "/" + "sigmas_compare.png");
} 
