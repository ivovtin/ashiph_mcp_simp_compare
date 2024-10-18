void npe_momentum_draw_parametrization2()
{
   float x, y1, y2, y3, y4, y5;
   vector< float > mom, curve1, curve2, curve3, curve4, curve5;
   vector< float > Ksigma_mu, Ksigma_pi, Ksigma_sipm;   

   std::ifstream f_in("graph.dat");
   while( !f_in.eof() )
   {
	  std::string s;
     while( getline(f_in, s) )
     {
         if( s.size()>=1 && s[0]!='#' )
         {
            std::istringstream i_str(s);
            i_str >> x >> y1 >> y2 >> y3 >> y4 >> y5;
            //i_str >> x >> y1 >> y2;
            mom.push_back(x);
            curve1.push_back(y1);
            if(y2<0.06) y2=0.06; 
            curve2.push_back(y2);
            curve3.push_back(y3);
            curve4.push_back(y4);
            curve5.push_back(y5);
            //cout<<x<<"\t"<<y1<<endl;

            //eff = 1-exp(-y3);
            //underthr_eff_mu = 1-exp(-y1);
            //underthr_eff_pi = 1-exp(-y2);

            //Ksigma = abs(sqrt(2.)*(TMath::ErfInverse(1-2*underthr_eff)+TMath::ErfInverse(1-2*(1-eff))));
            //Ksigma_mu.push_back(abs(sqrt(2.)*(TMath::ErfInverse(1-2*underthr_eff_mu)+TMath::ErfInverse(1-2*(1-eff)))));
            //Ksigma_pi.push_back(abs(sqrt(2.)*(TMath::ErfInverse(1-2*underthr_eff_pi)+TMath::ErfInverse(1-2*(1-eff)))));

            //if(x<445) Ksigma_sipm.push_back(4.5);  //p2 - u=53.5V, t=15C
            //if(x<445) Ksigma_sipm.push_back(4.0);  //p9 - u=53.5V, t=15C
            //if(x<445) Ksigma_sipm.push_back(4.0);  //p2 - u=57V, t=45C
            //if(x<445) Ksigma_sipm.push_back(3.0);  //p9 - u=57V, t=45C
         }
      }
   }
   f_in.close();
   
   gStyle->SetOptStat(0);

   TCanvas *c1 = new TCanvas();
   c1->cd();	

   TF1 *nphe11 = new TF1("nphe11","([2]+[1]*((x^2-(1/([0]^2-1)))/(x^2))*(1+TMath::Erf((x-abs(sqrt(1/([0]^2-1))))/0.0001))/2)+[1]/2*[3]*(1+1/(x^2))*(1/[4]-1/(x^2))*(1+TMath::Erf((x-abs(sqrt([4])))/0.0001))/2",0,20);
   nphe11->SetLineColor(4);
   nphe11->SetLineWidth(2);
   nphe11->SetLineStyle(10);

   nphe11->SetParameter(0,1.05);
   nphe11->SetParName(0,"n=1.05");
   nphe11->SetParName(1,"N_{pe}(#beta=1)");
   //nphe11->SetParameter(1,17.6);  // for n=1.05 6cm Integral SiPM OverVoltage=2V -> PDE(500)=25%
   //nphe11->SetParameter(1,2.40);  // mcp - snd 
   nphe11->SetParameter(1,7.14);    // u=53.5 t=15 n=1.05
   //nphe11->SetParameter(1,3.91);    // u=53.5 t=15 n=1.05
   //nphe11->SetParameter(1,5.13);    // u=53.5 t=15 n=1.051
   nphe11->SetParName(2,"Const.lev from SiPM noises");
   nphe11->SetParameter(2,0.351); // for 4.5MHz and window 78ns
   //nphe11->SetParameter(2,0.05);  // mcp - snd
   //nphe11->SetParameter(3,0.14);  
   nphe11->SetParameter(3,0.05);   // mcp - snd  
   nphe11->SetParName(3,"#ksi"); 
   nphe11->SetParameter(4,1.6); // for n=1.05
   //nphe11->SetParameter(4,0.447369); // for n=1.05
   nphe11->SetParName(4,"T_{Ch.thr}^{#delta}(n=1.05)");

   Float_t me=0.5;
   Float_t mmu=106.0; 
   Float_t mpi=139.0;
   Float_t mka=498.0;

   TProfile* npee=new TProfile("npee","",1000,-0.5,990.5); // npe vs momentum for electrons
   TProfile* npemu=new TProfile("npemu","",1000,-0.5,990.5); // npe vs momentum for muons
   TProfile* npepi=new TProfile("npepi","",1000,-0.5,990.5); // npe vs momentum for pions

   TProfile* Ksigma=new TProfile("Ksigma","",250,-0.5,990.5); // sigma vs momentum for pions

   //TH2F* null2=new TH2F("null2","",200,0,1000,55,-0.05,3.505); // for n=1.05 and t=6cm
   TH2F* null2=new TH2F("null2","",200,0,1000,55,-0.05,9.05); // for n=1.05 and t=6cm
   null2->SetXTitle("P, #frac{MeV}{c}");
   null2->SetYTitle("N_{pe}");

   TH2F* null3=new TH2F("null3","",200,0,1000,55,-0.05,7.5); 
   null3->SetXTitle("P, #frac{MeV}{c}");
   null3->SetYTitle("#sigma");

   for (int i=0; i<1000; i++)
   {
      const int ncounts = 30;
      int N0ee[ncounts]={0}, N0pi[ncounts]={0}, Ntotal=0;   
      float Npethr[ncounts] = {0.};

      for (int l=0; l<1000; l++) 
      { 
         float mu_e = gRandom->Poisson(nphe11->Eval((float)i/me));
         float mu_pi = gRandom->Poisson(nphe11->Eval((float)i/mpi));
         float mu_mu = gRandom->Poisson(nphe11->Eval((float)i/mmu));
         
         npee->Fill((float)i,mu_e);
         npemu->Fill((float)i,mu_mu);
         npepi->Fill((float)i,mu_pi);
         
         for( int k=0; k<ncounts; k++ )
         {
            Npethr[k] = k*0.2;
            if( mu_e<Npethr[k] ) N0ee[k]++;
            if( mu_pi<Npethr[k] ) N0pi[k]++;
         }
         Ntotal++;   
      }
      float underthr_eff_pi[ncounts];
      float eff[ncounts];   
   
      //for( int k=0; k<ncounts; k++ )
      for( int k=14; k<15; k++ )
      {
         underthr_eff_pi[k] = (float)N0pi[k]/(float)Ntotal;
         eff[k] = (float)N0ee[k]/(float)Ntotal;   
         if(i<300) cout<<i<<"\t"<<Ntotal<<"\t"<<Npethr[k]<<"\t"<<N0ee[k]<<"\t"<<N0pi[k]<<endl;
      }

      //if(i<300) cout<<i<<"\t"<<Ntotal<<"\t"<<N0ee<<"\t"<<N0pi<<endl;
      //if(i<300) cout<<i<<"\t"<<underthr_eff_pi<<"\t"<<eff<<endl;   

      //Ksigma->Fill((float)i,abs(sqrt(2.)*(TMath::ErfInverse(1-2*underthr_eff_pi)+TMath::ErfInverse(1-2*(1-eff)))));
      Ksigma->Fill((float)i,abs(sqrt(2.)*(TMath::ErfInverse(1-2*underthr_eff_pi[14])+TMath::ErfInverse(1-2*(1-eff[14])))));
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

   TGraphErrors* g1 = new TGraphErrors(mom.size(), &mom[0], &curve1[0], 0, 0);
   g1->SetMarkerStyle(21);
   g1->SetMarkerSize(0.8);
   g1->SetMarkerColor(8);
   g1->GetYaxis()->SetRangeUser(0, 3);
   g1->Draw("same");

   TLegend* leg1=new TLegend(100,1.2,200,1.9,"","");
   TLegendEntry *le1=leg1->AddEntry(npee,"e","elp");
   le1->SetTextColor(kRed);
   TLegendEntry *le2=leg1->AddEntry(npemu,"#mu","elp");
   le2->SetTextColor(kGreen);
   TLegendEntry *le3=leg1->AddEntry(npepi,"#pi","elp");
   le3->SetTextColor(kBlue);
   leg1->Draw("same");

   gPad->Update();

   c1->SaveAs("res.png");

   TCanvas *c2 = new TCanvas();
   c2->cd();   

   null3->Draw();

   Ksigma->SetLineColor(kRed);
   Ksigma->SetMarkerColor(kRed);
   Ksigma->Draw("same");

   gPad->Update();

   c2->SaveAs("resKsigma.png");
  
   /*
   TGraphErrors* g2 = new TGraphErrors(mom.size(), &mom[0], &curve2[0], 0, 0);
   g2->SetMarkerStyle(21);
   g2->SetMarkerSize(0.8);
   g2->SetMarkerColor(4);
   g2->GetYaxis()->SetRangeUser(0, 3);
   
   TGraphErrors* g3 = new TGraphErrors(mom.size(), &mom[0], &curve3[0], 0, 0);
   g3->SetMarkerStyle(21);
   g3->SetMarkerSize(0.8);
   g3->SetMarkerColor(2);
   g3->GetYaxis()->SetRangeUser(0, 3);

   TGraphErrors* g4 = new TGraphErrors(mom.size(), &mom[0], &curve4[0], 0, 0);
   g4->SetMarkerStyle(21);
   g4->SetMarkerSize(0.8);
   g4->SetMarkerColor(4);
   g4->GetYaxis()->SetRangeUser(0, 3);

   TGraphErrors* g5 = new TGraphErrors(mom.size(), &mom[0], &curve5[0], 0, 0);
   g5->SetMarkerStyle(21);
   g5->SetMarkerSize(0.8);
   g5->SetMarkerColor(4);
   g5->GetYaxis()->SetRangeUser(0, 3);
        
   //g1->Draw("AP");
   //g2->Draw("AP");
   //g3->Draw("AP");
   //g4->Draw("AP");
   //g5->Draw("AP");
  
   TMultiGraph *mg = new TMultiGraph();
   mg->Add(g1);
   mg->Add(g2);
   mg->Add(g3);
   //mg->Add(g4);
   //mg->Add(g5);
   mg->SetTitle("; P, MeV/c; A, ph.e.");
   mg->GetYaxis()->SetRangeUser(0, 3);
   mg->Draw("AP");
   //mg->Write("res");;
   c1->SaveAs("res.png");
    

   TCanvas *c2 = new TCanvas();
   c2->cd();   

   TGraphErrors* g6 = new TGraphErrors(mom.size(), &mom[0], &Ksigma_mu[0], 0, 0);
   g6->SetMarkerStyle(21);
   g6->SetMarkerSize(0.8);
   g6->SetMarkerColor(8);
   g6->GetYaxis()->SetRangeUser(0, 5);
   //g6->Draw("AP");

   TGraphErrors* g7 = new TGraphErrors(mom.size(), &mom[0], &Ksigma_pi[0], 0, 0);
   g7->SetMarkerStyle(21);
   g7->SetMarkerSize(0.8);
   g7->SetMarkerColor(4);
   g7->GetYaxis()->SetRangeUser(0, 5);

   TGraphErrors* g8 = new TGraphErrors(Ksigma_sipm.size(), &mom[0], &Ksigma_sipm[0], 0, 0);
   g8->SetMarkerStyle(24);
   g8->SetMarkerSize(0.8);
   g8->SetMarkerColor(4);
   g8->GetYaxis()->SetRangeUser(0, 5);

   TMultiGraph *mg2 = new TMultiGraph();
   mg2->Add(g6);
   mg2->Add(g7);
   mg2->Add(g8);

   mg2->SetTitle("; P, MeV/c; #sigma");
   mg2->GetYaxis()->SetRangeUser(0, 6);
   mg2->Draw("AP");
   c2->SaveAs("res_Ksigma.png");
   */
} 
