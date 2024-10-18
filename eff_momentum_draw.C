void eff_momentum_draw()
{
   float x1, x2, y1, y2;
   vector< float > mom, mom_sipm, mom1, mom2, curve1, curve2;
   vector< float > Ksigma, Ksigma_sipm;   
   float eff;
   float underthr_eff;

   std::ifstream f_in1("from_paper_effpi_n1.13.csv");
   while( !f_in1.eof() )
   {
	  std::string s1;
     while( getline(f_in1, s1) )
     {
         if( s1.size()>=1 && s1[0]!='#' )
         {
            std::istringstream i_str1(s1);
            i_str1 >> x1 >> y1;
            mom1.push_back(x1);
            curve1.push_back(y1);           
         }
      }
   }
   f_in1.close();


   std::ifstream f_in2("from_paper_effK_n1.13.csv");
   while( !f_in2.eof() )
   {
     std::string s2;
     while( getline(f_in2, s2) )
     {
         if( s2.size()>=1 && s2[0]!='#' )
         {
            std::istringstream i_str2(s2);
            i_str2 >> x2 >> y2;
            mom2.push_back(x2);
            curve2.push_back(y2);
         }
      }
   }
   f_in2.close();


   for(int i=25; i<900; i=i+50)
   {
      eff = 0;
      underthr_eff = 0;
      mom.push_back(i);
      float aver_curve1 = 0;
      int counts1 = 0;
      for (int k = 0; k < mom1.size(); k++) {
         if( mom1[k]>i-25 && mom1[k]<i+25 ){
            aver_curve1 += curve1[k];
            counts1++;   
         }   
      }
      eff = aver_curve1/counts1;

      float aver_curve2 = 0;
      int counts2 = 0;
      for (int k = 0; k < mom2.size(); k++) {
         if( mom2[k]>i-50 && mom2[k]<i+50 ){
            aver_curve2 += curve2[k];
            counts2++;   
         }   
      }
      underthr_eff = aver_curve2/counts2;

      cout<<i<<"\t"<<eff<<"\t"<<underthr_eff<<endl;

      Ksigma.push_back(abs(sqrt(2.)*(TMath::ErfInverse(1-2*underthr_eff)+TMath::ErfInverse(1-2*(1-eff)))));

      if(i>400) 
      {
         mom_sipm.push_back(i);   
         //Ksigma_sipm.push_back(6.4);  //p2 - u=53.5V, t=15C
         //Ksigma_sipm.push_back(5.1);  //p9 - u=53.5V, t=15C
         //Ksigma_sipm.push_back(5.4);  //p2 - u=57V, t=45C
         Ksigma_sipm.push_back(3.9);  //p9 - u=57V, t=45C
      }
   }

   
   TCanvas *c1 = new TCanvas();
   c1->cd();	
   
   TGraphErrors* g1 = new TGraphErrors(mom1.size(), &mom1[0], &curve1[0], 0, 0);
   g1->SetMarkerStyle(21);
   g1->SetMarkerSize(0.8);
   g1->SetMarkerColor(8);
   g1->GetYaxis()->SetRangeUser(0, 1.1);

   
   TGraphErrors* g2 = new TGraphErrors(mom2.size(), &mom2[0], &curve2[0], 0, 0);
   g2->SetMarkerStyle(21);
   g2->SetMarkerSize(0.8);
   g2->SetMarkerColor(4);
   g2->GetYaxis()->SetRangeUser(0, 1.1);
   
    
   //g1->Draw("AP");
   //g2->Draw("AP");
   //g2->Draw("SAME*");
   
   TMultiGraph *mg = new TMultiGraph();
   mg->Add(g1);
   mg->Add(g2);
   mg->SetTitle("; P, MeV/c; #varepsilon");
   mg->GetYaxis()->SetRangeUser(0, 1.1);
   mg->GetXaxis()->SetRangeUser(0, 900);
   mg->Draw("AP");
   
   c1->SaveAs("res.png");
    
   
   TCanvas *c2 = new TCanvas();
   c2->cd();   

   TGraphErrors* g6 = new TGraphErrors(mom.size(), &mom[0], &Ksigma[0], 0, 0);
   g6->SetMarkerStyle(21);
   g6->SetMarkerSize(0.8);
   g6->SetMarkerColor(4);
   g6->GetYaxis()->SetRangeUser(0, 5);
   //g6->Draw("AP");

   TGraphErrors* g7 = new TGraphErrors(mom_sipm.size(), &mom_sipm[0], &Ksigma_sipm[0], 0, 0);
   g7->SetMarkerStyle(24);
   g7->SetMarkerSize(0.8);
   g7->SetMarkerColor(4);
   g7->GetYaxis()->SetRangeUser(0, 5);
   
   TMultiGraph *mg2 = new TMultiGraph();
   mg2->Add(g6);
   mg2->Add(g7);
   
   mg2->SetTitle("; P, MeV/c; #sigma");
   mg2->GetYaxis()->SetRangeUser(0, 7.0);
   mg2->Draw("AP");
   c2->SaveAs("res_Ksigma.png");
    
} 
