void makeOtherInputFile(){


  int const N_bins=15;

  TString out_fileName="signal_WV_noRatio_mu.root";

  TFile * f_histo=new TFile("mu_boosted.root","read");
  f_histo->cd();
  TH1D * histo_SM=(TH1D*)(diboson->Clone("name_h1"));


  TFile * f_par=new TFile("ATGC_shape_coefficients.root","read");
  f_par->cd();
  // TProfile2D
  TProfile2D * h2d_lam_dk_par0=(TProfile2D*)(p0_lambda_dkg->Clone("name_p0_lam_dk"));
  TProfile2D * h2d_lam_dk_par1=(TProfile2D*)(p1_lambda_dkg->Clone("name_p1_lam_dk"));
  TProfile2D * h2d_lam_dk_par2=(TProfile2D*)(p2_lambda_dkg->Clone("name_p2_lam_dk"));
  TProfile2D * h2d_lam_dk_par3=(TProfile2D*)(p3_lambda_dkg->Clone("name_p3_lam_dk"));
  TProfile2D * h2d_lam_dk_par4=(TProfile2D*)(p4_lambda_dkg->Clone("name_p4_lam_dk"));
  TProfile2D * h2d_lam_dk_par5=(TProfile2D*)(p5_lambda_dkg->Clone("name_p5_lam_dk"));
  TProfile2D * h2d_lam_dk_par6=(TProfile2D*)(p6_lambda_dkg->Clone("name_p6_lam_dk"));
  
  TProfile2D * h2d_lam_dg_par0=(TProfile2D*)(p0_lambda_dg1->Clone("name_p0_lam_dg"));
  TProfile2D * h2d_lam_dg_par1=(TProfile2D*)(p1_lambda_dg1->Clone("name_p1_lam_dg"));
  TProfile2D * h2d_lam_dg_par2=(TProfile2D*)(p2_lambda_dg1->Clone("name_p2_lam_dg"));
  TProfile2D * h2d_lam_dg_par3=(TProfile2D*)(p3_lambda_dg1->Clone("name_p3_lam_dg"));
  TProfile2D * h2d_lam_dg_par4=(TProfile2D*)(p4_lambda_dg1->Clone("name_p4_lam_dg"));
  TProfile2D * h2d_lam_dg_par5=(TProfile2D*)(p5_lambda_dg1->Clone("name_p5_lam_dg"));
  TProfile2D * h2d_lam_dg_par6=(TProfile2D*)(p6_lambda_dg1->Clone("name_p6_lam_dg"));
  
  TProfile2D * h2d_dk_dg_par0=(TProfile2D*)(p0_dkg_dg1->Clone("name_p0_dk_dg"));
  TProfile2D * h2d_dk_dg_par1=(TProfile2D*)(p1_dkg_dg1->Clone("name_p1_dk_dg"));
  TProfile2D * h2d_dk_dg_par2=(TProfile2D*)(p2_dkg_dg1->Clone("name_p2_dk_dg"));
  TProfile2D * h2d_dk_dg_par3=(TProfile2D*)(p3_dkg_dg1->Clone("name_p3_dk_dg"));
  TProfile2D * h2d_dk_dg_par4=(TProfile2D*)(p4_dkg_dg1->Clone("name_p4_dk_dg"));
  TProfile2D * h2d_dk_dg_par5=(TProfile2D*)(p5_dkg_dg1->Clone("name_p5_dk_dg"));
  TProfile2D * h2d_dk_dg_par6=(TProfile2D*)(p6_dkg_dg1->Clone("name_p6_dk_dg"));


  TH2D *binContents_lam_dk[N_bins];
  TH2D *binContents_lam_dg[N_bins];
  TH2D *binContents_dk_dg[N_bins];
  //  TH2D *binContents_lam_dk[15];
  for (int i=0;i<N_bins;i++){

    ///////////////////////////////////////////////////////////////////////////////
    //
    //           2D grid: lambda, dkg
    //
    //

    TString test_lam_dk ="bin_content_lam_dk_";
    test_lam_dk+=i+1;
    test_lam_dk+="_";

    binContents_lam_dk[i]=new TH2D();
    binContents_lam_dk[i]= (TH2D)(h2d_lam_dk_par0->Clone(test_lam_dk));

    double pt=histo_SM->GetBinCenter(i+1);
    double pt_low=histo_SM->GetBinLowEdge(i+1);
    double pt_high=histo_SM->GetBinLowEdge(i+2);

    test_lam_dk+=pt_low;
    test_lam_dk+="to";
    test_lam_dk+=pt_high;
    test_lam_dk+="GeV";
    cout << "name: "<<test_lam_dk<< endl;


    double yield_SM=histo_SM->GetBinContent(i+1);

    int l_N=h2d_lam_dk_par0->GetXaxis()->GetNbins();
    int dk_N=h2d_lam_dk_par0->GetYaxis()->GetNbins();

    for (int dk_i=1;dk_i<dk_N+1 ;dk_i++){
      for (int l_i=1;l_i<l_N+1 ;l_i++){

	double l=h2d_lam_dk_par0->GetXaxis()->GetBinCenter(l_i);
	double dk=h2d_lam_dk_par0->GetYaxis()->GetBinCenter(dk_i);
	//    int l=0.;
	
	int binX=h2d_lam_dk_par0->GetXaxis()->FindBin(l) ;
	int binY=h2d_lam_dk_par0->GetYaxis()->FindBin(dk) ;
	double p0=h2d_lam_dk_par0->GetBinContent(binX,binY);
	double p1=h2d_lam_dk_par1->GetBinContent(binX,binY);
	double p2=h2d_lam_dk_par2->GetBinContent(binX,binY);
	double p3=h2d_lam_dk_par3->GetBinContent(binX,binY);
	double p4=h2d_lam_dk_par4->GetBinContent(binX,binY);
	double p5=h2d_lam_dk_par5->GetBinContent(binX,binY);
	double p6=h2d_lam_dk_par6->GetBinContent(binX,binY);
	
	double yield=yield_SM*(p0+p1*pt+p2*pt*pt+p3*pt*pt*pt+p4*pt*pt*pt*pt +p5*pt*pt*pt*pt*pt+p6*pt*pt*pt*pt*pt*pt);
	binContents_lam_dk[i]->SetBinContent(binX,binY,yield);
	//	cout <<"(dk="<<dk<<" , l="<<l<<" )" << " pt mean: "<< pt<<" SM yield: "<<yield_SM<<" yield: "<<yield << endl;
      }
    }

    binContents_lam_dk[i]->SetName(test_lam_dk);
    binContents_lam_dk[i]->SetTitle(test_lam_dk);
    binContents_lam_dk[i]->GetZaxis()->SetTitle("yield");


    ///////////////////////////////////////////////////////////////////////////////
    //
    //           2D grid: lambda, dg
    //
    //

    TString test_lam_dg ="bin_content_lam_dg_";
    test_lam_dg+=i+1;
    test_lam_dg+="_";

    binContents_lam_dg[i]=new TH2D();
    binContents_lam_dg[i]= (TH2D)(h2d_lam_dg_par0->Clone(test_lam_dg));

    double pt=histo_SM->GetBinCenter(i+1);
    double pt_low=histo_SM->GetBinLowEdge(i+1);
    double pt_high=histo_SM->GetBinLowEdge(i+2);

    test_lam_dg+=pt_low;
    test_lam_dg+="to";
    test_lam_dg+=pt_high;
    test_lam_dg+="GeV";
    cout << "name: "<<test_lam_dg<< endl;


    double yield_SM=histo_SM->GetBinContent(i+1);

    int l_N=h2d_lam_dg_par0->GetXaxis()->GetNbins();
    int dg_N=h2d_lam_dg_par0->GetYaxis()->GetNbins();

    for (int dg_i=1;dg_i<dg_N+1 ;dg_i++){
      for (int l_i=1;l_i<l_N+1 ;l_i++){

	double l=h2d_lam_dg_par0->GetXaxis()->GetBinCenter(l_i);
	double dg=h2d_lam_dg_par0->GetYaxis()->GetBinCenter(dg_i);
	//    int l=0.;
	
	int binX=h2d_lam_dg_par0->GetXaxis()->FindBin(l) ;
	int binY=h2d_lam_dg_par0->GetYaxis()->FindBin(dg) ;
	double p0=h2d_lam_dg_par0->GetBinContent(binX,binY);
	double p1=h2d_lam_dg_par1->GetBinContent(binX,binY);
	double p2=h2d_lam_dg_par2->GetBinContent(binX,binY);
	double p3=h2d_lam_dg_par3->GetBinContent(binX,binY);
	double p4=h2d_lam_dg_par4->GetBinContent(binX,binY);
	double p5=h2d_lam_dg_par5->GetBinContent(binX,binY);
	double p6=h2d_lam_dg_par6->GetBinContent(binX,binY);
	
	double yield=yield_SM*(p0+p1*pt+p2*pt*pt+p3*pt*pt*pt+p4*pt*pt*pt*pt +p5*pt*pt*pt*pt*pt+p6*pt*pt*pt*pt*pt*pt);
	binContents_lam_dg[i]->SetBinContent(binX,binY,yield);
	//	cout <<"(dg="<<dg<<" , l="<<l<<" )" << " pt mean: "<< pt<<" SM yield: "<<yield_SM<<" yield: "<<yield << endl;
      }
    }

    binContents_lam_dg[i]->SetName(test_lam_dg);
    binContents_lam_dg[i]->SetTitle(test_lam_dg);
    binContents_lam_dg[i]->GetZaxis()->SetTitle("yield");

    ///////////////////////////////////////////////////////////////////////////////
    //
    //           2D grid: dk, dg
    //
    //

    TString test_dk_dg ="bin_content_dk_dg_";
    test_dk_dg+=i+1;
    test_dk_dg+="_";

    binContents_dk_dg[i]=new TH2D();
    binContents_dk_dg[i]= (TH2D)(h2d_dk_dg_par0->Clone(test_dk_dg));

    double pt=histo_SM->GetBinCenter(i+1);
    double pt_low=histo_SM->GetBinLowEdge(i+1);
    double pt_high=histo_SM->GetBinLowEdge(i+2);

    test_dk_dg+=pt_low;
    test_dk_dg+="to";
    test_dk_dg+=pt_high;
    test_dk_dg+="GeV";
    cout << "name: "<<test_dk_dg<< endl;


    double yield_SM=histo_SM->GetBinContent(i+1);

    int dk_N=h2d_dk_dg_par0->GetXaxis()->GetNbins();
    int dg_N=h2d_dk_dg_par0->GetYaxis()->GetNbins();

    for (int dg_i=1;dg_i<dg_N+1 ;dg_i++){
      for (int dk_i=1;dk_i<l_N+1 ;dk_i++){

	double dk=h2d_dk_dg_par0->GetXaxis()->GetBinCenter(dk_i);
	double dg=h2d_dk_dg_par0->GetYaxis()->GetBinCenter(dg_i);
	//    int l=0.;
	
	int binX=h2d_dk_dg_par0->GetXaxis()->FindBin(dk) ;
	int binY=h2d_dk_dg_par0->GetYaxis()->FindBin(dg) ;
	double p0=h2d_dk_dg_par0->GetBinContent(binX,binY);
	double p1=h2d_dk_dg_par1->GetBinContent(binX,binY);
	double p2=h2d_dk_dg_par2->GetBinContent(binX,binY);
	double p3=h2d_dk_dg_par3->GetBinContent(binX,binY);
	double p4=h2d_dk_dg_par4->GetBinContent(binX,binY);
	double p5=h2d_dk_dg_par5->GetBinContent(binX,binY);
	double p6=h2d_dk_dg_par6->GetBinContent(binX,binY);
	
	double yield=yield_SM*(p0+p1*pt+p2*pt*pt+p3*pt*pt*pt+p4*pt*pt*pt*pt +p5*pt*pt*pt*pt*pt+p6*pt*pt*pt*pt*pt*pt);
	binContents_dk_dg[i]->SetBinContent(binX,binY,yield);
	//	cout <<"(dg="<<dg<<" , l="<<l<<" )" << " pt mean: "<< pt<<" SM yield: "<<yield_SM<<" yield: "<<yield << endl;
      }
    }

    binContents_dk_dg[i]->SetName(test_dk_dg);
    binContents_dk_dg[i]->SetTitle(test_dk_dg);
    binContents_dk_dg[i]->GetZaxis()->SetTitle("yield");

  }

  TCanvas * c1=new TCanvas();
  c1->Divide(1,2);
  c1->cd(1);
  binContents_lam_dk[1]->Draw("lego");
  c1->cd(2);
  binContents_lam_dk[14]->Draw("lego");

  cout <<"compare SM and calculated yield lam-dk:" << endl;
  for (int i=0;i<N_bins;i++){
    cout <<histo_SM->GetBinContent(i+1)<<" "<<binContents_lam_dk[i]->GetBinContent(binContents_lam_dk[i]->FindBin(0.,0.))  <<endl;
  }
  cout <<"compare SM and calculated yield lam-dg:" << endl;
  for (int i=0;i<N_bins;i++){
    cout <<histo_SM->GetBinContent(i+1)<<" "<<binContents_lam_dg[i]->GetBinContent(binContents_lam_dg[i]->FindBin(0.,0.))  <<endl;
  }
  cout <<"compare SM and calculated yield dk-dg:" << endl;
  for (int i=0;i<N_bins;i++){
    cout <<histo_SM->GetBinContent(i+1)<<" "<<binContents_dk_dg[i]->GetBinContent(binContents_dk_dg[i]->FindBin(0.,0.))  <<endl;
  }

  TFile *fileout=new TFile(out_fileName,"recreate");
  fileout->cd();
  for (int i=0;i<N_bins;i++){
    binContents_lam_dk[i]->Write();
    binContents_lam_dg[i]->Write();
    binContents_dk_dg[i]->Write();
  }


}
