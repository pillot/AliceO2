#include <fstream>
#include <vector>
#include <algorithm>
#include <functional>


#include "Validation.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHMappingInterface/CathodeSegmentation.h"
#include "MCHBase/Digit.h"
#include "MCHPreClustering/PreClusterBlock.h"
#include "MCHPreClustering/PreClusterFinder.h"
#include "DigitsFileReader.h"
#include "MCHClustering/ClusteringForTest.h"

#include "TMath.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TF2.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include <cmath>

using namespace o2::mch;
using namespace std;

namespace o2 {

namespace mch {

Validation::Validation(){};

Double_t myMathieson2D(Double_t *x, Double_t *par){
    Double_t pitch = 0.25;
    Float_t xx = x[0]-par[2];
    Float_t yy = x[1]-par[3];
    Double_t K2x = M_PI*0.5*(1-0.5*sqrt(par[0]));
    Double_t K1x = K2x * sqrt(par[0]) / 4 / atan(sqrt(par[0]));
    Double_t fx = K1x * ((1-pow(tanh(K2x*xx/pitch),2))/(1+ par[0]*pow(tanh(K2x*xx/pitch),2)));
    Double_t K2y = M_PI*0.5*(1-0.5*sqrt(par[1]));
    Double_t K1y = K2y * sqrt(par[1]) / 4 / atan(sqrt(par[1]));
    Double_t fy = K1y * ((1-pow(tanh(K2y*yy/pitch),2))/(1+ par[1]*pow(tanh(K2y*yy/pitch),2)));
    Double_t f = fx*fy;
    return f;
}

void myMath(){
//    Double_t ( Validation::*func ) ( Double_t*, Double_t* );
//    func = &Validation::myMathieson2D;
//    TF2 *f1 = new TF2("myMath",(this->*func),0,10,0,10,4);

    TF2 *f1 = new TF2("myMath",myMathieson2D,-10,10,-10,10,4);
    f1->SetParameters(0.5085, 0.5840, 0., 0.);
    f1->SetParNames("K3x", "K3y", "Mean x", "Mean y");
    //f1->Draw("colz");
}

void Validation::PlotMathieson2D(){
    
    digits.clear();
    
    TH2F* hb(NULL);
    TH2F* hnb(NULL);
    myMath();
    
    o2::mch::mapping::CathodeSegmentation catsegb(809, kTRUE);
    o2::mch::mapping::CathodeSegmentation catsegnb(809, kFALSE);

    int nopadsb = catsegb.nofPads();
    int nopadsnb = catsegnb.nofPads();
    int padid;
    
    
    //Conversion des vecteurs en arrays bending
    
    double xlowsb[lowxsb.size()];
    double ylowsb[lowysb.size()];
    std::copy(lowxsb.begin(), lowxsb.end(), xlowsb);
    std::copy(lowysb.begin(), lowysb.end(), ylowsb);
    
    //Conversion des vecteurs en arrays non-bending
    
    double xlowsnb[lowxsnb.size()];
    double ylowsnb[lowysnb.size()];
    std::copy(lowxsnb.begin(), lowxsnb.end(), xlowsnb);
    std::copy(lowysnb.begin(), lowysnb.end(), ylowsnb);
    
    //Création et remplissage histogrammes bending
    
    cout << "Generating histograms bending and non-bending..."<< endl;
    
    hb = new TH2F("hb","hb",lowxsb.size()-1,xlowsb,lowysb.size()-1,ylowsb);
    gRandom->SetSeed(0);
    hb->FillRandom("myMath", /*nsamples=*/2000);
    
    TCanvas *c2 = new TCanvas("c2", "The road to PlotDorado",
                              200,10,600,400);
    c2->Divide(1,2);
    c2->cd(1);
    
    hb->SetTitle("Mathieson 2D Distribution - Bending plane - DE809");
    hb->Draw("LEGO2");
    
    c2->Update();
    c2->Draw();
    
    TCanvas *c3 = new TCanvas("c3", "The road to PlotDorado",
    200,10,600,400);
    c3->Divide(1,2);
    c3->cd(1);
    
    hb->Draw("colz");

    c3->Update();
    c3->Draw();
    
    //Création et remplissage histogrammes non-bending
    
    hnb = new TH2F("hnb","hnb",lowxsnb.size()-1,xlowsnb,lowysnb.size()-1,ylowsnb);
    hnb->FillRandom("myMath", /*nsamples=*/2000);
    
//    TCanvas *c21 = new TCanvas("c21", "The road to PlotDorado",
//                              200,10,600,400);
    
    c2->cd(2);
    hnb->SetTitle("Mathieson 2D Distribution - Non-Bending plane - DE809");
    hnb->Draw("LEGO2");

    c2->Update();
    c2->Draw();
//
//    TCanvas *c31 = new TCanvas("c31", "The road to PlotDorado",
//    200,10,600,400);
    
    c3->cd(2);
    hnb->Draw("colz");

    c3->Update();
    c3->Draw();
    
    
    int nbbinsb = (lowxsb.size()-1)*(lowysb.size()-1);
    int nbbinsnb = (lowxsnb.size()-1)*(lowysnb.size()-1);
    double bincontent;
    double binxcent;
    double binycent;
    double charge[nopadsb + nopadsnb];
    
    for(int i=0; i<nopadsb + nopadsnb; i++ ){
        charge[i] = 0.;
    }
    
    cout << "Number of bending bins = " << lowxsb.size()-1 << " x " << lowysb.size()-1 << " = " << nbbinsb << endl;
    cout << "Number of non-bending bins = " << lowxsnb.size()-1 << " x " << lowysnb.size()-1 << " = " << nbbinsnb << endl;
    
    
    //Remplissage digits bending plane
    
    cout << "Filling digits from bending plane..." << endl;
    
    for(int binx = 1; binx<lowxsb.size(); binx++){
        for(int biny = 1; biny<lowysb.size(); biny++){
            bincontent = hb->GetBinContent(binx,biny);
            binxcent = (hb->GetXaxis())->GetBinCenter(binx);
            binycent = (hb->GetYaxis())->GetBinCenter(biny);
            
//            if((binx*biny)%4 == 0){
//                cout << "Content of bending bin (binx= " << binx << ", biny= " << biny << "): " << bincontent << endl;
//                cout << "Center of bending bin (binx= " << binx << ", biny= " << biny << "): (" << binxcent << "," << binycent << ")" << endl;
//            }
            
            padid = catsegb.findPadByPosition(binxcent,binycent);
            charge[padid] += bincontent;
            
//            if(charge[padid] > 2){
//               cout << "Content of bending bin (binx= " << binx << ", biny= " << biny << "): " << bincontent << endl;
//               cout << "Center of bending bin (binx= " << binx << ", biny= " << biny << "): (" << binxcent << "," << binycent << ")" << endl;
//               cout << "Padid: " << padid << "  Charge: " << charge[padid] <<endl;
//           }
        }
    }
    
    //Remplissage digits non-bending plane
    
    cout << "Filling digits from non-bending plane..." << endl;
        
        for(int binx = 1; binx<lowxsnb.size(); binx++){
            for(int biny = 1; biny<lowysnb.size(); biny++){
                bincontent = hnb->GetBinContent(binx,biny);
                binxcent = (hnb->GetXaxis())->GetBinCenter(binx);
                binycent = (hnb->GetYaxis())->GetBinCenter(biny);
                
                padid = catsegnb.findPadByPosition(binxcent,binycent) + nopadsb;
                charge[padid] += bincontent;
                
//                if(charge[padid] > 2){
//                    cout << "Content of non-bending bin (binx= " << binx << ", biny= " << biny << "): " << bincontent << endl;
//                    cout << "Center of non-bending bin (binx= " << binx << ", biny= " << biny << "): (" << binxcent << "," << binycent << ")" << endl;
//                    cout << "Padid: " << padid << "  Charge: " << charge[padid] <<endl;
//                }
            }
        }
    
    cout << "Nombre de bins bending : " << nbbinsb << " Nombre de pads bending : " << nopadsb << endl;
    cout << "Nombre de bins non-bending : " << nbbinsnb << " Nombre de pads non-bending : " << nopadsnb << endl;
    cout << "IL Y A UN PB SUR LE NON BENDING (MERGING DE BINS A FAIRE). On a mergé les bins artificiellement lors de la formation des digits." << endl;
    
    for(int i=0; i<nopadsb + nopadsnb; i++){
        if(charge[i] > 2){  //Couper le bruit type

          //  cout << "Un digit est en cours de création..." << endl;

            digits.push_back( std::make_unique<mch::Digit>() );
            mch::Digit* digit = digits.back().get();

               digit->setDetID(809);
               digit->setPadID(i);
               digit->setADC(charge[i]);

            cout << "Digit padid " << digit->getPadID() << " et ADC " << digit->getADC() << endl;
        }
    }
    
    
}

void Validation::InfoDE809b(){
    
    int detElemId =809;
    bool isbending = kTRUE;
    o2::mch::mapping::CathodeSegmentation catseg(detElemId, isbending);

    int nopads = catseg.nofPads();
    
    double padposx;
    double padposy;
    double padsizex;
    double padsizey;
    double lowerx;
    double lowery;
    lowxsb.clear();
    lowysb.clear();
    int valuesx = 0;
    int valuesy = 0;
    
    
    cout << "There are " << nopads << " pads on bending plane" << endl;
    
    for(int catPadindex = 0; catPadindex<nopads; catPadindex++){
        padposx = catseg.padPositionX(catPadindex);
        padposy = catseg.padPositionY(catPadindex);
        padsizex = catseg.padSizeX(catPadindex);
        padsizey = catseg.padSizeY(catPadindex);
        lowerx = padposx - (padsizex/2);
        lowery = padposy - (padsizey/2);
        
        if(!(find(lowxsb.begin(), lowxsb.end(), lowerx) != lowxsb.end())){
            lowxsb.push_back(lowerx);
            valuesx++;
        }
        if(!(find(lowysb.begin(), lowysb.end(), lowery) != lowysb.end())){
            lowysb.push_back(lowery);
            valuesy++;
        }
        
    //    if(catPadindex%10==0){
    //        cout << "catPadindex " << catPadindex << "\n"
    //        << "Position: (" << padposx << "," << padposy
    //        << ")\n"
    //        << "Size: (" << padsizex << "," << padsizey
    //        << ")\n";
    //        << "LowerLimits: (" << lowerx << "," << lowery
    //        << ")" << endl;
    //    }
        
    }
    
    cout << "Il doit y avoir " << valuesx << " valeurs dans le vecteur lowxsb"
    << "et " << valuesy << " valeurs dans le vecteur lowysb" << endl;
    cout << "Il y a " << lowxsb.size() << " valeurs dans le vecteur lowxsb"
    << "et " << lowysb.size() << " valeurs dans le vecteur lowysb" << endl;
    
    cout << "Rangement des bords de pads par ordre croissant ..." << endl;
    
    sort(lowxsb.begin(), lowxsb.end());
    sort(lowysb.begin(), lowysb.end());
    
    cout << "Suppression des doublons par erreur d'arrondi ..." << endl;
    
    auto lastxb = std::unique(lowxsb.begin(), lowxsb.end(),
    [](double l, double r) { return std::abs(l - r) < 0.001; });
    lowxsb.erase(lastxb, lowxsb.end());
    
    auto lastyb = std::unique(lowysb.begin(), lowysb.end(),
    [](double l, double r) { return std::abs(l - r) < 0.001; });
    lowysb.erase(lastyb, lowysb.end());
    
    cout << "Il y a maintenant " << lowxsb.size() << " valeurs dans le vecteur lowxsb"
    << " et " << lowysb.size() << " valeurs dans le vecteur lowysb" << endl;
    
    
    cout << "Ajout de la borne sup" << endl;
    lowxsb.push_back(2*lowxsb[lowxsb.size()-1]-lowxsb[lowxsb.size()-2]);
    lowysb.push_back(2*lowysb[lowysb.size()-1]-lowysb[lowysb.size()-2]);
    
    cout << "xlowsb:" << lowxsb[0] << "," << lowxsb[1] << "..." << lowxsb[valuesx-1] << "," << lowxsb[valuesx] << endl;
    cout << "ylowsb:" << lowysb[0] << "," << lowysb[1] << "..." << lowysb[valuesy-1] << "," << lowysb[valuesy] << endl;
    
//    for(int i=0; i<lowxsb.size(); i++){
//        cout << lowxsb[i] << endl;
//    }
//    for(int i=0; i<lowysb.size(); i++){
//        cout << lowysb[i] << endl;
//    }
    

//    int b, nb;
//    bool found = seg.findPadPairByPosition(x, y, b, nb);
//    int nopads = seg.nofPads();
//    int nods = seg.nofDualSampas();
//
//    if (seg.isValid(b)) {
//    std::cout << "There is a bending pad at position " << x << "," << y << "\n"
//    << " which belongs to dualSampa " << seg.padDualSampaId(b)
//    << " and has a PadId " << b
//    << " and has a x-size of " << seg.padSizeX(b) << " cm\n"
//    << " and has a y-size of " << seg.padSizeY(b) << " cm\n"
//    << " and has a x-position of " << seg.padPositionX(b) << " cm\n"
//    << " and has a y-position of " << seg.padPositionY(b) << " cm\n"
//    << "Il y a" << nopads << "pads au total sur les deux faces"
//        << "et" << nods << "DualSampas";
//    }

    //assert(b == seg.findPadByFEE(76, 9));
    
}

void Validation::InfoDE809nb(){
    
    int detElemId =809;
    bool isbending = kFALSE;
    o2::mch::mapping::CathodeSegmentation catseg(detElemId, isbending);

    int nopads = catseg.nofPads();
    
    double padposx;
    double padposy;
    double padsizex;
    double padsizey;
    double lowerx;
    double lowery;
    lowxsnb.clear();
    lowysnb.clear();
    int valuesx = 0;
    int valuesy = 0;
    
    
    cout << "There are " << nopads << " pads on non-bending plane" << endl;
    
    for(int catPadindex = 0; catPadindex<nopads; catPadindex++){
        padposx = catseg.padPositionX(catPadindex);
        padposy = catseg.padPositionY(catPadindex);
        padsizex = catseg.padSizeX(catPadindex);
        padsizey = catseg.padSizeY(catPadindex);
        lowerx = padposx - (padsizex/2);
        lowery = padposy - (padsizey/2);
        
        if(!(find(lowxsnb.begin(), lowxsnb.end(), lowerx) != lowxsnb.end())){
            lowxsnb.push_back(lowerx);
            valuesx++;
        }
        if(!(find(lowysnb.begin(), lowysnb.end(), lowery) != lowysnb.end())){
            lowysnb.push_back(lowery);
            valuesy++;
        }
        
//        if(catPadindex%10==0){
//            cout << "catPadindex " << catPadindex << "\n"
//            << "Position: (" << padposx << "," << padposy
//            << ")\n"
//            << "Size: (" << padsizex << "," << padsizey
//            << ")\n"
//            << "LowerLimits: (" << lowerx << "," << lowery
//            << ")" << endl;
//        }
        
    }
    
    cout << "Il doit y avoir " << valuesx << " valeurs dans le vecteur lowxsnb"
    << "et " << valuesy << " valeurs dans le vecteur lowysnb" << endl;
    cout << "Il y a " << lowxsnb.size() << " valeurs dans le vecteur lowxsnb"
    << "et " << lowysnb.size() << " valeurs dans le vecteur lowysnb" << endl;
    
    cout << "Rangement des bords de pads par ordre croissant ..." << endl;
    
    sort(lowxsnb.begin(), lowxsnb.end());
    sort(lowysnb.begin(), lowysnb.end());
    
    cout << "Suppression des doublons par erreur d'arrondi ..." << endl;
    
    auto lastxnb = std::unique(lowxsnb.begin(), lowxsnb.end(),
    [](double l, double r) { return std::abs(l - r) < 0.001; });
    lowxsnb.erase(lastxnb, lowxsnb.end());
    
    auto lastynb = std::unique(lowysnb.begin(), lowysnb.end(),
    [](double l, double r) { return std::abs(l - r) < 0.001; });
    lowysnb.erase(lastynb, lowysnb.end());
    
    
    cout << "Il y a maintenant " << lowxsnb.size() << " valeurs dans le vecteur lowxsnb"
    << " et " << lowysnb.size() << " valeurs dans le vecteur lowysnb" << endl;
    
    cout << "Ajout de la borne sup" << endl;
    lowxsnb.push_back(2*lowxsnb[lowxsnb.size()-1]-lowxsnb[lowxsnb.size()-2]);
    lowysnb.push_back(2*lowysnb[lowysnb.size()-1]-lowysnb[lowysnb.size()-2]);
    
//    cout << "xlowsnb:" << lowxsnb[0] << "," << lowxsnb[1] << "..." << lowxsnb[lowxsnb.size()-2] << "," << lowxsnb[lowxsnb.size()-1] << endl;
//    cout << "ylowsnb:" << lowysnb[0] << "," << lowysnb[1] << "..." << lowysnb[lowysnb.size()-2] << "," << lowysnb[lowysnb.size()-1] << endl;
    
//    for(int i=0; i<lowxsnb.size(); i++){
//        cout << lowxsnb[i] << endl;
//    }
//    for(int i=0; i<lowysnb.size(); i++){
//        cout << lowysnb[i] << endl;
//    }
    
}

void Validation::TestClustering(){
    
    cout << "Filling buffer of digits..." << endl;
    
    nDigits = getNumberOfDigits();
    digitsBuffer = (mch::Digit*)realloc(digitsBuffer, sizeof(mch::Digit) * nDigits);
    storeDigits(digitsBuffer);
    
      mch::PreClusterFinder preClusterFinder;
      mch::PreClusterBlock preClusterBlock;
      mch::Clustering clustering;
      
      std::string fname;
      preClusterFinder.init(fname);

      char* preClustersBuffer = NULL;
      std::vector<mch::Clustering::Cluster> clusters(0);

        // load the digits from the memory buffer and run the pre-clustering phase
        preClusterFinder.reset();
        preClusterFinder.loadDigits(digitsBuffer, nDigits);
        preClusterFinder.run();

        // get number of pre-clusters and store them into a memory buffer
        auto preClustersSize = preClusterBlock.getPreClustersBufferSize(preClusterFinder);
        printf("preClustersSize: %d\n", (int)preClustersSize);
        preClustersBuffer = (char*)realloc(preClustersBuffer, preClustersSize);
        preClusterBlock.storePreClusters(preClusterFinder, preClustersBuffer);

        //continue;
        printf("\n\n==========\nReading clusters\n\n");

        std::vector<mch::PreClusterStruct> preClusters;
        preClusterBlock.readPreClusters(preClusters, preClustersBuffer, preClustersSize);
          
          printf("\n\n==========\nRunning Clustering\n\n");
          
        //Runs the clustering of preClusters following a CenterOfGravity algorithm. Fills clusters.
        clustering.runFinderCOG(preClusters, clusters);
        printf("Number of clusters obtained and saved: %lu\n", clusters.size());
          
      // clustering.runFinderSimpleFit(preClusters, clusters);
          
     //     clustering.runFinderGaussianFit(preClusters, clusters);
}

ssize_t Validation::getNumberOfDigits()
{
  return digits.size();
}


void Validation::storeDigits(void* bufferPtr)
{
  mch::Digit* ptr = (mch::Digit*)bufferPtr;
  for(unsigned int di = 0; di < digits.size(); di++) {

    memcpy(ptr, digits[di].get(), sizeof(mch::Digit));
    ptr += 1;
  }
}

}
}
