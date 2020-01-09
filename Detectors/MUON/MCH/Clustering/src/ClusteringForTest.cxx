//
//  ClusteringForTest.cxx
//
//  Algo de clustring repris de AliRoot, ClusterCOG
//
//  Created by PERRIN Sébastien on 31/10/2019.
//

#include "MCHClustering/ClusteringForTest.h"
#include "MCHBase/Digit.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHPreClustering/PreClusterFinder.h"
#include "MCHPreClustering/PreClusterBlock.h"
//#include "TVirtualFitter.h"
//#include "TObjArray.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"

//#include "TVector2.h"
#include <chrono>
#include <memory>
#include <stdexcept>
#include <vector>
#include <TMath.h>
#include <stdio.h>
#include <math.h>

#include <fairmq/Tools.h>
#include <FairMQLogger.h>

namespace o2
{
namespace mch
{

using namespace std;

//_________________________________________________________________________________________________
void Clustering::runFinderCOG(std::vector<PreClusterStruct>& preClusters, std::vector<Cluster>& clusters)
{
    
    cout << "preClusters.size():" << preClusters.size() << endl;
    int sizedigit;
    int jete = 0;
    Digit digittmp;
    Cluster clustertmp;
    
    for (Int_t i = 0; i < preClusters.size(); i++)   //On boucle sur les preclusters identifiés
    {
        if(preClusters[i].nDigits == 1){         // Je vire les preclusters qui n'ont qu'un seul digit
            jete++;
            continue;
        }
        cout << "Precluster avec " << preClusters[i].nDigits << " digits." << endl;  //Nombre de digits dans le precluster aue l'on regarde
        
        std::vector<Digit> preclustertmp;
        const Digit* ptrdigit = preClusters[i].digits;
        
        for (int j=0; j<preClusters[i].nDigits; j++){          //Je récupère les digits et les mets dans un   vecteur de digits, preclustertmp
            digittmp = *ptrdigit;
            cout << "Digit numero:" << j << endl;
            cout << "digittmp.getADC():" << digittmp.getADC() << endl;
            cout << "digittmp.getDetID():" << digittmp.getDetID() << endl;
            cout << "digittmp.getPadID():" << digittmp.getPadID() << endl;
//            cout << "type digittmp:" << typeid(digittmp).name() << endl;
//            cout << "type ptrdigit:" << typeid(ptrdigit).name() << endl;
//            cout << "ptrdigit:" << ptrdigit << endl;
            preclustertmp.push_back(digittmp);
            sizedigit = sizeof(digittmp);
            ptrdigit++;
        }
        
        //Je clusterise mon vecteur de digits, preclustertmp et l'ajoute au vecteur de clusters, nommé clusters.
        
       // cout << "preclustertmp[0].getADC():" << preclustertmp[0].getADC() << endl;
        clustertmp = FinderCOG(preclustertmp);
      //  cout << "COORD X DU CLUSTER AJOUTÉ A CLUSTERS: " << clustertmp.getx() << endl;
        clusters.push_back(clustertmp);
    }
    
    cout << "Nombre de preclusters au départ:" << preClusters.size() << endl;
    cout << "Nombre de preclusters de taille 1, abandonnés:" << jete << endl;
    
}

//_________________________________________________________________________________________________
Clustering::Cluster Clustering::FinderCOG(std::vector<Digit> &precluster)
{
    cout << "\n\n==========\nRunning COG Algorithm\n\n" << endl;

    Double_t xmin = 1E9;
    Double_t ymin = 1E9;
    Double_t xmax = -1E9;
    Double_t ymax = -1E9;
    Double_t charge[2] = { 0.0, 0.0 };
    Double_t multiplicity[2] = { 0.0, 0.0 };
    Cluster cluster;
    
    Double_t x[] = { 0.0, 0.0 };
    Double_t y[] = { 0.0, 0.0 };
    
    Double_t xsize[] = { 0.0, 0.0 };
    Double_t ysize[] = { 0.0, 0.0 };
    
         
    for ( Int_t i = 0; i < precluster.size(); ++i ) //On boucle sur les digits de notre precluster
    {
       int detid = precluster[i].getDetID();
       int padid = precluster[i].getPadID();
       cout << "\nDetID:" << detid << " PadID:" << padid << endl;
       mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);
       
       
       
       //On définit le vecteur position et taille du pad en question
     std::vector<Double_t> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
     std::vector<Double_t> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};
       
       cout << "PadPosition: "<< padPosition[0] << " " << padPosition[1] << endl;
       cout << "PadSize: "<< padSize[0] << " " << padSize[1] << endl;
       
       //On met à jour les xminmax et yminmax
     xmin = TMath::Min(padPosition[0]-0.5*padSize[0],xmin);
     xmax = TMath::Max(padPosition[0]+0.5*padSize[0],xmax);
     ymin = TMath::Min(padPosition[1]-0.5*padSize[1],ymin);
     ymax = TMath::Max(padPosition[1]+0.5*padSize[1],ymax);
        
        for ( Int_t cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
         {
           
             if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
             {
                 cout << "Cathode:" << cathode << endl;
                 
                 //On ajoute un terme au numerateur du COG et de la précision
               x[cathode] += padPosition[0]*precluster[i].getADC();
               y[cathode] += padPosition[1]*precluster[i].getADC();
               xsize[cathode] += padSize[0];
               ysize[cathode] += padSize[1];
                 
                 //On ajoute un terme au dénominateur du COG et de la précision
               charge[cathode] += precluster[i].getADC();
               multiplicity[cathode] += 1;
             }
             
         }
                
    }
           
    for ( Int_t cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
    {
//        cout << "\nCathode:" << cathode << endl;
//        cout << "x:" << x[cathode] << endl;
//        cout << "y:" << y[cathode] << endl;
//        cout << "Charge:" << charge[cathode] << endl;
//        cout << "xsize:" << xsize[cathode] << endl;
//        cout << "ysize:" << ysize[cathode] << endl;
//        cout << "Multiplicity:" << multiplicity[cathode] << endl;
        
                 //On fait les calculs de COG et précision
               if ( charge[cathode] != 0 )
               {
                 x[cathode] /= charge[cathode];
                 y[cathode] /= charge[cathode];
               }
               if ( multiplicity[cathode] != 0 )
               {
                 xsize[cathode] /= multiplicity[cathode];
                 ysize[cathode] /= multiplicity[cathode];
               }
        cout << "\nPost-calcul cathode:" << cathode << endl;
        cout << "x:" << x[cathode] << endl;
        cout << "y:" << y[cathode] << endl;
        cout << "Charge:" << charge[cathode] << endl;
        cout << "xsize:" << xsize[cathode] << endl;
        cout << "ysize:" << ysize[cathode] << endl;
        cout << "Multiplicity:" << multiplicity[cathode] << endl;
    }
    

  Double_t xCOG = 0.0;
  Double_t yCOG = 0.0;
  Double_t ex = 0.0;
  Double_t ey = 0.0;

  // On regarde sur quelle cathode la précision est la meilleure pour définir x et y du COG.
  xCOG = ( xsize[0] < xsize[1] ) ? x[0] : x[1];
  yCOG = ( ysize[0] < ysize[1] ) ? y[0] : y[1];
    
    ex = ( xsize[0] < xsize[1] ) ? xsize[0] : xsize[1];
    ey = ( ysize[0] < ysize[1] ) ? ysize[0] : ysize[1];
    
//    cout << "ex:" << ex << endl;
//    cout << "ey:" << ey << endl;
    
    double timestamp = precluster[0].getTimeStamp();  // On récupère le timestamp du premier digit du precluster
//
//    // On print en console les infos importantes pour pouvoir checker les problèmes
//    printf("\n\nCluster multiplicity %ld (x,y)=(%e,%e) boundaries=(xmin,ymin,xmax,ymax)=(%e,%e,%e,%e)"
//                  " (x0,y0,x1,y1)=(%e,%e,%e,%e) timestamp = %lf \n\n",
//                  precluster.size(),xCOG,yCOG,xmin,ymin,xmax,ymax,
//                  x[0],y[0],x[1],y[1],timestamp);
//
    //On remplit un cluster avec les infos nécessaires
    cluster.setx(xCOG);
    cluster.sety(yCOG);
    cluster.setex(ex);
    cluster.setey(ey);
    cluster.settimestamp(timestamp);    // En première approx. on peut affecter le timestamp du pre;ier digit au cluster
    
    printf("\n\nCluster multiplicity %ld \n(x,y)=(%e,%e) \nboundaries=(xmin,ymin,xmax,ymax)=(%e,%e,%e,%e) \nprecision=(ex,ey)=(%e,%e) \ntimestamp = %lf \n\n",precluster.size(),cluster.getx(),cluster.gety(),xmin,ymin,xmax,ymax,cluster.getex(),cluster.getey(),cluster.gettimestamp());
    
    return cluster;
    
    
//     PRINCIPE GENERAL COG
//     Pour chaque cathode
//        Pour chaque pad dans le precluster
//            Mettre à jour les valeurs des limites du cluster
//            Ajouter les termes au numerateur de x,y COG et au calcul de la precision
//        Diviser par les denominateurs
//     Obtenir les parametres du cluster
     
    
}

//_________________________________________________________________________________________________

















//_________________________________________________________________________________________________
void Clustering::runFinderSimpleFit(std::vector<PreClusterStruct>& preClusters, std::vector<Cluster>& clusters)
{
    cout << "preClusters.size():" << preClusters.size() << endl;
    int sizedigit;
    int jete = 0;
    Digit digittmp;
    Cluster clustertmpCOG;
    Cluster clustertmp;
    
    for (Int_t i = 0; i < preClusters.size(); i++)   //On boucle sur les preclusters identifiés
    {
        if(preClusters[i].nDigits == 1){         // Je vire les preclusters qui n'ont qu'un seul digit
            jete++;
            continue;
        }
        cout << "Precluster avec " << preClusters[i].nDigits << " digits." << endl;  //Nombre de digits dans le precluster aue l'on regarde
        
        std::vector<Digit> preclustertmp;
        const Digit* ptrdigit = preClusters[i].digits;
        
        for (int j=0; j<preClusters[i].nDigits; j++){          //Je récupère les digits et les mets dans un   vecteur de digits, preclustertmp
            digittmp = *ptrdigit;
            cout << "Digit numero:" << j << endl;
            cout << "digittmp.getADC():" << digittmp.getADC() << endl;
            cout << "digittmp.getDetID():" << digittmp.getDetID() << endl;
            cout << "digittmp.getPadID():" << digittmp.getPadID() << endl;
//            cout << "type digittmp:" << typeid(digittmp).name() << endl;
//            cout << "type ptrdigit:" << typeid(ptrdigit).name() << endl;
//            cout << "ptrdigit:" << ptrdigit << endl;
            preclustertmp.push_back(digittmp);
            sizedigit = sizeof(digittmp);
            ptrdigit++;
        }
        
        //Je clusterise mon vecteur de digits, preclustertmp et l'ajoute au vecteur de clusters, nommé clusters.
        
       // cout << "preclustertmp[0].getADC():" << preclustertmp[0].getADC() << endl;
        clustertmpCOG = FinderCOG(preclustertmp);
      //  cout << "COORD X DU CLUSTER AJOUTÉ A CLUSTERS: " << clustertmp.getx() << endl;
        if(preClusters.size() < 3){
            clustertmp = clustertmpCOG;
        }
        else{
        clustertmp = FinderSimpleFit(preclustertmp, clustertmpCOG);
         //   clustertmp = ComputePosition(preclustertmp, clustertmpCOG);
        }
        
        clusters.push_back(clustertmp);
    }
    
    cout << "Nombre de preclusters au départ:" << preClusters.size() << endl;
    cout << "Nombre de preclusters de taille 1, abandonnés:" << jete << endl;
    
}


//_________________________________________________________________________________________________
Clustering::Cluster Clustering::FinderSimpleFit(std::vector<Digit> &precluster, Clustering::Cluster clustertmpCOG)
{
    cout << "\n\n==========\nRunning SimpleFit Algorithm\n\n" << endl;
    
    Cluster cluster;

    //Création des vecteurs position du hit sur chaque cqthode
    Double_t xhit[2] = { clustertmpCOG.getx(), clustertmpCOG.getx() };
    Double_t yhit[2] = { clustertmpCOG.gety(), clustertmpCOG.gety() };
    
  //  cout << "clustertmpCOG.getx()" << clustertmpCOG.getx() << endl;
    
    //Création du float porteur de la coordonnée la plus précise
    Double_t xhitfinal = 0;
    Double_t yhitfinal = 0;
    
    //Création des vecteurs de charge totale, multiplicité et tailles sur chaque cathode
    Double_t chargetot[2] = { 0.0, 0.0 };
    Double_t multiplicity[2] = { 0.0, 0.0 };
    Double_t xsize[2] = { 0.0, 0.0 };
    Double_t ysize[2] = { 0.0, 0.0 };
    
    //Détermination des vecteurs de charge totale, multiplicité et taille qui devient précision
    for ( Int_t i = 0; i < precluster.size(); ++i ){
        int detid = precluster[i].getDetID();
        int padid = precluster[i].getPadID();
        mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);
            
        for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
        {
                
            if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
                    {
                        std::vector<Double_t> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
                        std::vector<Double_t> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};

                        chargetot[cathode] += precluster[i].getADC();
                        multiplicity[cathode] += 1;
                        xsize[cathode] += padSize[0];
                        ysize[cathode] += padSize[1];
                    }
        }
    }
    for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
    {
        if ( multiplicity[cathode] != 0 )
        {
          xsize[cathode] /= multiplicity[cathode];
          ysize[cathode] /= multiplicity[cathode];
        }
    }
    
    //Minimisation du Chi2 en 2D sur chaque cathode
    for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
    {
    
        Double_t Kx3 = 0.5085;
        Double_t Ky3 = 0.5840;
        
        Double_t pasxhit = 0.2;
        Double_t pasyhit = 0.2;
        
        int compteurboucle = 0;
        
        Double_t ChiDeux = Chi2Mathieson(cathode, chargetot[cathode], precluster, clustertmpCOG, xhit[cathode], yhit[cathode], Kx3, Ky3);
        Double_t OldChiDeux = 1;
        Double_t NewChiDeux = 100;
        Double_t deltaxhit = 100;
        Double_t deltayhit = 100;
        Double_t limitstrict = 0.0001;
        Double_t limitlean = 0.005;
        Double_t limitx = 0.;
        Double_t limity = 0.;
        
        if(cathode == 0){
            limitx = limitstrict;
            limity = limitlean;
        }
        else{
            limity = limitstrict;
            limitx = limitlean;
        }
        
        while ((sqrt(pow((NewChiDeux-OldChiDeux)/OldChiDeux,2)) > 0.001)||(abs(deltaxhit)>limitx)||(abs(deltayhit)>limity)){
            compteurboucle++;
            
            //Iteration and improvement of estimate
            Double_t dChidxhit = (Chi2Mathieson(cathode, chargetot[cathode], precluster, clustertmpCOG, xhit[cathode]+pasxhit, yhit[cathode], Kx3, Ky3)-Chi2Mathieson(cathode, chargetot[cathode], precluster, clustertmpCOG, xhit[cathode], yhit[cathode], Kx3, Ky3))/pasxhit;
            Double_t dChidyhit = (Chi2Mathieson(cathode, chargetot[cathode], precluster, clustertmpCOG, xhit[cathode], yhit[cathode]+pasyhit, Kx3, Ky3)-Chi2Mathieson(cathode, chargetot[cathode], precluster, clustertmpCOG, xhit[cathode], yhit[cathode], Kx3, Ky3))/pasyhit;
            
//            float dChidKx3 = (Chi2Mathieson(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], Kx3+pasKx3, Ky3)-Chi2Mathieson(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], Kx3, Ky3))/pasKx3;


//            float dChidKy3 = (Chi2Mathieson(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], Kx3, Ky3+pasKy3)-Chi2Mathieson(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], Kx3, Ky3))/pasKy3;

            deltaxhit = pasxhit*dChidxhit*1000;
            deltayhit = pasyhit*dChidyhit*1000;
            xhit[cathode] = xhit[cathode] - pasxhit*dChidxhit*1000;
            yhit[cathode] = yhit[cathode] - pasyhit*dChidyhit*1000;
//            Kx3 = Kx3 - pasKx3*dChidKx3;
//            Ky3 = Ky3 - pasKy3*dChidKy3;
            
            
            // Update ChiDeux
            OldChiDeux = ChiDeux;
            ChiDeux = Chi2Mathieson(cathode, chargetot[cathode], precluster, clustertmpCOG, xhit[cathode], yhit[cathode], Kx3, Ky3);
            NewChiDeux = ChiDeux;
            cout << "\nBoucle While numero " << compteurboucle << endl;
            cout << "xhit: " << xhit[cathode] << endl;
            cout << "yhit: " << yhit[cathode] << endl;
            cout << "deltaxhit: " << deltaxhit << endl;
            cout << "deltayhit: " << deltayhit << endl;
            cout << "OldChiDeux: " << OldChiDeux << endl;
            cout << "NewChiDeux: " << NewChiDeux << endl;
            
            if(compteurboucle > 1000){
                cout << "Déjà 1.000 boucles pour trouver Chi2, utilisation résultats COG..." << endl;
                return clustertmpCOG;
            }
        }
        
        cout << "\n\ncathode " << cathode << endl;
        cout << "Boucle Finale numero " << compteurboucle << endl;
        cout << "xhit final: " << xhit[cathode] << endl;
        cout << "yhit final: " << yhit[cathode] << endl;
        cout << "Kx3: " << Kx3 << endl;
        cout << "Ky3: " << Ky3 << endl;
        cout << "Chi2 Final: " << NewChiDeux << endl;
    
    }
    
    
    //Détermination de la position finale du hit à l'aide de la résolution sur chaque cathode
    xhitfinal = ( xsize[0] < xsize[1] ) ? xhit[0] : xhit[1];
    yhitfinal = ( ysize[0] < ysize[1] ) ? yhit[0] : yhit[1];
    Double_t ex = ( xsize[0] < xsize[1] ) ? xsize[0] : xsize[1];
    Double_t ey = ( ysize[0] < ysize[1] ) ? ysize[0] : ysize[1];
    
    double timestamp = precluster[0].getTimeStamp();  // On récupère le timestamp du premier digit du precluster
    
    
    //On remplit un cluster avec les infos nécessaires
    cluster.setx(xhitfinal);
    cluster.sety(yhitfinal);
    cluster.setex(ex);
    cluster.setey(ey);
    cluster.settimestamp(timestamp);    // En première approx. on peut affecter le timestamp du pre;ier digit au cluster
    
            printf("\n\nCluster par simple fit:\nCluster multiplicity %ld \n(x,y)=(%e,%e) \nprecision=(ex,ey)=(%e,%e) \ntimestamp = %lf \n\n",precluster.size(),cluster.getx(),cluster.gety(),cluster.getex(),cluster.getey(),cluster.gettimestamp());
    
    return cluster;

}


//_________________________________________________________________________________________________
Double_t Clustering::Chi2Mathieson(int cathode, Double_t chargetot, std::vector<Digit> &precluster, Clustering::Cluster clustertmp, Double_t xhit, Double_t yhit, Double_t Kx3, Double_t Ky3)
{
    
    Double_t sum = 0;
    for(Int_t i = 0; i < precluster.size(); ++i){
        
        int detid = precluster[i].getDetID();
        int padid = precluster[i].getPadID();
        mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);
        
        if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
        {
              //On définit le vecteur position et taille du pad en question
            std::vector<Double_t> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
            std::vector<Double_t> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};
            
            std::vector<Double_t> lowerLeft = {xhit - padPosition[0] - 0.5*padSize[0], yhit - padPosition[1] - 0.5*padSize[1]};
            std::vector<Double_t> upperRight= {lowerLeft[0] + padSize[0], lowerLeft[1] + padSize[1]};
            
            Double_t qfit = IntMathiesonXY(lowerLeft[0], lowerLeft[1], upperRight[0], upperRight[1], Kx3, Ky3)/chargetot;

            Double_t qtrue = precluster[i].getADC()/chargetot;
            
            Double_t delta = qfit-qtrue;
            
            sum += pow(delta,2);
        }
    }
    
    return sum;
}


//Il faudra un visuel


//_________________________________________________________________________________________________
Double_t Clustering::IntMathiesonXY(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t Kx3, Double_t Ky3)
{
    Double_t pitch = 0.25;                                          //Adapted fron AliMuonMathieson
    Double_t Kx2 = M_PI / 2. * (1. - 0.5 * sqrt(Kx3));
    Double_t cx1 = Kx2 * sqrt(Kx3) / 4. / atan(sqrt(Kx3));
    Double_t Kx4 = cx1 / Kx2 / sqrt(Kx3);
    
    Double_t Ky2 = M_PI / 2. * (1. - 0.5 * sqrt(Ky3));
    Double_t cy1 = Ky2 * sqrt(Ky3) / 4. / atan(sqrt(Ky3));
    Double_t Ky4 = cy1 / Ky2 / sqrt(Ky3);
    
    x1 /= pitch;
    x2 /= pitch;
    y1 /= pitch;
    y2 /= pitch;

    Double_t ux1 = sqrt(Kx3)*tanh(Kx2*x1);
    Double_t ux2 = sqrt(Kx3)*tanh(Kx2*x2);
    
    Double_t uy1 = sqrt(Ky3)*tanh(Ky2*y1);
    Double_t uy2 = sqrt(Ky3)*tanh(Ky2*y2);
    
    return Double_t(4.*Kx4*(atan(ux2)-atan(ux1))*
                   Ky4*(atan(uy2)-atan(uy1)));
}

//_________________________________________________________________________________________________



















//_________________________________________________________________________________________________

double RosenBrock(const double *xx )
{
  const Double_t x = xx[0];
  const Double_t y = xx[1];
  const Double_t tmp1 = y-x*x;
  const Double_t tmp2 = 1-x;
  return 100*tmp1*tmp1+tmp2*tmp2;
}

int NumericalMinimization(const char * minName = "Minuit2",
                          const char *algoName = "" ,
                          int randomSeed = -1)
{
   // create minimizer giving a name and a name (optionally) for the specific
   // algorithm
   // possible choices are:
   //     minName                  algoName
   // Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
   //  Minuit2                     Fumili2
   //  Fumili
   //  GSLMultiMin                ConjugateFR, ConjugatePR, BFGS,
   //                              BFGS2, SteepestDescent
   //  GSLMultiFit
   //   GSLSimAn
   //   Genetic
   ROOT::Math::Minimizer* minimum =
      ROOT::Math::Factory::CreateMinimizer(minName, algoName);
   // set tolerance , etc...
   minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
   minimum->SetMaxIterations(10000);  // for GSL
   minimum->SetTolerance(0.001);
   minimum->SetPrintLevel(1);
   // create function wrapper for minimizer
   // a IMultiGenFunction type
   ROOT::Math::Functor f(&RosenBrock,2);
   double step[2] = {0.01,0.01};
   // starting point
   double variable[2] = { -1.,1.2};
   if (randomSeed >= 0) {
      TRandom2 r(randomSeed);
      variable[0] = r.Uniform(-20,20);
      variable[1] = r.Uniform(-20,20);
   }
   minimum->SetFunction(f);
   // Set the free variables to be minimized !
   minimum->SetVariable(0,"x",variable[0], step[0]);
   minimum->SetVariable(1,"y",variable[1], step[1]);
   // do the minimization
   minimum->Minimize();
   const double *xs = minimum->X();
   std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): "
             << minimum->MinValue()  << std::endl;
   // expected minimum is 0
   if ( minimum->MinValue()  < 1.E-4  && f(xs) < 1.E-4)
      std::cout << "Minimizer " << minName << " - " << algoName
                << "   converged to the right minimum" << std::endl;
   else {
      std::cout << "Minimizer " << minName << " - " << algoName
                << "   failed to converge !!!" << std::endl;
      Error("NumericalMinimization","fail to converge");
   }
   return 0;
}
//_________________________________________________________________________________________________
























//_________________________________________________________________________________________________
//void Clustering::FitFunction(Int_t& /*notused*/, Double_t* /*notused*/,
//                            Double_t& sum, Double_t* par,
//                             Int_t /*notused*/)
//  {
//    /// Chi2 Function to minimize: Mathieson charge distribution in 2 dimensions
//
//    TObjArray* userObjects = static_cast<TObjArray*>(TVirtualFitter::GetFitter()->GetObjectFit());
//
//    std::vector<Digit> precluster = static_cast<std::vector<Digit>>(userObjects->At(0));
//    int cathode = static_cast<int>(userObjects->At(1));
//    Double_t chargetot = static_cast<Double_t>(userObjects->At(2));
//    Double_t Kx3 = static_cast<Double_t>(userObjects->At(3));
//    Double_t Ky3 = static_cast<Double_t>(userObjects->At(4));
//
//    sum = 0.0;
//
//    for(Int_t i = 0; i < precluster.size(); ++i){
//
//    int detid = precluster[i].getDetID();
//    int padid = precluster[i].getPadID();
//    mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);
//
//    if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
//    {
//          //On définit le vecteur position et taille du pad en question
//        std::vector<Double_t> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
//        std::vector<Double_t> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};
//
//        std::vector<Double_t> lowerLeft = {par[0] - padPosition[0] - 0.5*padSize[0], par[1] - padPosition[1] - 0.5*padSize[1]};
//        std::vector<Double_t> upperRight= {lowerLeft[0] + padSize[0], lowerLeft[1] + padSize[1]};
//
//        Double_t qfit = IntMathiesonXY(lowerLeft[0], lowerLeft[1], upperRight[0], upperRight[1], Kx3, Ky3)/chargetot;
//
//        Double_t qtrue = precluster[i].getADC()/chargetot;
//
//        Double_t delta = qfit-qtrue;
//
//        sum += pow(delta,2);
//    }
//  }
//}
//
////_________________________________________________________________________________________________
//Clustering::Cluster Clustering::ComputePosition(std::vector<Digit> &precluster, Clustering::Cluster clustertmpCOG)
//{
//  /// Compute the position of the given cluster, by fitting a Mathieson
//  /// charge distribution to it
//
//  TVirtualFitter* fitter = TVirtualFitter::Fitter(0,2);
//  fitter->SetFCN(FitFunction);
//
//  cout << "\n\n==========\nRunning MinuitMathieson Algorithm\n\n" << endl;
//
//  Cluster cluster;
//
//  //Création des vecteurs position du hit sur chaque cqthode
//  Double_t xhit[2] = { clustertmpCOG.getx(), clustertmpCOG.getx() };
//  Double_t yhit[2] = { clustertmpCOG.gety(), clustertmpCOG.gety() };
//
//    //Création des vecteurs de charge totale, multiplicité et tailles sur chaque cathode
//       Double_t chargetot[2] = { 0.0, 0.0 };
//       Double_t multiplicity[2] = { 0.0, 0.0 };
//       Double_t xsize[2] = { 0.0, 0.0 };
//       Double_t ysize[2] = { 0.0, 0.0 };
//
//   //Détermination des vecteurs de charge totale, multiplicité et taille qui devient précision
//   for ( Int_t i = 0; i < precluster.size(); ++i ){
//       int detid = precluster[i].getDetID();
//       int padid = precluster[i].getPadID();
//       mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);
//
//       for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
//       {
//
//           if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
//                   {
//                       std::vector<Double_t> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
//                       std::vector<Double_t> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};
//
//                       chargetot[cathode] += precluster[i].getADC();
//                       multiplicity[cathode] += 1;
//                   }
//       }
//   }
//
//  Float_t stepX = 0.01; // cm
//  Float_t stepY = 0.01; // cm
//    Double_t Kx3 = 0.5085;
//    Double_t Ky3 = 0.5840;
//
//    Double_t arg(-1); // disable printout
//
//    for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
//    {
//
//        fitter->ExecuteCommand("SET PRINT",&arg,1);
//
//        fitter->SetParameter(0,"cluster X position",xhit[cathode],stepX,0,0);
//        fitter->SetParameter(1,"cluster Y position",yhit[cathode],stepY,0,0);
//
//        TObjArray userObjects;
//
//        userObjects.Add(precluster);
//        userObjects.Add(cathode);
//        userObjects.Add(chargetot[cathode]);
//        userObjects.Add(Kx3);
//        userObjects.Add(Ky3);
//
//        fitter->SetObjectFit(&userObjects);
//
//        Int_t val = fitter->ExecuteCommand("MIGRAD",0,0);
//        if ( val )
//        {
//        // fit failed. Using COG results, with big errors
//        cout << "Fit failed. Using COG results for cluster" << endl;
//            return clustertmpCOG;
//        }
//
//        xhit[cathode] = fitter->GetParameter(0);
//        yhit[cathode] = fitter->GetParameter(1);
//        xsize[cathode] = fitter->GetParError(0);
//        ysize[cathode] = fitter->GetParError(1);
//
//        Double_t amin, edm, errdef;
//        Int_t nvpar, nparx;
//
//        fitter->GetStats(amin, edm, errdef, nvpar, nparx);
//
//        Double_t chi2 = amin;
//
//        printf("Cluster fitted to (x,y)=(%e,%e) (xerr,yerr)=(%e,%e) \n chi2=%e ndf=%d",
//                        xhit[cathode],yhit[cathode],
//                        xsize[cathode],ysize[cathode],chi2,fitter->GetNumberFreeParameters());
//    }
//
//    //Détermination de la position finale du hit à l'aide de la résolution sur chaque cathode
//    Double_t xhitfinal = ( xsize[0] < xsize[1] ) ? xhit[0] : xhit[1];
//    Double_t yhitfinal = ( ysize[0] < ysize[1] ) ? yhit[0] : yhit[1];
//    Double_t ex = ( xsize[0] < xsize[1] ) ? xsize[0] : xsize[1];
//    Double_t ey = ( ysize[0] < ysize[1] ) ? ysize[0] : ysize[1];
//
//    double timestamp = precluster[0].getTimeStamp();  // On récupère le timestamp du premier digit du precluster
//
//
//    //On remplit un cluster avec les infos nécessaires
//    cluster.setx(xhitfinal);
//    cluster.sety(yhitfinal);
//    cluster.setex(ex);
//    cluster.setey(ey);
//    cluster.settimestamp(timestamp);    // En première approx. on peut affecter le timestamp du pre;ier digit au cluster
//
//            printf("\n\nCluster par simple fit:\nCluster multiplicity %ld \n(x,y)=(%e,%e) \nprecision=(ex,ey)=(%e,%e) \ntimestamp = %lf \n\n",precluster.size(),cluster.getx(),cluster.gety(),cluster.getex(),cluster.getey(),cluster.gettimestamp());
//
//    return cluster;
//}

//_________________________________________________________________________________________________











/*

//_________________________________________________________________________________________________
void Clustering::runFinderGaussianFit(std::vector<PreClusterStruct>& preClusters, std::vector<Cluster>& clusters)
{
    cout << "preClusters.size():" << preClusters.size() << endl;
    int sizedigit;
    int jete = 0;
    Digit digittmp;
    Cluster clustertmpCOG;
    Cluster clustertmp;
    
    for (Int_t i = 0; i < preClusters.size(); i++)   //On boucle sur les preclusters identifiés
    {
        if(preClusters[i].nDigits == 1){         // Je vire les preclusters qui n'ont qu'un seul digit
            jete++;
            continue;
        }
        cout << "Precluster avec " << preClusters[i].nDigits << " digits." << endl;  //Nombre de digits dans le precluster aue l'on regarde
        
        std::vector<Digit> preclustertmp;
        const Digit* ptrdigit = preClusters[i].digits;
        
        for (int j=0; j<preClusters[i].nDigits; j++){          //Je récupère les digits et les mets dans un   vecteur de digits, preclustertmp
            digittmp = *ptrdigit;
            cout << "Digit numero:" << j << endl;
            cout << "digittmp.getADC():" << digittmp.getADC() << endl;
            cout << "digittmp.getDetID():" << digittmp.getDetID() << endl;
            cout << "digittmp.getPadID():" << digittmp.getPadID() << endl;
//            cout << "type digittmp:" << typeid(digittmp).name() << endl;
//            cout << "type ptrdigit:" << typeid(ptrdigit).name() << endl;
//            cout << "ptrdigit:" << ptrdigit << endl;
            preclustertmp.push_back(digittmp);
            sizedigit = sizeof(digittmp);
            ptrdigit++;
        }
        
        //Je clusterise mon vecteur de digits, preclustertmp et l'ajoute au vecteur de clusters, nommé clusters.
        
       // cout << "preclustertmp[0].getADC():" << preclustertmp[0].getADC() << endl;
        clustertmpCOG = FinderCOG(preclustertmp);
      //  cout << "COORD X DU CLUSTER AJOUTÉ A CLUSTERS: " << clustertmp.getx() << endl;
        if(preClusters.size() < 3){
            clustertmp = clustertmpCOG;
        }
        else{
        clustertmp = FinderGaussianFit(preclustertmp, clustertmpCOG);
        }
        
        clusters.push_back(clustertmp);
    }
    
    cout << "Nombre de preclusters au départ:" << preClusters.size() << endl;
    cout << "Nombre de preclusters de taille 1, abandonnés:" << jete << endl;
    
}


//_________________________________________________________________________________________________
Clustering::Cluster Clustering::FinderGaussianFit(std::vector<Digit> &precluster, Clustering::Cluster clustertmpCOG)
{
    cout << "\n\n==========\nRunning GaussianFit Algorithm\n\n" << endl;
    
    Cluster cluster;

    //Création des vecteurs position du hit sur chaque cqthode
    Double_t xhit[2] = { clustertmpCOG.getx(), clustertmpCOG.getx() };
    Double_t yhit[2] = { clustertmpCOG.gety(), clustertmpCOG.gety() };
    
  //  cout << "clustertmpCOG.getx()" << clustertmpCOG.getx() << endl;
    
    //Création du float porteur de la coordonnée la plus précise
    Double_t xhitfinal = 0;
    Double_t yhitfinal = 0;
    
    //Création des vecteurs de charge totale, multiplicité et tailles sur chaque cathode
    Double_t chargetot[2] = { 0.0, 0.0 };
    Double_t multiplicity[2] = { 0.0, 0.0 };
    Double_t xsize[2] = { 0.0, 0.0 };
    Double_t ysize[2] = { 0.0, 0.0 };
    
    //Détermination des vecteurs de charge totale, multiplicité et taille qui devient précision
    for ( Int_t i = 0; i < precluster.size(); ++i ){
        int detid = precluster[i].getDetID();
        int padid = precluster[i].getPadID();
        mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);
            
        for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
        {
                
            if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
                    {
                        std::vector<Double_t> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
                        std::vector<Double_t> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};

                        chargetot[cathode] += precluster[i].getADC();
                        multiplicity[cathode] += 1;
                        xsize[cathode] += padSize[0];
                        ysize[cathode] += padSize[1];
                    }
        }
    }
    for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
    {
        if ( multiplicity[cathode] != 0 )
        {
          xsize[cathode] /= multiplicity[cathode];
          ysize[cathode] /= multiplicity[cathode];
        }
    }
    
    //Minimisation du Chi2 en 2D sur chaque cathode
    for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
    {
    
        Double_t sigx = 0.7604;
        Double_t sigy = 0.7808;
        
        Double_t pasxhit = 0.2;
        Double_t pasyhit = 0.2;
        
        int compteurboucle = 0;
        
        Double_t ChiDeux = Chi2Gauss(cathode, chargetot[cathode], precluster, clustertmpCOG, xhit[cathode], yhit[cathode], sigx, sigy);
        Double_t OldChiDeux = 1;
        Double_t NewChiDeux = 100;
        Double_t deltaxhit = 100;
        Double_t deltayhit = 100;
        Double_t limitstrict = 0.0001;
        Double_t limitlean = 0.005;
        Double_t limitx = 0.;
        Double_t limity = 0.;
        
        if(cathode == 0){
            limitx = limitstrict;
            limity = limitlean;
        }
        else{
            limity = limitstrict;
            limitx = limitlean;
        }
        
        while ((sqrt(pow((NewChiDeux-OldChiDeux)/OldChiDeux,2)) > 0.001)||(abs(deltaxhit)>limitx)||(abs(deltayhit)>limity)){
            compteurboucle++;
            
            //Iteration and improvement of estimate
            Double_t dChidxhit = (Chi2Gauss(cathode, chargetot[cathode], precluster, clustertmpCOG, xhit[cathode]+pasxhit, yhit[cathode], sigx, sigy)-Chi2Gauss(cathode, chargetot[cathode], precluster, clustertmpCOG, xhit[cathode], yhit[cathode], sigx, sigy))/pasxhit;
            Double_t dChidyhit = (Chi2Gauss(cathode, chargetot[cathode], precluster, clustertmpCOG, xhit[cathode], yhit[cathode]+pasyhit, sigx, sigy)-Chi2Gauss(cathode, chargetot[cathode], precluster, clustertmpCOG, xhit[cathode], yhit[cathode], sigx, sigy))/pasyhit;
            

            deltaxhit = pasxhit*dChidxhit*1000;
            deltayhit = pasyhit*dChidyhit*1000;
            xhit[cathode] = xhit[cathode] - pasxhit*dChidxhit*1000;
            yhit[cathode] = yhit[cathode] - pasyhit*dChidyhit*1000;
            
            
            // Update ChiDeux
            OldChiDeux = ChiDeux;
            ChiDeux = Chi2Gauss(cathode, chargetot[cathode], precluster, clustertmpCOG, xhit[cathode], yhit[cathode], sigx, sigy);
            NewChiDeux = ChiDeux;
            cout << "\nBoucle While numero " << compteurboucle << endl;
            cout << "xhit: " << xhit[cathode] << endl;
            cout << "yhit: " << yhit[cathode] << endl;
            cout << "deltaxhit: " << deltaxhit << endl;
            cout << "deltayhit: " << deltayhit << endl;
            cout << "OldChiDeux: " << OldChiDeux << endl;
            cout << "NewChiDeux: " << NewChiDeux << endl;
            
            if(compteurboucle > 1000){
                cout << "Déjà 1.000 boucles pour trouver Chi2, utilisation résultats COG..." << endl;
                return clustertmpCOG;
            }
        }
        
        cout << "\n\ncathode " << cathode << endl;
        cout << "Boucle Finale numero " << compteurboucle << endl;
        cout << "xhit final: " << xhit[cathode] << endl;
        cout << "yhit final: " << yhit[cathode] << endl;
        cout << "sigx: " << sigx << endl;
        cout << "sigy: " << sigy << endl;
        cout << "Chi2 Final: " << NewChiDeux << endl;
    
    }
    
    
    //Détermination de la position finale du hit à l'aide de la résolution sur chaque cathode
    xhitfinal = ( xsize[0] < xsize[1] ) ? xhit[0] : xhit[1];
    yhitfinal = ( ysize[0] < ysize[1] ) ? yhit[0] : yhit[1];
    Double_t ex = ( xsize[0] < xsize[1] ) ? xsize[0] : xsize[1];
    Double_t ey = ( ysize[0] < ysize[1] ) ? ysize[0] : ysize[1];
    
    double timestamp = precluster[0].getTimeStamp();  // On récupère le timestamp du premier digit du precluster
    
    
    //On remplit un cluster avec les infos nécessaires
    cluster.setx(xhitfinal);
    cluster.sety(yhitfinal);
    cluster.setex(ex);
    cluster.setey(ey);
    cluster.settimestamp(timestamp);    // En première approx. on peut affecter le timestamp du pre;ier digit au cluster
    
            printf("\n\nCluster par simple fit:\nCluster multiplicity %ld \n(x,y)=(%e,%e) \nprecision=(ex,ey)=(%e,%e) \ntimestamp = %lf \n\n",precluster.size(),cluster.getx(),cluster.gety(),cluster.getex(),cluster.getey(),cluster.gettimestamp());
    
    return cluster;

}

//_________________________________________________________________________________________________
Double_t Clustering::Chi2Gauss(int cathode, Double_t chargetot, std::vector<Digit> &precluster, Clustering::Cluster clustertmp, Double_t xhit, Double_t yhit, Double_t sigx, Double_t sigy)
{
    
    Double_t sum = 0;
    for(Int_t i = 0; i < precluster.size(); ++i){
        
        int detid = precluster[i].getDetID();
        int padid = precluster[i].getPadID();
        mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);
        
        if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
        {
              //On définit le vecteur position et taille du pad en question
            std::vector<double> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
            std::vector<double> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};
            
            std::vector<double> lowerLeft = {xhit - padPosition[0] - 0.5*padSize[0], yhit - padPosition[1] - 0.5*padSize[1]};
            std::vector<double> upperRight= {lowerLeft[0] + padSize[0], lowerLeft[1] + padSize[1]};
            
            Double_t qfit = IntGaussXY(lowerLeft[0], lowerLeft[1], upperRight[0], upperRight[1], sigx, sigy)/chargetot;

            Double_t qtrue = precluster[i].getADC()/chargetot;
            
            Double_t delta = qfit-qtrue;
            
            sum += pow(delta,2);
        }
    }
    
    return sum;
}


//_________________________________________________________________________________________________
Double_t Clustering::IntGaussXY(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t sigx, Double_t sigy)
{
    Double_t pitch = 0.25;
    x1 /= pitch;
    x2 /= pitch;
    y1 /= pitch;
    y2 /= pitch;
    
    Double_t erfx = 0.5*( erf( (x2)/(sigx*sqrt(2)) ) - erf( (x1)/(sigx*sqrt(2)) ) );
    Double_t erfy = 0.5*( erf( (y2)/(sigy*sqrt(2)) ) - erf( (y1)/(sigy*sqrt(2)) ) );
    Double_t integral = erfx*erfy;
    
    return integral;
}


//_________________________________________________________________________________________________
void Clustering::ComparisonMathiesonGauss(std::string& buffer)
{
    Cluster* clusterCOG;
    Cluster* clusterMathieson;
    Cluster* clusterGauss;
    
    clusterCOG = runFinderCOG(buffer);
    clusterMathieson = runFinderSimpleFit(buffer);
    clusterGauss = runFinderGaussianFit(buffer);

}

*/

} // namespace mch
} // namespace o2
