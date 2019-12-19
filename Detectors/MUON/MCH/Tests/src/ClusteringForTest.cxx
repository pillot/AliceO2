//
//  ClusteringForTest.cxx
//
//  Algo de clustring repris de AliRoot, ClusterCOG
//
//  Created by PERRIN Sébastien on 31/10/2019.
//

#include "ClusteringForTest.h"
#include "MCHBase/Digit.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHPreClustering/PreClusterFinder.h"
#include "MCHPreClustering/PreClusterBlock.h"

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
        cout << "preClusters[i].nDigits:" << preClusters[i].nDigits << endl;  //Nombre de digits dans le precluster aue l'on regarde
        
        std::vector<Digit> preclustertmp;
        const Digit* ptrdigit = preClusters[i].digits;
        
        for (int j=0; j<preClusters[i].nDigits; j++){          //Je récupère les digits et les mets dans un   vecteur de digits, preclustertmp
            digittmp = *ptrdigit;
            cout << "j:" << j << endl;
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
        cout << "COORD X DU CLUSTER AJOUTÉ A CLUSTERS: " << clustertmp.getx() << endl;
        clusters.push_back(clustertmp);
    }
    
    cout << "Nombre de preclusters de taille 1, abandonnés:" << jete << endl;
    
}

//_________________________________________________________________________________________________
Clustering::Cluster Clustering::FinderCOG(std::vector<Digit> &precluster)
{

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
    
    for ( Int_t cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
     {
       for ( Int_t i = 0; i < precluster.size(); ++i ) //On boucle sur les digits de notre precluster
       {
           int detid = precluster[i].getDetID();
           int padid = precluster[i].getPadID();
           mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);
           
           //On définit le vecteur position et taille du pad en question
         std::vector<Double_t> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
         std::vector<Double_t> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};
           
           //On met à jour les xminmax et yminmax
         xmin = TMath::Min(padPosition[0]-0.5*padSize[0],xmin);
         xmax = TMath::Max(padPosition[0]+0.5*padSize[0],xmax);
         ymin = TMath::Min(padPosition[1]-0.5*padSize[1],ymin);
         ymax = TMath::Max(padPosition[1]+0.5*padSize[1],ymax);
           
           
         if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
         {
             
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
         
     }
    

  Double_t xCOG = 0.0;
  Double_t yCOG = 0.0;

  // On regarde sur quelle cathode la précision est la meilleure pour définir x et y du COG.
  xCOG = ( xsize[0] < xsize[1] ) ? x[0] : x[1];
  yCOG = ( ysize[0] < ysize[1] ) ? y[0] : y[1];
    
    double timestamp = precluster[0].getTimeStamp();  // On récupère le timestamp du premier digit du precluster
    
    // On print en console les infos importantes pour pouvoir checker les problèmes
    printf("\n\nCluster multiplicity %ld (x,y)=(%e,%e) boundaries=(xmin,ymin,xmax,ymax)=(%e,%e,%e,%e)"
                  " (x0,y0,x1,y1)=(%e,%e,%e,%e) timestamp = %lf \n\n",
                  precluster.size(),xCOG,yCOG,xmin,ymin,xmax,ymax,
                  x[0],y[0],x[1],y[1],timestamp);
    
    //On remplit un cluster avec les infos nécessaires
    cluster.setx(xCOG);
    cluster.sety(yCOG);
    cluster.setex(( xsize[0] < xsize[1] ) ? xsize[0] : xsize[1]);
    cluster.setey(( ysize[0] < ysize[1] ) ? ysize[0] : ysize[1]);
    cluster.settimestamp(timestamp);    // En première approx. on peut affecter le timestamp du pre;ier digit au cluster
    
    printf("\n\nCluster multiplicity %ld (x,y)=(%e,%e) precision=(ex,ey)=(%e,%e)"
    "timestamp = %lf \n\n",
    precluster.size(),cluster.getx(),cluster.gety(),cluster.getex(),cluster.getey(),cluster.gettimestamp());
    
    return cluster;
    
    
//     PRINCIPE GENERAL COG
//     Pour chaque cathode
//        Pour chaque pad dans le precluster
//            Mettre à jour les valeurs des limites du cluster
//            Ajouter les termes au numerateur de x,y COG et au calcul de la precision
//        Diviser par les denominateurs
//     Obtenir les parametres du cluster
     
    
}

/*
//_________________________________________________________________________________________________
Clustering::Cluster* Clustering::runOldFinderSimpleFit(std::vector<Digit> &preclusterslist, std::vector<int> &preclusterssizes)
{
    int index = 0;
    Cluster* clustertmp;
    Cluster* cluster;
    
    for (Int_t i = 0; i < preclusterssizes.size(); i++)
    {
        if(preclusterssizes[i] < 3){         // Je vire les preclusters qui n'ont qu'un seul digit ou deux                                     digits, pas assez pour Mathieson Fit
            index += preclusterssizes[i];
            continue;
        }
        
        std::vector<Digit> preclustertmp;                       //J'extrais les preclusters
        for (int j=index; j<index+preclusterssizes[i]; j++){
            preclustertmp.push_back(preclusterslist[j]);
        }
        
        index += preclusterssizes[i];                           //Je clusterise chaque precluster
        clustertmp = FinderCOG(preclustertmp);
        
        cluster = FinderSimpleFit(preclustertmp, clustertmp);
        return cluster;
    }
    
}

//_________________________________________________________________________________________________
Clustering::Cluster* Clustering::runFinderSimpleFit(std::string& buffer)
{
    int index = 0;
    Cluster* clustertmp;
    Cluster* cluster;
    void* bufferInspectorPointer;                       //Je crée un ptr de type variable
    bufferInspectorPointer = buffer;
    
    while (index < buffer.size)                         //Tant que j'ai pqs fini le buffer
    {
        
        int* ptrint = (int*)bufferInspectorPointer;        //Je récupère la taille avec un ptr int
        int sizetmp = *ptrint;                             //et avance index et ptr
        index += sizeof(*ptrint);
        (int*)bufferInspectorPointer = ptrint++;
        
        Digit* ptrdigit = (Digit*)bufferInspectorPointer; //Je transforme mon ptr en ptr digit
        
        if(sizetmp < 3){                                   //Virer les preclusters de taille < 3
            index += sizetmp*sizeof(*ptrdigit);
            for(int a =0; a<sizetmp; a++){                  //J'avance index et ptr digit
                ptrdigit++;
            }
            (Digit*)bufferInspectorPointer = ptrdigit;     //Je ramène mon ptr originel
            continue;
        }
        
        std::vector<Digit> preclustertmp;
        
        for(int i=0, i < sizetmp, i++){                 //Je remplis le vecteur precluster avec les digits
            preclustertmp.push_back(*ptrdigit);
            index += sizeof(*ptrdigit);
            ptrdigit++;
        }
        (Digit*)bufferInspectorPointer = ptrdigit;
                                                                //Je clusterise chaque precluster
        clustertmp = FinderCOG(preclustertmp);
        
        cluster = FinderSimpleFit(preclustertmp, clustertmp);   //Je fais le fit
        return cluster;
    }
    
}

//_________________________________________________________________________________________________
Clustering::Cluster* Clustering::FinderSimpleFit(std::vector<Digit> &precluster, Clustering::Cluster* clustertmp)
{
    Cluster* cluster;

    //Création des vecteurs position du hit sur chaque cqthode
    float xhit[2] = { clustertmp.getx(), clustertmp.getx() };
    float yhit[2] = { clustertmp.gety(), clustertmp.gety() };
    
    //Création du float porteur de la coordonnée la plus précise
    float xhitfinal = 0;
    float yhitfinal = 0;
    
    //Création des vecteurs de charge totale, multiplicité et tailles sur chaque cathode
    Double_t chargetot[2] = { 0.0, 0.0 };
    Double_t multiplicity[2] = { 0.0, 0.0 };
    Double_t xsize[2] = { 0.0, 0.0 };
    Double_t ysize[2] = { 0.0, 0.0 };
    
    //Détermination des vecteurs de charge totale, multiplicité et taille qui devient précision
    for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
    {
        for ( Int_t i = 0; i < precluster.size(); ++i ){
            int detid = precluster[i].getDetID();
            int padid = precluster[i].getPadID();
            mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);
            
            if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
                    {
                        std::vector<double> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
                        std::vector<double> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};

                        chargetot[cathode] += precluster[i].getADC();
                        multiplicity[cathode] += 1;
                        xsize[cathode] += padSize[0];
                        ysize[cathode] += padSize[1];
                    }
        }
        
        if ( multiplicity[cathode] != 0 )
        {
          xsize[cathode] /= multiplicity[cathode];
          ysize[cathode] /= multiplicity[cathode];
        }
    }
    
    //Minimisation du Chi2 en 2D sur chaque cathode
    for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
    {
    
        float Kx3 = 1;
        float Ky3 = 1;
        
        float pasxhit = 0.001;
        float pasyhit = 0.001;
        float pasKx3 = 0.01;
        float pasKy3 = 0.01;
        
        int compteurboucle = 0;
        
        Double_t ChiDeux = Chi2Mathieson(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], Kx3, Ky3);
        Double_t OldChiDeux = 1;
        Double_t NewChiDeux = 100;
        
        while (sqrt(pow((NewChiDeux-OldChiDeux)/OldChiDeux,2)) > 0.0001){
            compteurboucle++;
            
            //Iteration and improvement of estimate
            float dChidxhit = (Chi2Mathieson(cathode, precluster, clustertmp, xhit[cathode]+pasxhit, yhit[cathode], Kx3, Ky3)-Chi2Mathieson(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], Kx3, Ky3))/pasxhit;
            float dChidyhit = (Chi2Mathieson(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode]+pasyhit, Kx3, Ky3)-Chi2Mathieson(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], Kx3, Ky3))/pasyhit;
            
//            float dChidKx1 = (Chi2Mathieson(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], Kx1+pasKx1, Kx2, Kx3, Ky1, Ky2, Ky3)-Chi2Mathieson(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], Kx1, Kx2, Kx3, Ky1, Ky2, Ky3))/pasKx1;
//            float dChidKx2 = (Chi2Mathieson(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], Kx1, Kx2+pasKx2, Kx3, Ky1, Ky2, Ky3)-Chi2Mathieson(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], Kx1, Kx2, Kx3, Ky1, Ky2, Ky3))/pasKx2;
            float dChidKx3 = (Chi2Mathieson(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], Kx3+pasKx3, Ky3)-Chi2Mathieson(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], Kx3, Ky3))/pasKx3;

//            float dChidKy1 = (Chi2Mathieson(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], Kx1, Kx2, Kx3, Ky1+pasKy1, Ky2, Ky3)-Chi2Mathieson(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], Kx1, Kx2, Kx3, Ky1, Ky2, Ky3))/pasKy1;
//            float dChidKy2 = (Chi2Mathieson(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], Kx1, Kx2, Kx3, Ky1, Ky2+pasKy2, Ky3)-Chi2Mathieson(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], Kx1, Kx2, Kx3, Ky1, Ky2, Ky3))/pasKy2;
            float dChidKy3 = (Chi2Mathieson(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], Kx3, Ky3+pasKy3)-Chi2Mathieson(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], Kx3, Ky3))/pasKy3;

            
            xhit[cathode] = xhit[cathode] - pasxhit*dChidxhit;
            yhit[cathode] = yhit[cathode] - pasyhit*dChidyhit;
            Kx3 = Kx3 - pasKx3*dChidKx3;
            Ky3 = Ky3 - pasKy3*dChidKy3;
            
            
            // Update ChiDeux
            OldChiDeux = ChiDeux;
            ChiDeux = Chi2Mathieson(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], Kx3, Ky3);
            NewChiDeux = ChiDeux;
            cout << "Boucle While numero " << compteurboucle << endl;
            cout << "xhit: " << xhit << endl;
            cout << "yhit: " << yhit << endl;
            cout << "Kx3: " << Kx3 << endl;
            cout << "Ky3: " << Ky3 << endl;
            cout << "OldChiDeux: " << OldChiDeux << endl;
            cout << "NewChiDeux: " << NewChiDeux << endl;
            
            if(compteurboucle > 1000){
                cout << "Déjà 1000 boucles pour trouver Chi2, utilisation résultats COG..." << endl;
                return clustertmp
            }
        }
        
        cout << "cathode " << cathode << endl;
        cout << "Boucle Finale numero " << compteurboucle << endl;
        cout << "xhit final: " << xhit[cathode] << endl;
        cout << "yhit final: " << yhit[cathode] << endl;
        cout << "Kx3 final: " << Kx3 << endl;
        cout << "Ky3 final: " << Ky3 << endl;
        cout << "Chi2 Final: " << NewChiDeux << endl;
    
    }
    
    
    //Détermination de la position finale du hit à l'aide de la résolution sur chaque cathode
    xhitfinal = ( xsize[0] < xsize[1] ) ? xhit[0] : xhit[1];
    yhitfinal = ( ysize[0] < ysize[1] ) ? yhit[0] : yhit[1];
    
    double timestamp = precluster[0].getTimeStamp();  // On récupère le timestamp du premier digit du precluster
    
    
    //On remplit un cluster avec les infos nécessaires
    cluster->setx(xhitfinal);
    cluster->sety(yhitfinal);
    cluster->setex(( xsize[0] < xsize[1] ) ? xsize[0] : xsize[1]);
    cluster->setey(( ysize[0] < ysize[1] ) ? ysize[0] : ysize[1]);
    cluster->settimestamp(timestamp);    // En première approx. on peut affecter le timestamp du pre;ier digit au cluster
    
    return cluster;

}

//_________________________________________________________________________________________________
float Clustering::Chi2Mathieson(int cathode, Double_t chargetot, std::vector<Digit> &precluster, Clustering::Cluster* clustertmp, float xhit, float yhit, float Kx3, float Ky3)
{
    
    float sum = 0;
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
            
            float qfit = IntMathiesonXY(lowerLeft[0], lowerLeft[1], upperRight[0], upperRight[1], Kx3, Ky3)/chargetot[cathode];

            float qtrue = precluster[i].getADC()/chargetot[cathode];
            
            float delta = qfit-qtrue;
            
            sum += pow(delta,2);
        }
    }
}


//Il faudra un visuel


//_________________________________________________________________________________________________
float Clustering::IntMathiesonXY(float x1, float y1, float x2, float y2, float Kx3, float Ky3)
{
    float pitch = 1.0;                                          //Adapted fron AliMuonMathieson
    float Kx2 = M_PI / 2. * (1. - 0.5 * sqrt(Kx3));
    float cx1 = Kx2 * sqrt(Kx3) / 4. / atan(sqrt(Kx3));
    float Kx4 = cx1 / Kx2 / sqrt(Kx3);
    
    float Ky2 = M_PI / 2. * (1. - 0.5 * sqrt(Ky3));
    float cy1 = Ky2 * sqrt(Ky3) / 4. / atan(sqrt(Ky3));
    float Ky4 = cy1 / Ky2 / sqrt(Ky3);
    
    x1 /= pitch;
    x2 /= pitch;
    y1 /= pitch;
    y2 /= pitch;

    Double_t ux1 = sqrt(Kx3)*tanh(Kx2*x1);
    Double_t ux2 = sqrt(Kx3)*tanh(Kx2*x2);
    
    Double_t uy1 = sqrt(Ky3)*tanh(Ky2*y1);
    Double_t uy2 = sqrt(Ky3)*tanh(Ky2*y2);
    
    return Float_t(4.*Kx4*(atan(ux2)-atan(ux1))*
                   Ky4*(atan(uy2)-atan(uy1)));
}



//_________________________________________________________________________________________________
Clustering::Cluster* Clustering::runFinderGaussianFit(std::string& buffer)
{
    int index = 0;
    Cluster* clustertmp;
    Cluster* cluster;
    void* bufferInspectorPointer;                       //Je crée un ptr de type variable
    bufferInspectorPointer = buffer;
    
    while (index < buffer.size)                         //Tant que j'ai pqs fini le buffer
    {
        
        int* ptrint = (int*)bufferInspectorPointer;        //Je récupère la taille avec un ptr int
        int sizetmp = *ptrint;                             //et avance index et ptr
        index += sizeof(*ptrint);
        (int*)bufferInspectorPointer = ptrint++;
        
        Digit* ptrdigit = (Digit*)bufferInspectorPointer; //Je transforme mon ptr en ptr digit
        
        if(sizetmp < 3){                                   //Virer les preclusters de taille < 3
            index += sizetmp*sizeof(*ptrdigit);
            for(int a =0; a<sizetmp; a++){                  //J'avance index et ptr digit
                ptrdigit++;
            }
            (Digit*)bufferInspectorPointer = ptrdigit;     //Je ramène mon ptr originel
            continue;
        }
        
        std::vector<Digit> preclustertmp;
        
        for(int i=0, i < sizetmp, i++){                 //Je remplis le vecteur precluster avec les digits
            preclustertmp.push_back(*ptrdigit);
            index += sizeof(*ptrdigit);
            ptrdigit++;
        }
        (Digit*)bufferInspectorPointer = ptrdigit;
                                                                //Je clusterise chaque precluster
        clustertmp = FinderCOG(preclustertmp);
        
        cluster = FinderGaussianFit(preclustertmp, clustertmp);             //Je fais le fit
        return cluster;
    }
    
}

//_________________________________________________________________________________________________
Clustering::Cluster* Clustering::FinderGaussianFit(std::vector<Digit> &precluster, Clustering::Cluster* clustertmp)
{
    Cluster* cluster;

    //Création des vecteurs position du hit sur chaque cqthode
    float xhit[2] = { clustertmp.getx(), clustertmp.getx() };
    float yhit[2] = { clustertmp.gety(), clustertmp.gety() };
    
    //Création du float porteur de la coordonnée la plus précise
    float xhitfinal = 0;
    float yhitfinal = 0;
    
    //Création des vecteurs de charge totale, multiplicité et tailles sur chaque cathode
    Double_t chargetot[2] = { 0.0, 0.0 };
    Double_t multiplicity[2] = { 0.0, 0.0 };
    Double_t xsize[2] = { 0.0, 0.0 };
    Double_t ysize[2] = { 0.0, 0.0 };
    
    //Détermination des vecteurs de charge totale, multiplicité et taille qui devient précision
    for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
    {
        for ( Int_t i = 0; i < precluster.size(); ++i ){
            int detid = precluster[i].getDetID();
            int padid = precluster[i].getPadID();
            mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);
            
            if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
                    {
                        std::vector<double> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
                        std::vector<double> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};

                        chargetot[cathode] += precluster[i].getADC();
                        multiplicity[cathode] += 1;
                        xsize[cathode] += padSize[0];
                        ysize[cathode] += padSize[1];
                    }
        }
        
        if ( multiplicity[cathode] != 0 )
        {
          xsize[cathode] /= multiplicity[cathode];
          ysize[cathode] /= multiplicity[cathode];
        }
    }
    
    //Minimisation du Chi2 en 2D sur chaque cathode
    for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
    {
    
        float sigx = 1;
        float sigy = 1;
        
        float pasxhit = 0.001;
        float pasyhit = 0.001;
        float passigx = 0.01;
        float passigy = 0.01;
        
        int compteurboucle = 0;
        
        Double_t ChiDeux = Chi2Gauss(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], sigx, sigy);
        Double_t OldChiDeux = 1;
        Double_t NewChiDeux = 100;
        
        while (sqrt(pow((NewChiDeux-OldChiDeux)/OldChiDeux,2)) > 0.0001){
            compteurboucle++;
            
            //Iteration and improvement of estimate
            float dChidxhit = (Chi2Gauss(cathode, precluster, clustertmp, xhit[cathode]+pasxhit, yhit[cathode], sigx, sigy)-Chi2Gauss(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], sigx, sigy))/pasxhit;
            float dChidyhit = (Chi2Gauss(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode]+pasyhit, sigx, sigy)-Chi2Gauss(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], sigx, sigy))/pasyhit;
            
            float dChidsigx = (Chi2Gauss(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], sigx+passigx, sigy)-Chi2Gauss(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], sigx, sigy))/passigx;
            
            float dChidsigy = (Chi2Gauss(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], sigx, sigy+passigy)-Chi2Gauss(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], sigx, sigy))/passigy;

            
            xhit[cathode] = xhit[cathode] - pasxhit*dChidxhit;
            yhit[cathode] = yhit[cathode] - pasyhit*dChidyhit;
            sigx = sigx - passigx*dChidsigx;
            sigy = sigy - passigy*dChidsigy;

            
            // Update ChiDeux
            OldChiDeux = ChiDeux;
            ChiDeux = Chi2Gauss(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], sigx, sigy);
            NewChiDeux = ChiDeux;
            cout << "Boucle While numero " << compteurboucle << endl;
            cout << "xhit: " << xhit << endl;
            cout << "yhit: " << yhit << endl;
            cout << "sigx: " << sigx << endl;
            cout << "sigy: " << sigy << endl;
            cout << "OldChiDeux: " << OldChiDeux << endl;
            cout << "NewChiDeux: " << NewChiDeux << endl;
            
            if(compteurboucle > 10000){
                cout << "Déjà 10000 boucles pour trouver Chi2, utilisation résultats COG..." << endl;
                return clustertmp
            }
        }
        
        cout << "cathode " << cathode << endl;
        cout << "Boucle Finale numero " << compteurboucle << endl;
        cout << "xhit final: " << xhit[cathode] << endl;
        cout << "yhit final: " << yhit[cathode] << endl;
        cout << "sigx final: " << sigx << endl;
        cout << "sigy final: " << sigy << endl;
        cout << "Chi2 Final: " << NewChiDeux << endl;
    
    }
    
    
    //Détermination de la position finale du hit à l'aide de la résolution sur chaque cathode
    xhitfinal = ( xsize[0] < xsize[1] ) ? xhit[0] : xhit[1];
    yhitfinal = ( ysize[0] < ysize[1] ) ? yhit[0] : yhit[1];
    
    double timestamp = precluster[0].getTimeStamp();  // On récupère le timestamp du premier digit du precluster
    
    
    //On remplit un cluster avec les infos nécessaires
    cluster->setx(xhitfinal);
    cluster->sety(yhitfinal);
    cluster->setex(( xsize[0] < xsize[1] ) ? xsize[0] : xsize[1]);
    cluster->setey(( ysize[0] < ysize[1] ) ? ysize[0] : ysize[1]);
    cluster->settimestamp(timestamp);    // En première approx. on peut affecter le timestamp du pre;ier digit au cluster
    
    return cluster;

}

//_________________________________________________________________________________________________
float Clustering::Chi2Gauss(int cathode, Double_t chargetot, std::vector<Digit> &precluster, Clustering::Cluster* clustertmp, float xhit, float yhit, float sigx, float sigy)
{
    
    float sum = 0;
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
            
            float qfit = IntGaussXY(lowerLeft[0], lowerLeft[1], upperRight[0], upperRight[1], xhit, sigx, yhit, sigy)/chargetot[cathode];

            float qtrue = precluster[i].getADC()/chargetot[cathode];
            
            float delta = qfit-qtrue;
            
            sum += pow(delta,2);
        }
    }
    
    return sum;
}


//_________________________________________________________________________________________________
float Clustering::IntGaussXY(float x1, float y1, float x2, float y2, float xhit, float sigx, float yhit, float sigy)
{
    float erfx = 0.5*( erf( (xhit-x1)/(sigx*sqrt(2)) ) - erf( (xhit-x2)/(sigx*sqrt(2)) ) );
    float erfy = 0.5*( erf( (yhit-y1)/(sigy*sqrt(2)) ) - erf( (yhit-y2)/(sigy*sqrt(2)) ) );
    float integral = erfx*erfy;
    
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
