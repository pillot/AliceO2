//
//  ClusteringForTest.h
//  
//
//  Created by PERRIN Sébastien on 31/10/2019.
//

#ifndef ClusteringForTest_h
#define ClusteringForTest_h

#include <stdio.h>
#include <cassert>
#include <cstdint>
#include <unordered_map>
#include <vector>

//include michael
//include segmentation.h

#include "MCHBase/Digit.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHPreClustering/PreClusterBlock.h"

namespace o2
{
namespace mch
{

// Classes Cluster

class Clustering
{
    public:
    struct Cluster {
    Double_t x;            // x position du cluster
    Double_t y;            // y position du cluster
    // float z;            // z position du cluster
    Double_t ex;           // resolution en x
    Double_t ey;           // resolution en y
    //  float ez;           // resolution en z
    Double_t timestamp;    // Timestamp associé au cluster
    //  int idPreCluster;   // Id du PreCluster associé à ce cluster
        
        
    Double_t getx() const { return mx; }
    void setx(Double_t x) { mx = x; }
        
    Double_t gety() const { return my; }
    void sety(Double_t y) { my = y; }
        
    Double_t getex() const { return mex; }
    void setex(Double_t x) { mex = ex; }
        
    Double_t getey() const { return mey; }
    void setey(Double_t ey) { mey = ey; }

    Double_t gettimestamp() const { return mtimestamp; }
    void settimestamp(Double_t timestamp) { mtimestamp = timestamp; }

    private:
    Double_t mx;
    Double_t my;
    Double_t mex;
    Double_t mey;
    Double_t mtimestamp;
    };
    
    //   struct PreClusterStruct {
    //       std::vector<Digit> digits;
    //   };
    //
    //   struct PreClusterSizes {
    //       std::vector<int> sizes;
    //   };

  Clustering() = default;
  ~Clustering() = default;

  Clustering(const Clustering&) = delete;
  Clustering& operator=(const Clustering&) = delete;
  Clustering(Clustering&&) = delete;
  Clustering& operator=(Clustering&&) = delete;

    
    void runFinderCOG(std::vector<PreClusterStruct>& preClusters, std::vector<Cluster>& clusters);
    Clustering::Cluster FinderCOG(std::vector<Digit> &precluster);
//    
//    Clustering::Cluster* runOldFinderSimpleFit(std::vector<Digit> &preclusterslist, std::vector<int> &preclusterssizes);
//    Clustering::Cluster* runFinderSimpleFit(std::string& buffer);
//    Clustering::Cluster* FinderSimpleFit(std::vector<Digit> &precluster, Clustering::Cluster* clustertmp);
//    float Chi2Mathieson(int cathode, Double_t chargetot, std::vector<Digit> &precluster, Clustering::Cluster* clustertmp, float xhit, float yhit, float Kx3, float Ky3);
//    
//    float IntMathiesonXY(float x1, float y1, float x2, float y2, float Kx3, float Ky3);
//    
//    Clustering::Cluster* runFinderGaussianFit(std::string& buffer);
//    Clustering::Cluster* FinderGaussianFit(std::vector<Digit> &precluster, Clustering::Cluster* clustertmp);
//    float Chi2Gauss(int cathode, Double_t chargetot, std::vector<Digit> &precluster, Clustering::Cluster* clustertmp, float xhit, float yhit, float sigx, float sigy);
//    
//    float IntGaussXY(float x1, float y1, float x2, float y2, float xhit, float sigx, float yhit, float sigy);
//
//    
//    void ComparisonMathiesonGauss(std::string& buffer);
};


} // namespace mch
} // namespace o2

#endif /* ClusteringForTest_h */
