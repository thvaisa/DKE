#include <vector>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <pcl/point_cloud.h>
#include <pcl/octree/octree_search.h>
#include "lodepng.h"
#include <Eigen/Dense>
#include <jsoncpp/json/json.h>
#include <string>

//Read points from the file
template<class VecType> 
void read_file(std::vector<VecType>& positions, const char* fname, double scale){
    std::ifstream infile(fname);
    double a,b,c;
    double d,e,f;
    bool first = true;
    while (infile >> a >> b >> c)
    {
        if(first){
            d = 0;//a;
            e = 0;//b;
            f = 0;//c;
            first = false;
        } 
        positions.push_back(VecType((a-d)*scale,(b-e)*scale,(c-f)*scale));
    }

}

template<class VecType> 
void read_file_bin(std::vector<VecType>& positions, const char* fname, double scale){
    FILE* file;
    file = fopen(fname, "rb");
    std::ifstream infile(fname);
    double a,b,c;
    double d,e,f;
    bool first = true;
    while (!feof(file))
    {
        if(first){
            d = 0;//a;
            e = 0;//b;
            f = 0;//c;
            first = false;
        } 
        double tmp[3];
        std::size_t val = fread(tmp,sizeof(double),3,file); 
        positions.push_back(VecType((tmp[0]-d)*scale,(tmp[1]-e)*scale,(tmp[2]-f)*scale));
    }
    fclose(file);
}


//Get dimensions
void get_dimensions(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, double* dimensions){
    dimensions[0] = std::numeric_limits<double>::max();
    dimensions[1] = -std::numeric_limits<double>::max();
    dimensions[2] = dimensions[0],
    dimensions[3] = dimensions[1];
    dimensions[4] = dimensions[0];
    dimensions[5] = dimensions[1];
    for (std::size_t i = 0; i < cloud->points.size(); ++i)
    {
        pcl::PointXYZ position = cloud->points[i];
        double pos[3] = {position.x,position.y,position.z};
        for(std::size_t j=0;j<3;++j){
            dimensions[j*2] = std::min(dimensions[j*2],pos[j]);
            dimensions[j*2+1] = std::max(dimensions[j*2+1],pos[j]);
        }
    }
}



//Create point cloud
pcl::PointCloud<pcl::PointXYZ>::Ptr createPointCloud(const char*  fname, double scale){
    std::vector<pcl::PointXYZ> positions = std::vector<pcl::PointXYZ>();
    read_file_bin<pcl::PointXYZ>(positions,fname,scale);

    // Generate pointcloud data
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);

    cloud->width = positions.size();
    cloud->height = 1;
    cloud->points.resize (cloud->width * cloud->height);

    for (std::size_t i = 0; i < positions.size (); ++i)
    {
        cloud->points[i].x = positions[i].x;
        cloud->points[i].y = positions[i].y;
        cloud->points[i].z = positions[i].z;
    }
    return cloud;
}


//Create octree
pcl::octree::OctreePointCloudSearch<pcl::PointXYZ>::Ptr createOctree(pcl::PointCloud<pcl::PointXYZ>::Ptr pointCloud, double resolution){
    pcl::octree::OctreePointCloudSearch<pcl::PointXYZ>::Ptr octree(new pcl::octree::OctreePointCloudSearch<pcl::PointXYZ>(resolution));

    octree->setInputCloud (pointCloud);
    octree->addPointsFromInputCloud ();
    return octree;

}





//Create point cloud
pcl::PointCloud<pcl::PointXYZ>::Ptr createEstimationGrid(std::size_t* nPoints, double* bbox){
    //std::vector<pcl::PointXYZ>::Ptr positions = new std::vector<pcl::PointXYZ>();

    // Generate pointcloud data
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);

    double dX = (bbox[1]-bbox[0])/(nPoints[0]-1);
    double dY = (bbox[3]-bbox[2])/(nPoints[1]-1);
    double dZ = (bbox[5]-bbox[4])/(nPoints[2]-1);

    double minDX = std::min(dX,std::min(dY,dZ));
    nPoints[0] = (int)std::ceil((bbox[1]-bbox[0])/minDX);
    nPoints[1] = (int)std::ceil((bbox[3]-bbox[2])/minDX);
    nPoints[2] = (int)std::ceil((bbox[5]-bbox[4])/minDX);

    dX = minDX;
    dY = minDX;
    dZ = minDX;

   
    cloud->width = nPoints[0]*nPoints[1]*nPoints[2];
    cloud->height = 1;
    cloud->points.resize (cloud->width * cloud->height);

    
    int indx = 0;
    for (std::size_t i = 0; i < nPoints[0]; ++i)
    {
        for (std::size_t j= 0; j < nPoints[1]; ++j)
        {
            for (std::size_t k = 0; k < nPoints[2]; ++k)
            {
                cloud->points[indx].x = bbox[0]+i*dX;
                cloud->points[indx].y = bbox[2]+j*dY;
                cloud->points[indx].z = bbox[4]+k*dZ;



                ++indx;
            }
        }
    }
    return cloud;
}





double kernelFunction(Eigen::Vector3d& x, Eigen::Vector3d & mean, Eigen::Matrix3d& covinv, double detCov){

    Eigen::Vector3d a = x-mean;
    double value = a.transpose()*(covinv*a);
    value = std::exp(-0.5*value)/(std::pow(2*M_PI,3)*std::sqrt(detCov));
    return value;
}


//Openvbd modified cookbook
int main(int argc, char* argv[]){

    //std for the normal distribution
    double std = 1.0;
    //how many points are used to discretize the space
    std::size_t nPs = 100;
    //scale
    double scale = 1.0e-4;

    double bbox[6] = {-10,10,-10,10,-10,10};
    if(argc<4){
        std::cout << "Give me some parameters: " << std::endl;
        return EXIT_FAILURE;
    }

    std::string str;
    std::string output_file;
    std::string stamp = std::string(argv[2]);
    std = std::stod(std::string(argv[3]));
    scale = std::stod(std::string(argv[4]));
    nPs = std::stod(std::string(argv[5]));
    output_file = std::string(argv[6]);


    std::cout << "output to" << std::endl;
    std::string pngFile = output_file+stamp+std::string(".png");
    std::string auxFile = output_file+stamp+std::string("aux.json");
    std::cout << pngFile << " " << auxFile << std::endl;

    std::cout << std << " " << scale << " " << nPs << " " << output_file << std::endl;
    std::cout << stamp << " " << argv[1] << std::endl;


    //Create point cloud
    auto cloud = createPointCloud(argv[1],scale);
    get_dimensions(cloud, bbox);

    double extra = (bbox[1]-bbox[0])*0.01;

    for(std::size_t i=0;i<3;++i){
        bbox[i*2] = bbox[i*2]-extra;
        bbox[i*2+1] = bbox[i*2+1]+extra;
    }


    float resolution = extra*0.1;
    
    //auto octree = createOctree(cloud, resolution);
    std::cout << resolution << std::endl;
    //Create Octree for the isosurface extraction
    std::size_t nPoints[3] = {nPs,nPs,nPs};
    auto cloud2 =createEstimationGrid(nPoints, bbox);
    auto octree2 = createOctree(cloud2, resolution);
    auto densityMap = std::vector<double>(nPoints[0]*nPoints[1]*nPoints[2]);
    
    //Iterate all the points from the simulation
    std::vector<int> pointIdxRadiusSearch;
    std::vector<float> pointRadiusSquaredDistance;
    pcl::PointXYZ searchPoint;
    double maxVal = 0;    
    double radius = 2*std;

    Eigen::Matrix3d A(3,3);

    A(0,0) = std;
    A(0,1) = 0;
    A(0,2) = 0;
    A(1,0) = 0;
    A(1,1) = std;
    A(1,2) = 0;
    A(2,0) = 0;
    A(2,1) = 0;
    A(2,2) = std;

    Eigen::Matrix3d invcov = A.inverse();
    double detCov = A.determinant();


    for(std::size_t i=0;i<cloud->points.size();++i){
        searchPoint.x = cloud->points[i].x;
        searchPoint.y = cloud->points[i].y;
        searchPoint.z = cloud->points[i].z;

        Eigen::Vector3d mean(searchPoint.x,searchPoint.y,searchPoint.z);
        if(std::sqrt(mean.dot(mean))<0.5) continue;

        if (octree2->radiusSearch (searchPoint, radius, pointIdxRadiusSearch, pointRadiusSquaredDistance) > 0)
        {
            for (std::size_t j = 0; j < pointIdxRadiusSearch.size (); ++j){
                
                Eigen::Vector3d x(cloud2->points[pointIdxRadiusSearch[j]].x,
                                cloud2->points[pointIdxRadiusSearch[j]].y,
                                cloud2->points[pointIdxRadiusSearch[j]].z);
                //std::cout << kernelFunction(mean,x,invcov,detCov) << std::endl;
                
                densityMap[pointIdxRadiusSearch[j]] = densityMap[pointIdxRadiusSearch[j]]
                                                        +kernelFunction(mean,x,invcov,detCov);
                maxVal = std::max(densityMap[pointIdxRadiusSearch[j]],maxVal);
            }
        }

    }

    std::cout << maxVal << std::endl;


    //std::vector<unsigned char> pngBuffer(nPoints[0]*nPoints[1]*nPoints[2]*4);
    std::vector<unsigned char> pngBuffer;
    pngBuffer.reserve(nPoints[0]*nPoints[1]*nPoints[2]*4);
    int indx=0;
    for (std::size_t i = 0; i < nPoints[0]; ++i)
    {
        for (std::size_t j= 0; j < nPoints[1]; ++j)
        {
            for (std::size_t k = 0; k < nPoints[2]; ++k)
            {
                pngBuffer.push_back((std::uint8_t)(densityMap[indx]*1.0/maxVal*255));
                pngBuffer.push_back((std::uint8_t)(densityMap[indx]*1.0/maxVal*255));
                pngBuffer.push_back((std::uint8_t)(densityMap[indx]*1.0/maxVal*255));
                pngBuffer.push_back((std::uint8_t)(255));
                //std::cout << "ahashd" << std::endl;
                indx++;
            }
        }
    }


    

    //std::vector<std::uint8_t> ImageBuffer;
    //lodepng::encode(ImageBuffer, pngBuffer, nPoints[1]*nPoints[2], nPoints[0]);
    //lodepng::save_file(ImageBuffer, "SomeImage.png");
    unsigned error =lodepng::encode(pngFile, pngBuffer,  nPoints[1]*nPoints[2], nPoints[0]);
    //if there's an error, display it
     if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;

    //print
    
    //for(std::size_t i=0;i<cloud2->points.size();++i){
        //if(densityMap[i]>=1){
        //    pcl::PointXYZ point = cloud2->points[i];
        //    std::cout << point.x*1.0<< " " << point.y << " " << point.z << std::endl;
            //std::cout << densityMap[i]*1.0/maxVal << std::endl;
        //}



    //    std::cout << densityMap[i]*1.0/maxVal << std::endl;
        //if(i>10000){
       //     break;
        //}
    //}


    /*
    std::ofstream myfile;
    myfile.open (auxFile);

    for(std::size_t i=0;i<3;++i){
        myfile <<  std::abs(bbox[i*2+1]-bbox[i*2]) << std::endl; 
    }
    for(std::size_t i=0;i<3;++i){
        myfile <<  (bbox[i*2+1]+bbox[i*2])/2.0 << std::endl; 
    }
    myfile <<   nPoints[0] << " " << nPoints[1] << " " << nPoints[2] << std::endl; 
    myfile.close();
    */

    std::ofstream file_id;
    file_id.open(auxFile);
    Json::Value dimension(Json::arrayValue);
    Json::Value origin(Json::arrayValue);
    Json::Value indices(Json::arrayValue);
    Json::Value event;
    for(std::size_t i=0;i<3;++i){
        dimension.append(Json::Value(std::abs(bbox[i*2+1]-bbox[i*2]))); 
    }

    for(std::size_t i=0;i<3;++i){
        origin.append(Json::Value((bbox[i*2+1]+bbox[i*2])/2.0)); 
    }

    for(std::size_t i=0;i<3;++i){
        indices.append(Json::Value(nPoints[i]*1.0)); 
    }

    event["1"]["dimension"]=dimension;
    event["1"]["origin"]=origin;
    event["1"]["indices"]=indices;
    
    Json::StyledWriter styledWriter;
    file_id << styledWriter.write(event);

    file_id.close();

    return EXIT_SUCCESS;
}



