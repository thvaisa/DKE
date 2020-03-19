#include <vector>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
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
template<class VecType> 
void get_dimensions(std::vector<VecType> positions, double* dimensions){
    dimensions[0] = std::numeric_limits<double>::max();
    dimensions[1] = -std::numeric_limits<double>::max();
    dimensions[2] = dimensions[0],
    dimensions[3] = dimensions[1];
    dimensions[4] = dimensions[0];
    dimensions[5] = dimensions[1];
    for (std::size_t i = 0; i < positions.size(); ++i)
    {
        VecType position = positions[i];
        double pos[3] = {position[0],position[1],position[2]};
        for(std::size_t j=0;j<3;++j){
            dimensions[j*2] = std::min(dimensions[j*2],pos[j]);
            dimensions[j*2+1] = std::max(dimensions[j*2+1],pos[j]);
        }
    }
}


double kernelFunction(Eigen::Vector3d& dist, Eigen::Matrix3d& covinv, double detCov){

    Eigen::Vector3d a = dist;
    double value = a.transpose()*(covinv*a);
    value = std::exp(-0.5*value)/(std::pow(2*M_PI,3)*std::sqrt(detCov));
    return value;
}



//Create point cloud
double* create_estimation_grid(std::size_t* nPoints, double* bbox){
    //std::vector<pcl::PointXYZ>::Ptr positions = new std::vector<pcl::PointXYZ>();

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


    double* grid = new double[nPoints[0]*nPoints[1]*nPoints[2]];

    return grid;
}



inline void get_points_within_radius(Eigen::Vector3d point, double radius, 
                            double* bbox, std::size_t* nPoints,
                            double* densityMap,Eigen::Matrix3d& invcov, double detCov){
    double dX = (bbox[1]-bbox[0])/nPoints[0];

    int sx = std::max((int)std::ceil((point[0]-bbox[0])/dX),0);
    int sy = std::max((int)std::ceil((point[1]-bbox[2])/dX),0);
    int sz = std::max((int)std::ceil((point[2]-bbox[4])/dX),0);

    int nDim = (int)std::floor(2*radius/dX); 
    Eigen::Vector3d pos = Eigen::Vector3d(0,0,0);
    for(int i=sx;i<=sx+nDim && i<nPoints[0];++i){
        for(int j=sy;j<=sy+nDim && j<nPoints[1];++j){
            for(int k=sz;i<=sz+nDim && k<nPoints[2];++k){
                pos[0] = bbox[0]+i*dX;
                pos[1] = bbox[2]+j*dX;
                pos[2] = bbox[4]+k*dX;
                Eigen::Vector3d dist = pos-point;
                if((dist).dot(dist)<radius*radius){
                    int indx = i*nPoints[1]*nPoints[2]+j*nPoints[2]+k;
                    densityMap[indx] = densityMap[indx]
                                    +kernelFunction(dist,invcov,detCov);
                }
            }
        }
    }

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
    std::vector<Eigen::Vector3d> points = std::vector<Eigen::Vector3d>();
    read_file_bin<Eigen::Vector3d>(points,argv[1],scale);
    get_dimensions<Eigen::Vector3d>(points, bbox);
    

    double extra = (bbox[1]-bbox[0])*0.01;

    for(std::size_t i=0;i<3;++i){
        bbox[i*2] = bbox[i*2]-extra;
        bbox[i*2+1] = bbox[i*2+1]+extra;
    }

    


    std::size_t nPoints[3] = {nPs,nPs,nPs};
    double * grid = create_estimation_grid(nPoints, bbox);
    double * densityMap = new double[nPoints[0]*nPoints[1]*nPoints[2]];
    
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

    Eigen::Vector3d searchPoint = Eigen::Vector3d();
    auto dist_vector = std::vector<Eigen::Vector3d>(); 
    auto indxs = std::vector<int>();

    for(std::size_t i=0;i<points.size();++i){
        //Eigen::Vector3d mean(searchPoint.x,searchPoint.y,searchPoint.z);
        //if(std::sqrt(mean.dot(mean))<0.5) continue;
        get_points_within_radius(points[i], radius, 
                            bbox, nPoints,
                            densityMap,invcov,detCov);
    }


    for (std::size_t j = 0; j < nPoints[0]*nPoints[1]*nPoints[2]; ++j){
        maxVal = std::max(densityMap[j],maxVal);
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
    delete[] grid;
    delete[] densityMap;
    return EXIT_SUCCESS;
}



