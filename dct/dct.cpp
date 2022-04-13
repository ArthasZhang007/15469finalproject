#include <bits/stdc++.h>
#include "matrix.h"
#include "../CImg/CImg.h"
using namespace std;

using namespace cimg_library; 
string image_name;// = "elden-ring-ranni.jpg";
typedef double data_t;
typedef unsigned char pixel_t; 

const int N = 8;
data_t Quality = 50;
matrix<data_t> T(N,N);
matrix<data_t> Q(N,N);
int reconstruct = 0;
void preprocess()
{
    // the cosine matrix
    data_t sqt2 = sqrt(2.0);
    data_t inv_sqtN = 1.0/sqrt(N);
    for(int j = 0; j < N; j++)
    {
        for(int i = 0; i < N; i++)
        {
            if(i==0)T.set(i,j, inv_sqtN);
            else T.set(i,j, sqt2 * inv_sqtN * cos(
                (2.0*j + 1.0)* i * M_PI/(2.0 * N)
            ));
        }
    }
    //std::cout<<T<<std::endl;
    //auto q = T.zigzag(1e-6);
    //std::cout<<q.size()<<std::endl;
    //for(auto &e : q)std::cout<<e<<' ';
    //std::cout<<std::endl;

    //std::cout<<T.T()<<std::endl;

    // the quant matrix
    ifstream qin("quant.txt");
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            data_t v;
            qin>>v;
            Q.set(i,j,v);
        }
    }
    //std::cout<<Q<<std::endl;
    if(Quality > 50) Q = Q * ((100.0 - Quality)/50.0);
    if(Quality < 50) Q = Q * (50.0 / Quality);
    //std::cout<<Q<<std::endl;
    /*matrix<data_t> D(N,N),M(N,N);
    D = T * M * T.T();*/
}
CImg<pixel_t> graying(CImg<pixel_t> src)
{
    int width = src.width();
    int height = src.height();
    int channel = src.spectrum();
    CImg<pixel_t> res(src);
    //pixel_t *data = src.data();
    for(int y = 0; y < height; y++)
    {
        for(int x = 0; x < width; x++)
        {
            int intermid = 0;
            for(int c = 0; c < channel; c++)
            {
                intermid += src(x,y,0,c);
                //src(x,y,0,2) = src(x,y,0,2)/2;//min(255,(int)src(x,y,0,c)*3/2);
            }
            intermid/=channel;
            for(int c = 0; c < channel; c++)res(x,y,0,c) = intermid;
        }
    }
    return res;
}
int main(int argc, char *argv[])
{
    int opt = 0;
    do {
        opt = getopt(argc, argv, "f:p:i:");
        switch (opt) {
        case 'f':
            image_name = optarg;
            break;

        case 'p':
            Quality = atof(optarg);
            break;

        case 'i':
            reconstruct = atoi(optarg);
            break;

        case -1:
            break;

        default:
            break;
        }
    } while (opt != -1);






    clock_t start = clock();
    preprocess();
    std::cout<<"preprocessing : "<<((double)(clock() - start))/1000<<" ms"<<std::endl; start = clock();
    CImg<pixel_t> src(image_name.c_str());
    CImg<pixel_t> res(src);
    CImg<pixel_t> gray;
    gray = graying(src);
    int width = src.width();
    int height = src.height();
    int channel = src.spectrum();

    std::cout<<"image processing : "<<((double)(clock() - start))/1000<<" ms"<<std::endl; start = clock();
    double total_rate = 0;
    for(int y = 0; y < height; y += N)
    {
        for(int x = 0; x < width; x += N)
        {
            matrix<double> M(N,N), C(N,N), F(N,N);
            for(int dy = 0; dy < N; dy ++)
            {
                for(int dx = 0; dx < N; dx ++)
                {
                    M.set(dx, dy, src(x+dx, y+dy, 0, 0));
                }
            }
            
            C = (T * (M - 128) * T.T())/Q;
            C.round();
            double compression_rate = (double)C.zigzag(1e-5).size()/(N*N);
            //if(compression_rate < 1.0)std::cout<<compression_rate<<std::endl;
            total_rate += compression_rate;

            F = (T.T() * (C ^ Q) * T) + 128; 


            /*if(++cnt < 5)
            {
                std::cout<<C<<std::endl;
                //std::cout<<F<<std::endl;
            }*/

            for(int dy = 0; dy < N; dy ++)
            {
                for(int dx = 0; dx < N; dx ++)
                {
                    for(int c = 0; c < channel; c++)
                        res(x+dx, y+dy, 0, c) = F.get(dx, dy);
                }
            }
        }
    }
    std::cout<<"dct : "<<((double)(clock() - start))/1000<<" ms"<<std::endl; start = clock();
    string final_compress = "ranni_compressed_";
    final_compress += std::to_string((int)Quality)+".jpg";

    if(!reconstruct)res.save(final_compress.c_str());
    else{
        for(int y = 0; y < height; y++)
        {
            for(int x = 0; x < width; x++)
            {
                for(int c = 0; c < channel; c++)
                {
                    //int val = 0;
                    //if(gray(x,y,0,c))val = std::min(255, (int)src(x,y,0,c) * res(x,y,0,c) / gray(x,y,0,c));
                    //val = std::max(val, 0);
                    //src(x,y,0,c) = val;

                    src(x,y,0,c) = res(x,y,0,c);

                    //src(x,y,0,2) = src(x,y,0,2)/2;//min(255,(int)src(x,y,0,c)*3/2);
                }
                //std::cout<<res(x,y,0,0) - gray(x,y,0,0c)<<std::endl;
            }
        }
        src.save(final_compress.c_str());
    }
    
    std::cout<<"save : "<<((double)(clock() - start))/1000<<" ms"<<std::endl; start = clock();
    std::cout<<"compression rate : "<<total_rate * (N*N)/(width * height) * 100<<"%"<<std::endl;
    return 0;
}