#include <bits/stdc++.h>
#include "matrix.h"
#include "time.h"
#include "../CImg/CImg.h"
#include <omp.h>
#include <unistd.h>
#define MAX_NUMTHREADS 16
#define timediff ((double)(clock() - start))/CLOCKS_PER_SEC * 1000
#define printTime(msg) {std::cout<<(msg)<<" : "<<timediff<<" ms "<<std::endl; start = clock();}

using namespace std;

using namespace cimg_library; 

#if defined _MSC_VER
#include <direct.h>
#elif defined __GNUC__
#include <sys/types.h>
#include <sys/stat.h>
#endif
typedef double data_t;
typedef unsigned char pixel_t; 
typedef char cp_t;
string remove_suf(string x)
{
    size_t s = x.size();
    for(size_t i=0;i<x.size();i++)
    {
        if(x[i]=='.')s = i;
    }
    return x.substr(0, s);
}
void createDir(string dir) {
#if defined _MSC_VER
    _mkdir(dir.data());
#elif defined __GNUC__
    mkdir(dir.data(), 0777);
#endif
}


string image_name;// = "ranni.jpg";
string image_pure;// = "ranni"


const int N = 8;

matrix<data_t> T(N,N);
matrix<data_t> Q(N,N);
int reconstruct = 0;
vector<double> dct_times;
vector<double> cps_rates;


long GetFileSize(std::string filename)
{
    struct stat stat_buf;
    int rc = stat(filename.c_str(), &stat_buf);
    return rc == 0 ? stat_buf.st_size : -1;
}

//helping file functions
template<typename T>
void flush_disk(vector<vector<T>> src, string dest)
{
    FILE* file = fopen(dest.c_str(), "wb");
    size_t siz = src.size();
    fwrite(&siz, 1, sizeof(size_t), file);
    for(auto &e : src)
    {
        char sub_siz = (char)e.size();
        fwrite(&sub_siz, 1, sizeof(char), file);
        fwrite(&e[0], 1, sub_siz * sizeof(T), file);
    }

    fclose(file);
}
template<typename T>
void fetch_disk(vector<vector<T>> &ans, string src) {
    FILE* file = fopen(src.c_str(), "rb");
    size_t siz;
    fread(&siz, 1, sizeof(size_t), file);
    for(size_t i=0;i<siz;i++)
    {
        vector<T> temp;
        char n;
        fread(&n, 1, sizeof(char), file);
        temp.resize(n);
        fread(&temp[0], 1, n * sizeof(T), file);
        ans.push_back(std::move(temp));
    }
    fclose(file);
}


template<typename T>
bool check_equal(vector<vector<T>> a, vector<vector<T>> b)
{
    if(a.size()!=b.size())return false;
    for(int i=0;i<a.size();i++)
    {
        if(a[i].size()!=b[i].size())return false;
        for(int j=0;j<a[i].size();j++)
        {
            if(a[i][j]!=b[i][j])return false;
        }
    }
    return true;
}



void preprocess(data_t Quality)
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

    if(Quality > 50) Q = Q * ((100.0 - Quality)/50.0);
    if(Quality < 50) Q = Q * (50.0 / Quality);

}

CImg<pixel_t> graying(CImg<pixel_t> src)
{
    int width = src.width();
    int height = src.height();
    int channel = src.spectrum();
    CImg<pixel_t> res(src);


    for(int y = 0; y < height; y++)
    {
        for(int x = 0; x < width; x++)
        {
            int intermid = 0;
            for(int c = 0; c < channel; c++)
            {
                intermid += src(x,y,0,c);
            }
            intermid/=channel;
            for(int c = 0; c < channel; c++)res(x,y,0,c) = intermid;
        }
    }
    return res;
}

uint64_t timeSinceEpochMillisec() {
  using namespace std::chrono;
  return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}


void run(data_t Quality, int nthreads)
{
    clock_t start = clock();
    preprocess(Quality);
    printTime("preprocessing"); 
    CImg<pixel_t> src(image_name.c_str());
    CImg<pixel_t> res(src);
    CImg<pixel_t> gray;
    gray = graying(src);
    string grayfile = image_pure + "/" + image_pure + "_grayscale.jpg";
    gray.save(grayfile.c_str());
    int width = src.width();
    int height = src.height();
    int channel = src.spectrum();

    printTime("image processing"); 

    double total_rate[MAX_NUMTHREADS];
    string final_compress = "";
    final_compress += image_pure + "/" + image_pure + "_compressed_";
    final_compress += std::to_string((int)Quality);

    vector<vector<cp_t>> cps_vecs[MAX_NUMTHREADS];
    
    int blocks = height * width/(N*N);
    int step = (blocks + nthreads - 1)/nthreads;  
    uint64_t cloc = timeSinceEpochMillisec();
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        total_rate[tid] = 0;

        for(int i = tid * step; i < blocks && i < (tid+1)*step; i++)
        {
            int y = i/(width/N);
            int x = i%(width/N);
            y*=N;
            x*=N;
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
            auto cps_vec = C.zigzag(1e-5);
            double compression_rate = (double)cps_vec.size()/(N*N);
            //matrix<data_t> C_(N,N);
            //C_.zagzig(cps_vec);
            cps_vecs[tid].push_back(std::move(cps_vec));
            
            total_rate[tid] += compression_rate;

            F = (T.T() * (C ^ Q) * T) + 128; 


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
    //cout<< timeSinceEpochMillisec() - cloc<<endl;
    double sum_rate = 0;
    vector<vector<cp_t>> cps_final;
    for(int i=0;i<nthreads;i++)
    {
        sum_rate += total_rate[i];
        for(auto &x : cps_vecs[i])
        {
            cps_final.push_back(std::move(x));
        }
    }
    dct_times.push_back(timeSinceEpochMillisec() - cloc); start = clock();
    cout<<"dct : "<<dct_times.back()<<" ms"<<endl;
    

    
    vector<vector<char>> d_com;
    flush_disk(cps_final, final_compress + ".cps");
    fetch_disk(d_com, final_compress + ".cps");
    //assert(check_equal(cps_final, d_com));


    printTime("flushing"); 

    if(!reconstruct)
    {
        res.save((final_compress + ".jpg").c_str());
    }
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
    printTime("save"); 
    auto rate = sum_rate * (N*N)/(width * height) * 100;
    
    cps_rates.push_back(rate);
    std::cout<<"compression rate : "<<rate<<"%"<<std::endl;
}
int main(int argc, char *argv[])
{

    data_t Quality = 50;
    int nthreads = 1;
    
    int opt = 0;
    do {
        opt = getopt(argc, argv, "f:p:n:i:");
        switch (opt) {
            case 'f':
                image_name = optarg;
                break;

            case 'p':
                Quality = atof(optarg);
                break;

            case 'n':
                nthreads = atoi(optarg);
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


    //std::cout<<reconstruct<<std::endl;
    omp_set_num_threads(nthreads);

    





    image_pure = remove_suf(image_name);

    createDir(image_pure);
    
    //run(Quality, nthreads);
    vector<double> actual_rates;
    for(data_t q = 5; q < 100; q+=5)
    {
        string dest = image_pure + "/" + image_pure+"_compressed_" + std::to_string((int)q) + ".log";
        freopen(dest.c_str(), "w+", stdout);
        run(q, nthreads);
        string cps = image_pure + "/" + image_pure+"_compressed_" + std::to_string((int)q) + ".cps";
        string ori = image_pure + "/" + image_pure+"_grayscale.jpg";
        cout<<GetFileSize(cps)<<' '<<GetFileSize(ori);
        double actual_rate = (double)GetFileSize(cps)/GetFileSize(ori)*100.0;
        actual_rates.push_back(actual_rate);
        
    }
    string dest = image_pure + "/" + image_pure+"_summary_" + std::to_string((int)nthreads) + ".log";
    freopen(dest.c_str(), "w+", stdout);
    uint64_t total_time = 0;
    for(size_t i=0;i<cps_rates.size();i++)
    {
        cout<<(i*5 + 5)<<' '<<cps_rates[i]<<"% "<<actual_rates[i]<<"% "<<dct_times[i]<<endl;
        total_time +=dct_times[i];
    }
    cout<<total_time<<endl;
    

    return 0;
}