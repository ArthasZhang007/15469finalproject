#pragma once
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>


template<typename data_t>
class matrix{
    public:
        data_t *data;
        int w,h;
        matrix(int x, int y) {

            w = x; 
            h = y;
            data = new data_t[w*h];
        }
        ~matrix()
        {
            //delete[] data;
        }
        data_t get(int x, int y)
        {
            int id = y * w + x;
            assert(id >= 0 && id < w * h);
            return data[y * w + x];
        }
        void set(int x, int y, data_t c)
        {
            int id = y * w + x;
            assert(id >= 0 && id < w * h);
            data[y * w + x] = c;
        }
        friend matrix<data_t> operator *(matrix<data_t> a, matrix<data_t> b)
        {
            assert(a.h == b.w);
            matrix<data_t> c(a.w, b.h);
            memset(c.data, 0, sizeof(data_t) * c.w * c.h);
            for(int j = 0; j < c.h; j++)
            {
                for(int i = 0; i < c.w; i++)
                {
                    for(int k = 0; k < a.h; k++)
                        c.set(i,j, c.get(i,j) + a.get(i,k) * b.get(k,j));
                }
            }
            return c;
        } 
        friend matrix<data_t> operator *(matrix<data_t> a, data_t k)
        {
            matrix<data_t> c(a.w, a.h);
            for(int j = 0; j < c.h; j++)
            {
                for(int i = 0; i < c.w; i++)
                {
                    c.set(i,j, a.get(i,j) * k);
                }
            }
            return c;
        } 
        friend matrix<data_t> operator /(matrix<data_t> a, data_t k)
        {
            matrix<data_t> c(a.w, a.h);
            for(int j = 0; j < c.h; j++)
            {
                for(int i = 0; i < c.w; i++)
                {
                    c.set(i,j, a.get(i,j) / k);
                }
            }
            return c;
        } 
        friend matrix<data_t> operator /(matrix<data_t> a, matrix<data_t> b)
        {
            assert(a.w == b.w && a.h == b.h);
            matrix<data_t> c(a.w, a.h);
            for(int j = 0; j < c.h; j++)
            {
                for(int i = 0; i < c.w; i++)
                {
                    c.set(i,j, a.get(i,j) / b.get(i,j));
                }
            }
            return c;
        }
        friend matrix<data_t> operator ^(matrix<data_t> a, matrix<data_t> b)
        {
            assert(a.w == b.w && a.h == b.h);
            matrix<data_t> c(a.w, a.h);
            for(int j = 0; j < c.h; j++)
            {
                for(int i = 0; i < c.w; i++)
                {
                    c.set(i,j, a.get(i,j) * b.get(i,j));
                }
            }
            return c;
        }
        friend matrix<data_t> operator +(matrix<data_t> a, data_t k)
        {
            matrix<data_t> c(a.w, a.h);
            for(int j = 0; j < c.h; j++)
            {
                for(int i = 0; i < c.w; i++)
                {
                    c.set(i,j, a.get(i,j) + k);
                }
            }
            return c;
        }  
        friend matrix<data_t> operator -(matrix<data_t> a, data_t k)
        {
            matrix<data_t> c(a.w, a.h);
            for(int j = 0; j < c.h; j++)
            {
                for(int i = 0; i < c.w; i++)
                {
                    c.set(i,j, a.get(i,j) - k);
                }
            }
            return c;
        } 
        

        bool inbound(int x, int y)
        {
            return x >= 0 && x < w && y >= 0 && y < h;
        }
        void round()
        {
            for(int j = 0; j < h; j++)
            {
                for(int i = 0; i < w; i++)
                {
                    set(i,j, int(get(i,j) + 0.5));
                }
            }
        }
        std::vector<data_t> zigzag(data_t eps)
        {
            std::vector<data_t> res;
            int x = 0, y = 0, lx,ly,mode = 0, offset = 0;
            res.push_back(get(0,0)); 
            do
            {
                
                lx = x; ly = y;
                if(mode == 0)
                {
                    x++;
                    mode = (1 + offset) % 4;
                }
                else if(mode == 1)
                {
                    if(inbound(x-1, y+1))
                    {
                        x--;
                        y++;
                    }
                    else mode = (2 + offset) % 4;
                }
                else if(mode == 2)
                {
                    y++;
                    mode = (3 + offset) % 4;
                    //else mode = 0;
                }
                else // mode = 3
                {
                    if(inbound(x+1, y-1))
                    {
                        x++;
                        y--;
                    }
                    else mode = (0 + offset) % 4;
                }
                
                if(lx!=x || ly!=y )
                {
                    
                    if(!inbound(x,y))
                    {
                        offset = 2;
                    }
                    else 
                    {
                        res.push_back(get(x,y));
                        //std::cout<<x<<' '<<y<<std::endl;
                    }
                }
            }while(x != w - 1 || y != h - 1);
            while(!res.empty() && std::abs(res.back()) < eps)res.pop_back();
            return res;
        }

        matrix<data_t> T()
        {
            matrix<data_t> c(h,w);
            for(int j = 0; j < c.h; j++)
            {
                for(int i = 0; i < c.w; i++)
                {
                    c.set(i,j, get(j,i));
                }
            }
            return c;
        }
        
};
template<typename data_t>
std::ostream& operator <<(std::ostream &os, matrix<data_t> A)
{
    os << std::setiosflags(std::ios::fixed)<<std::setprecision(4);
    os << A.w <<' '<<A.h<<std::endl;
    for(int i = 0; i < A.w; i++)
    {
        for(int j = 0; j < A.h; j++)
        {
            os<<A.get(i,j)<<' ';
        }
        os<<std::endl;
    }
    return os;
}

