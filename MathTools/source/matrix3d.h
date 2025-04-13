#pragma once
#include <math.h>
#include "vectorNd.h"
namespace GeomElements
{

    template<int ndim, typename T>
    class matrix3d
    {
        double components[ndim*ndim]={0};
    public:
        matrix3d(double* aoa)
        {
            if(ndim==2)
            {
                components[0] = cos(*aoa);
                components[1] = -sin(*aoa);
                components[2] = sin(*aoa);
                components[3] = cos(*aoa);
            }
        }

        vector3d<ndim,T> dot_product(vector3d<ndim,T>& column)
        {
            vector3d<ndim,T> result;
            for(int i = 0;i<ndim;i++)
            {
                result[i] = 0;
                for(int j = 0;j<ndim;j++)
                {
                    result[i] += components[ndim*i+j] * column[j];
                }
            }
            return result;
        }




    };
}