#include <Eigen/Dense>
#include <iostream>
#include "StencilHandler.h"
int main()
{
    StencilHandler* hdler = new StencilHandler[1];
    hdler->searchCentralStencil2D("naca0012.grd");
    std::cout << "Hello Eigen\n";
}