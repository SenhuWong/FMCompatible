#pragma once
#include <string>

class StencilHandler
{
    int nnodes = 0;
    int nedges = 0;
    int ncells = 0;

    int* nodeEachCell = nullptr;
    int* nodeOffsetEachCell = nullptr;
    
    //two 1-d array indicating cells tied to each vert
    int* cellEachNode = nullptr;
    int* cellOffsetEachNode = nullptr;

    std::string folder_name = "DefaultStencilFolder";
public:
    void searchCentralStencil2D(std::string mesh_path);

    void partitioning();

    

};