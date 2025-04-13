#pragma once
#include <string>
#include <Scene.h>
#include <Entity.h>
#include <tioga.h>
#include <FileManager.h>

template<int ndim>
class MeshBlockInfo
{
public:
	int ntypes = 0;
	int* eachNodeCount = NULL;//Unchanged
	int* eachCellCount = NULL;//Unchanged
	int nnodes = 0;//
	int nedges = 0;
	int ncells = 0;//
	double* rxyz = NULL;//

	int* ibl = NULL;//Will be modifies by tioga
	
    int* wbc = NULL;//Unchanged
	int* obc = NULL;//Unchanged

	int nwbc = 0;
	int nobc = 0;
	 
    
    int** vconn = NULL;//Unchanged
 
	MeshBlockInfo()
	{

    }
	~MeshBlockInfo()//This might not be that necessary cuz tioga will manage it
	{

    }

};
template<int ndim>
class MeshLoader
{
public:
    const int cellBase = 1;
    const int nodeBase = 1;
    int signWall = -1;
    int signOver = -2;
	int signSymm = -3;

	MeshBlockInfo<ndim> cur_bi;

	int meshtagCounter = 1;

	const std::string folder_name = "InterMesh";
	std::string InterFileName;

    MeshLoader()
    {
		myCreateDirectory(folder_name);
    }
    ~MeshLoader()
    {

    }

    void preprocessFiles(std::string mesh_path);

    void loadFiles(std::string mesh_path, Glamour::Entity entity);
    
    void selectReadForUnstruct(std::string mesh_path, Glamour::Entity entity);
	void selectReadForUniTioga(std::string mesh_path, Glamour::Entity entity);

	void readWallBoundary(std::string mesh_path, Glamour::Entity entity);
};



