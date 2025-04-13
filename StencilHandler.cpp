#include "StencilHandler.h"
#include <fstream>
#include <vector>
#include <set>
#include <math.h>
#include "CellBrick.h"
#include <hdf5.h>
#include <Eigen/Dense>
void StencilHandler::searchCentralStencil2D(std::string filename)
{
    const int ndim = 2;
    const int cellBase = 1;
    int to_throw_away;
	std::ifstream fin;
	fin.open(filename);
	if (!fin.is_open())
	{
		throw std::runtime_error("Open file Failure\n");
	}
	int temp_nnodes,temp_nedges,temp_ncells;
			
    fin >> to_throw_away >> to_throw_away >> to_throw_away;
	fin >> temp_nnodes >> temp_nedges >> temp_ncells >> to_throw_away >> to_throw_away;

	double *xyzs = new double[ndim * temp_nnodes];

	for (int i = 0; i < temp_nnodes; i++)
	{
		for (int j = 0; j < ndim; j++)
		{
			fin >> xyzs[ndim * i + j];
		}
	}
	int *edge2Cell = new int[2 * temp_nedges];
    int temp_counts;

	int *temp_nodes3D = NULL;
	int temp_nodes[4];
    Brick4Cobalt2D * temp_cells = nullptr;

    if (ndim == 2)
	{
		temp_cells = new Brick4Cobalt2D[temp_ncells];
		for (int i = 0; i < temp_nedges; i++)
		{
			fin >> temp_counts;
			for (int j = 0; j < temp_counts; j++)
			{
				fin >> temp_nodes[j];
			}
			fin >> edge2Cell[2 * i] >> edge2Cell[2 * i + 1];

			temp_cells[edge2Cell[2 * i] - cellBase].insert(temp_nodes, (double*)&edge2Cell[2*i], 2);
					
			if (edge2Cell[2 * i + 1] > 0)
			{
				temp_cells[edge2Cell[2 * i + 1] - cellBase].insert(temp_nodes,(double*)&edge2Cell[2*i+1], 2);
			}
		}
	}
    fin.close();

    //Allocation of this class' member

    nodeOffsetEachCell = new int[temp_ncells+1];
    nodeOffsetEachCell[0] = 0;
    for(int i = 0;i < temp_ncells;i++)
    {
        nodeOffsetEachCell[i+1] = nodeOffsetEachCell[i] + temp_cells[i].sizeIs();
    }
    
    cellOffsetEachNode = new int[temp_nnodes+1];
    for(int i = 0;i < temp_nnodes+1;i++)
    {
        cellOffsetEachNode[i] = 0;
    }
    //Now counting how many each cell and each node contains
    int *int_iter;
    for(int i = 0;i < temp_ncells;i++)
    {
        int_iter = temp_cells[i].begin();
        for(int j = 0;j < temp_cells[i].sizeIs();j++)
        {
            temp_counts = *int_iter - 1;//The vert ind starts from 1 in temp_cells

            cellOffsetEachNode[temp_counts + 1]++;
            int_iter++;
        }
    }
    for(int i = 0;i < temp_nnodes;i++)
    {
        cellOffsetEachNode[i+1] = cellOffsetEachNode[i] + cellOffsetEachNode[i+1];
    }

    //Now fill in each node of each cell 
    std::vector<int> nodeEachCellOffset(temp_ncells,0);
    std::vector<int> cellEachNodeOffset(temp_nnodes,0);
    nodeEachCell = new int[nodeOffsetEachCell[temp_ncells]];
    cellEachNode = new int[cellOffsetEachNode[temp_nnodes]];
    for(int i = 0;i < temp_ncells;i++)
    {
        int_iter = temp_cells[i].begin();
        for(int j = 0;j < temp_cells[i].sizeIs();j++)
        {
            temp_counts = *int_iter - 1;

            nodeEachCell[nodeOffsetEachCell[i]+nodeEachCellOffset[i]] = temp_counts;
            cellEachNode[cellOffsetEachNode[temp_counts]+cellEachNodeOffset[temp_counts]] = i;
            nodeEachCellOffset[i]++;
            cellEachNodeOffset[temp_counts]++;
            int_iter++;
        }
    }

    //Now that we have the node 2 cell and cell 2 node mapping, time to find stencil

    std::cout << "End of the game\n";

    for(int i = 0;i < temp_ncells;i++)
    {
        std::cout << "cell " << i << ":" << nodeOffsetEachCell[i] << " to " << nodeOffsetEachCell[i+1] << "\n";
        for(int j = nodeOffsetEachCell[i];j < nodeOffsetEachCell[i+1];j++)
        {
            std::cout << nodeEachCell[j] << "\n";

        }
        
    }

    for(int i = 0;i < temp_nnodes;i++)
    {
        std::cout << "node " << i << ":" << cellOffsetEachNode[i] << " to " << cellOffsetEachNode[i+1] << "\n";
        for(int j = cellOffsetEachNode[i];j < cellOffsetEachNode[i+1];j++)
        {
            std::cout << cellEachNode[j] << "\n";

        }
        
    }

    //Rerun fin for debug
    {
    std::ifstream fin;
	fin.open(filename);
	if (!fin.is_open())
	{
		throw std::runtime_error("Open file Failure\n");
	}
	int temp_nnodes,temp_nedges,temp_ncells;
			
    fin >> to_throw_away >> to_throw_away >> to_throw_away;
	fin >> temp_nnodes >> temp_nedges >> temp_ncells >> to_throw_away >> to_throw_away;

	double *xyzs = new double[ndim * temp_nnodes];

	for (int i = 0; i < temp_nnodes; i++)
	{
		for (int j = 0; j < ndim; j++)
		{
			fin >> xyzs[ndim * i + j];
		}
	}
	int *edge2Cell = new int[2 * temp_nedges];
    int temp_counts;

	int *temp_nodes3D = NULL;
	int temp_nodes[4];
    Brick4Cobalt2D * temp_cells = nullptr;

    if (ndim == 2)
	{
		temp_cells = new Brick4Cobalt2D[temp_ncells];
		for (int i = 0; i < temp_nedges; i++)
		{
			fin >> temp_counts;
			for (int j = 0; j < temp_counts; j++)
			{
				fin >> temp_nodes[j];
			}
			fin >> edge2Cell[2 * i] >> edge2Cell[2 * i + 1];

            

            for (int j = 0; j < temp_counts; j++)
            {
                int found_vert_in_left = 0;//false;
                int found_vert_in_right = 0;//false;
                int found_left_in_vert = 0;//false;
                int found_right_in_vert = 0;//false;
                
                int nodeId = temp_nodes[j] - 1;
                int leftid = edge2Cell[2 * i] - 1;
                int rightid = edge2Cell[2 * i + 1] - 1;

                for(int m = nodeOffsetEachCell[leftid];m < nodeOffsetEachCell[leftid+1];m++)
                {
                    found_vert_in_left += (nodeEachCell[m] == nodeId);
                }
                if(rightid >= 0)
                {
                    for(int m = nodeOffsetEachCell[rightid];m < nodeOffsetEachCell[rightid+1];m++)
                    {
                        found_vert_in_right += (nodeEachCell[m] == nodeId);
                    }
                }
                for(int m = cellOffsetEachNode[nodeId];m < cellOffsetEachNode[nodeId+1];m++)
                {
                    found_left_in_vert += (cellEachNode[m] == leftid);
                    if(rightid >= 0)
                    {
                        found_right_in_vert += (cellEachNode[m] == rightid);
                    }
                }

                if(found_vert_in_left!=1)
                {
                    throw std::runtime_error("not found vert in left\n");
                }
                if(found_vert_in_right!=1)
                {
                    if(rightid>=0)
                    {
                        throw std::runtime_error("not found vert in right\n");
                    }
                }
                if(found_left_in_vert!=1)
                {
                    throw std::runtime_error("not found left in vert\n");
                }
                if(found_right_in_vert!=1)
                {
                    if(rightid>=0)
                    {
                        throw std::runtime_error("not found right in vert\n");
                    }
                }


            }

			temp_cells[edge2Cell[2 * i] - cellBase].insert(temp_nodes, (double*)&edge2Cell[2*i], 2);
					
			if (edge2Cell[2 * i + 1] > 0)
			{
				temp_cells[edge2Cell[2 * i + 1] - cellBase].insert(temp_nodes,(double*)&edge2Cell[2*i+1], 2);
			}
		}
	}
    fin.close();
    delete[] xyzs;
    delete[] edge2Cell;
    delete[] temp_cells;
    }


    //Now I am pretty sure it works, time to find stencil for each node

    //Now let's say I will have first all neighbouring cells of each node
    //The stencil for a node is found by:
    //First identify all cells that are connected to it.
    //THen for each cell, tag its neighbour cells.

    int* edgeOffsetEachCell = new int[temp_ncells+1];
    for(int i = 0;i < temp_ncells+1;i++)
    {
        edgeOffsetEachCell[i] = 0;
    }
    for(int i = 0;i < temp_nedges;i++)
    {
        int lC = edge2Cell[2*i] - 1;
        int rC = edge2Cell[2*i+1] - 1;
        edgeOffsetEachCell[lC+1]++;
        if(rC>=0)
        {
            edgeOffsetEachCell[rC+1]++;
        }
    }
    for(int i = 0;i < temp_ncells;i++)
    {
        edgeOffsetEachCell[i+1] = edgeOffsetEachCell[i] + edgeOffsetEachCell[i+1];
    }
    
    int* edgeEachCell = new int[edgeOffsetEachCell[temp_ncells]];
    // std::cout << "debug 1\n";
    std::vector<int> edgeEachCellOffset(temp_ncells,0);
    for(int i = 0;i < temp_nedges;i++)
    {
        int lC = edge2Cell[2*i] - 1;
        int rC = edge2Cell[2*i+1] - 1;
        edgeEachCell[edgeOffsetEachCell[lC]+edgeEachCellOffset[lC]] = i;
        edgeEachCellOffset[lC]++;
        if(rC>=0)
        {
            edgeEachCell[edgeOffsetEachCell[rC]+edgeEachCellOffset[rC]] = i;
            edgeEachCellOffset[rC]++;
        }
    }
    // std::cout << "debug 2\n";
    //Rerun fin for debug
    {
    if (ndim == 2)
	{
        for (int i = 0; i < temp_nedges; i++)
        {
            int lC = edge2Cell[2*i] - 1;
            int rC = edge2Cell[2*i+1] - 1;
            int found_edge_in_left = 0;
            int found_edge_in_right = 0;
            for(int m = edgeOffsetEachCell[lC];m < edgeOffsetEachCell[lC+1];m++)
            {
                found_edge_in_left += (edgeEachCell[m] == i);
            }
            if(rC >= 0)
            {
                for(int m = edgeOffsetEachCell[rC];m < edgeOffsetEachCell[rC+1];m++)
                {
                    found_edge_in_right += (edgeEachCell[m] == i);
                }
            }
            if(found_edge_in_left!=1)
            {
                throw std::runtime_error("not found edge in left\n");
            }
            if(found_edge_in_right!=1)
            {
                if(rC>=0)
                {
                    throw std::runtime_error("not found edge in right\n");
                }
            }
        }
		
	}

    }

    //Here is the part where the stencil is found
    std::set<int,std::less<int>> neighbouring_cells;
    std::set<int,std::less<int>> donating_nodes;
    std::vector<int> stencil_cells_each_vert;
    std::vector<int> stencil_cells_offset_each_vert;
    stencil_cells_offset_each_vert.push_back(0);
    for(int i = 0;i < temp_nnodes;i++)
    {
        //For node i, find all its cells
        
        neighbouring_cells.clear();
        for(int j = cellOffsetEachNode[i];j < cellOffsetEachNode[i+1];j++)
        {
            int cellId = cellEachNode[j];
            neighbouring_cells.insert(cellId);
            //Now for each neighbouring cell, find its vertices.
        }
        for(auto iter = neighbouring_cells.begin();iter != neighbouring_cells.end();iter++)
        {
            stencil_cells_each_vert.push_back(*iter);
        }
        stencil_cells_offset_each_vert.push_back(stencil_cells_each_vert.size());
    }


    //Now it's time to write into the HDF5 Format file so we could find its 
    {
        hid_t file_id = H5Fcreate("Stencil", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        //These are some basic attributes describing this mesh
        {
            hid_t attr_id = H5Screate(H5S_SCALAR);
            hid_t attr = H5Acreate2(file_id, "nnodes", H5T_NATIVE_INT,attr_id, H5P_DEFAULT, H5P_DEFAULT);
            H5Awrite(attr, H5T_NATIVE_INT, &nnodes);

            attr_id = H5Screate(H5S_SCALAR);
            attr = H5Acreate2(file_id, "ncells", H5T_NATIVE_INT,attr_id, H5P_DEFAULT, H5P_DEFAULT);
            H5Awrite(attr, H5T_NATIVE_INT, &ncells);
            
            attr_id = H5Screate(H5S_SCALAR);
            attr = H5Acreate2(file_id, "nedges", H5T_NATIVE_INT,attr_id, H5P_DEFAULT, H5P_DEFAULT);
            H5Awrite(attr, H5T_NATIVE_INT, &nedges);
        }

        //Let's first store all verts' coordinates
        {
            hsize_t dim[2] = {temp_nnodes,ndim};
            hid_t dataspace_id = H5Screate_simple(2, dim, NULL);
            std::string dataset_name = "coordinates";
            hid_t dataset_id = H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_DOUBLE,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, NULL, dataspace_id, H5P_DEFAULT, xyzs);
            H5Sclose(dataspace_id);
            H5Dclose(dataset_id);
        }
        //cell-vert relations
        {
            {
                hsize_t dim[1] = {temp_ncells+1};
                hid_t dataspace_id = H5Screate_simple(1, dim, NULL);
                std::string dataset_name = "vertOffsetEachCell";
                hid_t dataset_id = H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_INT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Dwrite(dataset_id, H5T_NATIVE_INT, NULL, dataspace_id, H5P_DEFAULT, nodeOffsetEachCell);
                H5Sclose(dataspace_id);
                H5Dclose(dataset_id);
            }
            {
                hsize_t dim[1] = {nodeOffsetEachCell[temp_ncells]};
                hid_t dataspace_id = H5Screate_simple(1, dim, NULL);
                std::string dataset_name = "vertEachCell";
                hid_t dataset_id = H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_INT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Dwrite(dataset_id, H5T_NATIVE_INT, NULL, dataspace_id, H5P_DEFAULT, nodeEachCell);
                H5Sclose(dataspace_id);
                H5Dclose(dataset_id);

            }
        }
        //vert-cell relations
        {
            {
                hsize_t dim[1] = {temp_nnodes+1};
                hid_t dataspace_id = H5Screate_simple(1, dim, NULL);
                std::string dataset_name = "cellOffsetEachVert";
                hid_t dataset_id = H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_INT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Dwrite(dataset_id, H5T_NATIVE_INT, NULL, dataspace_id, H5P_DEFAULT, cellOffsetEachNode);
                H5Sclose(dataspace_id);
                H5Dclose(dataset_id);
            }
            {
                hsize_t dim[1] = {cellOffsetEachNode[temp_nnodes]};
                hid_t dataspace_id = H5Screate_simple(1, dim, NULL);
                std::string dataset_name = "cellEachVert";
                hid_t dataset_id = H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_INT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Dwrite(dataset_id, H5T_NATIVE_INT, NULL, dataspace_id, H5P_DEFAULT, cellEachNode);
                H5Sclose(dataspace_id);
                H5Dclose(dataset_id);
            }
        }
        //cell-edge relations
        {
            {
                hsize_t dim[1] = {temp_ncells+1};
                hid_t dataspace_id = H5Screate_simple(1, dim, NULL);
                std::string dataset_name = "edgeOffsetEachCell";
                hid_t dataset_id = H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_INT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Dwrite(dataset_id, H5T_NATIVE_INT, NULL, dataspace_id, H5P_DEFAULT, edgeOffsetEachCell);
                H5Sclose(dataspace_id);
                H5Dclose(dataset_id);
            }
            {
                hsize_t dim[1] = {edgeOffsetEachCell[temp_ncells]};
                hid_t dataspace_id = H5Screate_simple(1, dim, NULL);
                std::string dataset_name = "edgeEachCell";
                hid_t dataset_id = H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_INT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Dwrite(dataset_id, H5T_NATIVE_INT, NULL, dataspace_id, H5P_DEFAULT, edgeEachCell);
                H5Sclose(dataspace_id);
                H5Dclose(dataset_id);
            }
        }
        //cell-edge relations
        {
            {
                hsize_t dim[1] = {2*temp_nedges};
                hid_t dataspace_id = H5Screate_simple(1, dim, NULL);
                std::string dataset_name = "cellEachEdge";
                hid_t dataset_id = H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_INT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Dwrite(dataset_id, H5T_NATIVE_INT, NULL, dataspace_id, H5P_DEFAULT, edge2Cell);
                H5Sclose(dataspace_id);
                H5Dclose(dataset_id);
            }
        }
        //edge-vert relations
        {
            {
                int max_size = 1 + 2 + (ndim == 2 ? 2:4);
                hsize_t dim[2] = {temp_nedges,max_size};
                hid_t dataspace_id = H5Screate_simple(2, dim, NULL);
                std::string dataset_name = "vertEachEdge";
                hid_t dataset_id = H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_INT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Dwrite(dataset_id, H5T_NATIVE_INT, NULL, dataspace_id, H5P_DEFAULT, edge2Cell);


                hsize_t count[2] = {1, max_size};
                hsize_t offset[2] = {0,0};
                hsize_t stride[2] = {1,1};
                hsize_t block[2] = {1,1};

                hid_t memspace_id = H5Screate_simple(2, count, NULL);


                //Rerun fin for debug
                {
                std::ifstream fin;
	            fin.open(filename);
	            if (!fin.is_open())
	            {
	            	throw std::runtime_error("Open file Failure\n");
	            }
	            int temp_nnodes,temp_nedges,temp_ncells;

                fin >> to_throw_away >> to_throw_away >> to_throw_away;
	            fin >> temp_nnodes >> temp_nedges >> temp_ncells >> to_throw_away >> to_throw_away;

	            double to_throw_away_double;
	            for (int i = 0; i < temp_nnodes; i++)
	            {
	            	for (int j = 0; j < ndim; j++)
	            	{
	            		fin >> to_throw_away_double;
	            	}
	            }

		        for (int i = 0; i < temp_nedges; i++)
		        {
                    int temp_nodes[max_size] = {-1};
			        fin >> temp_nodes[0];
			        for (int j = 0; j < temp_nodes[0]; j++)
			        {
				        fin >> temp_nodes[j+1];
                        temp_nodes[j+1] -= 1;
			        }
			        fin >> temp_nodes[max_size-2] >> temp_nodes[max_size-1];
                    temp_nodes[max_size-2] -= 1;
                    temp_nodes[max_size-1] -= 1;
                    offset[0] = i;
                    hid_t status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, stride, count, block);

                    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT, temp_nodes);
                }
                H5Sclose(memspace_id);
                H5Sclose(dataspace_id);
                H5Dclose(dataset_id);
                fin.close();
                }
            }
        }
        // std::cin.get();
        //-vert-stencil relations
        {
            // std::vector<int> stencil_cells_each_vert;
            // std::vector<int> stencil_cells_offset_each_vert;
            {
                hsize_t dim[1] = {temp_nnodes+1};
                hid_t dataspace_id = H5Screate_simple(1, dim, NULL);
                std::string dataset_name = "stencilOffsetEachVert";
                hid_t dataset_id = H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_INT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Dwrite(dataset_id, H5T_NATIVE_INT, NULL, dataspace_id, H5P_DEFAULT, stencil_cells_offset_each_vert.data());
                H5Sclose(dataspace_id);
                H5Dclose(dataset_id);
            }
            // std::cin.get();
            {
                hsize_t dim[1] = {stencil_cells_offset_each_vert[temp_nnodes]};
                hid_t dataspace_id = H5Screate_simple(1, dim, NULL);
                std::string dataset_name = "stencilEachVert";
                hid_t dataset_id = H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_INT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Dwrite(dataset_id, H5T_NATIVE_INT, NULL, dataspace_id, H5P_DEFAULT, stencil_cells_each_vert.data());
                H5Sclose(dataspace_id);
                H5Dclose(dataset_id);
            }
        }
        std::cout << "Hello Eigen ended\n";

        //Okay, it's worth noting that now its time to compute the matrix 
        std::vector<int> stencil_constant_ptr(2*temp_nnodes+2,0);
        for(int i = 0;i < temp_nnodes;i++)
        {
            int nNeighbours = stencil_cells_offset_each_vert[i+1] - stencil_cells_offset_each_vert[i];
            stencil_constant_ptr[2*(i+1)] = nNeighbours;
            stencil_constant_ptr[2*(i+1)+1] = stencil_constant_ptr[2*i+1] + nNeighbours*nNeighbours+nNeighbours;


        }
        double* stencil_matrix_data = new double[stencil_constant_ptr[2*(temp_nnodes+1)-1]];
        for(int i = 0;i < temp_nnodes;i++)
        {
            int nNeighbours = stencil_cells_offset_each_vert[i+1] - stencil_cells_offset_each_vert[i];
            //We will store a [nxn] matrix and a [nx1] vector
            double* matrix = &stencil_matrix_data[stencil_constant_ptr[2*i]];
            //Now compute the distance between each 2 of the neighbours
            double max_distance = -10000;
            double min_distance = 10000;
            double mean_distance = 0;
            
            double* cell_center_coords = new double[ndim*nNeighbours];
            for(int j = 0;j < ndim*nNeighbours;j++)
            {
                cell_center_coords[j] = 0;
            }
            for(int j = 0;j < nNeighbours;j++)
            {
                int j_id = stencil_cells_each_vert[stencil_cells_offset_each_vert[i]+j];
                for(int l = nodeOffsetEachCell[j_id];l < nodeOffsetEachCell[j_id+1];l++)
                {
                    int vert_id = nodeEachCell[l];
                    for(int k = 0;k < ndim;k++)
                    {
                        cell_center_coords[j*ndim+k] += xyzs[vert_id*ndim+k];
                    }
                }
                for(int k = 0;k < ndim;k++)
                {
                    cell_center_coords[j*ndim+k] /= nodeOffsetEachCell[j_id+1] - nodeOffsetEachCell[j_id];
                }
                for(int l = 0;l <= j;l++)
                {
                    // int l_id = stencil_cells_each_vert[stencil_cells_offset_each_vert[i]+l];
                    double distance = 0;
                    for(int k = 0;k < ndim;k++)
                    {
                        double dis = (cell_center_coords[ndim*j+k]-cell_center_coords[ndim*l+k]);
                        distance += dis*dis;
                    }
                    max_distance = std::max(max_distance,distance);
                    min_distance = std::min(min_distance,distance);
                    // distance = sqrt(distance);
                    matrix[l*nNeighbours+j] = distance;
                    matrix[j*nNeighbours+l] = distance;

                    //Now it's time to decide which kernel to use, I will go with Gaussian for now.
                    // double kernel_val = exp(-distance/(epsilon*epsilon));
                }
                double distance = 0;
                for(int k = 0;k < ndim;k++)
                {
                    double dis = (cell_center_coords[ndim*j+k]-xyzs[ndim*i+k]);
                    distance += dis*dis;
                }
                mean_distance += distance;
                matrix[nNeighbours*nNeighbours+j] = distance;
                std::cout << "distance to target: " << j << ":" <<distance << "\n";
            }
            for(int j = 0;j < nNeighbours;j++)
            {
                std::cout << "cell " << j << "center:\n";
                std::cout << cell_center_coords[ndim*j] << " , " << cell_center_coords[ndim*j+1] << "\n";
            }
            delete[] cell_center_coords;
            mean_distance = mean_distance/nNeighbours;
            for(int j = 0;j < nNeighbours*nNeighbours+nNeighbours;j++)
            {
                matrix[j] = exp(-matrix[j]/mean_distance);
            }
            for(int j = nNeighbours*nNeighbours;j<nNeighbours*nNeighbours+nNeighbours;j++)
            {
                std::cout << matrix[j] << "\n";
            }
            std::cout << "node " << i << " out of " << temp_nnodes << "\n";
            std::cout << "before eigen\n";
            // for(int j = 0;j < nNeighbours;j++)
            // {
            //     int j_id = stencil_cells_each_vert[stencil_cells_offset_each_vert[i]+j];
            //     std::cout << j << "th id is " << j_id << "\n";
            //     std::cout << xyzs[ndim*j_id+0] << "\t" << xyzs[ndim*j_id+1] << "\n";
            // }
            // std::cin.get();
            Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> mat;
            mat.resize(nNeighbours,nNeighbours);
            for(int j = 0;j < nNeighbours;j++)
            {
                for(int k = 0;k < nNeighbours;k++)
                {
                    mat(j,k) = matrix[j*nNeighbours+k];
                }
            }
            std::cout << mat << "\n";
            std::cout << mat.determinant() << "\n";
            if(mat.determinant() == 0)
            {
                throw std::runtime_error("invertible which is impossible\n");
            }
            else
            {
                Eigen::MatrixXd mat_inv = mat.inverse();
                std::cout << mat_inv << "\n";
                for(int j = 0;j < nNeighbours;j++)
                {
                    for(int k = 0;k < nNeighbours;k++)
                    {
                        matrix[j*nNeighbours+k] = mat_inv(j,k);
                    }
                }

            }
            
            //Now it's time to do an inverse here.

        }

        {
            // std::vector<int> stencil_cells_each_vert;
            // std::vector<int> stencil_cells_offset_each_vert;
            {
                hsize_t dim[1] = {2*(temp_nnodes+1)};
                hid_t dataspace_id = H5Screate_simple(1, dim, NULL);
                std::string dataset_name = "stencil_constant_ptr";
                hid_t dataset_id = H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_INT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Dwrite(dataset_id, H5T_NATIVE_INT, NULL, dataspace_id, H5P_DEFAULT, stencil_constant_ptr.data());
                H5Sclose(dataspace_id);
                H5Dclose(dataset_id);
            }
            // std::cin.get();
            {
                hsize_t dim[1] = {stencil_constant_ptr[2*(temp_nnodes)+1]};
                hid_t dataspace_id = H5Screate_simple(1, dim, NULL);
                std::string dataset_name = "stencil_constant_data";
                hid_t dataset_id = H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_DOUBLE,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, NULL, dataspace_id, H5P_DEFAULT, stencil_matrix_data);
                H5Sclose(dataspace_id);
                H5Dclose(dataset_id);
            }
        }

    }
}


void StencilHandler::partitioning()
{
    int num_proc,cur_proc;
    MPI_Comm_rank(MPI_COMM_WORLD,&cur_proc);
    MPI_Comm_size(MPI_COMM_WORLD,&num_proc);

    MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Info info = MPI_INFO_NULL;
	hid_t plist_id=H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id,comm,info);
	hid_t file_id = H5Fopen("Stencil",H5P_DEFAULT,H5FD_MPIO_INDEPENDENT);
	hid_t status = H5Pclose(plist_id);

    int nnodes,ncells,nedges;
    
    {
        hid_t attr_id =  H5Aopen(file_id, "nnodes", H5P_DEFAULT);
        H5Aread(attr_id, H5T_NATIVE_INT, &nnodes);

        attr_id =  H5Aopen(file_id, "ncells", H5P_DEFAULT);
        H5Aread(attr_id, H5T_NATIVE_INT, &ncells);

        attr_id =  H5Aopen(file_id, "nedges", H5P_DEFAULT);
        H5Aread(attr_id, H5T_NATIVE_INT, &nedges);
    }

    //Now read from the HDF5 file the layout of cells to do METIS
    {
        std::string dataset_name = "coordinates";
        
    }

    {
        std::string dataset_name = "vertOffsetEachCell";
    }
    {
        std::string dataset_name = "vertEachCell";
    }
    {
        std::string dataset_name = "cellOffsetEachVert";
    }
    {
        std::string dataset_name = "cellEachVert";
    }





}