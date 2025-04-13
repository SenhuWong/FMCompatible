#include <hdf5.h>
#include <H5FDmpi.h>
#include <fstream>
#include <metis.h>
#include <numeric>
#include "MeshLoader.h"

#include <CellBrick.h>
#include <LocalGlobalMap.h>

#include <UnstructBlock.h>
#include <Components.h>
#include <tioga.h>
template<int ndim>
void MeshLoader<ndim>::preprocessFiles(std::string filename)
{
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
#ifdef DIVIDED_BY_THOUSAND
	for (int i = 0; i < temp_nnodes * ndim;i++)
	{
		xyzs[i] = xyzs[i] / 1000;
	}
#endif	
	int *edge2Cell = new int[2 * temp_nedges];
	int temp_counts;

	int *temp_nodes3D = NULL;
	int temp_nodes[4]; 
	Brick4Cobalt3D *temp_cells3D = NULL;
	Brick4Cobalt2D *temp_cells2D = NULL;

	if (ndim == 2)
	{
		temp_cells2D = new Brick4Cobalt2D[temp_ncells];
		for (int i = 0; i < temp_nedges; i++)
		{
			fin >> temp_counts;
			for (int j = 0; j < temp_counts; j++)
			{
				fin >> temp_nodes[j];
			}
			fin >> edge2Cell[2 * i] >> edge2Cell[2 * i + 1];

			temp_cells2D[edge2Cell[2 * i] - cellBase].insert(temp_nodes, (double*)&edge2Cell[2*i], 2);
					
			if (edge2Cell[2 * i + 1] > 0)
			{
				temp_cells2D[edge2Cell[2 * i + 1] - cellBase].insert(temp_nodes,(double*)&edge2Cell[2*i+1], 2);
			}
		}
	}
	else if (ndim == 3)
	{
		temp_nodes3D = new int[4 * temp_nedges];//4 times the num of edges
		temp_cells3D = new Brick4Cobalt3D[temp_ncells];
		for (int i = 0; i < temp_nedges; i++)
		{
			int indx = 4*i;

			fin >> temp_counts;
			for (int j = 0; j < temp_counts; j++)
			{
				fin >> temp_nodes3D[indx + j];
			}
			fin >> edge2Cell[2 * i] >> edge2Cell[2 * i + 1];
			temp_cells3D[edge2Cell[2 * i] - cellBase].insert(&temp_nodes3D[indx],1, xyzs, temp_counts);
			if (edge2Cell[2 * i + 1] > 0)
			{
				temp_cells3D[edge2Cell[2 * i + 1] - cellBase].insert(&temp_nodes3D[indx],-1, xyzs, temp_counts);
			}
        }
		for(int i = 0;i<temp_ncells;i++)
		{
			temp_cells3D[i].reOrdering(xyzs);
			temp_cells3D[i].checkRepeat();
		}
		delete[] temp_nodes3D;
	}
	fin.close();

    int num_proc,cur_proc;
    MPI_Comm_rank(MPI_COMM_WORLD,&cur_proc);
    MPI_Comm_size(MPI_COMM_WORLD,&num_proc);

	idx_t num_elements = temp_ncells;
	idx_t num_nodes = temp_nnodes;
	idx_t index_offset = 0;
	std::vector<idx_t> elemPtr;
	std::vector<idx_t> elemInd;
	int *int_iter;
	for (int i = 0; i < temp_ncells; i++)
	{
		elemPtr.push_back(index_offset);
		if (ndim == 2)
		{
			index_offset += static_cast<idx_t>(temp_cells2D[i].sizeIs());
			int_iter = temp_cells2D[i].begin();
			for (int j = 0; j < temp_cells2D[i].sizeIs(); j++)
			{
				temp_counts = *int_iter;
				// temp_counts should start from 1;
				elemInd.push_back(static_cast<idx_t>(temp_counts - 1)); // Splitting index start from 0
				int_iter++;
			}
		}
		else if (ndim == 3)
		{
			index_offset += static_cast<idx_t>(temp_cells3D[i].sizeIs());
			int_iter = temp_cells3D[i].begin();
			for (int j = 0; j < temp_cells3D[i].sizeIs(); j++)
			{
				temp_counts = *int_iter;
				// temp_counts should start from 1;
				elemInd.push_back(static_cast<idx_t>(temp_counts - 1)); // Splitting index start from 0
				int_iter++;
			}
		}
	}
	elemPtr.push_back(static_cast<idx_t>(index_offset));
	std::vector<idx_t> cellPartition(num_elements, -1);
	idx_t NCOMM = 2; // 2D case
	idx_t NPART = static_cast<idx_t>(num_proc);
	idx_t objVal = 0;
	std::vector<idx_t> nodePartition(num_nodes, 0);
	idx_t options[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options);
	options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
	options[METIS_OPTION_CONTIG] = 1;
	options[METIS_OPTION_NUMBERING] = 0; // Start from 0;
	std::cout << "Calling Mestis part\n";
	if (num_proc > 1)
	{
		METIS_PartMeshDual(&num_elements, &num_nodes, elemPtr.data(), elemInd.data(), NULL, NULL, &NCOMM, &NPART, NULL, options, &objVal, cellPartition.data(), nodePartition.data());
	}
	else
	{
		for (auto iter = cellPartition.begin(); iter != cellPartition.end(); iter++)
		{
			*iter = cur_proc;
		}
	}
	std::cout << "Called Mestis part\n";

	int *eachCount = new int[num_proc];
	for (int i = 0; i < num_proc; i++)
	{
		eachCount[i] = 0;
	}
			
	for (auto itera = cellPartition.begin(); itera != cellPartition.end(); itera++)
	{
				
		// std::cout << *itera << '\n';
		++eachCount[*(itera)];
	}
			
	MPI_Request *send_request = new MPI_Request[num_proc - 1];
	MPI_Status *send_status = new MPI_Status[num_proc - 1];
			
			
	delete[] xyzs;
	delete[] edge2Cell;

	delete[] eachCount;
			
	if(cur_proc!=0)
	{
		return;

	}
	// outout a cell file so that we don't need to reinsert to get all the cellsS

	hid_t file_id = H5Fcreate(InterFileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
				
	int eachCellSize = 1 + 1 + (ndim == 2 ? 4 : 8);
	hsize_t dim[2] = {temp_ncells, eachCellSize};
	hid_t dataspace_id = H5Screate_simple(2, dim, NULL);
	std::string dataset_name = "cell_distribution";
	int buffer[10];
	hid_t dset_id = H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	hsize_t count[2] = {1, eachCellSize};
	hsize_t offset[2] = {0, 0};
	hsize_t stride[2] = {1, 1};
	hsize_t block[2] = {1, 1};
				
	hid_t memspace_id = H5Screate_simple(2, count, NULL);

	if (ndim == 2)
	{
		for (int i = 0; i < temp_ncells; i++)
		{
			offset[0] = i;
			hid_t status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, stride, count, block);
			buffer[0] = temp_cells2D[i].sizeIs();
			buffer[1] = cellPartition[i];
			int_iter = temp_cells2D[i].begin();
			for (int j = 0; j < temp_cells2D[i].sizeIs(); j++)
			{
				buffer[2 + j] = (*int_iter);
				++int_iter;
			}
			for (int j = temp_cells2D[i].sizeIs(); j < 4; j++)
			{
				buffer[2 + j] = -1;
			}
			status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT, buffer);
		}
	}
	else if (ndim == 3)
	{
		for (int i = 0; i < temp_ncells; i++)
		{
			offset[0] = i;
			hid_t status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, stride, count, block);
			buffer[0] = temp_cells3D[i].sizeIs();
			buffer[1] = cellPartition[i];
			int_iter = temp_cells3D[i].begin();
			for (int j = 0; j < temp_cells3D[i].sizeIs(); j++)
			{
				buffer[2 + j] = (*int_iter);
				++int_iter;
			}
			for (int j = temp_cells3D[i].sizeIs(); j < 8; j++)
			{
				buffer[2 + j] = -1;
			}
			status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT, buffer);
		}
	}
	H5Sclose(memspace_id);
	H5Dclose(dset_id);
	H5Sclose(dataspace_id);
	H5Fclose(file_id);
	if (temp_cells2D)
		delete[] temp_cells2D;
	if (temp_cells3D)
		delete[] temp_cells3D;
	delete[] send_request;
	delete[] send_status;
}

template<int ndim>
void MeshLoader<ndim>::loadFiles(std::string mesh_path,Glamour::Entity entity)
{
    //First check if there is a h5 file.
	std::size_t botDirPos =  mesh_path.find_last_of("/");
	std::string mesh_name = mesh_path.substr(botDirPos+1,mesh_path.length());

	InterFileName = folder_name + "/_cell_" + std::to_string(UNSTRUCT_MPI::Get().Size()) + "_" + mesh_name + ".h5";
    if(!myFileExists(InterFileName))
    {
		std::cout<<"Doing Metis Partition\n";
        preprocessFiles(mesh_path);
    }
	else
	{
		std::cout<<"Skipping Metis Partition\n";
	}
    //Then read for Unstruct and Tioga, these 2 might need to be merged together
    selectReadForUnstruct(mesh_path,entity);
    selectReadForUniTioga(mesh_path,entity);
	readWallBoundary(mesh_path,entity);
    //Finally registering in tioga and initializing in meshComponent
    TIOGA::tioga::Get().registerGridData(meshtagCounter,
        cur_bi.nnodes,cur_bi.rxyz,cur_bi.ibl,cur_bi.nwbc,cur_bi.nobc,
        cur_bi.wbc,cur_bi.obc,cur_bi.ntypes,cur_bi.eachNodeCount,cur_bi.eachCellCount,
        cur_bi.vconn);
    
    ++meshtagCounter;

}

template<int ndim>
void MeshLoader<ndim>::selectReadForUnstruct(std::string mesh_name,Glamour::Entity entity)
{
    std::string filename = mesh_name;
	UnstructBlock* commBlocks = new UnstructBlock[1];
	
    std::ifstream finGrid;
	finGrid.open(filename, std::ios::in);
	if (!finGrid.is_open())
	{
		throw std::runtime_error("Open grid file failure\n");
	}

	int temp_nnodes, temp_nedges, temp_ncells;
	int to_throw_away;
	finGrid >> to_throw_away >> to_throw_away >> to_throw_away;
	finGrid >> temp_nnodes >> temp_nedges >> temp_ncells >> to_throw_away >> to_throw_away;

	double *temp_xyz = new double[temp_nnodes * ndim];
	int temp_count;
	for (temp_count = 0; temp_count < temp_nnodes; temp_count++)
	{
		for (int j = 0; j < ndim; j++)
		{
			finGrid >> temp_xyz[ndim * temp_count + j];
		}
	}
#ifdef DIVIDED_BY_THOUSAND
	for (int i = 0;i < temp_nnodes * ndim;i++)
	{
		temp_xyz[i] = temp_xyz[i]/1000;
	}
#endif
	int *edge2Cell = new int[2 * temp_nedges];
	std::vector<int> edgeNodePtr(temp_nedges + 1, -1);
	std::vector<int> edgeNodeInd;
	edgeNodeInd.reserve(2 * (ndim - 1) * temp_nedges);
	int temp_nodes[4];
	edgeNodePtr[0] = 0;
	for (int i = 0; i < temp_nedges; i++)
	{
		finGrid >> temp_count;
		edgeNodePtr[i + 1] = edgeNodePtr[i] + temp_count;
		for (int j = 0; j < temp_count; j++)
		{
			finGrid >> temp_nodes[j];
			edgeNodeInd.push_back(temp_nodes[j]);
		}
		finGrid >> edge2Cell[2 * i] >> edge2Cell[2 * i + 1];
	}
	finGrid.close();
	// Extend cells distributed on cur_proc by nfringe.
	int *buffered_flag = new int[temp_ncells];
	int *buffered_flagtmp = new int[temp_ncells];
	for (int i = 0; i < temp_ncells; i++)
	{
		buffered_flag[i] = buffered_flagtmp[i] = 0;
	}

    int num_proc,cur_proc;
    MPI_Comm_size(MPI_COMM_WORLD,&num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD,&cur_proc);
		
	// Use HDF5 to fill in buffered_flag to be 0 where the METIS original split is.
	int *cellRealDistribution = new int[temp_ncells];
	/*
		Here comes the parallel HDF5
	*/
		
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Info info = MPI_INFO_NULL;
	hid_t plist_id=H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id,comm,info);
	hid_t file_id = H5Fopen(InterFileName.c_str(),H5P_DEFAULT,H5FD_MPIO_INDEPENDENT);
	hid_t status = H5Pclose(plist_id);
	std::string curDatasetname = "cell_distribution";
	//H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
	hid_t dataset_id = H5Dopen(file_id,curDatasetname.c_str(),H5P_DEFAULT);


	hid_t file_dataspace = H5Dget_space(dataset_id);
	int eachCellSize = 1 + 1 + (ndim == 2 ? 4 : 8);
	int *buffer = new int[eachCellSize];
	hsize_t count[2] = {1, eachCellSize};
	hsize_t offset[2] = {0, 0};
	hsize_t stride[2] = {1, 1};
	hsize_t block[2] = {1, 1};
	hid_t memspace_id = H5Screate_simple(2, count, NULL);
	for (int i = 0; i < temp_ncells;i++ )
	{
		offset[0] = i;
		hid_t status = H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, offset, stride, count, block);
		H5Dread(dataset_id, H5T_NATIVE_INT, memspace_id, file_dataspace, H5P_DEFAULT, buffer);
		cellRealDistribution[i] = buffer[1];
		if(buffer[1]==cur_proc)
		{
			buffered_flag[i] = buffered_flagtmp[i] = BIT(0);
		}
	}
		
	int nfringe = 2;
	for (int i = 0; i < nfringe; i++)
	{
		int curLayer = BIT(i);
		int nextLayer = BIT(i+1);
		for (int j = 0; j < temp_nedges; j++)
		{
			if (edge2Cell[2 * j] < 0 or edge2Cell[2 * j + 1] < 0)
			{
				continue;
			}
			bool leftInside = buffered_flag[edge2Cell[2 * j] - cellBase] > 0;
			bool rightInside = buffered_flag[edge2Cell[2 * j + 1] - cellBase] > 0;
			if (leftInside and !rightInside)
			{
				buffered_flagtmp[edge2Cell[2 * j + 1] - cellBase] = nextLayer;
			}
			if (!leftInside and rightInside)
			{
				buffered_flagtmp[edge2Cell[2 * j] - cellBase] = nextLayer;
			}
		}
		for (int j = 0; j < temp_ncells; j++)
		{
			buffered_flag[j] = buffered_flagtmp[j];
		}
	}
	// Find all Edges and Nodes related to the extended localCells;
	std::set<int, std::less<int>> *hereNodes = new std::set<int, std::less<int>>[1];
	std::set<int, std::less<int>> *hereCells = new std::set<int, std::less<int>>[1];

	//Wall nodes needed by SST Model
	std::set<int, std::less<int>> *globalWbc = new std::set<int, std::less<int>>[1];
	for (int i = 0; i < temp_ncells; i++)
	{
			
		offset[0] = i;
		hid_t status = H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, offset, stride, count, block);
		H5Dread(dataset_id, H5T_NATIVE_INT, memspace_id, file_dataspace, H5P_DEFAULT, buffer);
		cellRealDistribution[i] = buffer[1];
		if(buffered_flag[i]>=BIT(0))
		{
			for (int j = 0; j < buffer[0]; j++)
			{
				hereNodes->insert(buffer[2+j]);
			}
			hereCells->insert(i + cellBase);	
		}
	}
		

	LocalGlobalMap hereNodeMapping(hereNodes->size(), hereNodes->begin(), hereNodes->end());
	LocalGlobalMap hereCellMapping(hereCells->size(), hereCells->begin(), hereCells->end());
	// Allocate Point3d.
	GeomElements::point3d<2> *pt2d;
	GeomElements::point3d<3> *pt3d;
	if (ndim == 2)
	{
		pt2d = new GeomElements::point3d<2>[hereNodes->size()];
		int tmp_count = 0;
		for (auto iter = hereNodes->begin(); iter != hereNodes->end(); iter++)
		{
			pt2d[tmp_count++] = GeomElements::vector3d<2, double>(&(temp_xyz[2 * ((*iter) - nodeBase)]));
		}
	}
	else if (ndim == 3)
	{
		pt3d = new GeomElements::point3d<3>[hereNodes->size()];
		int tmp_count = 0;
		for (auto iter = hereNodes->begin(); iter != hereNodes->end(); iter++)
		{
			pt3d[tmp_count++] = GeomElements::vector3d<3, double>(&(temp_xyz[3 * ((*iter) - nodeBase)]));
		}
	}
	// Allocate Cell3d
	GeomElements::cell3d<2> *cel2d;
	GeomElements::cell3d<3> *cel3d;
	// Sort out the distribution of communication cells on each Proc
	std::vector<int> nEachProcBufferCell(num_proc, 0);
	int nlocalCells = 0;
	// How many communication cells there is on each proc
	for (int i = 0; i < temp_ncells; i++)
	{
		if (buffered_flag[i] == BIT(0))
		{
			++nlocalCells;
		}
		else if (buffered_flag[i] > BIT(0))
		{
			++nEachProcBufferCell[cellRealDistribution[i]];
		}
	}
	int nbufferCells = std::accumulate(nEachProcBufferCell.begin(), nEachProcBufferCell.end(), 0);
	// How many related proc there is on this proc.

	temp_count = 0;

	std::map<int, int> *index2Proc = new std::map<int, int>[1];
	std::map<int, int> *proc2Index = new std::map<int, int>[1];

	for (int i = 0; i < num_proc; i++)
	{
		if (nEachProcBufferCell[i] != 0)
		{
			(*index2Proc)[temp_count] = i;
			(*proc2Index)[i] = temp_count;
			temp_count++;
		}
	}
	int nRelatedProcs = index2Proc->size();

	// Map is formed here, we should give it to BLock later.
		
	nlocalCells += nbufferCells;

	if (ndim == 2)
	{
		cel2d = new GeomElements::cell3d<2>[nlocalCells];
		int tmp_count = 0;
		for(int i = 0;i<temp_ncells;i++)
		{
			if(buffered_flag[i]>=BIT(0))
			{
				offset[0]=i;
				hid_t status = H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, offset, stride, count, block);
				H5Dread(dataset_id, H5T_NATIVE_INT, memspace_id, file_dataspace, H5P_DEFAULT, buffer);
				for (int j = 0; j < buffer[0]; j++)
				{
					cel2d[tmp_count].push_back(hereNodeMapping.global2local(buffer[2+j]) - 1);
				}
				if (buffered_flag[i] > BIT(0)) // Is a buffered cell,we could compute its local id and globalid(stored in remote id for now) and proc for sending
				{
					int localid = hereCellMapping.global2local(i + 1) - 1; // Starts from 0;
					commBlocks->pushRecvCell(CommunicationCell(localid, proc2Index->find(cellRealDistribution[i])->second, i + 1,buffered_flag[i]));
					if(buffered_flag[i]>BIT(2))
					{
						throw std::runtime_error("bufferFlag should not be greater than BIT(2)");
					}
				}
				tmp_count++;
			}
		}
	}
	else if(ndim==3)
	{
		cel3d = new GeomElements::cell3d<3>[nlocalCells];
		int tmp_count = 0;
		for(int i = 0;i<temp_ncells;i++)
		{
			if(buffered_flag[i]>=BIT(0))
			{
				offset[0]=i;
				hid_t status = H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, offset, stride, count, block);
				H5Dread(dataset_id, H5T_NATIVE_INT, memspace_id, file_dataspace, H5P_DEFAULT, buffer);
				for (int j = 0; j < buffer[0]; j++)
				{
					cel3d[tmp_count].push_back(hereNodeMapping.global2local(buffer[2+j]) - 1);
				}
				if (buffered_flag[i] > BIT(0)) // Is a buffered cell,we could compute its local id and globalid(stored in remote id for now) and proc for sending
				{
					int localid = hereCellMapping.global2local(i + 1) - 1; // Starts from 0;
					commBlocks->pushRecvCell(CommunicationCell(localid, proc2Index->find(cellRealDistribution[i])->second, i + 1,buffered_flag[i]));
				}
				tmp_count++;
			}
		}
	}
	status = H5Sclose(memspace_id);
	status = H5Sclose(file_dataspace);
	status = H5Dclose(dataset_id);
	status = H5Fclose(file_id);
	delete[] buffer;

	// Now for Cell communication,we need
	// 1.Proc associated with current proc and its mapping  (done)
	// 2.How many CommunicationCell there is on each proc
	// 3.Send the receive request and make those send request registered on remote side.

	// 1.make relevant proc mapping

	// 3.Send the receive request by 2 steps.
	// First let each other know how many cells there is.
	std::vector<int> nEachProcSendCells(nRelatedProcs, 0);
	std::vector<int> nEachProcRecvCells(nRelatedProcs, 0);
	for (auto iter = index2Proc->begin(); iter != index2Proc->end(); iter++)
	{
		nEachProcRecvCells[iter->first] = nEachProcBufferCell[iter->second];
	}

	std::vector<MPI_Request> EachPrcoRequestsSend(nRelatedProcs);
	std::vector<MPI_Request> EachPrcoRequestsRecv(nRelatedProcs);
	std::vector<MPI_Status> EachProcStatusRecv(nRelatedProcs);
	std::vector<MPI_Status> EachProcStatusSend(nRelatedProcs);

	for (int i = 0; i < nRelatedProcs; i++)
	{
		MPI_Irecv(&nEachProcSendCells[i], 1, MPI_INT, index2Proc->find(i)->second, 0, MPI_COMM_WORLD, &(EachPrcoRequestsRecv[i]));
		MPI_Isend(&nEachProcRecvCells[i], 1, MPI_INT, index2Proc->find(i)->second, 0, MPI_COMM_WORLD, &(EachPrcoRequestsSend[i]));
	}

	MPI_Waitall(nRelatedProcs, EachPrcoRequestsRecv.data(), EachProcStatusRecv.data());
	MPI_Waitall(nRelatedProcs, EachPrcoRequestsSend.data(), EachProcStatusSend.data());

	// Then Allocate space for sending and recving the communicationCells
	// First Alocate satially contiguous memory for send Recv at EachProc;
	int nlocal_recvCell_buffer = std::accumulate(nEachProcRecvCells.begin(), nEachProcRecvCells.end(), 0);
	int nlocal_sendCell_buffer = std::accumulate(nEachProcSendCells.begin(), nEachProcSendCells.end(), 0);
	int *local_recvCell_buffer = new int[3 * nlocal_recvCell_buffer];
	int *local_sendCell_buffer = new int[3 * nlocal_sendCell_buffer];
	// Calculate the offset for each Proc
	std::vector<int> ProcOffset_s(nRelatedProcs + 1, 0);
	std::vector<int> ProcOffset_r(nRelatedProcs + 1, 0);

	for (int i = 0; i < nRelatedProcs; i++)
	{
		ProcOffset_s[i + 1] = ProcOffset_s[i] + (nEachProcSendCells[i]);
		ProcOffset_r[i + 1] = ProcOffset_r[i] + (nEachProcRecvCells[i]);
	}


	temp_count = 0;
	std::vector<CommunicationCell> &cur_recvCell = commBlocks->getRecvCells();
	std::sort(cur_recvCell.begin(), cur_recvCell.end(), rankWithLocal); // Rank with the local id of the receive cells
	for (auto iter = cur_recvCell.begin(); iter != cur_recvCell.end(); iter++)
	{
		local_recvCell_buffer[temp_count++] = iter->d_localId;
		local_recvCell_buffer[temp_count++] = iter->d_remoteId;
		local_recvCell_buffer[temp_count++] = iter->d_layer_level;
	}

	// Sending receive cells to the sender side
	for (int i = 0; i < nRelatedProcs; i++)
	{
		MPI_Irecv(&local_sendCell_buffer[3 * ProcOffset_s[i]], 3 * nEachProcSendCells[i], MPI_INT, index2Proc->find(i)->second, 0, MPI_COMM_WORLD, &(EachPrcoRequestsRecv[i]));
		MPI_Isend(&local_recvCell_buffer[3 * ProcOffset_r[i]], 3 * nEachProcRecvCells[i], MPI_INT, index2Proc->find(i)->second, 0, MPI_COMM_WORLD, &(EachPrcoRequestsSend[i]));
	}

	MPI_Waitall(nRelatedProcs, EachPrcoRequestsRecv.data(), EachProcStatusRecv.data());
	MPI_Waitall(nRelatedProcs, EachPrcoRequestsSend.data(), EachProcStatusSend.data());
		
	temp_count = 0;
	for (int i = 0; i < nRelatedProcs; i++)
	{
		int proc_id = i;
		for (int j = ProcOffset_s[i]; j < ProcOffset_s[i + 1]; j++)
		{
			int remote_id = local_sendCell_buffer[3*temp_count];
			int global_id = local_sendCell_buffer[3*temp_count+1];
			int local_id = hereCellMapping.global2local(global_id) - 1;
			int layer_level = local_sendCell_buffer[3*temp_count+2];
			commBlocks->pushSendCell(CommunicationCell(local_id, proc_id, remote_id,layer_level));
			local_sendCell_buffer[temp_count] = local_id;
			++temp_count;
		}
	}
	for (int i = 0; i < nRelatedProcs; i++)
	{
		MPI_Irecv(&local_recvCell_buffer[ProcOffset_r[i]], nEachProcRecvCells[i], MPI_INT, index2Proc->find(i)->second, 1, MPI_COMM_WORLD, &(EachPrcoRequestsRecv[i]));
		MPI_Isend(&local_sendCell_buffer[ProcOffset_s[i]], nEachProcSendCells[i], MPI_INT, index2Proc->find(i)->second, 1, MPI_COMM_WORLD, &(EachPrcoRequestsSend[i]));
	}
	MPI_Waitall(nRelatedProcs, EachPrcoRequestsRecv.data(), EachProcStatusRecv.data());
	MPI_Waitall(nRelatedProcs, EachPrcoRequestsSend.data(), EachProcStatusSend.data());


	temp_count = 0;
	{
		std::vector<CommunicationCell> &cur_recvCell = commBlocks->getRecvCells();
		for (auto iter = cur_recvCell.begin(); iter != cur_recvCell.end(); iter++)
		{
			iter->d_remoteId = local_recvCell_buffer[temp_count++];
		}
		std::vector<CommunicationCell> &cur_sendCell = commBlocks->getSendCells();
		std::sort(cur_sendCell.begin(), cur_sendCell.end(), rankWithRemote); // Rank with the local id of the recv cells on the other side.
	}
	// Allocate Edge3d
	GeomElements::edge3d<2> *edg2d;
	GeomElements::edge3d<3> *edg3d;
	std::vector<int> local_edge_inds;
	local_edge_inds.reserve(temp_nedges / num_proc);
	
    std::vector<int> unFluxedEdgeInd;
	std::vector<int> fluxedEdgeInd;

	std::vector<int> wallBoundEdgeInd;
	std::vector<int> overBoundEdgeInd;
	std::vector<int> symmBoundEdgeInd;

	for (int i = 0; i < temp_nedges; i++)
	{
		// natural Boundary Type just find a minimal between left and right.
		// So if the minimum is a negative, then it is a natural boundary.
		int naturalBoundaryTypes = std::min(edge2Cell[2 * i], edge2Cell[2 * i + 1]);
		int naturalBoundaryOtherSide = std::max(edge2Cell[2 * i], edge2Cell[2 * i + 1]);
		bool naturalBoundaryExists = naturalBoundaryTypes < 0;
		bool withInLocalBufferedArea = false;
		// If a natural boundary exists, then check its bufferd_flag at the other side and if  
		if (naturalBoundaryExists) // Then check if the other side is in buffered area.
		{
			
			// If buffered_flag tag is 0, then it is inside the split of Metis
			// If is greater or equal to 1, then it is inside the buffer region.
			withInLocalBufferedArea = buffered_flag[naturalBoundaryOtherSide - cellBase] >= BIT(0);
			
			if (withInLocalBufferedArea) // The other side is in buffered area, then
			{
				// std::cout << "WithInLocalBufferedArea";
				if (naturalBoundaryTypes == signWall)
				{
					// std::cout << "SignWall found\n";
					wallBoundEdgeInd.push_back(local_edge_inds.size());
				}
				else if (naturalBoundaryTypes == signOver)
				{
					// std::cout << "SignOver found\n";
					overBoundEdgeInd.push_back(local_edge_inds.size());
				}
				else if(naturalBoundaryTypes == signSymm)
				{
					symmBoundEdgeInd.push_back(local_edge_inds.size());
				}
				else
				{
					throw std::runtime_error("Undefined Boundary Type\n");
				}
				local_edge_inds.push_back(i);
			}
			else // The other side is not in buffered area, leave it alone
			{
			}
		}
		else // Natural Boundary doesn't exists, then all we need to do is check if bothside is in buffered area.
		{
			withInLocalBufferedArea = buffered_flag[naturalBoundaryTypes - cellBase] >= BIT(0) and buffered_flag[naturalBoundaryOtherSide - cellBase] >= BIT(0);
			if (withInLocalBufferedArea) // Still need to check if it is a fluxed edge or unfluxed
			{
				bool isFluxedEdge = buffered_flag[naturalBoundaryTypes - cellBase] == BIT(0) or buffered_flag[naturalBoundaryOtherSide - cellBase] == BIT(0);
				if (isFluxedEdge) // std::cout<<"withInBUfferedLocal found\n";
				{
					local_edge_inds.push_back(i);
				}
				else
				{
					unFluxedEdgeInd.push_back(i);
				}
			}
			else
			{
			}
		}
	}
	int nFluxedEdge = local_edge_inds.size();
	for (int i = 0; i < unFluxedEdgeInd.size(); i++)
	{
		local_edge_inds.push_back(unFluxedEdgeInd[i]);
	}
		

	if (ndim == 2)
	{
		
		int tmp_count = 0;
		edg2d = new GeomElements::edge3d<2>[local_edge_inds.size()];
		for (auto iter = local_edge_inds.begin(); iter != local_edge_inds.end(); iter++)
		{
			// std::cout<<temp_count<<'\n';
			for (int j = edgeNodePtr[*iter]; j < edgeNodePtr[*iter + 1]; j++)
			{
				edg2d[tmp_count].push_back(hereNodeMapping.global2local(edgeNodeInd[j]) - 1);
			}
			// std::cout<<temp_count<<'\n';
			int lC_ind = edge2Cell[2 * (*iter)];
			int rC_ind = edge2Cell[2 * (*iter) + 1]; // hereCellMapping.global2local(edge2Cell[2 * (*iter)]); // LeftSide of
			if (lC_ind > 0)
			{
				int lC = hereCellMapping.global2local(lC_ind);
				if (lC > 0)
				{
					edg2d[tmp_count].setLeft(lC - cellBase);
					cel2d[lC - cellBase].push_edge(tmp_count);
					
				}
				else
				{
					throw std::runtime_error("Can't be the case where left is NULL\n");
				}
			}
			if (rC_ind > 0)
			{
				int rC = hereCellMapping.global2local(rC_ind);
				if (rC > 0)
				{
					edg2d[tmp_count].setRight(rC - cellBase);
					cel2d[rC-cellBase].push_edge(tmp_count);
				}
			}
			tmp_count++;	
		}
		for (auto iter = wallBoundEdgeInd.begin(); iter != wallBoundEdgeInd.end(); iter++)
		{
			if (edg2d[*iter].rCInd() > 0)
			{
				throw std::runtime_error("Right side of the wallBoundary is already set\n");
			}
			edg2d[*iter].setRight(GeomElements::edge::BoundaryType::WALL);
		}
		for (auto iter = overBoundEdgeInd.begin(); iter != overBoundEdgeInd.end(); iter++)
		{
			if (edg2d[*iter].rCInd() > 0)
			{
				throw std::runtime_error("Right side of the overBoundary is already set\n");
			}
			// std::cout<<"Over found\n";
			edg2d[*iter].setRight(GeomElements::edge::BoundaryType::FARFIELD);
		}
		for (auto iter = symmBoundEdgeInd.begin(); iter != symmBoundEdgeInd.end(); iter++)
		{
			if(edg2d[*iter].rCInd() > 0)
			{
				throw std::runtime_error("Right side of the symmBoundary is already set\n");
			}
			edg2d[*iter].setRight(GeomElements::edge::BoundaryType::SYMMETRY);
		}

		//Here we output the edge 2 cell map of 
	}
	else if (ndim == 3)
	{
		int tmp_count = 0;
		edg3d = new GeomElements::edge3d<3>[local_edge_inds.size()];
		for (auto iter = local_edge_inds.begin(); iter != local_edge_inds.end(); iter++)
		{
			for (int j = edgeNodePtr[*iter]; j < edgeNodePtr[*iter + 1]; j++)
			{
				edg3d[tmp_count].push_back(hereNodeMapping.global2local(edgeNodeInd[j]) - 1);
			}
			int lC_ind = edge2Cell[2 * (*iter)];
			int rC_ind = edge2Cell[2 * (*iter) + 1];
			if (lC_ind > 0)
			{
				int lC = hereCellMapping.global2local(lC_ind);
				if (lC > 0)
				{
					edg3d[tmp_count].setLeft(lC - cellBase);
					cel3d[lC - cellBase].push_edge(tmp_count);
				}
				else
				{
					throw std::runtime_error("Can't be the case where left is NULL\n");
				}
			}
			if (rC_ind > 0)
			{
				int rC = hereCellMapping.global2local(rC_ind);
				if (rC > 0)
				{
					edg3d[tmp_count].setRight(rC - cellBase);
					cel3d[rC-cellBase].push_edge(tmp_count);
				}
			}
			tmp_count++;	
		}
		for (auto iter = wallBoundEdgeInd.begin(); iter != wallBoundEdgeInd.end(); iter++)
		{
			if (edg3d[*iter].rCInd() > 0)
			{
				throw std::runtime_error("Right side of the wallBoundary is already set\n");
			}
			edg3d[*iter].setRight(GeomElements::edge::BoundaryType::WALL);
		}
		for (auto iter = overBoundEdgeInd.begin(); iter != overBoundEdgeInd.end(); iter++)
		{
			if (edg3d[*iter].rCInd() > 0)
			{
				throw std::runtime_error("Right side of the overBoundary is already set\n");
			}
			// std::cout<<"Over found\n";
			edg3d[*iter].setRight(GeomElements::edge::BoundaryType::FARFIELD);
		}
		for (auto iter = symmBoundEdgeInd.begin(); iter != symmBoundEdgeInd.end(); iter++)
		{
			if(edg3d[*iter].rCInd() > 0)
			{
				throw std::runtime_error("Right side of the symmBoundary is already set\n");
			}
			edg3d[*iter].setRight(GeomElements::edge::BoundaryType::SYMMETRY);
		}
	}
	if(ndim==2)
	{
		Glamour::MeshComponent<2>& curMeshCompo = entity.AddComponent<Glamour::MeshComponent<2>>();
		curMeshCompo.setMeshComponent(hereNodes->size(), hereCells->size(), local_edge_inds.size(),pt2d,cel2d,edg2d,nFluxedEdge);
		commBlocks->setLocalCommunication(proc2Index, index2Proc);
		curMeshCompo.setCommComponent(commBlocks);
		curMeshCompo.mesh_ind = meshtagCounter-1;
	}
	else if(ndim==3)
	{
		Glamour::MeshComponent<3>& curMeshCompo = entity.AddComponent<Glamour::MeshComponent<3>>();
		curMeshCompo.setMeshComponent(hereNodes->size(), hereCells->size(), local_edge_inds.size(),pt3d,cel3d,edg3d,nFluxedEdge);
		commBlocks->setLocalCommunication(proc2Index, index2Proc);
		curMeshCompo.setCommComponent(commBlocks);
		curMeshCompo.mesh_ind = meshtagCounter-1;
	}
}

template<int ndim>
void MeshLoader<ndim>::selectReadForUniTioga(std::string mesh_name, Glamour::Entity entity)
{

	std::string filename = mesh_name;
	std::ifstream finGrid;
	finGrid.open(filename, std::ios::in);
	if (!finGrid.is_open())
	{
		throw std::runtime_error("Open grid file failure\n");
	}
	int temp_nnodes, temp_nedges, temp_ncells;
	int to_throw_away;
	finGrid >> to_throw_away >> to_throw_away >> to_throw_away;
	finGrid >> temp_nnodes >> temp_nedges >> temp_ncells >> to_throw_away >> to_throw_away;
	int temp_count = 0;
	double *temp_xyz = new double[temp_nnodes * ndim];
	for (temp_count = 0; temp_count < temp_nnodes; temp_count++)
	{
		for (int j = 0; j < ndim; j++)
		{
			finGrid >> temp_xyz[ndim * temp_count + j];
		}
	}
#ifdef DIVIDED_BY_THOUSAND
	for (int i = 0;i < temp_nnodes * ndim;i++)
	{
		temp_xyz[i] = temp_xyz[i]/1000;
	}
#endif
	int *edge2Cell = new int[2*temp_nedges];
	int temp_nodes[4];
	for (int i = 0; i < temp_nedges; i++)
	{
		finGrid >> temp_count;
		for (int j = 0; j < temp_count; j++)
		{
			finGrid >> temp_nodes[j];
		}
		finGrid >> edge2Cell[2 * i] >> edge2Cell[2 * i + 1];
	}
	finGrid.close();
	int *buffered_flag = new int[temp_ncells];
	int *buffered_flagtmp = new int[temp_ncells];
	for(int i = 0; i < temp_ncells; i++)
	{
		buffered_flag[i] = buffered_flagtmp[i] = -1;
	}
	
    int num_proc,cur_proc;
    MPI_Comm_size(MPI_COMM_WORLD,&num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD,&cur_proc);
    /*
		Here comes the parallel HDF5 to read all the local Nodes and Cells
	*/	
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Info info = MPI_INFO_NULL;
	hid_t plist_id=H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id,comm,info);
	hid_t file_id = H5Fopen(InterFileName.c_str(),H5P_DEFAULT,H5FD_MPIO_INDEPENDENT);
	hid_t status = H5Pclose(plist_id);
	std::string curDatasetname = "cell_distribution";
	hid_t dataset_id = H5Dopen(file_id,curDatasetname.c_str(),H5P_DEFAULT);
	hid_t file_dataspace = H5Dget_space(dataset_id);
	int eachCellSize = 1 + 1 + (ndim == 2 ? 4 : 8);
	int *buffer = new int[eachCellSize];
	hsize_t count[2] = {1, eachCellSize};
	hsize_t offset[2] = {0, 0};
	hsize_t stride[2] = {1, 1};
	hsize_t block[2] = {1, 1};
	hid_t memspace_id = H5Screate_simple(2, count, NULL);

	for(int i = 0; i < temp_ncells; i++)
	{
		offset[0] = i;
		hid_t status = H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, offset, stride, count, block);
		H5Dread(dataset_id, H5T_NATIVE_INT, memspace_id, file_dataspace, H5P_DEFAULT, buffer);
		if(buffer[1]==cur_proc)
		{
			buffered_flag[i] = buffered_flagtmp[i] = 0;
		}
	}
	//Extend the already tagged buffered_flag by nfringe+1 layers.
	int nfringe = 2;
	for (int i = 0; i < nfringe + 1; i++)
	{
		for (int j = 0; j < temp_nedges; j++)
		{
			if (edge2Cell[2 * j] < 0 or edge2Cell[2 * j + 1] < 0)
			{
				continue;
			}
			bool leftInside = buffered_flag[edge2Cell[2 * j] - cellBase] >= 0;
			bool rightInside = buffered_flag[edge2Cell[2 * j + 1] - cellBase] >= 0;
			if (leftInside and !rightInside)
			{
				buffered_flagtmp[edge2Cell[2 * j + 1] - cellBase] = i+1;
			}
			if (!leftInside and rightInside)
			{
				buffered_flagtmp[edge2Cell[2 * j] - cellBase] = i+1;
			}
		}
		for (int j = 0; j < temp_ncells; j++)
		{
			buffered_flag[j] = buffered_flagtmp[j];
		}
	}
	std::set<int, std::less<int>> *hereCells = new std::set<int, std::less<int>>[1];
	std::set<int, std::less<int>> *unwantedNodes = new std::set<int, std::less<int>>[1];
	std::set<int, std::less<int>> *wantedNodes = new std::set<int, std::less<int>>[1];
	std::set<int, std::less<int>> *commNodes = new std::set<int, std::less<int>>[1];
	
	std::set<int, std::less<int>> *overNodes = new std::set<int, std::less<int>>[1];
	std::set<int, std::less<int>> *wallNodes = new std::set<int, std::less<int>>[1];
	std::set<int, std::less<int>> *hereNodes = new std::set<int, std::less<int>>[1];
	std::set<int, std::less<int>> *hereWbcNodes = new std::set<int, std::less<int>>[1];
	std::set<int, std::less<int>> *hereObcNodes = new std::set<int, std::less<int>>[1];
	for (int i = 0; i < temp_ncells; i++)
	{
		if(buffered_flag[i]>=0)
		{
		
			offset[0] = i;
			hid_t status = H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, offset, stride, count, block);
			H5Dread(dataset_id, H5T_NATIVE_INT, memspace_id, file_dataspace, H5P_DEFAULT, buffer);
		
			if(buffered_flag[i]==nfringe+1)
			{
				for(int j = 0; j < buffer[0]; j++)
				{
					unwantedNodes->insert(buffer[2+j]);
				}
			}
			else
			{
				for (int j = 0; j < buffer[0]; j++)
				{
					wantedNodes->insert(buffer[2+j]);
				}
				hereCells->insert(i + cellBase);
			}	
		}
	}
	std::insert_iterator<std::set<int, std::less<int>>> hereNodesInserter(*hereNodes, hereNodes->begin());
	std::insert_iterator<std::set<int, std::less<int>>> commNodesInserter(*commNodes, commNodes->begin());
	
	std::set_intersection(wantedNodes->begin(),wantedNodes->end(), unwantedNodes->begin(), unwantedNodes->end(), commNodesInserter);
	std::set_difference(wantedNodes->begin(), wantedNodes->end(), unwantedNodes->begin(), unwantedNodes->end(), hereNodesInserter);
	std::set_union(hereNodes->begin(),hereNodes->end(),commNodes->begin(),commNodes->end(),hereNodesInserter);
	
	delete[] wantedNodes;
	delete[] unwantedNodes;
	Brick4Cobalt2D *hereCells2D = NULL;
	Brick4Cobalt3D *hereCells3D = NULL;
	if(ndim==2)
	{
		hereCells2D = new Brick4Cobalt2D[hereCells->size()];
	}		
	else if(ndim==3)
	{
		hereCells3D = new Brick4Cobalt3D[hereCells->size()];
	}
	/*
		Here comes the parallel HDF5 to find vconn
	*/
	
	temp_count = 0;
	int temp_count_scratch = 0;
	if (ndim == 2)
	{
		for (int i = 0; i < temp_ncells; i++)
		{
			offset[0] = i;
			hid_t status = H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, offset, stride, count, block);
			
			H5Dread(dataset_id, H5T_NATIVE_INT, memspace_id, file_dataspace, H5P_DEFAULT, buffer);
			if (buffered_flag[i] >= 0 and buffered_flag[i] <= nfringe)
			{
				hereCells2D[temp_count].setSize(buffer[0]);
				for (int j = 0; j < buffer[0]; j++)
				{
					hereCells2D[temp_count].setNode(buffer[2 + j], j);
				}
				++temp_count;
			}
			
		}
	}
	else if (ndim == 3)
	{
		for (int i = 0; i < temp_ncells; i++)
		{
			offset[0] = i;
			hid_t status = H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, offset, stride, count, block);
			H5Dread(dataset_id, H5T_NATIVE_INT, memspace_id, file_dataspace,H5P_DEFAULT, buffer);
			if (buffered_flag[i] >= 0 and buffered_flag[i] <= nfringe)
			{
				hereCells3D[temp_count].setSize(buffer[0]);
				for (int j = 0; j < buffer[0]; j++)
				{
					hereCells3D[temp_count].setNode(buffer[2 + j], j);
				}
				++temp_count;
			}
			
		}
	}
	
	// continue;
	status = H5Sclose(memspace_id);
	status = H5Sclose(file_dataspace);
	status = H5Dclose(dataset_id);
	status = H5Fclose(file_id);
	delete[] buffer;
	// Making the mapping from global to local
	LocalGlobalMap hereMapping(hereNodes->size(), hereNodes->begin(), hereNodes->end());
	int nc[10] = {0};
	int size2Type[10] = {-1};
	int ntypes = 0;
	if (ndim == 2)
	{
		for (int i = 0; i < hereCells->size(); i++)
		{
			nc[hereCells2D[i].sizeIs() - 1]++;
		}
	}
	else if (ndim == 3)
	{
		for (int i = 0; i < hereCells->size(); i++)
		{
			nc[hereCells3D[i].sizeIs() - 1]++;
		}
	}
	for (int i = 0; i < 10; i++)
	{
		if (nc[i] > 0)
		{
			// std::cout<<i <<" "<<nc[i]<<'\n';
			ntypes++;
		}
	}
	cur_bi.ntypes = ntypes;
	cur_bi.vconn = new int *[ntypes];
	cur_bi.ncells = hereCells->size();
	cur_bi.eachNodeCount = new int[ntypes];
	cur_bi.eachCellCount = new int[ntypes];
	temp_count = 0;
	for (int i = 0; i < 10; i++)
	{
		if (nc[i] > 0)
		{
			size2Type[i] = temp_count;
			cur_bi.eachNodeCount[temp_count] = i + 1;
			cur_bi.eachCellCount[temp_count] = nc[i];
			temp_count++;
		}
	}
	// std::cout << "ntypes are " << ntypes << '\n';
	for (int i = 0; i < ntypes; i++)
	{
		// std::cout<<cur_bi.eachCellCount[i]<<" "<<cur_bi.eachNodeCount[i]<<'\n';
		// std::cout<<cur_bi.vconn[i]<<'\n';
		int size = cur_bi.eachCellCount[i] * cur_bi.eachNodeCount[i];
		// std::cout << "size is " << size << '\n';
		cur_bi.vconn[i] = new int[size];
	}
	std::vector<int> eachTypeCount(ntypes, 0);
	if (ndim == 2)
	{
		for (int i = 0; i < hereCells->size(); i++)
		{
			// std::cout<<hereCells2D[i].sizeIs();
			int curtype = size2Type[hereCells2D[i].sizeIs() - 1];
			// std::cout<<curtype<<" "<<cur_bi.eachNodeCount[0]<<'\n';
			for (int j = 0; j < cur_bi.eachNodeCount[curtype]; j++)
			{
				cur_bi.vconn[curtype][eachTypeCount[curtype] * cur_bi.eachNodeCount[curtype] + j] = hereMapping.global2local(hereCells2D[i].nodeIs(j));
				if(hereMapping.global2local(hereCells2D[i].nodeIs(j))<0)
				{
					std::cout<<hereCells2D[i].nodeIs(j)<<" is not found in local"<<'\n';
					cur_bi.vconn[curtype][eachTypeCount[curtype] * cur_bi.eachNodeCount[curtype] + j] = 1;
				}
			}
			eachTypeCount[curtype] += 1;
		}
	}
	else if (ndim == 3)
	{
		for (int i = 0; i < hereCells->size(); i++)
		{
			int curtype = size2Type[hereCells3D[i].sizeIs() - 1];
			for (int j = 0; j < cur_bi.eachNodeCount[curtype]; j++)
			{
				cur_bi.vconn[curtype][eachTypeCount[curtype] * cur_bi.eachNodeCount[curtype] + j] = hereMapping.global2local(hereCells3D[i].nodeIs(j));
			}
			eachTypeCount[curtype] += 1;
		}
	}
	if (hereCells2D)
		delete[] hereCells2D;
	if (hereCells3D)
		delete[] hereCells3D;
	// Now we have got vconn and nc and nv and ntypes all set, it's time to get node informations.
	
	int edgeNodes[10];
	int edgeCells[2];
	int temp_int = 0; // Used to record ptr
	
	finGrid.open(filename, std::ios::in);
	if (!finGrid.is_open())
	{
		throw std::runtime_error("Open grid file failure\n");
	}
	finGrid >> to_throw_away >> to_throw_away >> to_throw_away;
	finGrid >> to_throw_away >> to_throw_away >> to_throw_away >> to_throw_away >> to_throw_away;
	double to_throw_awayf64;
	for (temp_count = 0; temp_count < temp_nnodes; temp_count++)
	{
		for (int j = 0; j < ndim; j++)
		{
			finGrid >> to_throw_awayf64;
		}
	}
#ifdef DIVIDED_BY_THOUSAND
		
#endif
	for (int i = 0; i < temp_nedges; i++)
	{
		finGrid >> temp_count;
		for (int j = 0; j < temp_count; j++)
		{
			finGrid >> edgeNodes[j];
		}
		finGrid >> edgeCells[0] >> edgeCells[1];
		if (edgeCells[1] == signWall)
		{
			for (int j = 0; j < temp_count; j++)
			{
				wallNodes->insert(edgeNodes[j]);
			}
		}
		else if (edgeCells[1] == signOver)
		{
			for (int j = 0; j < temp_count; j++)
			{
				overNodes->insert(edgeNodes[j]);
			}
		}
		else if (edgeCells[1] == signSymm)
		{
		}
		else if (edgeCells[1] > 0)
		{
		}
		else
		{
			std::cout<< "edgeCell[1] is "<<edgeCells[1]<<'\n';
			std::cout << "Pre is " << edgeNodes[0] << " " << edgeNodes[1] << '\n';
			std::cout << "Something wrong,stop\n";
		}
	}
	finGrid.close();
	int local_nnodes = hereNodes->size();
	std::cout<<"Local nnode is "<<local_nnodes<<'\n';
	double *local_xyz = new double[local_nnodes * ndim];
	std::set<int, std::less<int>>::iterator itera = hereNodes->begin();
	// Now find all the wbc obc and local nodes coordinates
	temp_count = 0;
	for (auto iter = hereNodes->begin(); iter != hereNodes->end(); iter++)
	{
		for (int j = 0; j < ndim; j++)
		{
			local_xyz[temp_count++] = temp_xyz[ndim * ((*iter) - 1) + j];
		}
	}

	cur_bi.rxyz = local_xyz;
	cur_bi.nnodes = local_nnodes;
	std::insert_iterator<std::set<int, std::less<int>>> hereOBCInserter(*hereObcNodes, hereObcNodes->begin());
	std::insert_iterator<std::set<int, std::less<int>>> hereWBCInserter(*hereWbcNodes, hereWbcNodes->begin());
	std::set_intersection(hereNodes->begin(), hereNodes->end(), overNodes->begin(), overNodes->end(), hereOBCInserter);
	std::set_intersection(hereNodes->begin(), hereNodes->end(), wallNodes->begin(), wallNodes->end(), hereWBCInserter);
	delete[] hereNodes;
	delete[] overNodes;
	delete[] wallNodes;
	int *obc = NULL;
	int *wbc = NULL;
	if (hereObcNodes->size() > 0)
	{
		obc = new int[hereObcNodes->size()];
	}
	if (hereWbcNodes->size() > 0)
	{
		wbc = new int[hereWbcNodes->size()];
	}
	temp_count = 0;
	for (auto iter = hereObcNodes->begin(); iter != hereObcNodes->end(); iter++)
	{
		obc[temp_count] = hereMapping.global2local(*iter);
		++temp_count;
	}
	temp_count = 0;
	for (auto iter = hereWbcNodes->begin(); iter != hereWbcNodes->end(); iter++)
	{
		wbc[temp_count] = hereMapping.global2local(*iter);
		++temp_count;
	}
	temp_count = 0;
	cur_bi.nwbc = hereWbcNodes->size();
	cur_bi.nobc = hereObcNodes->size();
	delete[] hereWbcNodes;
	delete[] hereObcNodes;
	cur_bi.wbc = wbc;
	cur_bi.obc = obc;
	delete[] temp_xyz;
	
}


template<int ndim>
void MeshLoader<ndim>::readWallBoundary(std::string filename, Glamour::Entity entity)
{
	Glamour::MeshComponent<ndim>& curMeshComponent = entity.GetComponent<Glamour::MeshComponent<ndim>>();
	std::ifstream finGrid;
	// First read cellfile to fill the nodes
	// We need to first read how many cells there are.
	
	finGrid.open(filename, std::ios::in);
	if (!finGrid.is_open())
	{
		throw std::runtime_error("Open grid file failure\n");
	}
	int temp_nnodes, temp_nedges, temp_ncells;
	int to_throw_away;
	finGrid >> to_throw_away >> to_throw_away >> to_throw_away;
	finGrid >> temp_nnodes >> temp_nedges >> temp_ncells >> to_throw_away >> to_throw_away;
	// Now we have got vconn and nc and nv and ntypes all set, it's time to get node informations.
	double *temp_xyz = new double[temp_nnodes * ndim];
	int temp_count = 0;
	for (temp_count = 0; temp_count < temp_nnodes; temp_count++)
	{
		for (int j = 0; j < ndim; j++)
		{
			finGrid >> temp_xyz[ndim * temp_count + j];
		}
	}
#ifdef DIVIDED_BY_THOUSAND
	for (int i = 0;i < temp_nnodes * ndim;i++)
	{
		temp_xyz[i] = temp_xyz[i]/1000;
	}
#endif
	int edgeNodes[10];
	int edgeCells[2];
	int temp_int = 0; // Used to record ptr
	int wall_temp_int = 0;
	std::set<int,std::less<int>> wallNodeInds;
	std::vector<int> wallEdgeInds;
	std::vector<int> edgeNodePtr(temp_nedges + 1, -1);
	std::vector<int> edgeNodeInd;
	edgeNodeInd.reserve(2 * (ndim - 1) * temp_nedges);
	edgeNodePtr[0] = 0;
	for (int i = 0; i < temp_nedges; i++)
	{
		finGrid >> temp_count;
		edgeNodePtr[i + 1] = edgeNodePtr[i] + temp_count;
		
		for (int j = 0; j < temp_count; j++)
		{
			finGrid >> edgeNodes[j];
			edgeNodeInd.push_back(edgeNodes[j]);
		}
		finGrid >> edgeCells[0] >> edgeCells[1];
		if (edgeCells[1] == signWall)
		{
			wallEdgeInds.push_back(i);
			for (int j = 0; j < temp_count; j++)
			{
				wallNodeInds.insert(edgeNodes[j]);
			}
		}
	}
	finGrid.close();
	LocalGlobalMap hereNodeMapping(wallNodeInds.size(), wallNodeInds.begin(), wallNodeInds.end());
	std::cout<<"Global local mapping done\n";
	
	curMeshComponent.d_wallEdges = new GeomElements::edge3d<2>[wallEdgeInds.size()];
	curMeshComponent.d_wallNodes = new GeomElements::point3d<ndim>[wallNodeInds.size()];
	int tmp_count = 0;
	for (auto iter = wallNodeInds.begin(); iter != wallNodeInds.end(); iter++)
	{
		curMeshComponent.d_wallNodes[tmp_count++] = GeomElements::vector3d<2, double>(&(temp_xyz[ndim * ((*iter) - nodeBase)]));
	}
	tmp_count = 0;
	for (auto iter = wallEdgeInds.begin(); iter != wallEdgeInds.end(); iter++)
	{
		for (int j = edgeNodePtr[*iter]; j < edgeNodePtr[*iter + 1]; j++)
		{
			curMeshComponent.d_wallEdges[tmp_count].push_back(hereNodeMapping.global2local(edgeNodeInd[j]) - 1);
		}
		tmp_count++;
	}
	curMeshComponent.nwallEdge = wallEdgeInds.size();
	curMeshComponent.nwallNode = wallNodeInds.size();
	delete[] temp_xyz;

}