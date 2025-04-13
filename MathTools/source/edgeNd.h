#pragma once
#include "cellNd.h"
#include <vector>
#include <memory>
#define UNSTRUCT_DEBUG
namespace GeomElements
{
	class edge
	{
	protected:
		edge()
		{
		}

	public:
		enum BoundaryType
		{
			UNSET = -1,
			WALL = -2,
			FARFIELD = -3,
			OVERSET = -4,
			SYMMETRY = -5
		};
		
		virtual void clear() = 0;
		virtual int pointInd(int ind) const = 0;
		virtual int size() const = 0;
		virtual int lCInd() const = 0;
		virtual int rCInd() const = 0;
		virtual void setLeft(int ind) = 0;
		virtual void setRight(int ind) = 0;
		virtual double area() const = 0;
	};

	template <int ndim>
	class edge3d:public edge
	{

		static double Tolerance;
		static point3d<ndim> *d_boundingPoints;
		// static point3d<ndim>* d_nboundingPoints;

		int d_vertsInd[2 * (ndim - 1)]; // 2 in 2D, 4 in 3D
		unsigned int d_size = 0;
		int leftInd = BoundaryType::UNSET;
		int rightInd = BoundaryType::UNSET;
		vector3d<ndim, double> d_normal;
		vector3d<ndim, double> d_center;
		double d_area = 0;
	public:
		
		static void bindPoints(point3d<ndim> *toBind);
		edge3d()
			:edge()
		{
		}
		void push_back(int vertInd)
		{
			d_vertsInd[d_size++] = vertInd;
		}
		void clear()
		{
			d_size = 0;
		}
		int pointInd(int ind) const
		{
			return d_vertsInd[ind];
		}
		int size() const
		{
			return d_size;
		}
		int lCInd() const
		{
			return leftInd;
		}
		int rCInd() const
		{
			return rightInd;
		}
		void setLeft(int ind)
		{
			leftInd = ind;
		}
		void setRight(int ind)
		{
			rightInd = ind;
		}
		vector3d<ndim, double> center() const
		{
			return d_center;
		}
		vector3d<ndim, double> normal_vector() const
		{
			return d_normal;
		}
		double area() const
		{
			return d_area;
		}
		void inverse_normal_vector()
		{
			d_normal = d_normal*(-1);
		}
		void reverseOrient()
		{
			//In 2d case, swap 0 and 1
			//In 3d case, swap 0 and 2 
			int temp = d_vertsInd[ndim-1];
			d_vertsInd[ndim-1] = d_vertsInd[0];
			d_vertsInd[0] = temp;
			inverse_normal_vector();
		}
		void computeMetaData()
		{
			find_normal();
			find_area();
			find_center();
		}

		void newComputeMetaData()
		{
			find_center();
			if(d_size==3)
			{
				vector3d<ndim,double> A = this->point(0);
				vector3d<ndim,double> B = this->point(1);
				vector3d<ndim,double> C = this->point(2);

				double delta_xyA = (A[0] - B[0]) * (A[1] + B[1]);
				double delta_yzA = (A[1] - B[1]) * (A[2] + B[2]);
				double delta_xyB = (B[0] - C[0]) * (B[1] + C[1]);
				double delta_yzB = (B[1] - C[1]) * (B[2] + C[2]);
				double delta_xyC = (C[0] - A[0]) * (C[1] + A[1]);
				double delta_yzC = (C[1] - A[1]) * (C[2] + A[2]);

				double delta_zxA = (A[2] - B[2]) * (A[0] + B[0]);
				double delta_zxB = (B[2] - C[2]) * (B[0] + C[0]);
				double delta_zxC = (C[2] - A[2]) * (C[0] + A[0]);

				vector3d<ndim,double> S = vector3d<ndim,double>(
				delta_yzA + delta_yzB + delta_yzC,
				delta_zxA + delta_zxB + delta_zxC,
				delta_xyA + delta_xyB + delta_xyC) * 0.5;
				d_area = S.normalize();
				d_normal = S;
			}
			else if(d_size==4)
			{
				vector3d<ndim,double> A = this->point(0);
				vector3d<ndim,double> B = this->point(1);
				vector3d<ndim,double> C = this->point(2);
				vector3d<ndim,double> D = this->point(3);

				double delta_xA = D[0] - B[0];
				double delta_xB = C[0] - A[0];
				double delta_yA = D[1] - B[1];
				double delta_yB = C[1] - A[1];
				double delta_zA = D[2] - B[2];
				double delta_zB = C[2] - A[2];
				vector3d<ndim,double> S = vector3d<ndim,double>(
					delta_yA*delta_zB - delta_zA*delta_yB,
					delta_zA*delta_xB - delta_xA*delta_zB,
					delta_xA*delta_yB - delta_yA*delta_xB) * 0.5;
				d_area = S.normalize();
				d_normal = S;
			}

		}
		vector3d<ndim, double> &point(int ind,GeomElements::point3d<ndim>* edgeNodes)
		{
			return edgeNodes[d_vertsInd[ind]].d_pos;
		}

	public:
		vector3d<ndim, double> &point(int ind) const
		{
			return edge3d<ndim>::d_boundingPoints[d_vertsInd[ind]].d_pos;
		}
		
		// void make_clockwise();
		void find_center()
		{
			vector3d<ndim, double> result(0, 0, 0);
			for (int i = 0; i < d_size; i++)
			{
				result = result + point(i);
			}
			d_center = result / d_size;
			
		}
		void find_normal();
		void find_area();
		bool is_plane();
		bool is_clockwise();
	};
	template <int ndim>
	point3d<ndim> *edge3d<ndim>::d_boundingPoints = NULL;
	template <int ndim>
	// This function must be called only after static binding is done
	bool edge3d<ndim>::is_plane()
	{
		if (ndim == 2)
		{
			return true;
		}
		else if (ndim == 3)
		{
			if (d_size == 3)
				return true;
			else
			{
				GeomElements::vector3d<ndim, double> AB_cross_AC = (this->point(1) - this->point(0)).cross_product(this->point(2) - this->point(0));
				GeomElements::vector3d<ndim, double> AB_cross_AD = (this->point(1) - this->point(0)).cross_product(this->point(3) - this->point(0));
				GeomElements::vector3d<ndim, double> result = AB_cross_AC.cross_product(AB_cross_AD);
				for (int i = 0; i < ndim; i++)
				{
					if (result[i] < Tolerance)
					{
						return false;
					}
				}
				return true;
			}
		}
	}
	template <int ndim>
	double edge3d<ndim>::Tolerance = 1.0e-10;
	template <int ndim>
	void edge3d<ndim>::bindPoints(point3d<ndim> *toBind)
	{
		edge3d<ndim>::d_boundingPoints = toBind;
	}
	// This function must be called only after static binding is done
	template <int ndim>
	bool edge3d<ndim>::is_clockwise()
	{
		if (ndim == 2)
		{
			return true;
		}
		else if (ndim == 3)
		{
			if (d_size == 3)
				return true;
			else if (d_size == 4)
			{
				GeomElements::vector3d<ndim, double> AB_cross_AC = (this->point(1) - this->point(0)).cross_product(this->point(2) - this->point(0));
				GeomElements::vector3d<ndim, double> AB_cross_AD = (this->point(1) - this->point(0)).cross_product(this->point(3) - this->point(0));
				double CD_sameSide = AB_cross_AC.dot_product(AB_cross_AD);
				if (CD_sameSide > 0) // C and D are on the same side of AB
				{
					auto BC = this->point(2) - this->point(1);
					GeomElements::vector3d<ndim, double> BC_cross_BA = BC.cross_product(this->point(0) - this->point(1));
					GeomElements::vector3d<ndim, double> BC_cross_BD = BC.cross_product(this->point(3) - this->point(1));
					double AD_sameSide = BC_cross_BA.dot_product(BC_cross_BD);
					if (AD_sameSide > 0)
					{
						return true;
					}
					else
					{
						return false;
					}
				}
				else // AB is a diagonal
				{
					return false;
				}
				return true;
			}
		}
	}
	// This function must be called only after static binding is done
	template <int ndim>
	void edge3d<ndim>::find_normal()
	{
		find_area();
		if (ndim == 2)
		{
			// In 2d case, the normal vector would be the right hand side along the direction of the edge
			vector3d<ndim, double> direction = this->point(1) - this->point(0);
			d_normal[0] = direction[1] / d_area;
			d_normal[1] = -direction[0] / d_area;
		}
		else if (ndim == 3)
		{
			// In 3d case, the normal vector would be AB.cross(AC) after normalization.
			d_normal = (this->point(1) - this->point(0)).cross_product(this->point(2) - this->point(0));
			d_normal.normalize();
		}
	}
	// This function must be called only after static binding is done
	template <int ndim>
	void edge3d<ndim>::find_area()
	{
		if (ndim == 2)
		{
			vector3d<ndim, double> direction = this->point(1) - this->point(0);
			d_area = direction.normalize();
		}
		else if (ndim == 3)
		{
			d_area = 0;
			if (d_size == 3) // Triangle
			{
				d_area += ((this->point(1) - this->point(0)).cross_product(this->point(2) - this->point(0))).normalize() / 2;
			}
			else if (d_size == 4) // Quadrilateral but only in clockwise or counter clockwise
			{
				d_area += ((this->point(1) - this->point(0)).cross_product(this->point(2) - this->point(0))).normalize() / 2;
				d_area += ((this->point(2) - this->point(0)).cross_product(this->point(3) - this->point(0))).normalize() / 2;
			}
		}
	}
}
