/*
	Copyright (C) 2010 Samuele Salti, Federico Tombari, all rights reserved.

	This file is part of SHOT. SHOT has been developed by the 
	Computer Vision Laboratory of the University of Bologna
	(http://www.vision.deis.unibo.it)
	
	SHOT is an implementation of the work described in
	F. Tombari, S. Salti and L. Di Stefano 
	"Unique Signatures of Histograms for Local Surface Description"
	The 11th IEEE European Conference on Computer Vision (ECCV) 2010

	Contacts:
	Samuele Salti mailto:samuele.salti@unibo.it
	Federico Tombari mailto:federico.tombari@unibo.it


    SHOT is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SHOT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SHOT.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "utils.h"
#include "vtkPolyDataNormals.h"
#include "vtkOBJReader.h"
#include "vtkPLYReader.h"
#include "vtkErrorCode.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkTransform.h"
#include "vtkTransformFilter.h"
#include "vtkCleanPolyData.h"

#include <algorithm>
#include <cmath>
#include "randutil.h"

vtkPolyData* LoadPolyData(std::string filename)
{
	std::string ext = filename.substr(filename.length() - 3, filename.length());
	std::transform(ext.begin(), ext.end(), ext.begin(), tolower);

	if (ext == "ply")
		return LoadPly(filename.c_str());
	else if (ext == "off")
		return LoadOff(filename.c_str());
	else if (ext == "obj")
		return LoadOBJ(filename.c_str());
	
	printf("\nERROR LoadPolyData: %s\n", filename.c_str());
	printf("Type not supported.\n");
	return NULL;
}

vtkPolyData* LoadPly(const char* filename)
{
	vtkPLYReader  *plyReader = vtkPLYReader::New();
	plyReader->SetFileName(filename);
	plyReader->Update();
	int errCod = vtkErrorCode::GetLastSystemError();
//	if(errCod != vtkErrorCode::NoError)
//	{
//		const char* errString = vtkErrorCode::GetStringFromErrorCode(errCod);
//		printf("\nERROR LoadPly: %s\n", errString);
//		plyReader->Delete();
//		return NULL;
//	}

	vtkPolyData* polyData = vtkPolyData::New();
	polyData->ShallowCopy(plyReader->GetOutput());
	plyReader->Delete();
	return polyData;
}


vtkPolyData* LoadOff(const char* filename)
{
	FILE *in = fopen(filename, "r");
	if(!in)
	{
		perror("Problems opening OFF file");
		return NULL;
	}

	int n_pts, n_polys, nAux;
	
	
	fscanf(in, "OFF\n");


	char lineStart;
	fpos_t position;

	fgetpos (in, &position);
	fscanf(in, "%c", &lineStart);
	
  
	
	while (lineStart == '#')
	{
		while (lineStart != '\n')
			fscanf(in, "%c", &lineStart);

		fgetpos (in, &position);
		fscanf(in, "%c", &lineStart);
		
	}
	

	
  
  //fsetpos (in, &position);
	//fseek(in, -1, SEEK_CUR);

	char n_pts_str[10];
	char n_pts_start[2] = {lineStart, '\0'};

	fscanf(in, "%s %d %d\n", n_pts_str, &n_polys, &nAux);

	std::string n_pts_std_str = n_pts_start;
	n_pts_std_str += n_pts_str;
	n_pts = atoi(n_pts_std_str.c_str());

	vtkPolyData *poly = vtkPolyData::New();
	vtkPoints *points = vtkPoints::New();
	vtkCellArray *polys = vtkCellArray::New();

	// Read the points from the file
	float p[3];
	int i,j;
	for(i = 0; i < n_pts; i++)
	{
		fscanf(in, "%f %f %f\n", &p[0], &p[1], &p[2]);
		points->InsertPoint(i,p);
	}

	// Read the triangles from the file
	for(i = 0; i < n_polys; i++)
	{
		vtkIdType n_vertices;
		fscanf(in, "%d ", &n_vertices);
		polys->InsertNextCell(n_vertices);
		vtkIdType vertex_id;
		for(j = 0; j < n_vertices; j++)
		{
		   fscanf(in, "%d ", &vertex_id);
		   polys->InsertCellPoint(vertex_id);
		}
	}
	poly->SetPoints(points);
	points->Delete();
	poly->SetPolys(polys);
	polys->Delete();
	fclose(in);
	return poly;
}

vtkPolyData* LoadOBJ(const char* filename)
{
	vtkOBJReader *OBJReader = vtkOBJReader::New();
	OBJReader->SetFileName(filename);
	OBJReader->Update();
	int errCod = vtkErrorCode::GetLastSystemError();
	if(errCod != vtkErrorCode::NoError)
	{
		OBJReader->Delete();
		return NULL;
	}

	vtkPolyData* polyData = vtkPolyData::New();
	polyData->ShallowCopy(OBJReader->GetOutput());
	OBJReader->Delete();
	return polyData;
}

double computeMeshResolution(vtkPolyData* cloud)
{
	double meshResolution = 0;
	int numdistances = 0;

	vtkIdType* ptrCellIds = cloud->GetPolys()->GetPointer();
	int nEdges;
	int firstVer;
	int secondVer;
	double firstPoint[3], secondPoint[3];
	for(int ce=0; ce<cloud->GetPolys()->GetNumberOfCells(); ce++)
	{
		nEdges = *ptrCellIds;
		ptrCellIds++;
		for(int ed=0; ed<nEdges; ed++)
		{
			firstVer = *(ptrCellIds+ed);
			secondVer = *(ptrCellIds+((ed+1)%nEdges));
			cloud->GetPoint(firstVer, firstPoint);
			cloud->GetPoint(secondVer, secondPoint);
			meshResolution += EuclideanDistance(firstPoint, secondPoint );
			numdistances++;
		}
		ptrCellIds+=nEdges;
	}
	if(numdistances!=0)
	{
		meshResolution/=numdistances;
	}

	return meshResolution;
}

void rotate(vtkPolyData *polyData, const float angle, const float axis[3])
{
	vtkTransform* transform = vtkTransform::New();
	transform->RotateWXYZ(angle, axis[0], axis[1], axis[2]);
	vtkTransformFilter* transformFilter = vtkTransformFilter::New();
	transformFilter->SetInputData(polyData);
	transformFilter->SetTransform(transform);
	transformFilter->Update();
	polyData->ShallowCopy(transformFilter->GetPolyDataOutput());
	transformFilter->Delete();
	transform->Delete();
}

void addGaussianNoise(vtkPolyData* polyData, double noiseSigma, bool random)
{
	
	float* ptrPoints = (float*)polyData->GetPoints()->GetVoidPointer(0);

	for(int po = 0; po < polyData->GetNumberOfPoints(); po++)
	{
		
		*ptrPoints = *ptrPoints + box_muller(0.0, noiseSigma); 
		ptrPoints++;
		*ptrPoints = *ptrPoints + box_muller(0.0, noiseSigma); 
		ptrPoints++;
		*ptrPoints = *ptrPoints + box_muller(0.0, noiseSigma); 
		ptrPoints++;
	}
}

void cleanPolyData(vtkPolyData *polyData, const bool mergePoints, const double tolerance, bool removeNotPolysCells)
{
	if(polyData == NULL)
	{
		return;
	}

	if(removeNotPolysCells)
	{
		polyData->GetLines()->Reset();
		polyData->GetStrips()->Reset();
		polyData->GetVerts()->Reset();
	}

	vtkCleanPolyData* cleanPolyData = vtkCleanPolyData::New();
	cleanPolyData->SetInputData(polyData);
	cleanPolyData->SetAbsoluteTolerance(tolerance);
	cleanPolyData->SetPointMerging(mergePoints);
	cleanPolyData->ToleranceIsAbsoluteOn();
	//cleanPolyData->ConvertLinesToPointsOff();
	//cleanPolyData->ConvertPolysToLinesOff();
	//cleanPolyData->ConvertStripsToPolysOff();
	//cleanPolyData->PieceInvariantOff();

	cleanPolyData->Update();
	polyData->ShallowCopy( cleanPolyData->GetOutput() );
	cleanPolyData->Delete();

	if(removeNotPolysCells)
	{
		polyData->GetLines()->Reset();
		polyData->GetStrips()->Reset();
		polyData->GetVerts()->Reset();
	}
}

int GetBoundaryPoints(vtkPolyData *polydata, bool* &boundaryPointsIds)
{
	boundaryPointsIds = new bool[polydata->GetNumberOfPoints()];
	for(int po=0; po<polydata->GetNumberOfPoints(); po++)
	{
		boundaryPointsIds[po] = false;
	}
	vtkIdType* ptrCells = polydata->GetPolys()->GetPointer();
	int nVertex, ve1, ve2;
	vtkIdList* idList = vtkIdList::New();
	polydata->BuildLinks();
	//for every cell in polydata
	for(int ce=0; ce<polydata->GetNumberOfPolys(); ce++)
	{
		nVertex = *ptrCells;
		ptrCells++;
		//for every edge in polydata
		for(int ve=0; ve<nVertex; ve++)
		{
			ve1 = *(ptrCells + ve);
			ve2 = *(ptrCells + ((ve+1)%nVertex) );
			polydata->GetCellEdgeNeighbors(ce, ve1, ve2,idList);
			if(idList->GetNumberOfIds()<1)
			{
				boundaryPointsIds[ve1] = true;
				boundaryPointsIds[ve2] = true;
			}
		}
		ptrCells += nVertex;
	}
	idList->Delete();
	return polydata->GetNumberOfPoints();
}

void computeNormals(vtkPolyData *polyData)
{
	if(polyData == NULL)
	{
		return;
	}
	vtkPolyDataNormals* polydataNormals = vtkPolyDataNormals::New();
	

	polydataNormals->SplittingOff();
	
	polydataNormals->SetInputData(polyData);
	polydataNormals->ComputePointNormalsOn();
	polydataNormals->AutoOrientNormalsOff();
	polydataNormals->ConsistencyOff();
	polydataNormals->ComputeCellNormalsOff();
	
	polydataNormals->Update();

	polyData->ShallowCopy( polydataNormals->GetOutput() );
	polydataNormals->Delete();
}

/* boxmuller.c           Implements the Polar form of the Box-Muller
                         Transformation

                      (c) Copyright 1994, Everett F. Carter Jr.
                          Permission is granted by the author to use
			  this software for any application provided this
			  copyright notice is preserved.

*/


float frand() { return rand() * 1.0 / RAND_MAX; }


float box_muller(float m, float s)	/* normal random variate generator */
{				        /* mean m, standard deviation s */
	float x1, x2, w, y1;
	static float y2;
	static int use_last = 0;

	if (use_last)		        /* use value from previous call */
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {
			x1 = 2.0 * drand() - 1.0;
			x2 = 2.0 * drand() - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );

		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}

	return( m + y1 * s );
}
