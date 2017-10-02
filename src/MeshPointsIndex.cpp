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

#include "shot.h"


MeshPointsIndex::MeshPointsIndex(vtkPolyData* polyData)
:m_polyData(polyData)
{
	m_index = vtkPointLocator::New();
	m_index->SetDataSet(polyData);
	m_index->Update();
	m_pointsList = vtkIdList::New();
}

void MeshPointsIndex::SetPolyData(vtkPolyData* polyData)
{
	m_index->SetDataSet(polyData);
	m_index->Update();
}

MeshPointsIndex::~MeshPointsIndex()
{
	m_index->Delete();
	m_pointsList->Delete();
	
}


vtkIdList* MeshPointsIndex::findPointsWithinRadius(const double* point, double radius )
{
	m_index->FindPointsWithinRadius(radius, point, m_pointsList);
	return m_pointsList;
}

int MeshPointsIndex::findNearestPoint(double* point)
{
	return m_index->FindClosestPoint(point);
}


int MeshPointsIndex::findNearestPointWithinRadius(double* point, double radius, double & dist)
{
	return m_index->FindClosestPointWithinRadius(radius, point, dist);
}

int MeshPointsIndex::findNearestPointWithinRadius(float* point, float radius, double & dist)
{
	double doublePoint[] = {point[0],point[1],point[2]};
	return m_index->FindClosestPointWithinRadius(radius, doublePoint, dist);
}

double MeshPointsIndex::findNearestPoint(double* point, int &nearestPointId)
{
	nearestPointId = m_index->FindClosestPoint(point);
	double nearestPoint[3];
	m_polyData->GetPoint(nearestPointId, nearestPoint);
	return sqrt((nearestPoint[0]-point[0])*(nearestPoint[0]-point[0]) + (nearestPoint[1]-point[1])*(nearestPoint[1]-point[1]) + (nearestPoint[2]-point[2])*(nearestPoint[2]-point[2]));
}

vtkIdList* MeshPointsIndex::findNearestNPoints(int n, const double* point)
{
	m_index->FindClosestNPoints(n, point, m_pointsList);
	return m_pointsList;
}

