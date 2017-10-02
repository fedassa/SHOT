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

#ifndef SHOT_UTILS_H_
#define SHOT_UTILS_H_


#include "vtkPolyData.h"
#include <cmath>

inline float EuclideanDistance2(float* point1, float* point2){return (point1[0]-point2[0])*(point1[0]-point2[0]) + (point1[1]-point2[1])*(point1[1]-point2[1]) + (point1[2]-point2[2])*(point1[2]-point2[2]);  };

inline double EuclideanDistance2(double* point1, double* point2){return (point1[0]-point2[0])*(point1[0]-point2[0]) + (point1[1]-point2[1])*(point1[1]-point2[1]) + (point1[2]-point2[2])*(point1[2]-point2[2]);  };

inline float EuclideanDistance(float* point1, float* point2){ return sqrt((point1[0]-point2[0])*(point1[0]-point2[0]) + (point1[1]-point2[1])*(point1[1]-point2[1]) + (point1[2]-point2[2])*(point1[2]-point2[2])); };

inline double EuclideanDistance(double* point1, double* point2){ return sqrt((point1[0]-point2[0])*(point1[0]-point2[0]) + (point1[1]-point2[1])*(point1[1]-point2[1]) + (point1[2]-point2[2])*(point1[2]-point2[2])); };


vtkPolyData* LoadPolyData(std::string filename);

vtkPolyData* LoadPly(const char* filename);

vtkPolyData* LoadOff(const char* filename);

vtkPolyData* LoadOBJ(const char* filename);



void cleanPolyData(vtkPolyData *polyData, const bool mergePoints = true, const double tolerance = 0.0, bool removeNotPolysCells = true);

double computeMeshResolution(vtkPolyData* cloud);

int GetBoundaryPoints(vtkPolyData *polydata, bool* &boundaryPointsIds);

void addGaussianNoise(vtkPolyData* polyData, double noiseSigma = 0.01, bool random = true);

void rotate(vtkPolyData *polyData, const float angle, const float axis[3]);

void computeNormals(vtkPolyData *polyData);


/* boxmuller.c           Implements the Polar form of the Box-Muller
                         Transformation

                      (c) Copyright 1994, Everett F. Carter Jr.
                          Permission is granted by the author to use
			  this software for any application provided this
			  copyright notice is preserved.

*/


float frand();


float box_muller(float m, float s);	/* normal random variate generator */



#endif



			
