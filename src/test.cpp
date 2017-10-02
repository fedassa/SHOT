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
#include "utils.h"
#include <iostream>
#include <vector>

#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include "opencv2/opencv.hpp"
#include <fstream>

struct CmdLineParams
{
    float matchTh;
    double radiusMR;
    double sigmaNoiseMR;
    int nFeat;
    int minNeighbors;
    float rotationAngle;
    float rotationAxis[3];

    int nThreads;

    bool describeColor;
    bool describeShape;

    std::string datapath;
    std::string outputFile;

    int shotShapeBins;
    int shotColorBins;

    CmdLineParams() //default
    {
        shotShapeBins = 10;
        shotColorBins = 30;
        matchTh = 0.75f;
        radiusMR = 20;
        sigmaNoiseMR = 0.3;
        nFeat = 1000;
        minNeighbors = 5;
        describeShape = true;
        describeColor = false;
        rotationAngle = 60.0f;
        rotationAxis[0] = 0.75f;
        rotationAxis[1] = 0.1f;
        rotationAxis[2] = 1-0.75f*0.75f-0.1f*0.1f;
        datapath = "..\\data\\mario.ply";
        outputFile = "..\\data\\shot.txt";
        nThreads = 0;
    }
};

void drawBin(cv::Mat img, int i, int j, char val );

int main(int argc, char** argv)
{
    CmdLineParams params;

    std::ofstream outfile(params.outputFile.c_str());
    if (!outfile.is_open())
        std::cout << "\nWARNING\n It is not possible to write on the requested output file\n";

    double meshRes = 1/params.radiusMR;
    SHOTParams shotParams;
    shotParams.radius = meshRes * params.radiusMR;
    shotParams.localRFradius = meshRes * params.radiusMR;
    shotParams.minNeighbors = 0;
    shotParams.shapeBins = params.shotShapeBins;
    shotParams.colorBins = params.shotColorBins;
    shotParams.describeColor = params.describeColor;
    shotParams.describeShape = params.describeShape;
    shotParams.nThreads = params.nThreads;

    //set normal
    vtkSmartPointer<vtkDoubleArray> pointNormalsArray = vtkSmartPointer<vtkDoubleArray>::New();
    pointNormalsArray->SetNumberOfComponents(3); //3d normals (ie x,y,z)
    pointNormalsArray->SetNumberOfTuples(2);

    // Construct the normal vectors
    double pN1[3] = {0.0, 0.0, 1.0};
    double pN2[3] = {0.0, 0.0, 1.0};

    cv::RNG rng;
    int num = 0;
    for(int ii = 0; ii<1000 ; ii++ )
    {
        //set mesh
        vtkPolyData* mesh = vtkPolyData::New();

        //set point
        double point[3];
        point[0] = rng.uniform(-1.0,1.0);
        point[1] = rng.uniform(-1.0,1.0);
        point[2] = rng.uniform(-1.0,1.0);
        cout<<point[0]<<","<<point[1]<<","<<point[2]<<endl;
        if( point[0]*point[0]+point[1]*point[1]+point[2]*point[2] > 1 ) continue;
        num++;
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        points->InsertNextPoint(0.0, 0.0, 0.0);
        points->InsertNextPoint(point[0], point[1], point[2]);

        mesh->SetPoints(points);

        // Add the data to the normals array
        pointNormalsArray->SetTuple(0, pN1) ;
        pointNormalsArray->SetTuple(1, pN2) ;

        // Add the normals to the points in the polydata
        mesh->GetPointData()->SetNormals(pointNormalsArray);

        if (mesh)
        {
            Feature3D* feat = new Feature3D[1];

            int nActualFeat = 1;/*detector.extract(mesh, feat);*/
            feat->index = 0;
            feat->scale = -1;
            feat->score = 1;
            feat->x = 0;
            feat->y = 0;
            feat->z = 0;

            SHOTDescriptor descriptor(shotParams);

            double** desc;
            descriptor.describe(mesh, feat, desc, nActualFeat);


            //output as text
            cv::Mat features(params.nFeat, descriptor.getDescriptorLength(), CV_32FC1);

            if (outfile.is_open())
                outfile << feat[0].index << " " << feat[0].x << " " << feat[0].y << " " << feat[0].z;

            for (int j = 0; j < descriptor.getDescriptorLength(); j++)
            {
                features.at<float>(0,j) = desc[0][j];
                if (outfile.is_open())
                    outfile << " " << desc[0][j];
            }
            if (outfile.is_open())
                outfile << endl;

            //output as hist
            cv::Mat hist(4,8, CV_8UC1);
            for(int i = 0; i<32; i++)
            {
                double val = 0;
                for(int j = 0; j<11; j++)
                {
                    val += desc[0][i*11+j];
                }
                hist.at<char>(i%4, i/4) = (val==0)? 0: 127+128*val;
            }
            cout<<hist<<endl;
            cv::Mat histImg;
            cv::Size dsize;
            dsize.width = 800; dsize.height = 400;
            cv::resize(hist,histImg, dsize, 0,0, cv::INTER_NEAREST);
            imshow("hist", histImg);
//            cv::waitKey();

            //output as real hist
            cv::Mat histShow( 500, 1000, CV_8UC3, cv::Scalar(0,0,0) );

            //draw the bin
            for(int i = 0; i<4; i++)
            {
                for(int j = 0; j<8; j++)
                {
                    if( hist.at<char>(i,j))
                    {
                        drawBin(histShow, i, j, hist.at<char>(i,j) );
                    }
                }
            }

            //draw the point
            cv::Point2f  centerOfPlot;
            double distRatio = sqrt( (point[0]*point[0]+point[1]*point[1]+point[2]*point[2])
                    /(point[0]*point[0]+point[1]*point[1]) );
            centerOfPlot.x = 200*point[0]*distRatio + 250;
            centerOfPlot.y = -200*point[1]*distRatio + 250;
            if(point[2]>0) centerOfPlot.x += 500;
            circle(histShow, centerOfPlot , 5, cv::Scalar(0,0,255) );

            //draw the edges
            centerOfPlot.x = 250; centerOfPlot.y = 250;
            circle(histShow, centerOfPlot , 200, cv::Scalar(255,0,0) );
            circle(histShow, centerOfPlot , 100, cv::Scalar(255,0,0) );
            centerOfPlot.x = 750; centerOfPlot.y = 250;
            circle(histShow, centerOfPlot , 200, cv::Scalar(255,0,0) );
            circle(histShow, centerOfPlot , 100, cv::Scalar(255,0,0) );
            for(int i = 0; i <8; i++ )
            {
                cv::Point2f  endOfPlot;
                endOfPlot.x = 200*cos( i*DEG_45_TO_RAD ) + 250;
                endOfPlot.y = 200*sin( i*DEG_45_TO_RAD ) + 250;
                centerOfPlot.x = 250; centerOfPlot.y = 250;
                line(histShow, centerOfPlot, endOfPlot, cv::Scalar(255,0,0));
                endOfPlot.x = 200*cos( i*DEG_45_TO_RAD ) + 750;
                endOfPlot.y = 200*sin( i*DEG_45_TO_RAD ) + 250;
                centerOfPlot.x = 750; centerOfPlot.y = 250;
                line(histShow, centerOfPlot, endOfPlot, cv::Scalar(255,0,0));
            }


            imshow("hist2", histShow);
            cv::waitKey();

        }






        mesh->Delete();

    }

    getchar();
}

void drawBin(cv::Mat img, int i, int j, char val )
{
    int centX = 250, centY = 250, limtX = 200, limtY = 200;
    if( i%2 ) centX += 500;
    int ri = 100, ro = 200;
    if( i < 2 ) {limtX = 100; limtY = 100; ri = 0; ro = 100;}
    if( j > 3) limtY = -limtY;
    if( j<2 || j>5 ) limtX = -limtX;

    for(int jj = std::min( centX+limtX, centX ); jj< std::max( centX+limtX, centX ); jj++ )
    {
        for(int ii = std::min( centY+limtY, centY ); ii< std::max( centY+limtY, centY ); ii++ )
        {
            int localx = abs(jj-centX), localy = abs(ii-centY);
            int r = localx*localx+localy*localy;
            if( r>ro*ro || r<ri*ri ) continue;
            if( ((j%2)^((j/2)%2))^(localx<localy) ) continue;

            img.at<cv::Vec3b>(ii,jj)[0] = 000;
            img.at<cv::Vec3b>(ii,jj)[1] = val;
            img.at<cv::Vec3b>(ii,jj)[2] = 000;
        }
    }
}
