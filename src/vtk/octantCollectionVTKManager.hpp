#ifndef OCTANT_COLLECTION_VTK_MANAGER_HPP
#define OCTANT_COLLECTION_VTK_MANAGER_HPP

#include "octantCollectionVTKManager.h"

#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkVoxel.h>
#include <vtkCellArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>

#include <filesystem>

namespace ALBRT {

template<typename T, class Lattice, int nx, int ny, int nz, int nghost>
OctantCollectionVTKManager<T, Lattice, nx, ny, nz, nghost>::OctantCollectionVTKManager(
                                    OctantCollection<T, Lattice, nx, ny, nz, nghost>& octCollectionIn,
                                    std::string outputDirectory)
:   octCollection(octCollectionIn),
    outputDirectory(outputDirectory)
{
    std::filesystem::create_directory(outputDirectory);
}

template <typename T, class Lattice, int nx, int ny, int nz, int nghost>
void OctantCollectionVTKManager<T, Lattice, nx, ny, nz, nghost>::writeForest(long timestep, std::string filename)
{
    vtkNew<vtkPoints> points;
    vtkNew<vtkCellArray> octantCellArray;

    vtkNew<vtkFloatArray> octantMPIProc;
    vtkNew<vtkFloatArray> octantIndex;
    vtkNew<vtkFloatArray> octantLevel;

    octantMPIProc->SetName("MPIProcess");
    octantMPIProc->SetNumberOfComponents(1);
    octantMPIProc->SetNumberOfTuples(octCollection.localOctantInfos.size());

    octantIndex->SetName("LocalOctantIndex");
    octantIndex->SetNumberOfComponents(1);
    octantIndex->SetNumberOfTuples(octCollection.localOctantInfos.size());

    octantLevel->SetName("OctantLevel");
    octantLevel->SetNumberOfComponents(1);
    octantLevel->SetNumberOfTuples(octCollection.localOctantInfos.size());

    int currentOctantIndex = 0;
    for (auto octantInfo : octCollection.localOctantInfos) {
        
        //Insert the corners of the octant; store the first index;
        //VTK goes from low to high coordinate, in x-y-z order.
        auto idStart = points->InsertNextPoint( octantInfo.origin[0], 
                                                octantInfo.origin[1],  
                                                octantInfo.origin[2]);

        points->InsertNextPoint(octantInfo.origin[0] + octantInfo.physLength[0], 
                                octantInfo.origin[1], 
                                octantInfo.origin[2]);

        points->InsertNextPoint(octantInfo.origin[0], 
                                octantInfo.origin[1] + octantInfo.physLength[1], 
                                octantInfo.origin[2]);

        points->InsertNextPoint(octantInfo.origin[0] + octantInfo.physLength[0], 
                                octantInfo.origin[1] + octantInfo.physLength[1], 
                                octantInfo.origin[2]);

        points->InsertNextPoint(octantInfo.origin[0], 
                                octantInfo.origin[1],  
                                octantInfo.origin[2] + octantInfo.physLength[2]);

        points->InsertNextPoint(octantInfo.origin[0] + octantInfo.physLength[0], 
                                octantInfo.origin[1], 
                                octantInfo.origin[2] + octantInfo.physLength[2]);
                                
        points->InsertNextPoint(octantInfo.origin[0], 
                                octantInfo.origin[1] + octantInfo.physLength[1], 
                                octantInfo.origin[2] + octantInfo.physLength[2]);

        points->InsertNextPoint(octantInfo.origin[0] + octantInfo.physLength[0], 
                                octantInfo.origin[1] + octantInfo.physLength[1], 
                                octantInfo.origin[2] + octantInfo.physLength[2]);

        //Create the voxel and assign the relevant points.
        vtkNew<vtkVoxel> octant;

        for (int iPoint = 0; iPoint < 8; iPoint++)        
            octant->GetPointIds()->SetId(iPoint, idStart + iPoint);

        octantCellArray->InsertNextCell(octant);

        octantMPIProc->SetTuple1(currentOctantIndex, octantInfo.proc);
        octantIndex->SetTuple1(currentOctantIndex, octantInfo.index);
        octantLevel->SetTuple1(currentOctantIndex, octantInfo.level);

        currentOctantIndex++;
    }

    vtkNew<vtkUnstructuredGrid> octantRepresentation;
    octantRepresentation->SetPoints(points);
    octantRepresentation->SetCells(VTK_VOXEL, octantCellArray);
    
    octantRepresentation->GetCellData()->AddArray(octantMPIProc);
    octantRepresentation->GetCellData()->AddArray(octantIndex);
    octantRepresentation->GetCellData()->AddArray(octantLevel);

    std::string filenameProcessed = outputDirectory + filename + "_iT" + std::to_string(timestep) + ".pvtu";

    vtkNew<vtkXMLUnstructuredGridWriter> writer;
    writer->SetFileName(filenameProcessed.c_str());
    writer->SetInputData(octantRepresentation);
    writer->Write();
}

} //end namespace

#endif