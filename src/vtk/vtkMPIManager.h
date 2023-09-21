#ifndef VTK_MPI_MANAGER_H
#define VTK_MPI_MANAGER_H

#include <vtkNew.h>
#include <vtkMPIController.h>

namespace ALBRT {

struct VTKMPIManager {

    static vtkNew<vtkMPIController> controller;

    /**
     * Initialize the VTK MPI controller. Call this
     * after MPI has been initialized!
    */
    static void initialize();

    /**
     * Finalize any VTK resources.
    */
    static void finalize();
};

}

#endif