#ifndef ALBRT_H
#define ALBRT_H

#include "raja/RAJAInterface.h"
#include "raja/memoryManager.h"

#include "lbm/octant.h"
#include "lbm/octant.hpp"
#include "lbm/octantCollection.h"
#include "lbm/octantCollection.hpp"

#include "lbm/lattice/D3Q19.h"

#include "mpi/mpiManager.h"

#include "vtk/vtkMPIManager.h"
#include "vtk/octantCollectionVTKManager.h"
#include "vtk/octantCollectionVTKManager.hpp"

namespace ALBRT {
struct ALBRTMain {
    static void initialize();

    static void finalize();
};

} //end namespace

#endif