#include "vtkMPIManager.h"

namespace ALBRT {

vtkNew<vtkMPIController> VTKMPIManager::controller;

void VTKMPIManager::initialize()
{
    controller->Initialize();
    vtkMultiProcessController::SetGlobalController(controller);
}

void VTKMPIManager::finalize()
{
}

}