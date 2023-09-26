#include "ALBRT.h"
#include <p4est_base.h>

namespace ALBRT {

void ALBRTMain::initialize()
{
    MPIManager::initialize();
    p4est_init(NULL, SC_LP_PRODUCTION);
    VTKMPIManager::initialize();

}


void ALBRTMain::finalize()
{
    // Finalize in reverse order of initialization
    
    VTKMPIManager::finalize();
    // p4est has no finalize()
    MPIManager::finalize();    
}

} //end namespace