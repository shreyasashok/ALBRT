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
    VTKMPIManager::finalize();
    MPIManager::finalize();    
}

} //end namespace