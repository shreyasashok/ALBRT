#ifndef MEMORY_MANAGER_H
#define MEMORY_MANAGER_H

#include "RAJA/RAJA.hpp"
#include <memory>

namespace ALBRT {

/**
 * A lightweight RAJA wrapper for management of host and device memory.
*/
struct MemoryManager {

    /** Get host resource for memory management */
    static RAJA::resources::Resource getHostResource() { return RAJA::resources::Host();}

    /** Get device resource for memory management */
    static RAJA::resources::Resource getDeviceResource() { 
        #ifdef ALBRT_USE_CUDA
        return RAJA::resources::Cuda();
        #else
        return RAJA::resources::Host();
        #endif
    }

    /**
     * Returns a shared pointer to the host memory allocation. 
    */
    template <typename T>
    static std::shared_ptr<T> allocHost(size_t N) {
        return std::shared_ptr<T>(getHostResource().allocate<T>(N), [=] (T* p) {getHostResource().deallocate(p);} );
    }

    /**
     * Returns a shared pointer to the device memory allocation. 
    */
    template <typename T>
    static std::shared_ptr<T> allocDevice(size_t N) {
        return std::shared_ptr<T>(getDeviceResource().allocate<T>(N), [=] (T* p) {getDeviceResource().deallocate(p);});
    }

    /**
     * Returns a shared pointer to the device managed memory allocation. 
    */
    template <typename T>
    static std::shared_ptr<T> allocDeviceManaged(size_t N) {
        return std::shared_ptr<T>(getDeviceResource().allocate<T>(N, RAJA::resources::v1::MemoryAccess::Managed), 
                                                                        [=] (T* p) {getDeviceResource().deallocate(p);});
    }

    /**
     * Returns a raw (unmanaged) pointer to the host memory allocation. 
     * The user is responsible for freeing this memory later!
    */
    template <typename T>
    static T* rawAllocHost(size_t N) {
        return getHostResource().allocate<T>(N);
    }

    /**
     * Returns a raw (unmanaged) pointer to the host memory allocation. 
     * The user is responsible for freeing this memory later!
    */
    template <typename T>
    static T* rawAllocDevice(size_t N) {
        return getDeviceResource().allocate<T>(N);
    }
};

} //end namespace

#endif
