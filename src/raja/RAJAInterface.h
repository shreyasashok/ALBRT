#ifndef RAJA_INTERFACE_H
#define RAJA_INTERFACE_H

#include "RAJA/RAJA.hpp"

namespace ALBRT {

using host_exec = RAJA::seq_exec; //used for simple forall loops

#if defined(ALBRT_USE_CUDA)

using device_exec = RAJA::cuda_exec<256>; //used for simple forall loops
using device_exec_async = RAJA::cuda_exec_async<256>; //simple forall loops, async

/* Allocator for workgroup */
template < typename T >
struct device_allocator
{
  using value_type = T;

  device_allocator() = default;

  template < typename U >
  constexpr device_allocator(device_allocator<U> const&) noexcept
  { }

  /*[[nodiscard]]*/
  value_type* allocate(size_t num)
  {
    if (num > std::numeric_limits<size_t>::max() / sizeof(value_type)) {
      throw std::bad_alloc();
    }

    value_type *ptr = nullptr;
    cudaErrchk(cudaMallocHost((void **)&ptr, num*sizeof(value_type)));

    if (!ptr) {
      throw std::bad_alloc();
    }

    return ptr;
  }

  void deallocate(value_type* ptr, size_t) noexcept
  {
    cudaErrchk(cudaFreeHost(ptr));
  }
};

using device_workgroup_policy = RAJA::WorkGroupPolicy <
                                 RAJA::cuda_work_async<1024>,
                                 RAJA::unordered_cuda_loop_y_block_iter_x_threadblock_average,
                                 RAJA::constant_stride_array_of_objects,
                                 RAJA::indirect_function_call_dispatch >;

using device_workpool = RAJA::WorkPool< device_workgroup_policy,
                                     int,
                                     RAJA::xargs<>,
                                     device_allocator<char> >;

using device_workgroup = RAJA::WorkGroup< device_workgroup_policy,
                                       int,
                                       RAJA::xargs<>,
                                       device_allocator<char> >;

using device_worksite = RAJA::WorkSite< device_workgroup_policy,
                                     int,
                                     RAJA::xargs<>,
                                     device_allocator<char> >;

#else

using device_exec = RAJA::seq_exec; //used for simple forall loops
using device_exec_async = RAJA::seq_exec; //simple forall loops, async

/* Allocator for workgroup */
template < typename T >
struct device_allocator
{
  using value_type = T;

  device_allocator() = default;

  template < typename U >
  constexpr device_allocator(device_allocator<U> const&) noexcept
  { }

  /*[[nodiscard]]*/
  value_type* allocate(size_t num)
  {
    if (num > std::numeric_limits<size_t>::max() / sizeof(value_type)) {
      throw std::bad_alloc();
    }

    value_type *ptr = (value_type *) malloc(num*sizeof(value_type));

    if (!ptr) {
      throw std::bad_alloc();
    }

    return ptr;
  }

  void deallocate(value_type* ptr, size_t) noexcept
  {
    free(ptr);
  }
};

using device_workgroup_policy = RAJA::WorkGroupPolicy <
                                 RAJA::seq_work,
                                 RAJA::ordered,
                                 RAJA::ragged_array_of_objects,
                                 RAJA::indirect_function_call_dispatch >;

using device_workpool = RAJA::WorkPool< device_workgroup_policy,
                                     int,
                                     RAJA::xargs<>,
                                     device_allocator<char> >;

using device_workgroup = RAJA::WorkGroup< device_workgroup_policy,
                                       int,
                                       RAJA::xargs<>,
                                       device_allocator<char> >;

using device_worksite = RAJA::WorkSite< device_workgroup_policy,
                                     int,
                                     RAJA::xargs<>,
                                     device_allocator<char> >;

#endif

} //end namespace

#endif //end include guard
