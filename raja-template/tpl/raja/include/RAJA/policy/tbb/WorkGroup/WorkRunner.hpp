/*!
 ******************************************************************************
 *
 * \file
 *
 * \brief   Header file containing RAJA WorkRunner class specializations.
 *
 ******************************************************************************
 */

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Copyright (c) 2016-21, Lawrence Livermore National Security, LLC
// and RAJA project contributors. See the RAJA/COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#ifndef RAJA_tbb_WorkGroup_WorkRunner_HPP
#define RAJA_tbb_WorkGroup_WorkRunner_HPP

#include "RAJA/config.hpp"

#include "RAJA/policy/tbb/policy.hpp"

#include "RAJA/pattern/WorkGroup/WorkRunner.hpp"


namespace RAJA
{

namespace detail
{

/*!
 * Runs work in a storage container in order
 * and returns any per run resources
 */
template <typename ALLOCATOR_T,
          typename INDEX_T,
          typename ... Args>
struct WorkRunner<
        RAJA::tbb_work,
        RAJA::ordered,
        ALLOCATOR_T,
        INDEX_T,
        Args...>
    : WorkRunnerForallOrdered<
        RAJA::tbb_for_exec,
        RAJA::tbb_work,
        RAJA::ordered,
        ALLOCATOR_T,
        INDEX_T,
        Args...>
{ };

/*!
 * Runs work in a storage container in reverse order
 * and returns any per run resources
 */
template <typename ALLOCATOR_T,
          typename INDEX_T,
          typename ... Args>
struct WorkRunner<
        RAJA::tbb_work,
        RAJA::reverse_ordered,
        ALLOCATOR_T,
        INDEX_T,
        Args...>
    : WorkRunnerForallReverse<
        RAJA::tbb_for_exec,
        RAJA::tbb_work,
        RAJA::reverse_ordered,
        ALLOCATOR_T,
        INDEX_T,
        Args...>
{ };

}  // namespace detail

}  // namespace RAJA

#endif  // closing endif for header file include guard