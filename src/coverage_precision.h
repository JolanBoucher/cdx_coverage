/**
 * @file coverage_precision.h
 * @brief Shared enum selecting the granularity of the coverage computation.
 *
 * This tiny header exists purely so that `cli.hpp` (command-line surface)
 * and `gam_io.h` (the parallel GAM-processing engine) can agree on the same
 * vocabulary type without either one having to include the other. `cli.hpp`
 * has no business knowing about the GAM engine's internals, and `gam_io.h`
 * has no business knowing about CLI11/argument parsing, so the shared type
 * lives here instead, at global scope, next to the other CLI-facing enums
 * declared directly in `cli.hpp` (see `ComponentType`).
 */

#ifndef CDX_COVERAGE_COVERAGE_PRECISION_H
#define CDX_COVERAGE_COVERAGE_PRECISION_H

/**
 * @brief Selects how finely GAM alignments are turned into coverage values.
 *
 * - `Node`: the original, cheap behaviour. Every read that touches a node
 *   increments that node's coverage by exactly one, uniformly across the
 *   node's whole base-pair span once expanded. `process_gam` never looks at
 *   `Edit` submessages in this mode, so the cost of processing a GAM file is
 *   unaffected by how many indels/mismatches it contains, and no extra
 *   memory proportional to read/edit count is allocated. This is the right
 *   choice on memory-constrained machines or for very large/deep GAM files
 *   where the read/edit count itself (not the graph size) would otherwise
 *   dominate memory usage.
 *
 * - `Base`: refines the node-level result down to base-pair resolution by
 *   also tracking, per read, which portions of a touched node were *not*
 *   actually traversed at the base level: positions removed by a deletion
 *   edit, and positions outside a read's alignment start/end when that read
 *   begins or ends partway through a node (soft-clip-like boundary
 *   mappings). Those "coverage gaps" are subtracted from the node-level
 *   result after it has been expanded into base-pair space. See
 *   `coverage_gaps.h` for the geometry and `cov_projection.h` for the
 *   application step.
 */
enum class CoveragePrecision {
    Node,
    Base
};

#endif // CDX_COVERAGE_COVERAGE_PRECISION_H
