/**
 * @file coverage_gaps.h
 * @brief Pure geometry for locating base-pair-precision "coverage gaps".
 *
 * Background
 * ----------
 * The node-level coverage pipeline credits a whole node with +1 the moment
 * *any* mapping in a read's path touches it, then (elsewhere, in
 * cov_projection.cpp) paints that single value uniformly across every
 * base pair of the node. That uniform fill is exact for the common case -
 * a read that walks all the way through a node, matching it base for base -
 * but it over-counts in two situations that are visible only once a
 * `vg::Mapping`'s `Edit` list and `Position::offset`/`is_reverse` fields are
 * inspected:
 *
 *   1. A deletion edit (`from_length > 0, to_length == 0`) means some bases
 *      *inside* the node were skipped by the read entirely: they were never
 *      sequenced against, so they should not count as covered by that read.
 *
 *   2. The first and last mapping of a read's path may start or end
 *      partway through their node (`Position::offset() != 0`, or the sum of
 *      the mapping's edit `from_length`s falling short of the node's full
 *      length) - the read simply begins/ends its alignment there. Node-level
 *      coverage still credits the *entire* node, even though the portion
 *      before the start (or after the end) of the walk was never reached.
 *
 * This header captures only the coordinate geometry needed to describe
 * those two situations as "gaps": node-local base-pair ranges that received
 * a node-level +1 credit from a specific read but were not actually
 * traversed by it. It has **no** dependency on the vg Protobuf types
 * (`vg::Mapping`, `vg::Edit`) on purpose, so it can be exercised by fast,
 * standalone unit tests; `gam_io.cpp` is the only caller that bridges these
 * pure functions to real Protobuf-parsed alignments.
 *
 * Strand handling
 * ---------------
 * `vg::Position::offset` and the walk performed by a mapping's edits are
 * always expressed relative to the *strand being walked*: for a forward
 * mapping (`is_reverse == false`) that is the node's own forward sequence,
 * walked left to right starting at `offset`; for a reverse mapping
 * (`is_reverse == true`) it is the node's reverse complement, which is
 * equivalent to walking the node's *forward* coordinates right to left,
 * starting `offset` bases in from the node's forward end. Every function
 * below takes `node_length`, `offset` and `is_reverse` and returns ranges
 * already converted to the node's canonical forward-oriented coordinate
 * system (the same one `idx2bp` uses), so callers never have to reason
 * about strand again once they have called into this header.
 */

#ifndef CDX_COVERAGE_COVERAGE_GAPS_H
#define CDX_COVERAGE_COVERAGE_GAPS_H

#include "cdx_types.h"

#include <optional>
#include <vector>

/**
 * @brief A half-open [start, end) range of base pairs, in a node's
 *        canonical forward orientation (the same coordinate system used by
 *        `idx2bp`, independent of any individual read's strand).
 */
struct ForwardRange {
    cdx::SeqLen start; ///< Inclusive start offset, in node-local base pairs.
    cdx::SeqLen end;   ///< Exclusive end offset, in node-local base pairs.

    /// A range is only meaningful (worth emitting as a gap) if non-empty.
    [[nodiscard]] bool empty() const noexcept { return start >= end; }
};

/**
 * @brief A node-local coverage gap contributed by a single read.
 *
 * Represents a range of base pairs, within one node, that node-level
 * coverage credited to a specific read (because the read's path touched
 * that node) but which the read did not actually traverse at base-pair
 * resolution. Applying a gap means subtracting 1 from every base-pair
 * coverage position it spans, after node-level coverage has been expanded
 * into base-pair space.
 */
struct BpGap {
    /// Node identifier expressed as (node_id - nid_min): the same raw
    /// coordinate space used by process_gam's `target` coverage vector, not
    /// yet translated into a topological index or genomic position.
    cdx::Nid nid_offset;

    /// Node-local base-pair range covered by the gap, in the node's
    /// canonical forward orientation.
    ForwardRange range;
};

/**
 * @brief Converts a walk-order sub-range of a mapping's reference-consuming
 *        ("from") span into the node's canonical forward-oriented
 *        coordinates.
 *
 * A mapping walks a node starting at `offset`, consuming `from_length` bases
 * per edit as it goes (in the direction dictated by `is_reverse`). Given a
 * sub-range `[walk_start, walk_end)` expressed in that walk's own order (0 =
 * the first base consumed by the mapping, regardless of strand), this
 * returns the equivalent range in the node's fixed forward orientation.
 *
 * This is the single piece of orientation math every other computation in
 * this header builds on:
 *   - forward walk (`is_reverse == false`): the walk simply advances through
 *     the node's forward coordinates starting at `offset`, so walk position
 *     `w` sits at forward position `offset + w`.
 *   - reverse walk (`is_reverse == true`): the walk instead advances
 *     *backwards* through the node's forward coordinates, starting just
 *     before forward position `node_length - offset`, so walk position `w`
 *     sits at forward position `node_length - offset - 1 - w`.
 *
 * @param node_length Length of the node, in base pairs.
 * @param offset `Position::offset()` of the mapping - how far into the
 *        walked strand the mapping's first edit begins.
 * @param is_reverse `Position::is_reverse()` of the mapping.
 * @param walk_start Start of the sub-range, in walk order (0-based).
 * @param walk_end End of the sub-range (exclusive), in walk order.
 * @return The equivalent [start, end) range in the node's forward
 *         orientation. Clamped to `[0, node_length]` so that malformed or
 *         inconsistent input (e.g. edits that would overrun the node) can
 *         never produce an out-of-bounds or wrapped-around range; callers
 *         should treat a clamped-to-empty result as "nothing to record".
 */
[[nodiscard]] ForwardRange walkSpanToForwardRange(
    cdx::SeqLen node_length,
    cdx::SeqLen offset,
    bool is_reverse,
    cdx::SeqLen walk_start,
    cdx::SeqLen walk_end
);

/**
 * @brief Node-local range never reached by a mapping's walk because the
 *        mapping is the *first* one in its read's path and starts partway
 *        through the node.
 *
 * Only meaningful for the first mapping of a path: internal mappings always
 * start exactly where the previous one left off (offset 0 relative to their
 * own node), so this range is empty for them by construction.
 *
 * @param node_length Length of the node, in base pairs.
 * @param offset `Position::offset()` of the (first) mapping.
 * @param is_reverse `Position::is_reverse()` of the mapping.
 * @return The forward-oriented leading gap, or std::nullopt if the mapping
 *         starts exactly at the edge of the node it walks from (no gap).
 */
[[nodiscard]] std::optional<ForwardRange> leadingUncoveredRange(
    cdx::SeqLen node_length,
    cdx::SeqLen offset,
    bool is_reverse
);

/**
 * @brief Node-local range never reached by a mapping's walk because the
 *        mapping is the *last* one in its read's path and its edits stop
 *        short of the node's far end.
 *
 * Only meaningful for the last mapping of a path, symmetric to
 * `leadingUncoveredRange`.
 *
 * @param node_length Length of the node, in base pairs.
 * @param offset `Position::offset()` of the (last) mapping.
 * @param is_reverse `Position::is_reverse()` of the mapping.
 * @param total_from_consumed Sum of `from_length` across every edit in the
 *        mapping (deletions included, insertions contribute 0) - i.e. how
 *        many reference/node bases the mapping's walk actually consumed.
 * @return The forward-oriented trailing gap, or std::nullopt if the walk
 *         reaches exactly the far end of the node (no gap).
 */
[[nodiscard]] std::optional<ForwardRange> trailingUncoveredRange(
    cdx::SeqLen node_length,
    cdx::SeqLen offset,
    bool is_reverse,
    cdx::SeqLen total_from_consumed
);

#endif // CDX_COVERAGE_COVERAGE_GAPS_H
