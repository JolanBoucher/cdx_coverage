/**
 * @file coverage_gaps.cpp
 * @brief Implementation of the pure coverage-gap geometry declared in
 *        coverage_gaps.h.
 *
 * All three functions here work in a signed 64-bit intermediate type
 * (`Signed`) rather than directly in `cdx::SeqLen` (an unsigned 32-bit
 * type). The reverse-strand formulas below involve subtractions like
 * `node_length - offset - walk_end`, and on malformed/inconsistent input
 * (edits that sum to more than the node actually has left, an `offset`
 * larger than the node itself, etc.) those subtractions can go negative
 * before they are clamped back into range. Doing the subtraction in
 * unsigned arithmetic would silently wrap around to a huge positive value
 * instead of going negative, which `std::clamp` could not then distinguish
 * from a genuinely large, valid range - the signed intermediate is what
 * makes the defensive clamping below actually work.
 */

#include "coverage_gaps.h"

#include <algorithm>
#include <cstdint>

namespace {
    using Signed = std::int64_t;
}

ForwardRange walkSpanToForwardRange(
    const cdx::SeqLen node_length,
    const cdx::SeqLen offset,
    const bool is_reverse,
    const cdx::SeqLen walk_start,
    const cdx::SeqLen walk_end
) {
    const auto len = static_cast<Signed>(node_length);
    const auto off = static_cast<Signed>(offset);
    const auto walk_from = static_cast<Signed>(walk_start);
    const auto walk_to = static_cast<Signed>(walk_end);

    Signed start;
    Signed end;

    if (!is_reverse) {
        // Forward walk: walk position w sits at forward position offset + w
        // (the mapping simply advances left to right from its offset).
        start = off + walk_from;
        end = off + walk_to;
    } else {
        // Reverse walk: walk position w sits at forward position
        // node_length - offset - 1 - w (the mapping advances right to left,
        // starting `offset` bases in from the node's forward end). A
        // walk-order range [walk_from, walk_to) therefore lands, in
        // increasing forward-position order, on
        // [len - offset - walk_to, len - offset - walk_from).
        start = len - off - walk_to;
        end = len - off - walk_from;
    }

    // Defensive clamp: only ever engages on malformed/inconsistent input
    // (e.g. edits whose from_length sums past the node's actual length).
    // Well-formed alignments never hit these bounds.
    start = std::clamp<Signed>(start, 0, len);
    end = std::clamp<Signed>(end, 0, len);

    if (start >= end) {
        return ForwardRange{0, 0};
    }
    return ForwardRange{static_cast<cdx::SeqLen>(start), static_cast<cdx::SeqLen>(end)};
}

std::optional<ForwardRange> leadingUncoveredRange(
    const cdx::SeqLen node_length,
    const cdx::SeqLen offset,
    const bool is_reverse
) {
    const auto len = static_cast<Signed>(node_length);
    const auto off = std::clamp<Signed>(static_cast<Signed>(offset), 0, len);

    // Forward: everything strictly before the walk's own starting offset
    // was never reached by this mapping.
    // Reverse: the walk starts `offset` bases in from the node's forward
    // end, so it is everything *at or after* (node_length - offset) -
    // i.e. towards the forward end - that was never reached instead.
    const Signed start = is_reverse ? (len - off) : Signed{0};
    const Signed end = is_reverse ? len : off;

    if (start >= end) {
        return std::nullopt;
    }
    return ForwardRange{static_cast<cdx::SeqLen>(start), static_cast<cdx::SeqLen>(end)};
}

std::optional<ForwardRange> trailingUncoveredRange(
    const cdx::SeqLen node_length,
    const cdx::SeqLen offset,
    const bool is_reverse,
    const cdx::SeqLen total_from_consumed
) {
    const auto len = static_cast<Signed>(node_length);
    const auto off = std::clamp<Signed>(static_cast<Signed>(offset), 0, len);

    // Clamp the consumed span to what actually remains between the walk's
    // starting offset and the node's far end. Purely defensive: malformed
    // edits could otherwise claim to consume more than the node has left,
    // which would make the formulas below produce a negative (pre-clamp)
    // range on the wrong side of the node.
    const Signed remaining = len - off;
    const auto consumed = std::clamp<Signed>(static_cast<Signed>(total_from_consumed), 0, remaining);

    // Forward walk covers [offset, offset + consumed); anything after that,
    // up to the node's forward end, was never reached.
    // Reverse walk covers [len - offset - consumed, len - offset); anything
    // before that, down to forward position 0, was never reached.
    const Signed start = is_reverse ? Signed{0} : (off + consumed);
    const Signed end = is_reverse ? (len - off - consumed) : len;

    if (start >= end) {
        return std::nullopt;
    }
    return ForwardRange{static_cast<cdx::SeqLen>(start), static_cast<cdx::SeqLen>(end)};
}
