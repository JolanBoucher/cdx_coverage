/**
 * @file query_resolver_test.cpp
 * @brief Unit tests for src/query_resolver.cpp (ComponentResolver).
 *
 * One GoogleTest suite per public responsibility of ComponentResolver.
 * Each suite is preceded by a short block comment describing what is being
 * tested and any invariants assumed; each individual TEST() is preceded by
 * a one-line comment describing that specific case.
 *
 * Two invariants are assumed to be guaranteed upstream (by the CDX index
 * itself, built from a GBZ graph) and are therefore *not* exercised here as
 * "should still work" cases:
 *   - component IDs (cid) are assigned contiguously, 0..n-1;
 *   - each cid is registered exactly once (no re-registration).
 * Two distinct cids sharing the same name, however, is not guaranteed by
 * the CDX format and is actively rejected by register_component() as a
 * data-integrity guard (see DuplicateNameCaseInsensitiveThrows below) -
 * that rejection is itself the behavior under test, not an assumption we
 * work around.
 *
 * Only the public interface (register_component, resolve, get_name, size)
 * is reachable from tests; get_suggestions()'s content is only observable
 * indirectly through the "Did you mean" section of resolve()'s thrown
 * std::runtime_error message.
 */

#include "../src/query_resolver.h"

#include <gtest/gtest.h>

#include <stdexcept>
#include <string>

// =============================================================================
// register_component
//
// Registers a (cid, name) pair for later resolution. Names are compared
// case-insensitively for uniqueness; a duplicate name across two different
// cids is rejected, since it would make name-based resolution ambiguous.
// =============================================================================

// A single registered component is retrievable by cid afterward.
TEST(ComponentResolverRegisterTest, RegisterSingleComponent) {
    ComponentResolver resolver;
    resolver.register_component(0, "chr1");

    EXPECT_EQ(resolver.size(), 1u);
    EXPECT_EQ(resolver.get_name(0), "chr1");
}

// Several contiguous registrations all remain independently retrievable.
TEST(ComponentResolverRegisterTest, RegisterMultipleContiguousComponents) {
    ComponentResolver resolver;
    resolver.register_component(0, "chr1");
    resolver.register_component(1, "chr2");
    resolver.register_component(2, "chr3");

    EXPECT_EQ(resolver.size(), 3u);
    EXPECT_EQ(resolver.get_name(0), "chr1");
    EXPECT_EQ(resolver.get_name(1), "chr2");
    EXPECT_EQ(resolver.get_name(2), "chr3");
}

// Two different cids may not share the same name (case-insensitively) -
// this is a deliberate data-integrity guard, not an incidental limitation.
TEST(ComponentResolverRegisterTest, DuplicateNameCaseInsensitiveThrows) {
    ComponentResolver resolver;
    resolver.register_component(0, "chr1");

    EXPECT_THROW(resolver.register_component(1, "CHR1"), std::runtime_error);
}

// =============================================================================
// resolve() - numeric component ID path
//
// A purely-digit query is interpreted as a numeric cid. Bounds and
// overflow are both validated, always surfacing std::runtime_error (never
// an unrelated exception type leaking from the underlying std::stoull
// call).
// =============================================================================

// A well-formed numeric query resolves to the matching, previously
// registered component.
TEST(ComponentResolverResolveByIdTest, ResolveByExactNumericId) {
    ComponentResolver resolver;
    resolver.register_component(0, "chr1");
    resolver.register_component(1, "chr2");
    resolver.register_component(2, "chr3");

    const auto resolved = resolver.resolve("1");

    EXPECT_EQ(resolved.cid, 1u);
    EXPECT_EQ(resolved.name, "chr2");
}

// Leading zeros are tolerated: "007" is still parsed as the integer 7.
TEST(ComponentResolverResolveByIdTest, ResolveByNumericIdWithLeadingZeros) {
    ComponentResolver resolver;
    for (std::size_t cid = 0; cid <= 7; ++cid) {
        resolver.register_component(cid, "chr" + std::to_string(cid));
    }

    const auto resolved = resolver.resolve("007");

    EXPECT_EQ(resolved.cid, 7u);
    EXPECT_EQ(resolved.name, "chr7");
}

// A numeric query beyond the registered range raises, mentioning the
// range in the message.
TEST(ComponentResolverResolveByIdTest, ResolveNumericIdOutOfRangeThrows) {
    ComponentResolver resolver;
    resolver.register_component(0, "chr1");
    resolver.register_component(1, "chr2");
    resolver.register_component(2, "chr3");

    try {
        resolver.resolve("5");
        FAIL() << "expected std::runtime_error";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("out of range"), std::string::npos);
    }
}

// A digit-only query too large to fit in std::size_t must still surface
// as std::runtime_error (not std::out_of_range leaking from std::stoull),
// per the fix discussed: is_unsigned_integer() only guarantees "all
// digits", not that the value fits.
TEST(ComponentResolverResolveByIdTest, ResolveExcessivelyLargeNumericIdThrowsRuntimeError) {
    ComponentResolver resolver;
    resolver.register_component(0, "chr1");

    const std::string huge_query(30, '9'); // far beyond UINT64_MAX (~1.8e19)

    EXPECT_THROW(resolver.resolve(huge_query), std::runtime_error);
}

// A leading '-' makes the query fail the "purely digits" check, so it is
// treated as a (non-existent) name rather than risking a negative-to-cid
// conversion; must fail cleanly, not underflow or crash.
TEST(ComponentResolverResolveByIdTest, ResolveNegativeNumberIsTreatedAsName) {
    ComponentResolver resolver;
    resolver.register_component(0, "chr1");

    try {
        resolver.resolve("-1");
        FAIL() << "expected std::runtime_error";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("'-1' was not found"), std::string::npos);
    }
}

// =============================================================================
// resolve() - component name path
//
// A non-numeric query is resolved by case-insensitive name lookup.
// =============================================================================

// The exact registered name resolves directly.
TEST(ComponentResolverResolveByNameTest, ResolveByExactName) {
    ComponentResolver resolver;
    resolver.register_component(0, "chr1");

    const auto resolved = resolver.resolve("chr1");

    EXPECT_EQ(resolved.cid, 0u);
    EXPECT_EQ(resolved.name, "chr1");
}

// Name matching ignores case.
TEST(ComponentResolverResolveByNameTest, ResolveByNameCaseInsensitive) {
    ComponentResolver resolver;
    resolver.register_component(0, "chr1");

    const auto resolved = resolver.resolve("ChR1");

    EXPECT_EQ(resolved.cid, 0u);
    EXPECT_EQ(resolved.name, "chr1");
}

// =============================================================================
// resolve() - unresolved query error reporting and suggestions
//
// When neither the numeric-ID nor the name path resolves, resolve() builds
// an error message that includes ranked suggestions when close matches
// exist (prefix matches first, then Levenshtein-distance matches), or a
// generic "--inspect" hint when nothing is close enough.
// =============================================================================

// No registered name is remotely close: the message falls back to the
// generic hint, without a "Did you mean" section.
TEST(ComponentResolverSuggestionsTest, UnknownNameWithNoCloseMatchGetsGenericHint) {
    ComponentResolver resolver;
    resolver.register_component(0, "chr1");
    resolver.register_component(1, "chr2");

    try {
        resolver.resolve("zzzzzzzzzz");
        FAIL() << "expected std::runtime_error";
    } catch (const std::runtime_error& e) {
        const std::string message = e.what();
        EXPECT_NE(message.find("--inspect"), std::string::npos);
        EXPECT_EQ(message.find("Did you mean"), std::string::npos);
    }
}

// A one-character typo of a registered name is close enough (within the
// dynamic Levenshtein threshold) to be suggested.
TEST(ComponentResolverSuggestionsTest, TypoNameSuggestsClosestMatch) {
    ComponentResolver resolver;
    resolver.register_component(0, "apple");
    resolver.register_component(1, "banana");
    resolver.register_component(2, "cherry");

    try {
        resolver.resolve("aple"); // one deletion away from "apple"
        FAIL() << "expected std::runtime_error";
    } catch (const std::runtime_error& e) {
        const std::string message = e.what();
        EXPECT_NE(message.find("Did you mean"), std::string::npos);
        EXPECT_NE(message.find("apple"), std::string::npos);
    }
}

// A query that is a common prefix of several names surfaces all of them
// (prefix matches are preferred over, and returned ahead of, edit-distance
// matches).
TEST(ComponentResolverSuggestionsTest, PrefixQuerySuggestsMatchingNames) {
    ComponentResolver resolver;
    resolver.register_component(0, "chr1");
    resolver.register_component(1, "chr2");
    resolver.register_component(2, "chr10");
    resolver.register_component(3, "other");

    try {
        resolver.resolve("chr"); // prefix of chr1/chr2/chr10, not itself a name
        FAIL() << "expected std::runtime_error";
    } catch (const std::runtime_error& e) {
        const std::string message = e.what();
        EXPECT_NE(message.find("Did you mean"), std::string::npos);
        EXPECT_NE(message.find("chr1"), std::string::npos);
        EXPECT_NE(message.find("chr2"), std::string::npos);
        EXPECT_NE(message.find("chr10"), std::string::npos);
        EXPECT_EQ(message.find("other"), std::string::npos);
    }
}

// An empty query is a prefix of every name (rfind("", 0) == 0 always), so
// it surfaces the first registered names, up to the suggestion cap - a
// direct consequence of the prefix-match rule, not a special case for
// empty input.
TEST(ComponentResolverSuggestionsTest, EmptyQuerySuggestsFirstComponentsUpToLimit) {
    ComponentResolver resolver;
    resolver.register_component(0, "chr1");
    resolver.register_component(1, "chr2");
    resolver.register_component(2, "chr3");
    resolver.register_component(3, "chr4");

    try {
        resolver.resolve("");
        FAIL() << "expected std::runtime_error";
    } catch (const std::runtime_error& e) {
        const std::string message = e.what();
        EXPECT_NE(message.find("Did you mean"), std::string::npos);
        EXPECT_NE(message.find("chr1"), std::string::npos);
        EXPECT_NE(message.find("chr2"), std::string::npos);
        EXPECT_NE(message.find("chr3"), std::string::npos);
        EXPECT_EQ(message.find("chr4"), std::string::npos); // beyond the cap of 3
    }
}

// =============================================================================
// get_name
// =============================================================================

// A registered cid returns its associated name.
TEST(ComponentResolverGetNameTest, GetNameForValidCid) {
    ComponentResolver resolver;
    resolver.register_component(0, "chr1");
    resolver.register_component(1, "chr2");

    EXPECT_EQ(resolver.get_name(1), "chr2");
}

// A cid beyond the registered range raises std::out_of_range.
TEST(ComponentResolverGetNameTest, GetNameOutOfRangeThrows) {
    ComponentResolver resolver;
    resolver.register_component(0, "chr1");

    EXPECT_THROW(resolver.get_name(5), std::out_of_range);
}

// =============================================================================
// size()
// =============================================================================

// A freshly constructed resolver has no registered components.
TEST(ComponentResolverSizeTest, SizeIsZeroWhenEmpty) {
    const ComponentResolver resolver;

    EXPECT_EQ(resolver.size(), 0u);
}

// size() reflects the number of registered components (contiguous case).
TEST(ComponentResolverSizeTest, SizeReflectsRegisteredCount) {
    ComponentResolver resolver;
    resolver.register_component(0, "chr1");
    resolver.register_component(1, "chr2");
    resolver.register_component(2, "chr3");

    EXPECT_EQ(resolver.size(), 3u);
}
