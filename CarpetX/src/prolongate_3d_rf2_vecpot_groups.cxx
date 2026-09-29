// Parsing and validation of the `vector_potential="Ax Ay Az"` group tag.
//
// Kept out of driver.cxx so that it can be exercised by the standalone test in
// CarpetX/test/unit/, which stubs the handful of Cactus group queries below.
// Everything here runs at startup, once per group per level, and its only job
// is to turn a malformed declaration into a clear error rather than a wrong
// prolongation much later.

#ifndef CARPETX_VECPOT_STANDALONE_TEST
#include "prolongate_3d_rf2_vecpot.hxx"

#include <cctk.h>
#include <util_Table.h>
#endif

#include <array>
#include <cassert>
#include <string>
#include <vector>

namespace CarpetX {

// Defined in driver.cxx.
std::array<int, dim> get_group_indextype(int gi);
std::array<int, dim> get_group_nghostzones(int gi);
std::string get_group_prolongation_type(int gi);

// Raw parse of the `vector_potential="Ax Ay Az"` group tag into three group
// indices, or all -1 when the tag is absent.  Modelled on get_group_fluxes.
// No cross-checks; get_group_vector_potential validates.
static std::array<int, dim> parse_vector_potential_tag(const int gi) {
  assert(gi >= 0);
  const int tags = CCTK_GroupTagsTableI(gi);
  assert(tags >= 0);
  std::vector<char> buf(1000);
  const int iret =
      Util_TableGetString(tags, buf.size(), buf.data(), "vector_potential");
  if (iret == UTIL_ERROR_TABLE_NO_SUCH_KEY)
    buf[0] = '\0'; // default: not a vector potential
  else if (iret < 0)
    assert(0);

  const std::string str(buf.data());
  std::vector<std::string> strs;
  std::size_t end = 0;
  while (end < str.size()) {
    const std::size_t begin = str.find_first_not_of(' ', end);
    if (begin == std::string::npos)
      break;
    end = str.find(' ', begin);
    strs.push_back(str.substr(begin, end - begin));
  }

  std::array<int, dim> groups;
  groups.fill(-1);
  if (strs.empty())
    return groups;

  if (int(strs.size()) != dim)
    CCTK_VERROR("Group %s has a tag vector_potential=\"%s\" naming %d groups; "
                "exactly %d are required, in x, y, z order",
                CCTK_FullGroupName(gi), str.c_str(), int(strs.size()), dim);

  for (int d = 0; d < dim; ++d) {
    std::string str1 = strs[d];
    if (str1.find(':') == std::string::npos) {
      const char *const impl = CCTK_GroupImplementationI(gi);
      str1 = std::string(impl) + "::" + str1;
    }
    const int gi1 = CCTK_GroupIndex(str1.c_str());
    if (gi1 < 0)
      CCTK_VERROR("Group %s has a tag vector_potential=\"%s\" naming unknown "
                  "group %s",
                  CCTK_FullGroupName(gi), str.c_str(), str1.c_str());
    groups[d] = gi1;
  }
  return groups;
}

// The validated triple: the three group indices in x, y, z order, or all -1
// when this group is not part of one.  Every member carries the same tag,
// naming all three, so the triple is discoverable from any one of them.
//
// This is deliberately strict.  The three components are prolonged in one
// pass, so a malformed declaration would otherwise surface much later as a
// wrong prolongation rather than as an error.
std::array<int, dim> get_group_vector_potential(const int gi) {
  const std::array<int, dim> groups = parse_vector_potential_tag(gi);
  if (groups[0] < 0)
    return groups; // no tag

  for (int d = 0; d < dim; ++d)
    for (int d1 = d + 1; d1 < dim; ++d1)
      if (groups[d] == groups[d1])
        CCTK_VERROR("Group %s has a vector_potential tag naming the same group "
                    "more than once",
                    CCTK_FullGroupName(gi));

  // The caller must be one of the three, or the relation is not symmetric.
  bool found = false;
  for (int d = 0; d < dim; ++d)
    found = found || groups[d] == gi;
  if (!found)
    CCTK_VERROR("Group %s has a vector_potential tag that does not name the "
                "group itself",
                CCTK_FullGroupName(gi));

  cGroup gdata;
  int ierr = CCTK_GroupData(gi, &gdata);
  assert(!ierr);

  for (int d = 0; d < dim; ++d) {
    const int gi1 = groups[d];

    if (CCTK_GroupTypeI(gi1) != CCTK_GF)
      CCTK_VERROR("Group %s names %s in its vector_potential tag, but that is "
                  "not a grid function group",
                  CCTK_FullGroupName(gi), CCTK_FullGroupName(gi1));

    // Slot d must carry the edge centering for direction d: cell centred along
    // d, vertex centred across it.  (CarpetX indextype: 1 == cell.)
    const std::array<int, dim> it = get_group_indextype(gi1);
    for (int e = 0; e < dim; ++e)
      if (it[e] != int(e == d))
        CCTK_VERROR("Group %s is named in slot %d of a vector_potential tag and "
                    "must therefore have CENTERING={%s}, but its centering is "
                    "{%c%c%c}",
                    CCTK_FullGroupName(gi1), d,
                    d == 0 ? "cvv" : d == 1 ? "vcv" : "vvc",
                    it[0] ? 'c' : 'v', it[1] ? 'c' : 'v', it[2] ? 'c' : 'v');

    // Every member must name the same three groups, so that each of them
    // reaches the same dispatch decision.  Compare against the raw parse: the
    // full validator would recurse back into this group.
    if (gi1 != gi) {
      const std::array<int, dim> other = parse_vector_potential_tag(gi1);
      if (other != groups)
        CCTK_VERROR("Groups %s and %s declare different vector_potential tags; "
                    "every member of a triple must name all three, in the same "
                    "order",
                    CCTK_FullGroupName(gi), CCTK_FullGroupName(gi1));
    }

    // Prolongation is dispatched once for the triple, so a per-group override
    // on one member could not be honoured.
    if (get_group_prolongation_type(gi1) != "vecpot")
      CCTK_VERROR("Group %s is part of a vector_potential triple and must "
                  "therefore have prolongation_type=\"vecpot\"",
                  CCTK_FullGroupName(gi1));

    // The three MultiFabs are prolonged in one pass and must match in shape.
    cGroup gdata1;
    ierr = CCTK_GroupData(gi1, &gdata1);
    assert(!ierr);
    if (gdata1.numvars != gdata.numvars)
      CCTK_VERROR("Groups %s and %s are in one vector_potential triple but have "
                  "%d and %d variables; they must agree",
                  CCTK_FullGroupName(gi), CCTK_FullGroupName(gi1),
                  gdata.numvars, gdata1.numvars);
    if (get_group_nghostzones(gi1) != get_group_nghostzones(gi))
      CCTK_VERROR("Groups %s and %s are in one vector_potential triple but have "
                  "different numbers of ghost zones; they must agree",
                  CCTK_FullGroupName(gi), CCTK_FullGroupName(gi1));
  }

  return groups;
}
} // namespace CarpetX
