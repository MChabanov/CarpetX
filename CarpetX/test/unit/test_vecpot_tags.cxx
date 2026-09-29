// Standalone test of the vector_potential group-tag parser and validator in
// CarpetX/src/prolongate_3d_rf2_vecpot_groups.cxx.
//
// That code exists to turn a malformed declaration into a clear error instead
// of a wrong prolongation much later, so the error paths are the point: each
// case below is a declaration a user could plausibly write.
//
// The handful of Cactus group queries it uses are stubbed against a synthetic
// group registry.  CCTK_VERROR throws here instead of aborting, so failures
// can be asserted on.

#define CARPETX_VECPOT_STANDALONE_TEST 1

#include <array>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

namespace CarpetX {
constexpr int dim = 3;
}

// ---- Cactus stubs ----------------------------------------------------------

#define CCTK_GF 1
#define UTIL_ERROR_TABLE_NO_SUCH_KEY (-2)

struct cGroup {
  int numvars;
};

struct FakeGroup {
  std::string impl, name;
  std::array<int, 3> centering; // CarpetX convention: 1 == cell
  std::array<int, 3> nghost{{3, 3, 3}};
  int numvars = 1;
  int grouptype = CCTK_GF;
  std::string ptype = "vecpot";
  std::string tag; // the vector_potential tag, empty for none
  std::string full() const { return impl + "::" + name; }
};

static std::vector<FakeGroup> registry;
static std::string last_full_name;

static int CCTK_GroupIndex(const char *n) {
  for (size_t i = 0; i < registry.size(); ++i)
    if (registry[i].full() == n)
      return int(i);
  return -1;
}
static const char *CCTK_FullGroupName(int gi) {
  last_full_name = registry.at(gi).full();
  return last_full_name.c_str();
}
static const char *CCTK_GroupImplementationI(int gi) {
  return registry.at(gi).impl.c_str();
}
static int CCTK_GroupTypeI(int gi) { return registry.at(gi).grouptype; }
static int CCTK_GroupData(int gi, cGroup *g) {
  g->numvars = registry.at(gi).numvars;
  return 0;
}
static int CCTK_GroupTagsTableI(int gi) { return gi; }
static int Util_TableGetString(int table, int len, char *buf, const char *key) {
  if (std::strcmp(key, "vector_potential") != 0)
    return UTIL_ERROR_TABLE_NO_SUCH_KEY;
  const std::string &t = registry.at(table).tag;
  if (t.empty())
    return UTIL_ERROR_TABLE_NO_SUCH_KEY;
  if (int(t.size()) + 1 > len)
    return -1;
  std::strcpy(buf, t.c_str());
  return int(t.size());
}

#define CCTK_VERROR(...)                                                       \
  do {                                                                         \
    char _m[2048];                                                             \
    std::snprintf(_m, sizeof _m, __VA_ARGS__);                                 \
    throw std::runtime_error(_m);                                              \
  } while (0)

namespace CarpetX {
std::array<int, dim> get_group_indextype(int gi) {
  return registry.at(gi).centering;
}
std::array<int, dim> get_group_nghostzones(int gi) {
  return registry.at(gi).nghost;
}
std::string get_group_prolongation_type(int gi) { return registry.at(gi).ptype; }
} // namespace CarpetX

#include "../../src/prolongate_3d_rf2_vecpot_groups.cxx"

// ---- harness ---------------------------------------------------------------

static int failures = 0;

static void ok(const char *what, bool cond) {
  std::printf("  %-62s %s\n", what, cond ? "ok" : "FAIL");
  if (!cond)
    ++failures;
}

// Build a registry of three edge groups with the given tag on each.
static void make_triple(const std::string &tag = "Avec_x Avec_y Avec_z") {
  registry.clear();
  const char *nm[3] = {"Avec_x", "Avec_y", "Avec_z"};
  for (int d = 0; d < 3; ++d) {
    FakeGroup g;
    g.impl = "TestVP";
    g.name = nm[d];
    g.centering = {{d == 0, d == 1, d == 2}};
    g.tag = tag;
    registry.push_back(g);
  }
}

static std::string rejects(int gi) {
  try {
    CarpetX::get_group_vector_potential(gi);
  } catch (const std::exception &e) {
    return e.what();
  }
  return "";
}

static void expect_reject(const char *what, int gi, const char *needle) {
  const std::string msg = rejects(gi);
  const bool got = !msg.empty() && msg.find(needle) != std::string::npos;
  std::printf("  %-62s %s\n", what, got ? "ok" : "FAIL");
  if (!got) {
    ++failures;
    std::printf("      message was: %s\n", msg.empty() ? "(accepted)" : msg.c_str());
  }
}

int main() {
  std::printf("vector_potential group tag\n\n[1] accepted declarations\n");

  {
    make_triple();
    bool all = true;
    for (int d = 0; d < 3; ++d) {
      const auto g = CarpetX::get_group_vector_potential(d);
      all = all && g[0] == 0 && g[1] == 1 && g[2] == 2;
    }
    ok("a well-formed triple is accepted from every member", all);
  }
  {
    // Fully qualified names must work too.
    make_triple("TestVP::Avec_x TestVP::Avec_y TestVP::Avec_z");
    const auto g = CarpetX::get_group_vector_potential(1);
    ok("fully qualified group names resolve", g[0] == 0 && g[2] == 2);
  }
  {
    registry.clear();
    FakeGroup g;
    g.impl = "T";
    g.name = "plain";
    g.centering = {{0, 0, 0}};
    g.ptype = "ddf";
    registry.push_back(g);
    const auto r = CarpetX::get_group_vector_potential(0);
    ok("a group with no tag returns no triple", r[0] < 0);
  }

  std::printf("\n[2] rejected declarations\n");

  make_triple("Avec_x Avec_y");
  expect_reject("wrong number of names", 0, "exactly 3 are required");

  make_triple("Avec_x Avec_y Nope");
  expect_reject("unknown group named", 0, "unknown group");

  make_triple("Avec_x Avec_y Avec_y");
  expect_reject("same group named twice", 0, "more than once");

  make_triple();
  registry[2].tag = "Avec_x Avec_z Avec_y";
  expect_reject("members disagree about the triple", 0, "different vector_potential");

  make_triple();
  registry[1].centering = {{0, 0, 1}}; // {vvc}, belongs in slot z not y
  expect_reject("wrong centering for the slot", 0, "CENTERING={vcv}");

  make_triple();
  registry[2].ptype = "ddf";
  expect_reject("a member without prolongation_type=vecpot", 0,
                "must therefore have prolongation_type");

  make_triple();
  registry[1].numvars = 2;
  expect_reject("members with different variable counts", 0, "they must agree");

  make_triple();
  registry[2].nghost = {{2, 3, 3}};
  expect_reject("members with different ghost widths", 0, "ghost zones");

  make_triple();
  registry[0].grouptype = 0; // not CCTK_GF
  expect_reject("a member that is not a grid function", 0, "not a grid function");

  {
    // A group that names three others but not itself.
    make_triple();
    FakeGroup g;
    g.impl = "TestVP";
    g.name = "stray";
    g.centering = {{1, 0, 0}};
    g.tag = "Avec_x Avec_y Avec_z";
    registry.push_back(g);
    expect_reject("a tag that does not name the group itself", 3,
                  "does not name the");
  }

  std::printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "PASSED", failures,
              failures == 1 ? "" : "s");
  return failures ? 1 : 0;
}
