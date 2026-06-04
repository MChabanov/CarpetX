"""Open every PlanesX plane database under a directory with VisIt and plot it.

Drives the running VisIt over the plane output produced by the TestPlanesX
parfiles (e.g. the CI "plane-output" artifact, or a local test-planes-run/):

  visit -nowin -cli -s scripts/visit-check-planes.py <output-root> [<imgdir>]

For every Silo metafile index (*.silo_planes.visit), one leaf .silo file per
plane directory, and every openPMD HDF5 plane file (*.it*.h5), it opens the
database, makes a Pseudocolor plot of every scalar (capped per leaf file),
draws it, and queries MinMax -- forcing VisIt through metadata, mesh
construction (including the AMR mrgtree), and data reads. If <imgdir> is
given, one PNG per metafile is saved there. Exit code 0 iff every database
opened and every attempted plot drew.
"""

import glob
import os
import sys

LEAF_SCALAR_CAP = 16  # leaf files have one quadvar per (variable, box)

argv = [a for a in sys.argv[1:] if not a.endswith(".py")]
root = os.path.abspath(argv[0]) if argv else os.getcwd()
imgdir = os.path.abspath(argv[1]) if len(argv) > 1 else None
if imgdir and not os.path.isdir(imgdir):
    os.makedirs(imgdir)

failures = []
nchecked = 0


def open_db(db):
    """OpenDatabase, trying the explicit OpenPMD plugin for .h5 files (the
    reader does not always claim the generic .h5 extension)."""
    if OpenDatabase(db):
        return True
    if db.endswith(".h5"):
        for plugin in ("OpenPMD_1.0", "openPMD_1.0", "OpenPMD"):
            try:
                if OpenDatabase(db, 0, plugin):
                    return True
            except Exception:  # noqa: BLE001
                pass
    return False


def check_db(db, cap, save_image):
    global nchecked
    nchecked += 1
    rel = os.path.relpath(db, root)
    try:
        if not open_db(db):
            failures.append(rel)
            print("FAIL %-72s OpenDatabase failed" % rel)
            return
        md = GetMetaData(db)  # noqa: F821 (VisIt CLI builtin)
        scalars = [md.GetScalars(i).name for i in range(md.GetNumScalars())]
        plot = scalars[:cap] if cap else scalars
        ndrawn = 0
        for s in plot:
            DeleteAllPlots()  # noqa: F821
            if AddPlot("Pseudocolor", s) and DrawPlots():  # noqa: F821
                Query("MinMax")  # noqa: F821
                ndrawn += 1
        if save_image and imgdir and ndrawn:
            sw = SaveWindowAttributes()  # noqa: F821
            sw.outputToCurrentDirectory = 0
            sw.outputDirectory = imgdir
            sw.fileName = rel.replace(os.sep, "_").replace(".", "_")
            sw.format = sw.PNG
            SetSaveWindowAttributes(sw)  # noqa: F821
            SaveWindow()  # noqa: F821
        DeleteAllPlots()  # noqa: F821
        CloseDatabase(db)  # noqa: F821
        ok = ndrawn == len(plot) and ndrawn > 0
        print("%s %-72s scalars=%d plotted=%d" % (
            "OK  " if ok else "FAIL", rel, len(scalars), ndrawn))
        if not ok:
            failures.append(rel)
    except Exception as exc:  # noqa: BLE001
        failures.append(rel)
        print("FAIL %-72s %s" % (rel, exc))


# Silo metafile time-series indexes: the full multimesh/multivar (+ mrgtree
# AMR shadowing) path -- what a user actually opens.
for db in sorted(glob.glob(os.path.join(root, "*", "*.silo_planes.visit"))):
    check_db(db, cap=0, save_image=True)

# One leaf .silo per plane directory: leaf files must also open standalone
# (this is what plane_readers/silo_reader.py reads).
for d in sorted(glob.glob(os.path.join(root, "*", "*.silo_planes.dir"))):
    leaves = sorted(glob.glob(os.path.join(d, "*.silo")))
    if leaves:
        check_db(leaves[0], cap=LEAF_SCALAR_CAP, save_image=False)

# openPMD HDF5 plane files (the default .bp5 is not readable by VisIt's
# HDF5-only openPMD reader; planes-options writes .h5).
for db in sorted(glob.glob(os.path.join(root, "*", "*.it*.h5"))):
    check_db(db, cap=0, save_image=True)

print("checked %d databases, %d failure(s)" % (nchecked, len(failures)))
for f in failures:
    print("  FAILED: %s" % f)
sys.exit(1 if (failures or nchecked == 0) else 0)
