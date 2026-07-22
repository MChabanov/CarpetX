# How to generate reference output for a new test

Tests are auto-discovered from `<Thorn>/test/*.par` and compared against a checked-in reference directory of the same basename (`<Thorn>/test/<test>/`). To create that reference for a new `.par`, run the built executable on the par file once and copy the resulting output dir into the test tree.

Everything runs on the host: `$CACTUSX` (pre-set in the host shell) is the Cactus tree, this repo is symlinked into it at `$CACTUSX/arrangements/CarpetX`, and the executable is `$CACTUSX/exe/cactus_carpetx`.

```bash
# Run the par file; with IO::out_dir = $parfile the output lands in ./<test>/
cd <repo>/<Thorn>/test &&
  "$CACTUSX/exe/cactus_carpetx" <test>.par
```

Then curate the generated `<Thorn>/test/<test>/` directory to mirror an existing reference of the same family before checking it in:

- Match the **file set** of a sibling reference dir (e.g. `recover-openpmd/`). Recovery references typically omit the recovered iteration's output (`it000000`) and the `it00000000.bp5` / `performance.yaml` files — keep only the iterations the sibling keeps.
- For physics-invariant changes (e.g. a re-decomposition), `diff` the new `*.tsv` against the sibling reference and confirm they are bit-for-bit identical.

Finally run `./agent_scripts/test.sh` to confirm the new test is discovered and passes (`Number failed -> 0`).
