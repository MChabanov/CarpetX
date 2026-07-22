# How to access test logs and build artifacts

Both scripts build and test directly on the host in the Cactus tree at `$CACTUSX` (pre-set in the host shell), so logs and artifacts are ordinary host files. Inside an sbx sandbox (see `agent_scripts/sandbox/README.md`) everything below applies unchanged at the same `$CACTUSX`-relative locations, with `CACTUSX=/home/agent/cactus/Cactus` baked into the sandbox environment.

Per-test logs: `$CACTUSX/TEST/carpetx/<Thorn>/<test>.log`, with value-level differences in `<test>.diffs` (testsuite summary: `$CACTUSX/TEST/carpetx/summary.log`; when tests fail it ends with a `Tests failed:` section).
Full run logs from the scripts: `$CACTUSX/last-build.log` and `$CACTUSX/last-test.log` (the scripts print bounded excerpts on failure; Grep/Read these for anything beyond the excerpt).
Build artifacts: `$CACTUSX/configs/carpetx/`.
Executable: `$CACTUSX/exe/cactus_carpetx`.
