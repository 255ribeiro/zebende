# Per-platform PyPI wheels for zebende — options reference

This is a release-engineering reference, not a decided plan. It exists so you can read
through every viable way to stop shipping one universal wheel and start shipping
one properly-tagged wheel per platform, and pick the path you're comfortable with.

## 1. Current state

`zebende` publishes a single wheel via `poetry build` / `poetry publish`, tagged
`py3-none-any`. It bundles **all 6** prebuilt native binaries under
`zebende/zig_libs/<arch>-<os>/`, and at call time `zebende/get_platform.py` picks the
right one for the running machine and loads it with `ctypes.cdll.LoadLibrary`.

Why that's not great:

- Every install downloads/stores ~8 MB of binaries for platforms it will never run on.
- The wheel carries no manylinux/musllinux compliance signal, since it isn't really
  platform-specific — it's a purelib wheel with binaries hidden inside.
- Adding, dropping, or fixing one platform's binary means republishing the entire
  package for every platform.

## 2. The target model

One version number, **multiple wheel files** attached to that same release — the
same pattern numpy, pandas, and most compiled packages use. Each wheel contains only
its own platform's binary and carries a real platform tag, so `pip install zebende`
downloads only the ~0.5–3 MB relevant to the installer's machine and PyPI/pip do the
platform matching for you.

This is not "one version per platform" — it's one version, six wheel files.

## 3. Platform tag reference

| `build.zig` target | Wheel platform tag | Notes |
|---|---|---|
| `aarch64-macos` | `macosx_11_0_arm64` | pick a minimum macOS version you're willing to support |
| `x86_64-macos` | `macosx_10_9_x86_64` | ditto |
| `aarch64-linux` | `manylinux_2_17_aarch64` **or** `musllinux_1_1_aarch64` | see glibc note below |
| `x86_64-linux` | `manylinux_2_17_x86_64` **or** `musllinux_1_1_x86_64` | see glibc note below |
| `aarch64-windows` | `win_arm64` | |
| `x86_64-windows` | `win_amd64` | |

Two things worth knowing before picking tags:

- **PyPI rejects plain `linux_x86_64` / `linux_aarch64` wheels outright.** Linux
  wheels must claim a `manylinux_*` or `musllinux_*` tag or the upload is refused.
  This is the single easiest mistake to make here.
- Zig can target a specific glibc floor directly in the triple —
  `-Dtarget=x86_64-linux-gnu.2.17` — which is exactly what's needed to legitimately
  claim `manylinux_2_17`. The alternative is building against musl
  (`-Dtarget=x86_64-linux-musl`) and using a `musllinux_*` tag instead, which sidesteps
  glibc-version bookkeeping entirely at the cost of a statically-linked musl binary.
- Because the code only uses ctypes (no CPython C-API), the Python/ABI portion of
  every tag can stay `py3-none-...` — there's no `cp3xx`/`abi3` matrix to worry about,
  only the platform portion changes per wheel.

## 4. Build-backend options

How to actually produce a wheel whose tag you control, independent of the machine
building it (since all 6 targets get cross-compiled from one host today).

### setuptools + `bdist_wheel --plat-name`

The long-standing, most common mechanism for exactly this scenario.

**Steps**

1. Migrate `[tool.poetry]` metadata to PEP 621 `[project]`.
2. Add a minimal `setup.py` with a `bdist_wheel` subclass that forces
   `root_is_pure = False`.
3. For each of the 6 targets: stage only that platform's binary into the package
   directory, then run `python -m build --wheel` (or
   `setup.py bdist_wheel --plat-name=<tag>` directly), cleaning between runs.
4. Replace the `pbld`/`ppub` taskipy tasks: `poetry build` → `python -m build`,
   `poetry publish` → `twine upload`.

**Limitations**

- Loses Poetry's dependency/lockfile workflow for the *release build* specifically
  (Poetry can still be kept for day-to-day dev dependency management).
- `setup.py` is legacy-flavored and easy to get subtly wrong — forgetting
  `--plat-name` on one of the 6 runs silently ships a mistagged wheel.
- Needs a small helper script to swap which single binary is staged before each
  build.

### hatchling + custom build hook setting `build_data['tag']`

Modern, explicit, and verified cross-compile-safe: hatch's wheel builder exposes
`build_data['tag']` as a plain string you set yourself — unlike `infer_tag`, which
reads the actual host machine, `tag` is not affected by what machine is running the
build.

**Steps**

1. Migrate to `[project]` + `[tool.hatch.build...]`.
2. Write one `hatch_build.py` implementing `BuildHookInterface` that reads an env
   var (e.g. `ZEBENDE_TARGET`), stages only that platform's binary, and sets
   `build_data['tag']` plus `build_data['pure_python'] = False`.
3. Loop `ZEBENDE_TARGET=<x> hatch build -t wheel` over the 6 targets.
4. Publish with `twine upload` (or `hatch publish`, which exists as a built-in
   alternative to `poetry publish`).

**Limitations**

- Less battle-tested for this exact binary-bundling pattern than setuptools, though
  the mechanism itself is well documented.
- Still requires writing and maintaining a hook script.
- Same "forget one target in the loop" risk as any script-driven approach.

### Poetry + build script / env-var trick

Stay on the current backend and try to make it cooperate.

**Steps**

1. Add `[tool.poetry.build] script = "build.py"` that manipulates the
   `Distribution` to mark it non-pure.
2. Attempt to influence the resulting platform tag via `_PYTHON_HOST_PLATFORM` or
   a similar environment variable before each `poetry build` invocation.

**Limitations**

- **Unverified.** Poetry-core's wheel builder may or may not honor these overrides
  the way setuptools does — there's no documented, supported mechanism for this in
  poetry-core today. Before trusting this for a real release you'd need an
  empirical spike: build one wheel this way and inspect the actual tag it produced.
- Highest risk of the four real options of silently producing the wrong tag, or
  falling back to `any`.

### Hand-rolled / minimal PEP 517 backend

Full manual control, no framework in the loop at all.

**Steps**

1. Implement `build_wheel()` (and optionally `build_sdist()`) directly per PEP 517
   — roughly 50–100 lines.
2. Zip the package directory plus your own hand-written `METADATA`, `WHEEL`, and
   `RECORD` (with correct file hashes) into a `.whl` with whatever filename/tag you
   choose.

**Limitations**

- You own every edge case a framework normally handles for free — RECORD hashing,
  METADATA field formatting, future entry-point support if ever needed.
- Highest maintenance burden of any option here. Only worth it if the framework
  options above prove too fragile in practice.

### cibuildwheel

The standard industry tool for "build wheels for many platforms" — but a poor fit
*specifically for this project*.

**Steps**

- Would mean wrapping the build in cibuildwheel's config and running per-OS matrix
  jobs in CI.

**Limitations**

- cibuildwheel exists to solve "I need native compilers on N different OSes,"
  usually via per-OS CI runners or manylinux Docker containers. Zig already
  cross-compiles all 6 targets from a single machine today — cibuildwheel would add
  orchestration overhead without solving a problem this project actually has.
  Listed for completeness, not recommended here.

### flit

Ruled out for the same reason as plain Poetry: no supported mechanism to force a
platform-specific tag. Mentioned only so it's not a dangling question.

## 5. Automation options

Independent of which backend you pick above — how the multi-wheel build and publish
actually gets run.

### GitHub Actions, tag-triggered

**Steps**: one Linux runner job — no matrix needed, since Zig cross-compiles every
target from one host — that runs `zig build -Dball`, loops the chosen backend's
per-target wheel build, then `twine upload dist/*`, triggered on a version-tag push.

**Limitations**: publishing becomes irreversible-by-default on every tag push, so
this needs deliberate tagging discipline (only tag when genuinely ready) or an added
manual-approval gate in the workflow.

### GitHub Actions, manual `workflow_dispatch`

**Steps**: identical pipeline to the above, but triggered by a human clicking "Run
workflow" instead of firing on every tag push.

**Limitations**: still requires a human to remember to trigger it — not fully
hands-off — but removes the "tagged too early" failure mode.

### Local taskipy-only

**Steps**: restructure the existing `zba` / `libcp` / `pbld` / `ppub` tasks into a
loop over the 6 targets that produces 6 correctly-tagged wheels, then
`twine upload dist/*`.

**Limitations**: no new infrastructure to maintain, but your machine stays the
single point of failure for every release, and nothing prevents a partial/
inconsistent publish (e.g. the network drops after 4 of 6 wheels have uploaded).

## 6. sdist

A source distribution could still be published alongside the wheels for
completeness, but making `pip install` actually *build* from that sdist would
require the installing environment to invoke the Zig toolchain itself at install
time — a separate, larger piece of work than what's covered here. Recommendation:
treat sdist as optional/deprioritized rather than building that capability now.

## 7. At a glance

| Option | Effort | Risk | Cross-compile-safe? |
|---|---|---|---|
| setuptools + `--plat-name` | Medium | Low (proven pattern) | Yes |
| hatchling + build hook | Medium | Low–Medium (less prior art) | Yes |
| Poetry + env-var trick | Low (if it works) | High (unverified) | Unknown |
| Hand-rolled PEP 517 | High | Medium (you own every bug) | Yes |
| cibuildwheel | Medium–High | Low, but wrong tool for this job | Yes, but redundant |
| flit | — | Not viable | No |
| GitHub Actions (tag-triggered) | Medium setup, then free | Needs tagging discipline | — |
| GitHub Actions (`workflow_dispatch`) | Medium setup, then free | Lower (human gate) | — |
| Local taskipy-only | Low setup | Manual, no safety net | — |