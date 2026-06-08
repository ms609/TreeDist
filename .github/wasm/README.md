# webR/WebAssembly package patches

This directory holds WebAssembly (`wasm32-unknown-emscripten`) R package
binaries that have been **rebuilt from patched source** because the binaries
served by `repo.r-wasm.org` fail to load in webR.

## `fastmatch_1.1-8.tgz`

`fastmatch`'s `src/dummy.c` declares

```c
extern void R_registerRoutines(void);
extern void R_useDynamicSymbols(void);
void dummy(void) { R_registerRoutines(); R_useDynamicSymbols(); }
```

purely to silence a CRAN check — `dummy()` is never called, and the package
really registers its routines via `useDynLib` in `NAMESPACE`. On native R this
is harmless, but compiled to WebAssembly it emits a function **import** for
`R_registerRoutines` with type `() -> ()`, which fails to link against webR's
real `R_registerRoutines` (`(i32,i32,i32,i32,i32) -> i32`):

```
LinkError: WebAssembly.Instance(): Import #NN "env" "R_registerRoutines":
imported function does not match the expected type
```

`TreeTools` depends on `fastmatch`, so this breaks the whole `treespace`
Shinylive app at load. The committed binary is built from source with
`src/dummy.c` neutralised.

### Regenerating

Built by `.github/workflows/wasm-build.yml` (dispatch it). It must be compiled
with the **exact** webR toolchain that the bundled shinylive runtime uses —
currently webR **0.6.0** / Emscripten **4.0.8**. If shinylive's bundled webR
changes, re-run the workflow with matching `webr_version` / `emscripten_version`
inputs and replace this file. The `pkgdown.yml` splice is version-guarded and
will warn (rather than ship a mismatched binary) if `fastmatch`'s version
changes upstream.

### Removing

This whole workaround disappears once `TreeTools` moves `fastmatch` to
`Suggests` with a base-`match()` fallback — then webR never bundles `fastmatch`
at all.
