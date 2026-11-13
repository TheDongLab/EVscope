# EVscope smoke tests

These tests validate repository installability and selected script-level behavior without requiring the full human reference bundle. They are intended for CI and release checks. Full end-to-end examples require the Zenodo reference bundle documented in the main README.

When seaborn and ARIC are installed, the smoke test also loads the Step 20-22 deconvolution entry point with `--help` to check its imports. CI installs both dependencies; minimal local environments report a skip when either is absent.

Run:

```bash
bash tests/smoke/run_smoke.sh
```
