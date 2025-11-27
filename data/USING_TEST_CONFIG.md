# Using the test configuration

Use the `data/samples.test-new-analysis.json` config to point rarexsec at the combined Run 1 FHC ROOT file without editing the default samples list.

1. Build the project if you haven't already:
   ```bash
   make -j"$(nproc)"
   ```
2. From the repository root, export the test config so the wrapper uses it:
   ```bash
   export RAREXSEC_CFG="$PWD/data/samples.test-new-analysis.json"
   ```
3. Run any macro through the provided wrapper (this example inspects the samples):
   ```bash
   ./scripts/rarexsec-root.sh -b -q macros/inspect_simulation_samples.C
   ```

`RAREXSEC_BEAMLINE` and `RAREXSEC_PERIODS` default to `numi-fhc` and `run1`, which match the entries in the test config. Override them only if you add more beamlines or periods.
