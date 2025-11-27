# EVENT DISPLAY MACRO GUIDE

## How to run
Execute from the repository root:

- `scripts/root.sh analysis/event_display.C event_display_detector`
- `scripts/root.sh analysis/event_display.C event_display_semantic`

These produce per-plane PNGs plus a JSON manifest in `plots/event_display_detector` or `plots/event_display_semantic`.

## Things you may want to tweak
- If your run/subrun/event branch names differ, adjust `opt.cols.run`, `opt.cols.sub`, and `opt.cols.evt`.
- To show only recognised signal in MC semantic plots, uncomment `base_sel = "(" + base_sel + ") && is_signal && recognised_signal";` in `analysis/event_display.C`.
- Tune `det_min`, `det_max`, and `use_log_z` to match the ADC dynamic range of your detector images.
