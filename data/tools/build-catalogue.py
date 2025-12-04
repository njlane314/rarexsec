from __future__ import annotations
import argparse
import concurrent.futures
import json
import logging
import os
import sys
import sqlite3
from functools import lru_cache
from pathlib import Path
import xml.etree.ElementTree as ET
import re
import subprocess
import shutil
import uproot
from typing import Set, Tuple, List, Dict, Optional

# --- Paths & defaults ---------------------------------------------------------
DEFAULT_RUN_DB  = "/exp/uboone/data/uboonebeam/beamdb/run.db"
# Choose the recommended v3 beam-quality cut DB for NuMI (EA9*, tortgt*, ...):
DEFAULT_NUMI_DB = "/exp/uboone/data/uboonebeam/beamdb/numi_v3.db"

# Optional prescale DB helper
try:
    sys.path.append("/exp/uboone/data/uboonebeam/beamdb")
    import confDB  # provides confDB.confDB().getAllPrescaleFactors(run)
    _CONFDB = confDB.confDB()
except Exception:
    _CONFDB = None

HADD_TMPDIR = Path("/pnfs/uboone/scratch/users/nlane/tmp/")
MIN_FREE_GB = 5.0
DEFAULT_JOBS = min(8, os.cpu_count() or 1)
CATALOGUE_SUBDIR = Path("data") / "catalogues"
logging.basicConfig(level=logging.INFO, format="%(levelname)s | %(message)s")

# --- Small helpers ------------------------------------------------------------

def norm_run(run: str) -> str:
    digits = "".join(ch for ch in run if ch.isdigit())
    return f"r{digits}" if digits else run

def split_beam_key(beam_key: str) -> tuple[str, str]:
    for sep in ("_", "-"):
        if sep in beam_key:
            b, m = beam_key.split(sep, 1)
            return b, m
    return beam_key, "data"

def runset_token(runs: list[str]) -> str:
    if not runs:
        return "all"
    if len(runs) == 1:
        return norm_run(runs[0])
    try:
        nums = sorted(int("".join(ch for ch in r if ch.isdigit())) for r in runs)
        if nums == list(range(nums[0], nums[-1] + 1)):
            return f"r{nums[0]}-r{nums[-1]}"
    except Exception:
        pass
    return ",".join(norm_run(r) for r in runs)

def summarise_beams_for_name(beam_keys: list[str]) -> str:
    beams = set(beam_keys)
    has_bnb = any(k.startswith("bnb") for k in beams)
    has_nu  = any(k.startswith("numi") for k in beams)
    if has_bnb and has_nu:
        return "nu+bnb"
    if has_bnb:
        return "bnb"
    if has_nu:
        pols = sorted({
            mode.lower()
            for k in beams
            for mode in [split_beam_key(k)[1]]
            if mode and mode.lower() not in {"data", "ext"}
        })
        pol_str = f"[{','.join(pols)}]" if pols else ""
        return f"nu{pol_str}"
    return "multi"

def get_xml_entities(xml_path: Path, content: str | None = None) -> dict[str, str]:
    if content is None:
        content = xml_path.read_text()
    entity_regex = re.compile(r"<!ENTITY\s+([^\s]+)\s+\"([^\"]+)\">")
    return {match.group(1): match.group(2) for match in entity_regex.finditer(content)}

def run_command(command: list[str], execute: bool) -> bool:
    print(f"[COMMAND] {' '.join(command)}")
    if not execute:
        print("[INFO] Dry run mode. HADD command not executed.")
        return True
    if shutil.which(command[0]) is None:
        print(
            f"[ERROR] Command '{command[0]}' not found. Ensure ROOT is set up by running:\n"
            "             source /cvmfs/larsoft.opensciencegrid.org/products/common/etc/setups\n"
            "             setup root",
            file=sys.stderr,
        )
        return False
    try:
        subprocess.run(command, check=True)
        print("[STATUS] HADD Execution successful.")
        return True
    except subprocess.CalledProcessError as exc:
        print(f"[ERROR] HADD Execution failed: {exc}", file=sys.stderr)
        return False

@lru_cache(maxsize=256)
def list_root_files(input_dir: str) -> tuple[str, ...]:
    return tuple(str(p) for p in Path(input_dir).rglob("*.root"))

_POT_CACHE: dict[tuple[str, float], float] = {}

def get_total_pot_from_single_file(file_path: Path) -> float:
    if not file_path or not file_path.is_file():
        return 0.0
    try:
        st = file_path.stat()
        key = (str(file_path), st.st_mtime)
        if key in _POT_CACHE:
            return _POT_CACHE[key]
        tree_path = f"{file_path}:nuselection/SubRun"
        with uproot.open(tree_path) as tree:
            arr = tree["pot"].array(library="np")
            pot = float(arr.sum())
        _POT_CACHE[key] = pot
        return pot
    except Exception as exc:
        print(f"    Warning: Could not read POT from {file_path}: {exc}", file=sys.stderr)
        return 0.0

def get_total_pot_from_files_parallel(file_paths: list[str], max_workers: int) -> float:
    if not file_paths:
        return 0.0
    paths = [Path(p) for p in file_paths]
    mw = min(len(paths), max_workers)
    with concurrent.futures.ThreadPoolExecutor(max_workers=mw) as ex:
        return sum(ex.map(get_total_pot_from_single_file, paths))

def resolve_input_dir(stage_name: str | None, stage_outdirs: dict) -> str | None:
    return stage_outdirs.get(stage_name) if stage_name else None

def pot_sum_via_iterate(input_dir: str) -> float:
    pot_sum = 0.0
    expr = f"{input_dir}/**/*.root:nuselection/SubRun"
    for chunk in uproot.iterate(expr, filter_name=["pot"], library="np", step_size="50 MB"):
        pot_sum += float(chunk["pot"].sum())
    return pot_sum

# --- (run,subrun) extraction from ROOT ---------------------------------------

def _extract_pairs_from_file(f: Path) -> Set[Tuple[int, int]]:
    pairs: Set[Tuple[int, int]] = set()
    try:
        with uproot.open(f) as rf:
            candidates = [
                "nuselection/SubRun", "nuselection/SubRuns",
                "SubRun", "SubRuns",
                "subrun", "subruns",
                "nuselection/Events", "Events",
            ]
            def read_pairs(t) -> Set[Tuple[int, int]]:
                try:
                    bmap = {k.lower(): k for k in t.keys()}
                    if "run" in bmap and "subrun" in bmap:
                        run = t[bmap["run"]].array(library="np")
                        sub = t[bmap["subrun"]].array(library="np")
                        return set(zip(run.astype(int).tolist(), sub.astype(int).tolist()))
                except Exception:
                    pass
                return set()
            for name in candidates:
                try:
                    t = rf[name]
                except Exception:
                    continue
                got = read_pairs(t)
                if got:
                    return got
            try:
                for path, cls in (rf.classnames(recursive=True) or {}).items():
                    if cls == "TTree":
                        try:
                            got = read_pairs(rf[path])
                        except Exception:
                            continue
                        if got:
                            pairs |= got
                return pairs
            except Exception:
                return pairs
    except Exception as e:
        print(f"    Warning: {f}: failed extracting (run,subrun): {e}", file=sys.stderr)
    return pairs

def _collect_pairs_from_files(files: List[str]) -> Set[Tuple[int, int]]:
    pairs: Set[Tuple[int, int]] = set()
    for p in files:
        pairs |= _extract_pairs_from_file(Path(p))
    return pairs

# --- DB helpers: EXT triggers & NuMI metrics ----------------------------------

def _sum_ext_triggers_from_pairs(run_db: str, pairs: Set[Tuple[int, int]]
                                ) -> tuple[int, List[Tuple[int, int]], Dict[int, int]]:
    """
    Sum raw EXTTrig from runinfo for the provided (run,subrun) pairs.
    Returns total, missing pairs (not in DB), and a run->sum map.
    """
    if not pairs:
        return 0, [], {}
    conn = sqlite3.connect(run_db)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()
    cur.execute("PRAGMA temp_store=MEMORY;")
    cur.execute("CREATE TEMP TABLE pairs(run INTEGER, subrun INTEGER);")
    cur.executemany("INSERT INTO pairs(run, subrun) VALUES (?, ?);", list(pairs))
    total = cur.execute("""
        SELECT IFNULL(SUM(r.EXTTrig), 0)
        FROM runinfo r
        JOIN pairs p ON r.run = p.run AND r.subrun = p.subrun;
    """).fetchone()[0]
    by_run_rows = cur.execute("""
        SELECT r.run AS run, IFNULL(SUM(r.EXTTrig), 0) AS ext_sum
        FROM runinfo r
        JOIN pairs p ON r.run = p.run AND r.subrun = p.subrun
        GROUP BY r.run
        ORDER BY r.run;
    """).fetchall()
    missing_rows = cur.execute("""
        SELECT p.run, p.subrun
        FROM pairs p
        LEFT JOIN runinfo r ON r.run = p.run AND r.subrun = p.subrun
        WHERE r.run IS NULL;
    """).fetchall()
    conn.close()
    by_run = {int(r["run"]): int(r["ext_sum"]) for r in by_run_rows}
    missing_pairs = [(int(x["run"]), int(x["subrun"])) for x in missing_rows]
    return int(total), missing_pairs, by_run

def _sum_numi_metrics_from_pairs(numi_db: str, pairs: Set[Tuple[int, int]]
                                ) -> tuple[float, float, Dict[int, Dict[str, float]]]:
    """
    Sum EA9CNT_wcut and tortgt_wcut from the NuMI DB for the given (run,subrun) pairs.
    Returns (EA9CNT_wcut_total, tortgt_wcut_total, by_run_dict).
    """
    ea9_total, tortgt_total = 0.0, 0.0
    if not pairs:
        return ea9_total, tortgt_total, {}
    conn = sqlite3.connect(numi_db)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()
    cur.execute("PRAGMA temp_store=MEMORY;")
    cur.execute("CREATE TEMP TABLE pairs(run INTEGER, subrun INTEGER);")
    cur.executemany("INSERT INTO pairs(run, subrun) VALUES (?, ?);", list(pairs))
    row = cur.execute("""
        SELECT IFNULL(SUM(n.EA9CNT_wcut), 0.0) AS ea9,
               IFNULL(SUM(n.tortgt_wcut), 0.0) AS tortgt
        FROM numi n
        JOIN pairs p ON n.run = p.run AND n.subrun = p.subrun;
    """).fetchone()
    ea9_total   = float(row["ea9"] or 0.0)
    tortgt_total= float(row["tortgt"] or 0.0)
    by_run_rows = cur.execute("""
        SELECT n.run AS run,
               IFNULL(SUM(n.EA9CNT_wcut), 0.0) AS ea9,
               IFNULL(SUM(n.tortgt_wcut), 0.0) AS tortgt
        FROM numi n
        JOIN pairs p ON n.run = p.run AND n.subrun = p.subrun
        GROUP BY n.run
        ORDER BY n.run;
    """).fetchall()
    conn.close()
    by_run = {int(r["run"]): {"ea9_wcut": float(r["ea9"] or 0.0),
                              "tortgt_wcut": float(r["tortgt"] or 0.0)}
              for r in by_run_rows}
    return ea9_total, tortgt_total, by_run

def _prescale_for_run(run: int) -> Optional[float]:
    """
    Fetch the applicable prescale for EXT NUMI window for this run from confDB.
    Preference: post-2018May key, else legacy key, else any EXT_ key.
    Returns None if prescale is unavailable.
    """
    if _CONFDB is None:
        return None
    try:
        pf = _CONFDB.getAllPrescaleFactors(int(run))
        if not pf:
            return None
        for key in ("EXT_NUMIwin_2018May_FEMBeamTriggerAlgo",
                    "EXT_NUMIwin_FEMBeamTriggerAlgo"):
            if key in pf and pf[key] is not None:
                return float(pf[key])
        # Fallback: any EXT_ prescale key
        for k, v in pf.items():
            if k.startswith("EXT_") and v is not None:
                return float(v)
    except Exception:
        return None
    return None

def _ext_triggers_prescaled(by_run_ext: Dict[int, int]) -> tuple[float, Dict[int, Dict[str, float]]]:
    """
    Multiply EXTTrig(run) by the run-specific prescale factor (if available).
    Returns (total_prescaled, runwise_details).
    """
    total = 0.0
    details: Dict[int, Dict[str, float]] = {}
    for run, raw_count in by_run_ext.items():
        ps = _prescale_for_run(run)
        ps_used = float(ps) if (ps is not None and ps > 0) else 1.0
        total += ps_used * float(raw_count)
        details[run] = {"raw_exttrig": float(raw_count), "prescale": ps_used,
                        "exttrig_prescaled": ps_used * float(raw_count)}
    return total, details

# --- Classification -----------------------------------------------------------

def classify_kind(entry: dict) -> str:
    st = (entry.get("sample_type") or "").lower()
    key = (entry.get("sample_key") or "").lower()
    truth = (entry.get("truth_filter") or "").lower()
    if st == "ext":
        return "ext"
    if st == "dirt":
        return "dirt"
    if st == "data":
        return "data"
    if "strange" in key or "strange" in truth or "mc_n_strange" in truth:
        return "strangeness"
    return "beam"

# --- Core per-sample processing ----------------------------------------------

def process_sample_entry(
    entry: dict,
    processed_analysis_path: Path,
    stage_outdirs: dict,
    run_pot: float,
    ext_triggers_nominal: int,
    run_db: str,
    numi_db: str,
    jobs: int,
    is_detvar: bool = False,
) -> bool:
    """
    HADDs the sample, records bookkeeping, and (for data/ext) extracts pairs and DB metrics.
    Side-effects: embeds private keys (starting with "__") in 'entry' for period-level scaling.
    """
    if not entry.get("active", True):
        sample_key = entry.get("sample_key", "UNKNOWN")
        print(f"  Skipping {'detector variation' if is_detvar else 'sample'}: {sample_key} (marked as inactive)")
        return False

    stage_name = entry.get("stage_name")
    sample_key = entry.get("sample_key")
    sample_type = (entry.get("sample_type", "mc") or "mc").lower()

    print(f"  Processing {'detector variation' if is_detvar else 'sample'}: {sample_key} (from stage: {stage_name})")
    print(f"    HADD execution for this {'sample' if not is_detvar else 'detector variation'}: Enabled")

    input_dir = resolve_input_dir(stage_name, stage_outdirs)
    if not input_dir:
        print(f"    Warning: Stage '{stage_name}' not found in XML outdirs. Skipping '{sample_key}'.", file=sys.stderr)
        return False

    output_file = processed_analysis_path / f"{sample_key}.root"
    output_dir = output_file.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    if not os.access(output_dir, os.W_OK):
        print(f"    Error: Output directory '{output_dir}' is not writable.", file=sys.stderr)
        return False
    if output_file.exists():
        try:
            output_file.unlink()
        except OSError as exc:
            print(f"    Error: Cannot remove existing file '{output_file}': {exc}", file=sys.stderr)
            return False

    entry["relative_path"] = output_file.name
    entry["file"] = str(output_file)

    root_files = list(list_root_files(input_dir))
    if not root_files:
        print(f"    Warning: No ROOT files found in {input_dir}. HADD will be skipped.", file=sys.stderr)
        print("    Note: Proceeding to record metadata (if applicable).")
        source_files: List[str] = []
    else:
        use_parallel = jobs > 1
        chosen_tmp = HADD_TMPDIR
        if use_parallel:
            try:
                chosen_tmp.mkdir(parents=True, exist_ok=True)
                free_gb = shutil.disk_usage(chosen_tmp).free / (1024 ** 3)
            except Exception as e:
                print(f"    Warning: Could not evaluate free space in '{chosen_tmp}': {e}. Falling back to single-process hadd.")
                free_gb = 0.0
            if free_gb < MIN_FREE_GB:
                print(f"    Note: Only {free_gb:.1f} GB free in '{chosen_tmp}'. Falling back to single-process hadd.")
                use_parallel = False
        cmd = ["hadd", "-f"]
        if use_parallel:
            cmd += ["-j", str(jobs), "-d", str(chosen_tmp)]
        cmd += [str(output_file), *root_files]
        if not run_command(cmd, True):
            print(f"    Error: HADD failed for {sample_key}. Skipping further processing.", file=sys.stderr)
            return False
        source_files = [str(output_file)]

    # POT (from ntuple) for data/mc (pot_eff)
    pot_eff = 0.0
    try:
        if sample_type in {"mc"} or is_detvar or sample_type == "data":
            pot_eff = get_total_pot_from_files_parallel(source_files, jobs) if source_files else pot_sum_via_iterate(input_dir)
    except Exception as e:
        print(f"    Warning: POT evaluation failed for {sample_key}: {e}", file=sys.stderr)
        pot_eff = 0.0

    # Prepare defaults
    trig_eff = 0
    entry["__pairs"] = []
    entry["__by_run_ext"] = {}
    entry["__ext_trig_prescaled"] = 0.0
    entry["__ext_prescale_details"] = {}
    entry["__ea9_wcut"] = 0.0
    entry["__tortgt_wcut"] = 0.0
    entry["__numi_by_run"] = {}

    # For EXT/DATA: collect (run,subrun) and query DBs
    files_for_pairs = source_files if source_files else list(list_root_files(input_dir))
    pairs = _collect_pairs_from_files(list(files_for_pairs)) if files_for_pairs else set()
    entry["__pairs"] = sorted(list(pairs))

    if sample_type == "ext":
        # Sum EXTTrig (raw) and then apply prescales run-by-run
        total_ext, missing_pairs, by_run = _sum_ext_triggers_from_pairs(run_db, pairs)
        trig_eff = int(total_ext)
        entry["__by_run_ext"] = by_run
        prescaled_total, ps_details = _ext_triggers_prescaled(by_run)
        entry["__ext_trig_prescaled"] = float(prescaled_total)
        entry["__ext_prescale_details"] = ps_details

        print(f"    EXT triggers (raw, from run.db): {trig_eff}")
        if _CONFDB is None:
            print("    Note: confDB not available -> EXT prescale not applied (using prescale=1.0).", file=sys.stderr)
        print(f"    EXT triggers (prescaled): {entry['__ext_trig_prescaled']:.3f}")
        if missing_pairs:
            print(f"    Note: {len(missing_pairs)} (run,subrun) pairs in EXT files not found in run.db (showing up to 5): {missing_pairs[:5]}")

    if sample_type == "data":
        # Sum EA9CNT_wcut and tortgt_wcut for these (run,subrun)
        ea9, tortgt, by_run_numi = _sum_numi_metrics_from_pairs(numi_db, pairs)
        entry["__ea9_wcut"] = float(ea9)
        entry["__tortgt_wcut"] = float(tortgt)
        entry["__numi_by_run"] = by_run_numi
        print(f"    DATA NuMI metrics: EA9CNT_wcut={ea9:.3f}, tortgt_wcut={tortgt:.3f}")

    # Record public bookkeeping
    if sample_type == "mc" or is_detvar:
        entry["pot"] = float(run_pot)
        entry["pot_eff"] = float(pot_eff)
        entry["trig"] = 0
        entry["trig_eff"] = 0
    elif sample_type == "ext":
        entry["pot"] = 0.0
        entry["pot_eff"] = 0.0
        entry["trig"] = int(ext_triggers_nominal)  # from recipe (nominal, if provided)
        entry["trig_eff"] = int(trig_eff)          # from DB (raw)
    elif sample_type == "data":
        entry["pot"] = float(run_pot)
        entry["pot_eff"] = float(pot_eff or 0.0)
        entry["trig"] = 0
        entry["trig_eff"] = 0
    else:
        entry["pot"] = float(run_pot)
        entry["pot_eff"] = float(pot_eff)
        entry["trig"] = 0
        entry["trig_eff"] = 0

    # Clean out internal knobs from the downstream record
    entry.pop("triggers", None)
    entry.pop("stage_name", None)
    return True

# --- XML context --------------------------------------------------------------

def load_xml_context(xml_paths: list[Path]) -> tuple[dict[str, str], dict[str, str]]:
    entities: dict[str, str] = {}
    stage_outdirs: dict[str, str] = {}
    for xml in xml_paths:
        text = xml.read_text()
        entities.update(get_xml_entities(xml, text))
        root = ET.fromstring(text)
        proj = root.find("project")
        if proj is None:
            logging.error("Could not find <project> in XML '%s'", xml)
            continue
        for s in proj.findall("stage"):
            outdir_text = s.findtext("outdir") or ""
            for name, val in entities.items():
                outdir_text = outdir_text.replace(f"&{name};", val)
            stage_outdirs[s.get("name")] = outdir_text
    return entities, stage_outdirs

def default_xmls() -> list[Path]:
    return [
        Path("/exp/uboone/app/users/nlane/production/strangeness_mcc9/srcs/ubana/ubana/searchingforstrangeness/xml/numi_fhc_workflow_core.xml"),
        Path("/exp/uboone/app/users/nlane/production/strangeness_mcc9/srcs/ubana/ubana/searchingforstrangeness/xml/numi_fhc_workflow_detvar.xml"),
        Path("/exp/uboone/app/users/nlane/production/strangeness_mcc9/srcs/ubana/ubana/searchingforstrangeness/xml/numi_fhc_run1_dev.xml"),
    ]

# --- Main ---------------------------------------------------------------------

def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description="Aggregate ROOT samples from a recipe into a Hub-ready catalogue with EXT POT-equivalent scaling.")
    ap.add_argument("--recipe", type=Path, required=True, help="Path to recipe JSON (instance).")
    ap.add_argument("--run-db", type=str, default=DEFAULT_RUN_DB, help="Path to run.db (runinfo table).")
    ap.add_argument("--numi-db", type=str, default=DEFAULT_NUMI_DB, help="Path to NuMI DB (e.g. numi_v3.db).")
    ap.add_argument("--jobs", type=int, default=DEFAULT_JOBS, help="Parallel jobs for hadd / file scan.")
    args = ap.parse_args()

    jobs   = int(args.jobs)
    run_db = str(args.run_db)
    numi_db= str(args.numi_db)

    outdir = repo_root / CATALOGUE_SUBDIR

    with open(args.recipe) as f:
        cfg = json.load(f)

    if cfg.get("role") != "recipe":
        sys.exit(f"Expected role='recipe', found '{cfg.get('role')}'.")
    if cfg.get("recipe_kind", "instance") == "template":
        sys.exit("Refusing to run on a template. Copy it and set recipe_kind='instance'.")

    xml_paths = default_xmls()
    _entities, stage_outdirs = load_xml_context(xml_paths)

    ntuple_dir = Path(cfg["ntuple_base_directory"])
    ntuple_dir.mkdir(parents=True, exist_ok=True)

    beams_in: dict = cfg.get("beamlines", cfg.get("run_configurations", {}))
    beamlines_out: dict = {}

    for beam_key, run_block in beams_in.items():
        beam_active = bool(run_block.get("active", True))
        if not beam_active:
            logging.info("Skipping beam '%s' (active=false).", beam_key)
            continue

        beamline, mode = split_beam_key(beam_key)
        mode = mode.lower()

        for period, run_details in run_block.items():
            if period == "active":
                continue
            logging.info("Processing %s:%s", beam_key, period)
            is_ext_period = (mode == "ext")

            # These come from the recipe block; nominal_pot is used for data/mc, ext_triggers for ext
            run_pot = float(run_details.get("nominal_pot", run_details.get("pot", 0.0))) if not is_ext_period else 0.0
            ext_trig_nominal = int(run_details.get("ext_triggers", 0)) if is_ext_period else 0

            if not is_ext_period and run_pot == 0.0:
                logging.warning("No nominal POT provided for %s:%s (on-beam).", beam_key, period)

            samples_in = run_details.get("samples", []) or []
            if not samples_in:
                logging.info("Skipping %s:%s (no samples).", beam_key, period)
                continue

            samples_out: list[dict] = []
            # We retain references to the processed entries for period-level EXT scaling
            period_data_entries: List[dict] = []
            period_ext_entries:  List[dict] = []

            for sample in samples_in:
                s = dict(sample)
                ok = process_sample_entry(
                    s,
                    ntuple_dir,
                    stage_outdirs,
                    run_pot,
                    ext_trig_nominal,
                    run_db,
                    numi_db,
                    jobs,
                    is_detvar=False,
                )
                if ok and "detector_variations" in s:
                    new_vars = []
                    for dv in s["detector_variations"]:
                        dv2 = dict(dv)
                        if "sample_type" not in dv2:
                            dv2["sample_type"] = "mc"
                        process_sample_entry(
                            dv2,
                            ntuple_dir,
                            stage_outdirs,
                            run_pot,
                            ext_trig_nominal,
                            run_db,
                            numi_db,
                            jobs,
                            is_detvar=True,
                        )
                        new_vars.append(dv2)
                    s["detector_variations"] = new_vars

                # collect for later scaling
                kind = classify_kind(s)
                if kind == "data":
                    period_data_entries.append(s)
                if kind == "ext":
                    period_ext_entries.append(s)

                samples_out.append(s)

            # ----- Period-level EXT effective POT computation -------------------
            # Use the (first) DATA entry in this period as the reference for EA9 & tortgt (wcut)
            ref_data = period_data_entries[0] if period_data_entries else None
            ref_ea9  = float(ref_data["__ea9_wcut"]) if ref_data else 0.0
            ref_tort = float(ref_data["__tortgt_wcut"]) if ref_data else 0.0

            for s in period_ext_entries:
                denom_prescaled = float(s.get("__ext_trig_prescaled", 0.0))
                note_bits = []
                if denom_prescaled <= 0.0:
                    # Fallback to raw ext trig if prescale could not be applied
                    denom_prescaled = float(s.get("trig_eff", 0) or 0.0)
                    note_bits.append("prescale_missing_fallback_to_raw_EXTTrig")
                pot_equiv = 0.0
                if ref_ea9 > 0.0 and denom_prescaled > 0.0 and ref_tort > 0.0:
                    pot_equiv = ref_tort * (ref_ea9 / denom_prescaled)

                # Store all ingredients in the JSON for traceability
                s["pot_equiv"] = float(pot_equiv)
                s["pot_equiv_components"] = {
                    "ea9cnt_wcut_ref": ref_ea9,
                    "tortgt_wcut_ref": ref_tort,
                    "ext_trig_prescaled": denom_prescaled,
                    "ext_prescale_applied": (_CONFDB is not None),
                    "notes": ";".join(note_bits) if note_bits else ""
                }
                # Keep the per-run diagnostics too (small dicts)
                s["ext_prescale_details"] = s.get("__ext_prescale_details", {})
                s["ext_by_run_raw"]       = s.get("__by_run_ext", {})
                s["numi_by_run_ref"]      = ref_data.get("__numi_by_run", {}) if ref_data else {}

            # ----- Pack output for this period ---------------------------------
            for s in samples_out:
                # If "file" not present, rehydrate it from relative_path
                if "file" not in s:
                    rp = s.get("relative_path")
                    if rp:
                        s["file"] = str(ntuple_dir / rp)
                s["kind"] = classify_kind(s)
                dv_list = s.pop("detector_variations", []) or []
                detvars = {}
                for dv in dv_list:
                    dv_file = dv.get("file") or (str(ntuple_dir / dv.get("relative_path")) if dv.get("relative_path") else None)
                    tag = str(dv.get("variation_type") or dv.get("name") or dv.get("sample_key") or f"dv{len(detvars)+1}")
                    dv_desc = {"file": dv_file}
                    for meta_key in ("pot", "pot_eff", "trig", "trig_eff", "pot_equiv"):
                        if meta_key in dv:
                            dv_desc[meta_key] = dv[meta_key]
                    detvars[tag] = dv_desc
                if detvars:
                    s["detvars"] = detvars
                # Drop private keys from output
                for k in list(s.keys()):
                    if k.startswith("__"):
                        s.pop(k, None)
                # Drop fields not needed downstream
                for k in ("sample_key", "sample_type", "truth_filter", "exclusion_truth_filters",
                          "relative_path", "variation_type", "stage_name"):
                    s.pop(k, None)

            run_copy = dict(run_details)
            run_copy["samples"] = samples_out
            beamlines_out.setdefault(beam_key, {})[period] = run_copy

    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / "samples.json"
    catalogue = {"beamlines": beamlines_out}
    with open(out_path, "w") as f:
        json.dump(catalogue, f, indent=4)
    logging.info("Wrote catalogue: %s", out_path)

if __name__ == "__main__":
    main()
