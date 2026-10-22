"""Run the three saved GWAS-QTL examples once."""

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parent
CONFIG_PATH = ROOT / "config.json"


def load_json(path):
    return json.loads(path.read_text())


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def read_table(path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or len(rows) != len({row["snp"] for row in rows}):
        raise ValueError(f"Empty or duplicate SNPs: {path}")
    return rows


def check_tools(config):
    rscript = os.environ.get("RSCRIPT", "Rscript")
    smr = os.environ.get("SMR_BIN", "smr")
    rscript = shutil.which(rscript)
    smr = shutil.which(smr)
    if not rscript or not smr:
        raise RuntimeError("Rscript and SMR_BIN must both be installed")
    version = subprocess.run([smr, "--version"], capture_output=True,
                             text=True, timeout=20, check=False)
    if f"Version {config['software']['SMR']}" not in version.stdout + version.stderr:
        raise RuntimeError("SMR binary version does not match config.json")
    return rscript, smr


def build_pair(spec, config, gwas_rows, ld_hash, output):
    task_dir = ROOT / "data" / spec["slug"]
    meta = load_json(task_dir / "metadata.json")
    qtl_rows = read_table(task_dir / "qtl.tsv")
    qtl = {row["snp"]: row for row in qtl_rows}
    if set(qtl) != {row["snp"] for row in gwas_rows}:
        raise ValueError("GWAS and QTL SNP sets differ")
    if meta["identity"]["gwas_id"] != config["gwas_id"]:
        raise ValueError("GWAS identity differs")
    if meta["identity"]["feature_id"] != spec["feature"]:
        raise ValueError("Feature identity differs")
    if meta["identity"]["qtl_context_id"] != spec["context_id"]:
        raise ValueError("QTL context differs")
    if meta["window"] != config["window"]:
        raise ValueError("GWAS window differs")

    snps, positions, gb, gs, qb, qs, qm = [], [], [], [], [], [], []
    flips = 0
    for row in gwas_rows:
        snp = row["snp"]
        other = qtl[snp]
        chrom, pos, a1, a2 = snp.split(":")
        if (chrom != meta["chromosome"] or int(pos) != int(row["position"])
                or not meta["window_start"] <= int(pos) <= meta["window_end"]
                or (row["a1"], row["a2"]) != (a1, a2)
                or other["position"] != row["position"]):
            raise ValueError(f"Invalid position or allele: {snp}")
        if (other["a1"], other["a2"]) == (a1, a2):
            sign = 1
        elif (other["a1"], other["a2"]) == (a2, a1):
            sign = -1
            flips += 1
        else:
            raise ValueError(f"Unmatched QTL alleles: {snp}")
        snps.append(snp)
        positions.append(int(pos))
        gb.append(float(row["beta"]))
        gs.append(float(row["se"]))
        qb.append(sign * float(other["beta"]))
        qs.append(float(other["se"]))
        qm.append(float(other["MAF"]))

    paired = {
        "gwas": dict(meta["gwas_fields"], beta=gb, se=gs),
        "identity": meta["identity"],
        "position": positions,
        "qtl": dict(meta["qtl_fields"], MAF=qm, beta=qb, se=qs),
        "snp": snps,
        "task_id": meta["source_task_id"],
        "window": meta["window"],
    }
    source_hash = hashlib.sha256(json.dumps(
        paired, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    if source_hash != spec["input_sha256"]:
        raise ValueError("Reconstructed pair differs from the saved source input")
    for side in ("gwas", "qtl"):
        values = paired[side]
        if (not values["scale_certified"] or values["type"] not in ("cc", "quant")
                or any(not math.isfinite(b) or not math.isfinite(s) or s <= 0
                       for b, s in zip(values["beta"], values["se"]))):
            raise ValueError(f"Invalid {side} statistics")
        if values["type"] == "cc" and values["effect_scale"] != "log_odds_ratio":
            raise ValueError("GWAS effect scale differs")
        if values["type"] == "quant" and (values["N"] <= 1 or
                any(not 0 < maf <= 0.5 for maf in values["MAF"])):
            raise ValueError("Invalid quantitative sample size or MAF")

    susie = load_json(task_dir / "susie.json")
    if susie["task_id"] != paired["task_id"] or not set(susie["snp"]) <= set(snps):
        raise ValueError("SuSiE SNPs or task identity differ")
    for side in ("gwas", "qtl"):
        values = susie[side]
        if (values["snp"] != susie["snp"] or values["N"] <= 1
                or not values["scale_certified"]
                or Path(values["ld_binary"]).name != ld_hash + ".ld.bin"):
            raise ValueError(f"Invalid SuSiE {side} inputs")
    for suffix in ("besd", "epi", "esi"):
        if not (task_dir / "smr" / f"eqtl.{suffix}").is_file():
            raise FileNotFoundError(f"Missing SMR QTL {suffix}")

    save_json(output / "abf_input.json", paired)
    return {"snps": len(snps), "susie_snps": len(susie["snp"]),
            "allele_flips": flips, "source_input_sha256": source_hash}


def run_command(name, command, output, environment):
    started = time.monotonic()
    log = output / f"{name}.log"
    with log.open("x") as handle:
        try:
            result = subprocess.run(command, cwd=ROOT, env=environment,
                                    stdout=handle, stderr=subprocess.STDOUT,
                                    timeout=7200, check=False)
            code = result.returncode
        except subprocess.TimeoutExpired:
            handle.write("TIMEOUT; no automatic retry\n")
            code = 124
    return {"status": "PASS" if code == 0 else "FAIL",
            "exit_code": code, "seconds": time.monotonic() - started,
            "log": str(log)}


def run_smr(spec, config, smr, output, environment):
    shared = ROOT / "data/shared/smr"
    eqtl = ROOT / "data" / spec["slug"] / "smr/eqtl"
    settings = config["smr"]
    base = [smr, "--bfile", str(shared / "reference"),
            "--gwas-summary", str(shared / "gwas.ma"),
            "--beqtl-summary", str(eqtl),
            "--peqtl-smr", settings["peqtl_smr"],
            "--peqtl-heidi", settings["peqtl_heidi"],
            "--ld-lower-limit", settings["ld_lower_limit"],
            "--ld-upper-limit", settings["ld_upper_limit"],
            "--heidi-min-m", settings["heidi_min_m"],
            "--heidi-max-m", settings["heidi_max_m"],
            "--cis-wind", settings["cis_wind_kb"],
            "--thread-num", settings["thread_num"]]
    rows, checks = [], []
    for index in range(2):
        prefix = output / f"smr_{index}"
        check = run_command(f"SMR_{index}", base + ["--out", str(prefix)],
                            output, environment)
        checks.append(check)
        if check["exit_code"]:
            return checks, None
        with prefix.with_suffix(".smr").open(newline="") as handle:
            result = list(csv.DictReader(handle, delimiter="\t"))
        if len(result) != 1 or result[0]["probeID"] != spec["feature"]:
            raise ValueError("SMR returned the wrong feature or row count")
        rows.append(result[0])
    if rows[0] != rows[1]:
        raise ValueError("SMR repeatability check failed")
    value = {"p_smr": float(rows[0]["p_SMR"]),
             "p_heidi": float(rows[0]["p_HEIDI"]),
             "heidi_snps": int(rows[0]["nsnp_HEIDI"]),
             "repeatability": "PASS"}
    save_json(output / "smr.json", value)
    return checks, value


def main():
    parser = argparse.ArgumentParser(description="Run three saved GWAS-QTL examples")
    parser.add_argument("--output", default="results", help="New output directory")
    args = parser.parse_args()
    config = load_json(CONFIG_PATH)
    if (config["susie"]["maxit"] != 2000 or config["software"]["coloc"] != "5.2.3"
            or config["software"]["susieR"] != "0.14.2"):
        raise ValueError("Method profile differs from the documented test")
    rscript, smr = check_tools(config)
    gwas = read_table(ROOT / "data/gwas.tsv")
    ld_file = ROOT / "data/shared/ld.bin"
    ld_hash = hashlib.sha256(ld_file.read_bytes()).hexdigest()
    if ld_hash != "a1676766fcd058810de0b6f5a68067aca11fb601814bcda735e3aecf62840815":
        raise ValueError("Signed LD file differs from the saved test input")
    output_root = Path(args.output).resolve()
    output_root.mkdir(parents=True, exist_ok=False)
    save_json(output_root / "RUN_INTENT.json", {
        "config_sha256": hashlib.sha256(CONFIG_PATH.read_bytes()).hexdigest(),
        "ld_sha256": ld_hash, "started_epoch": time.time(),
        "tasks": [task["slug"] for task in config["tasks"]]})
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                       MKL_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")
    summary = []
    for spec in config["tasks"]:
        task_output = output_root / spec["slug"]
        task_output.mkdir()
        state = {"slug": spec["slug"], "status": "IN_PROGRESS", "steps": []}
        started = time.monotonic()
        try:
            state["preparation"] = build_pair(spec, config, gwas, ld_hash,
                                               task_output)
            state["preparation_seconds"] = time.monotonic() - started
        except Exception as exc:
            state.update(status="INPUT_HOLD", error=str(exc))
            summary.append(state)
            save_json(task_output / "status.json", state)
            continue

        commands = [
            ("ABF", [rscript, str(ROOT / "scripts/abf.R"), str(CONFIG_PATH),
                     str(task_output / "abf_input.json"), str(task_output / "abf.json")]),
            ("SuSiE", [rscript, str(ROOT / "scripts/susie.R"), str(CONFIG_PATH),
                       str(ROOT / "data" / spec["slug"] / "susie.json"),
                       str(ld_file), str(task_output / "susie.json")]),
        ]
        for name, command in commands:
            state["steps"].append({"name": name, **run_command(name, command,
                                                                 task_output, environment)})
        try:
            smr_checks, _ = run_smr(spec, config, smr, task_output, environment)
            state["steps"].extend({"name": f"SMR_{i}", **check}
                                  for i, check in enumerate(smr_checks))
        except Exception as exc:
            state["steps"].append({"name": "SMR", "status": "FAIL", "error": str(exc)})

        if all(step["status"] == "PASS" for step in state["steps"]):
            try:
                abf = load_json(task_output / "abf.json")
                susie = load_json(task_output / "susie.json")
                smr_result = load_json(task_output / "smr.json")
                actual = {"abf_pp_h4": abf["pp_h4"],
                          "susie_pp_h4": susie["pp_h4"],
                          "p_smr": smr_result["p_smr"],
                          "p_heidi": smr_result["p_heidi"]}
                checks = {key: math.isclose(value, spec["expected"][key],
                                            rel_tol=1e-9, abs_tol=1e-12)
                          for key, value in actual.items()}
                state.update(actual=actual, expected=spec["expected"],
                             comparison=checks,
                             status="PASS" if all(checks.values()) else "MISMATCH")
            except Exception as exc:
                state.update(status="INCOMPLETE", error=str(exc))
        else:
            state["status"] = "NATIVE_PARTIAL"
        save_json(task_output / "status.json", state)
        summary.append(state)
        save_json(output_root / "progress.json", summary)

    with (output_root / "results.tsv").open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["window", "context", "feature", "status", "ABF_PP_H4",
                         "SuSiE_PP_H4", "SMR_P", "HEIDI_P"])
        for spec, state in zip(config["tasks"], summary):
            values = state.get("actual", {})
            writer.writerow([config["window_id"], spec["context"], spec["feature"],
                             state["status"], values.get("abf_pp_h4", ""),
                             values.get("susie_pp_h4", ""),
                             values.get("p_smr", ""), values.get("p_heidi", "")])
    passed = all(state["status"] == "PASS" for state in summary)
    save_json(output_root / "QA.json", {"status": "PASS" if passed else "INCOMPLETE",
                                          "tasks": summary})
    print(f"{'PASS' if passed else 'INCOMPLETE'}: {output_root / 'results.tsv'}")
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
