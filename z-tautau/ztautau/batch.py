"""Run one per-file job over many files, robustly and resumably (adapted from z-mumu/zmumu/batch.py).

dCache NFS reads occasionally stall in uninterruptible I/O under many parallel readers, and EOS
reads over xrootd can fail transiently. A stalled multiprocessing worker cannot be killed, so every
file runs in its own subprocess with a wall-clock limit, and the per-file output marks it as done:
rerunning skips finished files, a failed or stalled file is retried (optionally from a fallback
source, e.g. EOS instead of dCache).

A script that uses this module must accept  --one-file SOURCE --key KEY --out OUTPUT  and write
OUTPUT atomically (write to a temporary name, then rename).
"""

from __future__ import annotations

import pickle
import subprocess
import sys
import time
from pathlib import Path


def write_pickle(obj, out: Path) -> None:
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + ".tmp")
    with open(tmp, "wb") as fh:
        pickle.dump(obj, fh)
    tmp.rename(out)


def run_files(tasks, script: Path, out_dir: Path, workers: int, extra_args=(), stall_s: float = 2400,
              retries: int = 2, fallback=None, log=print, suffix: str = ".root", python=None):
    """Run `python script --one-file SOURCE --key KEY --out out_dir/KEY+suffix` for each (source, key).

    Returns a list of notes (retries, replaced sources). Raises if a file fails `retries + 1` times.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    queue = [(str(src), key, 0) for src, key in tasks if not (out_dir / f"{key}{suffix}").exists()]
    log(f"[batch] {len(tasks) - len(queue)} of {len(tasks)} already done, {len(queue)} to run "
        f"with {workers} workers")
    notes, running, t0, done, failed = [], {}, time.time(), 0, []
    while queue or running:
        while queue and len(running) < workers:
            src, key, attempt = queue.pop(0)
            cmd = [python or sys.executable, str(script), "--one-file", src, "--key", key,
                   "--out", str(out_dir / f"{key}{suffix}"), *extra_args]
            err = open(out_dir / f"{key}.stderr", "w")
            proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=err, text=True)
            running[proc] = (src, key, attempt, time.time(), err)
        time.sleep(1)
        for proc, (src, key, attempt, started, err) in list(running.items()):
            alive = proc.poll() is None
            stalled = alive and time.time() - started > stall_s
            if alive and not stalled:
                continue
            if stalled:
                proc.kill()
                proc.wait()
            err.close()
            del running[proc]
            errfile = out_dir / f"{key}.stderr"
            if not stalled and proc.returncode == 0:
                done += 1
                errfile.unlink(missing_ok=True)
                if done % 5 == 0 or not (queue or running):
                    log(f"        {done} done, {len(running)} running, {len(queue)} queued ({time.time()-t0:.0f}s)")
                continue
            if stalled:
                reason = f"no result after {stall_s:.0f}s"
            else:
                lines = errfile.read_text().strip().splitlines() if errfile.exists() else []
                reason = lines[-1][:240] if lines else f"exit code {proc.returncode}"
            if attempt >= retries:
                failed.append(f"{key}: {reason}")
                log(f"        FAILED {key}: {reason}")
                continue
            new_src = fallback(src, key) if fallback else None
            note = f"{key}: {reason}; retry {attempt + 1}" + (f" from {new_src}" if new_src else "")
            log("        " + note)
            notes.append(note)
            queue.append((str(new_src or src), key, attempt + 1))
    if failed:
        raise RuntimeError(f"{len(failed)} files failed:\n  " + "\n  ".join(failed))
    return notes
