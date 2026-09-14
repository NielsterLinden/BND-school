"""Running one per-file function over many files, robustly and resumably.

Steps 1 and 2 read the local skim with a multiprocessing pool, which is fine.
The cross-checks read the unskimmed NanoAOD from dCache and the simulation from
CERN EOS, and both can misbehave: dCache NFS reads occasionally stall in
uninterruptible I/O under many parallel readers, and remote reads can fail
transiently. A stalled pool worker cannot be killed, so here every file runs in
its own subprocess with a wall-clock limit, and each result is written to its
own pickle -- rerunning skips files that are already done.

A script that uses this module must accept

    --one-file SOURCE --key KEY --out PART.pkl

and write its per-file result to PART.pkl with `write_part`.
"""

from __future__ import annotations

import pickle
import subprocess
import sys
import time
import zlib
from pathlib import Path

LCG_SETUP = "/cvmfs/sft.cern.ch/lcg/views/LCG_110/x86_64-el9-gcc14-opt/setup.sh"


# --------------------------------------------------------------------------
# Per-file results
# --------------------------------------------------------------------------
def write_part(obj, out: Path) -> None:
    """Atomic pickle write, so a killed job never leaves a half-written part."""
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".tmp")
    with open(tmp, "wb") as fh:
        pickle.dump(obj, fh)
    tmp.rename(out)


def merge(a, b):
    """Add b into a: numbers, arrays and hist objects add, lists extend, dicts recurse."""
    for k, v in b.items():
        if k not in a:
            a[k] = v
        elif isinstance(v, list):
            a[k] = a[k] + v
        elif isinstance(v, dict):
            merge(a[k], v)
        else:
            a[k] = a[k] + v
    return a


def load_parts(parts_dir: Path, start: dict = None) -> dict:
    total = start if start is not None else {}
    for part in sorted(Path(parts_dir).glob("*.pkl")):
        with open(part, "rb") as fh:
            total = merge(total, pickle.load(fh))
    return total


# --------------------------------------------------------------------------
# CERN Open Data catalogue and verified local copies
# --------------------------------------------------------------------------
def read_catalogue(tsv: Path) -> dict:
    """{file name: (xrootd URL, size, checksum)} from a
    `cernopendata-client get-file-locations --protocol xrootd --verbose` listing."""
    cat = {}
    for line in Path(tsv).read_text().splitlines():
        if line.strip():
            url, size, checksum = line.split("\t")
            cat[Path(url).name] = (url, int(size), checksum)
    return cat


def adler32(path: Path) -> str:
    value = 1
    with open(path, "rb") as fh:
        while True:
            block = fh.read(64 << 20)
            if not block:
                break
            value = zlib.adler32(block, value)
    return f"adler32:{value:08x}"


def cached_copy(name: str, subdir: str, catalogue: dict, cache_dir: Path, fetch: bool):
    """Path of a verified EOS copy in `cache_dir/subdir`, fetching it with xrdcp
    (from the LCG view) when `fetch` is set. None if absent and not fetched."""
    if name not in catalogue:
        return None
    url, size, checksum = catalogue[name]
    dest = Path(cache_dir) / subdir / name
    if dest.exists() and dest.stat().st_size == size:
        return dest
    if not fetch:
        return None
    dest.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["bash", "-c", f"source {LCG_SETUP} && xrdcp -f -s {url} {dest}"], check=True)
    if dest.stat().st_size != size or adler32(dest) != checksum:
        dest.unlink()
        raise IOError(f"copy of {name} does not match the CERN catalogue")
    return dest


# --------------------------------------------------------------------------
# The runner
# --------------------------------------------------------------------------
def run_files(tasks, script: Path, parts_dir: Path, workers: int, extra_args=(),
              stall_s: float = 1200, retries: int = 1, fallback=None, log=print,
              suffix: str = ".pkl"):
    """Run `script --one-file SOURCE --key KEY --out PART` for each (source, key).

    A task that fails or exceeds `stall_s` is retried `retries` times; before a
    retry, `fallback(source, key)` may return a replacement source (for example
    a verified local copy). Returns the list of notes about replaced sources.
    `suffix` is the extension of the per-file output that marks a task as done
    (".pkl" for the review scripts, ".root" for the v2 skims).
    """
    parts_dir = Path(parts_dir)
    parts_dir.mkdir(parents=True, exist_ok=True)
    queue = [(str(src), key, 0) for src, key in tasks
             if not (parts_dir / f"{key}{suffix}").exists()]
    log(f"[batch] {len(tasks) - len(queue)} of {len(tasks)} files already done, "
        f"{len(queue)} to run with {workers} workers")

    notes, running, t0, done = [], {}, time.time(), 0
    while queue or running:
        while queue and len(running) < workers:
            src, key, attempt = queue.pop(0)
            cmd = [sys.executable, str(script), "--one-file", src, "--key", key,
                   "--out", str(parts_dir / f"{key}{suffix}"), *extra_args]
            proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                                    text=True)
            running[proc] = (src, key, attempt, time.time())
        time.sleep(1)
        for proc, (src, key, attempt, started) in list(running.items()):
            still_running = proc.poll() is None
            stalled = still_running and time.time() - started > stall_s
            if still_running and not stalled:
                continue
            if stalled:
                proc.kill()
            del running[proc]
            if not stalled and proc.returncode == 0:
                done += 1
                if done % 5 == 0 or not (queue or running):
                    log(f"        {done} done, {len(running)} running, {len(queue)} queued"
                        f"  ({time.time()-t0:.0f}s)")
                continue
            if stalled:
                reason = f"no result after {stall_s:.0f}s"
            else:
                lines = (proc.stderr.read() or "").strip().splitlines()
                reason = lines[-1][:200] if lines else f"exit code {proc.returncode}"
            if attempt >= retries:
                raise RuntimeError(f"{key}: failed {attempt + 1} times, last: {reason}")
            new_src = fallback(src, key) if fallback else None
            note = f"{key}: {reason}; retrying" + (f" from {new_src}" if new_src else "")
            log("        " + note)
            notes.append(note)
            queue.append((str(new_src or src), key, attempt + 1))
    return notes
