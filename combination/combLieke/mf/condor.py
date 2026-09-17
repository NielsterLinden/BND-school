"""The fits as one HTCondor DAG (NIKHEF: schedd taai-007, Debian workers wrapped in the ATLAS almalinux9 image).

    lk_common ──┬── ch_<channel> (standalone fit + stat-only workspace, one node each) ── statonly
                ├── scan
                └── rank   (one job per nuisance parameter; the list is written by the common job)
    lk_<name>      (split, var_*, check_*: one node per likelihood, independent of the rest)
    final          (DAGMan FINAL node, runs whatever happened before: mergeranking + results --interim)

One node per likelihood and per channel, not one cluster: DAGMan removes every other job of a cluster as soon as one
of them fails. The ranking is one cluster of ~100 jobs, so `run.py rank` reports a failed refit without failing.

Every job is `condor/run_step.sh <combLieke> <run.py arguments>`. Submit files, parameter files, the DAG and all
HTCondor logs are generated into work/condor/ (git-ignored); only the two shell scripts in condor/ are sources.
"""

from __future__ import annotations

import subprocess

from .paths import HERE, WORK

IMAGE = "/cvmfs/atlas.cern.ch/repo/containers/fs/singularity/x86_64-almalinux9"
BATCH = "zcomb"

SUB = """executable = {here}/condor/run_step.sh
arguments  = {here} $(args)

output = {dir}/out/{node}_$(Cluster).$(Process).log
error  = {dir}/err/{node}_$(Cluster).$(Process).log
log    = {dir}/logs/{node}_$(Cluster).log

stream_output = True
stream_error  = True

accounting_group = atlas

request_cpus   = {cpus}
request_memory = {memory}

# express (10 min), short (4 h), medium (24 h), long (96 h)
+JobCategory = "{category}"

# the workers are Debian: the el9 LCG_110 view of setup.sh needs the almalinux9 image
+SingularityImage    = "{image}"
+SingularityBindCVMFS = True

JobBatchName = "{batch}_{node}"

queue args from {dir}/params_{node}.txt
"""


def write(likelihoods: dict, channels: list[str], submit: bool = False) -> None:
    cdir = WORK / "condor"
    for sub in ("out", "err", "logs"):
        (cdir / sub).mkdir(parents=True, exist_ok=True)
    if not (WORK / "common" / "multifit.config").exists():
        raise SystemExit("run `python run.py prepare` first")
    nodes = {
        **{f"lk_{n}": ([f"likelihood {n}"], "medium") for n in likelihoods},
        **{f"ch_{k}": ([f"channelfit {k}"], "short") for k in channels},
        "statonly": (["statonly"], "short"),
        "scan": (["scan"], "medium"),
        "rank": (None, "medium"),              # params_rank.txt is written by the common job (run.py write_rank_params)
        "final": (["final"], "short"),
    }
    for node, (params, category) in nodes.items():
        (cdir / f"{node}.sub").write_text(SUB.format(here=HERE, dir=cdir, node=node, cpus=4, memory="8GB" if node != "final" else "2GB",
                                                     category=category, image=IMAGE, batch=BATCH))
        if params is not None:
            (cdir / f"params_{node}.txt").write_text("".join(p + "\n" for p in params))
    (cdir / "params_rank.txt").unlink(missing_ok=True)
    for old in (WORK / "common" / "combination" / "Fits").glob("NPRanking*"):     # never merge the refits of an earlier fit
        old.unlink()
    dag = cdir / "comb.dag"
    for old in cdir.glob("comb.dag.*"):       # rescue files and DAGMan logs of an earlier submission
        old.unlink()
    dag.write_text("".join(f"JOB {n} {cdir}/{n}.sub\n" for n in nodes if n != "final")
                   + f"FINAL final {cdir}/final.sub\n"
                   + f"SCRIPT PRE rank {HERE}/condor/wait_for_file.sh {cdir}/params_rank.txt\n"
                   + f"PARENT lk_common CHILD scan rank {' '.join(f'ch_{k}' for k in channels)}\n"
                   + f"PARENT {' '.join(f'ch_{k}' for k in channels)} CHILD statonly\n")
    n_jobs = sum(len(p) for p, _ in nodes.values() if p)
    print(f"{dag}: {n_jobs} jobs + one ranking job per fitted parameter")
    if submit:
        subprocess.run(["condor_submit_dag", "-batch-name", BATCH, "-f", str(dag)], cwd=cdir, check=True)
        print(f"watch: /user/sjankovy/.claude/skills/nikhef-condor/condor_watch.sh --prefix {BATCH}")
