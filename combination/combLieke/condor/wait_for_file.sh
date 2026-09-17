#!/bin/bash
# DAGMan PRE script of the ranking node (runs on the schedd): the parameter list is written by the `common`
# job on a worker, so wait until the schedd's NFS view shows it (up to 10 min) before the node is submitted.
for _ in $(seq 120); do
    [ -s "$1" ] && exit 0
    sleep 5
done
echo "$1 did not appear" >&2
exit 1
