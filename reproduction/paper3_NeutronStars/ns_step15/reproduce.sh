#!/bin/sh
# Fresh-directory, sequential reproduction. No parallel scientific jobs.
set -eu
cd "$(dirname "$0")/../.."
output=${1:?Usage: reproduce.sh /absolute/new/output-directory}
case "$output" in /*) ;; *) echo 'Use an absolute output path' >&2; exit 2;; esac
[ ! -e "$output" ] || { echo 'Output already exists' >&2; exit 2; }
mkdir -p "$output/resource_evidence"
step=paper3_NeutronStars/ns_step15
launch() {
    stage=$1
    shift
    nice -n 19 python3 "$step/resource_control.py" --log "$output/resource_evidence/$stage.json" "$@"
}
launch probe --probe
launch scan "$step/scan.py" --output "$output/scan" --workers 1
launch refine "$step/refine.py" --scan "$output/scan" --output "$output/refinement" --workers 1
launch analyze "$step/analyze.py" --scan "$output/scan" --refinement "$output/refinement" --output "$output/analysis"
launch checks "$step/checks.py" --scan "$output/scan" --refinement "$output/refinement" --output "$output/checks" --workers 1
launch plot "$step/plot.py" --analysis "$output/analysis" --output "$output/figures"
