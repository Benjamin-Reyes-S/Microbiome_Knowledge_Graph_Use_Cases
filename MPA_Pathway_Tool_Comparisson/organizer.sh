#!/usr/bin/env bash
# Run the MPA_Pathway_Tool vs. KG comparison against one KG version.
# Usage: ./organizer.sh <neo4j_http_port>      e.g. ./organizer.sh 8034
# Inputs stay in data/; results go to KG_port<port>/{results,figures}.
set -euo pipefail
cd "$(dirname "$0")"

if [[ $# -ne 1 || ! "$1" =~ ^[0-9]+$ ]]; then
    echo "Usage: $0 <neo4j_http_port>   (e.g. $0 8034)" >&2
    exit 1
fi

source ../Human-GEM_Taxonomy_Audit/.env/bin/activate   # shared environment (not duplicated)
export NEO4J_URL="http://localhost:$1/db/neo4j/tx/commit"
echo "KG endpoint: $NEO4J_URL"

for s in 01_kg_module_definitions.py 02_calibrate_replication.py 03_kg_mapping.py \
         04_compare_and_plot.py 05_explain_R04347_difference.py; do
    echo "=== $s ==="
    python "$s"
done
