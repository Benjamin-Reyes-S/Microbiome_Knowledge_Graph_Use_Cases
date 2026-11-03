#!/usr/bin/env bash
# Run the Human-GEM taxonomy audit against one KG version.
# Usage: ./organizer.sh <neo4j_http_port>      e.g. ./organizer.sh 8004
# Results go to KG_port<port>/{data,results,figures}.
set -euo pipefail
cd "$(dirname "$0")"

if [[ $# -ne 1 || ! "$1" =~ ^[0-9]+$ ]]; then
    echo "Usage: $0 <neo4j_http_port>   (e.g. $0 8004)" >&2
    exit 1
fi

source .env/bin/activate
export NEO4J_URL="http://localhost:$1/db/neo4j/tx/commit"
echo "KG endpoint: $NEO4J_URL"

# 02 first: 01 looks up only the model's gene ids in the KG
for s in 02_extract_human_gem.py 01_extract_kg_human.py 03_match_model_to_kg.py \
         04_rescue_unmatched_reactions.py 05_xref_profile_unmatched.py; do
    echo "=== $s ==="
    python "$s"
done
