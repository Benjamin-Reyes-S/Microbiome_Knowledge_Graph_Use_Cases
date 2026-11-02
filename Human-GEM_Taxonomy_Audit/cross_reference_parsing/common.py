"""Shared helpers for the cross-reference parsing comparison (taxonomy of Human-GEM entities via public web APIs, no KG).

`fetch()` is the ONLY network entry point: GET/POST with per-host throttling, retries with back-off (UniProt answers 503 now
and then), and a persistent SQLite cache in data/api_cache.sqlite so every script can be re-run without repeating requests.
The parent audit folder is only READ (Human-GEM.xml, its match tables); everything is written inside this folder.
"""
import json, sqlite3, threading, time, urllib.error, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

BASE = Path(__file__).resolve().parent
PARENT = BASE.parent
DATA, RESULTS, FIGURES = BASE / "data", BASE / "results", BASE / "figures"
for p in (DATA, RESULTS, FIGURES):
    p.mkdir(exist_ok=True)
MODEL_XML = PARENT / "Human-GEM.xml"                 # read-only input
PARENT_RESULTS = PARENT / "results"                  # read-only (reaction_class, KG match tables)
HUMAN = 9606

# minimum seconds between two requests to the same host (politeness / documented rate limits) and parallel workers
MIN_GAP = {"rest.kegg.jp": 0.35, "rest.uniprot.org": 0.08, "www.rhea-db.org": 0.1, "bigg.ucsd.edu": 0.08,
           "rest.ensembl.org": 0.1, "reactome.org": 0.1, "rest.genenames.org": 0.15}
WORKERS = 6
_lock, _last = threading.Lock(), {}
_db_lock = threading.Lock()
_db = sqlite3.connect(DATA / "api_cache.sqlite", check_same_thread=False)
_db.execute("CREATE TABLE IF NOT EXISTS cache (k TEXT PRIMARY KEY, status INT, headers TEXT, body TEXT)")
_db.commit()
stats = {"cache_hit": 0, "request": 0, "failed": 0}


def fetch(url, data=None, headers=None, retries=6):
    """-> (status, headers dict (lower-case keys), body str). status 0 = gave up. 404 is a valid answer (id unknown to that DB)."""
    key = url + ("|" + data if data else "")
    with _db_lock:
        row = _db.execute("SELECT status,headers,body FROM cache WHERE k=?", (key,)).fetchone()
    if row:
        stats["cache_hit"] += 1
        return row[0], json.loads(row[1]), row[2]
    host = urllib.parse.urlsplit(url).netloc
    for attempt in range(retries):
        with _lock:                                   # one request slot per host every MIN_GAP seconds
            wait = _last.get(host, 0) + MIN_GAP.get(host, 0.1) - time.time()
            _last[host] = time.time() + max(wait, 0)
        if wait > 0:
            time.sleep(wait)
        req = urllib.request.Request(url, data.encode() if data else None,
                                     {"User-Agent": "humangem-xref-audit/1.0", "Accept": "*/*", **(headers or {})})
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                status, hd, body = r.status, {k.lower(): v for k, v in r.headers.items()}, r.read().decode("utf-8", "replace")
        except urllib.error.HTTPError as e:
            if e.code in (404, 400, 410):             # definitive "not found / bad id": cache it
                status, hd, body = e.code, {}, ""
            else:                                     # 429/5xx: back off and retry
                time.sleep(min(2 ** attempt, 30)); continue
        except Exception:
            time.sleep(min(2 ** attempt, 30)); continue
        stats["request"] += 1
        with _db_lock:
            _db.execute("INSERT OR REPLACE INTO cache VALUES (?,?,?,?)", (key, status, json.dumps(hd), body)); _db.commit()
        return status, hd, body
    stats["failed"] += 1
    return 0, {}, ""


def run_parallel(fn, items, workers=WORKERS, label=""):
    """Apply fn to every item with a thread pool (host throttling is inside fetch); progress printed every 10 %."""
    items, out, t0 = list(items), {}, time.time()
    step = max(1, len(items) // 10)
    with ThreadPoolExecutor(workers) as ex:
        for i, (it, res) in enumerate(zip(items, ex.map(fn, items)), 1):
            out[it] = res
            if i % step == 0:
                print(f"  {label} {i}/{len(items)} ({time.time() - t0:.0f}s, requests {stats['request']}, cache {stats['cache_hit']}, failed {stats['failed']})", flush=True)
    return out


def split(s):
    return [x for x in str(s).split(";") if x]


def log_calls(script, routes):
    """Append this script's API usage to results/api_call_log.csv (requests actually sent, cache hits, give-ups)."""
    import pandas as pd
    p = RESULTS / "api_call_log.csv"
    old = pd.read_csv(p) if p.exists() else pd.DataFrame()
    row = pd.DataFrame([dict(script=script, requests_sent=stats["request"], cache_hits=stats["cache_hit"], failed=stats["failed"],
                             routes=routes)])
    pd.concat([old[old.script != script] if len(old) else old, row]).to_csv(p, index=False)
