#!/usr/bin/env python3

import argparse
import csv
import re
import time
from pathlib import Path
from urllib.parse import urljoin

import requests
from bs4 import BeautifulSoup

"""
USAGE:
    python3 01_scrape_genome_size.py --pages 10 --outdir ./reptile_data --sleep 1.0

DESCRIPTION:
    A specialized scraper designed to extract C-value (genome size) records 
    from the Animal Genome Size Database (genomesize.com). 

    The script targets reptile search results, crawls individual species detail pages, 
    and parses C-value measurements (in picograms) using regex and table extraction. 
    Outputs include a summary CSV and raw table rows for downstream data cleaning 
    and integration into genomic pipelines.

ARGUMENTS:
    -o, --outdir : Directory to save CSVs and cached HTML files (default: genomes/records/...)
    --pages      : Total number of results pages to iterate through (default: 5)
    --sleep      : Rate-limiting delay between GET requests (default: 0.75s)
"""

BASE = "https://www.genomesize.com"


def normalize_species_name(organism_name: str) -> str:
    organism_name = organism_name.strip()
    if "_" in organism_name:
        parts = [p for p in organism_name.split("_") if p]
    else:
        parts = organism_name.split()

    if len(parts) >= 2:
        return f"{parts[0]}_{parts[1]}"
    elif len(parts) == 1:
        return parts[0]
    return "unknown_species"


def clean_text(x: str) -> str:
    return re.sub(r"\s+", " ", x or "").strip()


def get_soup(session, url, out_html=None, sleep=0.5):
    print(f"[GET] {url}")
    r = session.get(url, timeout=30)
    r.raise_for_status()

    if out_html:
        Path(out_html).write_text(r.text, encoding="utf-8")

    time.sleep(sleep)
    return BeautifulSoup(r.text, "html.parser")


def extract_species_links(soup):
    links = []
    for a in soup.find_all("a", href=True):
        href = a["href"]
        if "result_species.php?id=" in href:
            url = urljoin(BASE, href)
            species_text = clean_text(a.get_text(" "))
            links.append((species_text, url))

    # de-duplicate while preserving order
    seen = set()
    unique = []
    for text, url in links:
        if url not in seen:
            seen.add(url)
            unique.append((text, url))

    return unique


def extract_all_table_rows(soup):
    rows = []
    for table_i, table in enumerate(soup.find_all("table"), start=1):
        for tr_i, tr in enumerate(table.find_all("tr"), start=1):
            cells = [clean_text(td.get_text(" ")) for td in tr.find_all(["th", "td"])]
            cells = [c for c in cells if c]
            if cells:
                rows.append({
                    "table_index": table_i,
                    "row_index": tr_i,
                    "cells_joined": " | ".join(cells),
                    "cells": cells,
                })
    return rows


def find_c_values_from_text(text):
    """
    Flexible C-value regex. This catches patterns like:
    C-value: 2.30
    1C: 2.30 pg
    2.30 pg
    """
    text = clean_text(text)
    vals = []

    patterns = [
        r"C[-\s]?value[^0-9]{0,20}([0-9]+(?:\.[0-9]+)?)",
        r"\b1C[^0-9]{0,20}([0-9]+(?:\.[0-9]+)?)",
        r"\b([0-9]+(?:\.[0-9]+)?)\s*pg\b",
    ]

    for pat in patterns:
        for m in re.finditer(pat, text, flags=re.IGNORECASE):
            vals.append(m.group(1))

    # de-duplicate
    out = []
    seen = set()
    for v in vals:
        if v not in seen:
            seen.add(v)
            out.append(v)
    return out


def parse_species_page(session, species_name_from_link, species_url, html_dir, sleep=0.5):
    species_id = species_url.split("id=")[-1]
    html_path = html_dir / f"species_{species_id}.html"
    soup = get_soup(session, species_url, out_html=html_path, sleep=sleep)

    page_text = clean_text(soup.get_text(" "))
    table_rows = extract_all_table_rows(soup)
    c_values = find_c_values_from_text(page_text)

    # Try to infer scientific name from page text or link text
    species_name = species_name_from_link
    if not species_name or species_name.lower() in {"details", "more"}:
        # crude fallback: use title/header if present
        h_candidates = [clean_text(h.get_text(" ")) for h in soup.find_all(["h1", "h2", "h3"])]
        h_candidates = [h for h in h_candidates if h]
        if h_candidates:
            species_name = h_candidates[0]

    return {
        "species_id": species_id,
        "species_url": species_url,
        "species_name_from_link": species_name_from_link,
        "species_name_inferred": species_name,
        "species_key": normalize_species_name(species_name),
        "c_values_pg_detected": ";".join(c_values),
        "n_c_values_detected": len(c_values),
        "page_text": page_text,
        "table_rows": table_rows,
    }


def main():
    ap = argparse.ArgumentParser(
        description="Download reptile C-value records from Animal Genome Size Database result pages."
    )
    ap.add_argument(
        "-o", "--outdir",
        default="genomes/records/animal_genome_size_database/reptiles",
        help="Output directory."
    )
    ap.add_argument(
        "--pages",
        type=int,
        default=5,
        help="Number of reptile result pages to download. Default: 5."
    )
    ap.add_argument(
        "--sleep",
        type=float,
        default=0.75,
        help="Seconds to wait between requests."
    )
    args = ap.parse_args()

    outdir = Path(args.outdir)
    html_dir = outdir / "html"
    html_dir.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    session.headers.update({
        "User-Agent": "Mozilla/5.0 academic personal-use scraper; contact: local research use"
    })

    all_species_links = []

    # Start with search URL to establish whatever state/cookies the site needs.
    search_url = f"{BASE}/search.php?display=100&search=type&value=Reptiles"
    get_soup(session, search_url, out_html=html_dir / "search_reptiles.html", sleep=args.sleep)

    for page in range(1, args.pages + 1):
        url = f"{BASE}/results.php?page={page}"
        soup = get_soup(session, url, out_html=html_dir / f"results_page_{page}.html", sleep=args.sleep)

        links = extract_species_links(soup)
        print(f"[INFO] page {page}: found {len(links)} species/detail links")
        all_species_links.extend(links)

    # de-duplicate species links
    seen = set()
    species_links = []
    for text, url in all_species_links:
        if url not in seen:
            seen.add(url)
            species_links.append((text, url))

    print(f"[INFO] total unique species/detail links: {len(species_links)}")

    # Save link index
    links_csv = outdir / "reptile_species_links.csv"
    with links_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["species_name_from_link", "species_url"])
        w.writeheader()
        for text, url in species_links:
            w.writerow({"species_name_from_link": text, "species_url": url})

    records = []
    raw_rows = []

    for i, (species_text, species_url) in enumerate(species_links, start=1):
        print(f"[SPECIES {i}/{len(species_links)}] {species_text} {species_url}")
        parsed = parse_species_page(session, species_text, species_url, html_dir, sleep=args.sleep)

        records.append({
            "species_id": parsed["species_id"],
            "species_url": parsed["species_url"],
            "species_name_from_link": parsed["species_name_from_link"],
            "species_name_inferred": parsed["species_name_inferred"],
            "species_key": parsed["species_key"],
            "c_values_pg_detected": parsed["c_values_pg_detected"],
            "n_c_values_detected": parsed["n_c_values_detected"],
        })

        for row in parsed["table_rows"]:
            raw_rows.append({
                "species_id": parsed["species_id"],
                "species_url": parsed["species_url"],
                "species_name_inferred": parsed["species_name_inferred"],
                "species_key": parsed["species_key"],
                "table_index": row["table_index"],
                "row_index": row["row_index"],
                "cells_joined": row["cells_joined"],
            })

    summary_csv = outdir / "reptile_c_values_detected_summary.csv"
    with summary_csv.open("w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "species_id",
            "species_url",
            "species_name_from_link",
            "species_name_inferred",
            "species_key",
            "c_values_pg_detected",
            "n_c_values_detected",
        ]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(records)

    raw_rows_csv = outdir / "reptile_species_pages_raw_table_rows.csv"
    with raw_rows_csv.open("w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "species_id",
            "species_url",
            "species_name_inferred",
            "species_key",
            "table_index",
            "row_index",
            "cells_joined",
        ]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(raw_rows)

    print("\n[DONE]")
    print(f"Saved links:   {links_csv}")
    print(f"Saved summary: {summary_csv}")
    print(f"Saved raw rows:{raw_rows_csv}")
    print(f"Saved HTML:    {html_dir}")


if __name__ == "__main__":
    main()
