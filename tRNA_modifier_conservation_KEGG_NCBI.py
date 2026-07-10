# -*- coding: utf-8 -*-
"""Quantify conservation of selected tRNA-modification enzymes across bacteria.

For each bacterial phylum represented among KEGG bacterial genome records, the
script calculates the fraction of genome representatives with annotated KEGG
orthologues of MiaA, MnmE and MnmG. NCBI Taxonomy is used to assign phylum-level
classification from KEGG-linked taxonomy identifiers.

The analysis can retrieve KEGG records online or read previously downloaded
KEGG tables from a user-selected directory. Genome FASTA files are not required.

Requirements
------------
Python 3.9 or later and the following packages:
    requests, pandas, matplotlib, openpyxl, numpy

Outputs
-------
A user-selected output directory containing a timestamped results folder with
an Excel workbook, source tables, plotted values and publication-quality figure
files.
"""

# =============================================================================
# USER SETTINGS
# =============================================================================

DATA_SOURCE_MODE = "online"      # "online" or "manual"

TARGET_KOS = {
    "MiaA": "K00791",
    "MnmE": "K03650",
    "MnmG": "K03495",
}

OUTPUT_BASE_DIR = None  # Selected interactively at runtime.
MANUAL_INPUT_DIR = None   # Selected interactively when manual mode is used.
OUTPUT_FOLDER_PREFIX = "kegg_tRNA_modifier_conservation"
CREATE_TIMESTAMPED_OUTPUT_FOLDER = True
CACHE_FOLDER_NAME = "_cache"

COUNT_APPROXIMATE_SPECIES_INSTEAD_OF_KEGG_GENOMES = False

MIN_GENOMES_PER_PHYLUM_FOR_FIGURE = 20
MAX_PHYLUMS_IN_FIGURE = None
INCLUDE_ALL_THREE_IN_FIGURE = False

REQUEST_TIMEOUT_SECONDS = 120
POLITE_DELAY_SECONDS = 0.7
USE_CACHED_DOWNLOADS = True
NCBI_EUTILS_EMAIL = ""
NCBI_BATCH_SIZE = 100
NCBI_DELAY_SECONDS = 0.35

USE_MODERN_BACTERIAL_PHYLUM_NAMES = False

FIGURE_WIDTH_CM = 28
FIGURE_HEIGHT_CM = 21
DPI = 300
EXPORT_FORMATS = ["png", "svg", "pdf"]
SHOW_FIGURE = True

TITLE = "Conservation of MiaA, MnmE and MnmG across bacterial phyla"
X_AXIS_LABEL = "Bacterial genome representatives\nwith annotated ortholog (%)"
Y_AXIS_LABEL = ""

TITLE_FONT_SIZE = 13.2
AXIS_LABEL_FONT_SIZE = 17
TICK_FONT_SIZE = 15
LEGEND_FONT_SIZE = 17
N_LABEL_FONT_SIZE = 12

INTER_PHYLUM_SPACING = 2.75
BAR_GROUP_HEIGHT = 1.75
BAR_EDGE_WIDTH = 0.35

SHOW_N_BELOW_PHYLUM_LABEL = True
N_LABEL_PREFIX = "n="

TITLE_Y_POSITION = 0.988
LEGEND_Y_POSITION = 0.948
AXES_TOP_FOR_TITLE_AND_LEGEND = 0.935
LEFT_MARGIN = 0.245
RIGHT_MARGIN = 0.985
BOTTOM_MARGIN = 0.090
LEGEND_CENTER_MODE = "axes"  # "axes" or "figure"
Y_TICK_LABEL_PAD = 16

HIDE_TOP_AND_RIGHT_SPINES = True

GENE_COLORS = {
    "MiaA": "#E69F00",
    "MnmE": "#0072B2",
    "MnmG": "#009E73",
    "All three": "#444444",
}

# =============================================================================
# SCRIPT
# =============================================================================

import json
import tkinter as tk
from tkinter import filedialog
import re
import sys
import time
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path
from urllib.parse import urlencode

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests


KEGG_ENDPOINTS = {
    "genomes": "list/genome/bacteria",
    "taxonomy": "link/taxonomy/genome",
}

HEADERS = {
    "User-Agent": "tRNA-modifier-conservation/1.0",
    "Accept": "text/plain,*/*;q=0.8",
}

MODERN_PHYLUM_RENAME = {
    "Proteobacteria": "Pseudomonadota (Proteobacteria)",
    "Firmicutes": "Bacillota (Firmicutes)",
    "Actinobacteria": "Actinomycetota (Actinobacteria)",
    "Bacteroidetes": "Bacteroidota (Bacteroidetes)",
    "Cyanobacteria": "Cyanobacteriota (Cyanobacteria)",
    "Chloroflexi": "Chloroflexota (Chloroflexi)",
    "Spirochaetes": "Spirochaetota (Spirochaetes)",
    "Chlamydiae": "Chlamydiota (Chlamydiae)",
    "Nitrospirae": "Nitrospirota (Nitrospirae)",
    "Verrucomicrobia": "Verrucomicrobiota (Verrucomicrobia)",
    "Deinococcus-Thermus": "Deinococcota (Deinococcus-Thermus)",
    "Aquificae": "Aquificota (Aquificae)",
    "Thermotogae": "Thermotogota (Thermotogae)",
    "Fusobacteria": "Fusobacteriota (Fusobacteria)",
}


def select_directory(title: str) -> Path:
    """Open a directory picker and return the selected directory."""
    root = tk.Tk()
    root.withdraw()
    root.update()
    selected = filedialog.askdirectory(title=title)
    root.destroy()
    if not selected:
        raise RuntimeError(f"No directory was selected for: {title}")
    return Path(selected).expanduser().resolve()


def output_base_folder() -> Path:
    """Return the output base directory selected during program startup."""
    if OUTPUT_BASE_DIR is None:
        raise RuntimeError("The output directory has not been initialized.")
    base = Path(OUTPUT_BASE_DIR)
    base.mkdir(parents=True, exist_ok=True)
    return base


def cache_folder() -> Path:
    cache = output_base_folder() / CACHE_FOLDER_NAME
    cache.mkdir(parents=True, exist_ok=True)
    return cache


def make_output_folder() -> Path:
    base = output_base_folder()
    if CREATE_TIMESTAMPED_OUTPUT_FOLDER:
        stamp = datetime.now().strftime("%Y_%m_%d_%Hh%Mm%Ss")
        out = base / f"{OUTPUT_FOLDER_PREFIX}_{stamp}"
    else:
        out = base / OUTPUT_FOLDER_PREFIX
    out.mkdir(parents=True, exist_ok=True)
    return out


def manual_input_folder() -> Path:
    """Return the manually selected directory containing KEGG source tables."""
    if MANUAL_INPUT_DIR is None:
        raise RuntimeError("The manual KEGG input directory has not been initialized.")
    return Path(MANUAL_INPUT_DIR)


def cm_to_inch(x):
    return x / 2.54


def endpoint_to_file_name(endpoint: str) -> str:
    return endpoint.replace("/", "_") + ".tsv"


def endpoint_url(endpoint: str) -> str:
    return "https://rest.kegg.jp/" + endpoint.lstrip("/")


def read_manual_kegg_file(endpoint: str) -> str:
    path = manual_input_folder() / endpoint_to_file_name(endpoint)
    if not path.exists():
        raise FileNotFoundError(
            f"Manual KEGG file not found: {path}\n"
            "Expected files: list_genome_bacteria.tsv, link_taxonomy_genome.tsv, "
            "link_genes_K00791.tsv, link_genes_K03650.tsv, link_genes_K03495.tsv"
        )
    print(f"Reading manual KEGG file: {path}")
    return path.read_text(encoding="utf-8")


def download_text_from_kegg(endpoint: str) -> str:
    cache_file = cache_folder() / endpoint_to_file_name(endpoint)
    if USE_CACHED_DOWNLOADS and cache_file.exists() and cache_file.stat().st_size > 0:
        print(f"Using cached KEGG file: {cache_file}")
        return cache_file.read_text(encoding="utf-8")

    url = endpoint_url(endpoint)
    print(f"Downloading: {url}")
    response = requests.get(url, headers=HEADERS, timeout=REQUEST_TIMEOUT_SECONDS)
    if response.status_code != 200 or not response.text.strip():
        preview = response.text[:500].replace("\n", " ").replace("\r", " ")
        raise RuntimeError(
            f"KEGG download failed: {url}\n"
            f"HTTP status: {response.status_code}\n"
            f"Response preview: {preview!r}\n"
            "If online retrieval is unavailable, set DATA_SOURCE_MODE to 'manual'."
        )
    cache_file.write_text(response.text, encoding="utf-8")
    time.sleep(POLITE_DELAY_SECONDS)
    return response.text


def get_kegg_text(endpoint: str) -> str:
    mode = DATA_SOURCE_MODE.strip().lower()
    if mode == "manual":
        return read_manual_kegg_file(endpoint)
    if mode == "online":
        return download_text_from_kegg(endpoint)
    raise ValueError('DATA_SOURCE_MODE must be "online" or "manual".')


def copy_input_files_to_output(out: Path):
    endpoints = list(KEGG_ENDPOINTS.values()) + [f"link/genes/{ko}" for ko in TARGET_KOS.values()]
    for endpoint in endpoints:
        name = endpoint_to_file_name(endpoint)
        src = cache_folder() / name
        if DATA_SOURCE_MODE.strip().lower() == "manual":
            src = manual_input_folder() / name
        dst = out / name
        if src.exists() and not dst.exists():
            dst.write_text(src.read_text(encoding="utf-8"), encoding="utf-8")


def parse_bacterial_genome_list(text: str) -> pd.DataFrame:
    rows = []
    for line in text.splitlines():
        if not line.strip() or "\t" not in line:
            continue
        left, right = line.split("\t", 1)
        t_match = re.search(r"T\d{5}", left)
        if not t_match or ";" not in right:
            continue
        org_code, organism_name = right.split(";", 1)
        rows.append({
            "t_number": t_match.group(0),
            "org_code": org_code.strip(),
            "organism_name": organism_name.strip(),
        })
    df = pd.DataFrame(rows).drop_duplicates("t_number") if rows else pd.DataFrame()
    if df.empty:
        raise RuntimeError("Could not parse list_genome_bacteria.tsv")
    return df.reset_index(drop=True)


def parse_genome_taxonomy_links(text: str) -> pd.DataFrame:
    rows = []
    for line in text.splitlines():
        parts = line.strip().split("\t") if "\t" in line else line.strip().split()
        if len(parts) < 2:
            continue
        org_match = re.search(r"(?:gn:)?([A-Za-z0-9_]+)$", parts[0].strip())
        tax_match = re.search(r"(?:taxid|tax|taxonomy):(\d+)$", parts[1].strip())
        if org_match and tax_match:
            rows.append({"org_code": org_match.group(1), "taxonomy_taxid": tax_match.group(1)})
    df = pd.DataFrame(rows).drop_duplicates("org_code") if rows else pd.DataFrame()
    if df.empty:
        preview = text[:1000].replace("\n", "\\n")
        raise RuntimeError(f"Could not parse link_taxonomy_genome.tsv. Preview:\n{preview}")
    return df.reset_index(drop=True)


def parse_ko_gene_links(text: str) -> pd.DataFrame:
    rows = []
    for line in text.splitlines():
        parts = line.strip().split("\t") if "\t" in line else line.strip().split()
        if len(parts) < 2 or ":" not in parts[1]:
            continue
        gene_prefix, gene_id = parts[1].split(":", 1)
        rows.append({
            "ko_entry": parts[0].strip(),
            "gene_entry": parts[1].strip(),
            "gene_prefix": gene_prefix,
            "gene_id": gene_id,
        })
    return pd.DataFrame(rows)


def approx_species_name(name: str) -> str:
    name = re.sub(r"\s*\([^)]*\)", "", str(name)).strip()
    words = name.split()
    if len(words) >= 3 and words[0].lower() == "candidatus":
        return " ".join(words[:3])
    if len(words) >= 2:
        return " ".join(words[:2])
    return name or "Unknown species"


def load_cached_taxid_to_phylum() -> dict:
    path = cache_folder() / "ncbi_taxid_to_phylum.json"
    if USE_CACHED_DOWNLOADS and path.exists() and path.stat().st_size > 0:
        try:
            print(f"Using cached NCBI taxid-to-phylum map: {path}")
            return json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            return {}
    return {}


def save_cached_taxid_to_phylum(mapping: dict):
    path = cache_folder() / "ncbi_taxid_to_phylum.json"
    path.write_text(json.dumps(mapping, indent=2, sort_keys=True), encoding="utf-8")


def fetch_ncbi_taxonomy_batch(taxids) -> str:
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "taxonomy",
        "id": ",".join(taxids),
        "retmode": "xml",
        "tool": "kegg_tRNA_modifier_conservation",
    }
    if NCBI_EUTILS_EMAIL.strip():
        params["email"] = NCBI_EUTILS_EMAIL.strip()
    url = base + "?" + urlencode(params)
    response = requests.get(url, timeout=REQUEST_TIMEOUT_SECONDS)
    if response.status_code != 200 or not response.text.strip():
        preview = response.text[:500].replace("\n", " ").replace("\r", " ")
        raise RuntimeError(f"NCBI taxonomy download failed; status {response.status_code}; {preview!r}")
    return response.text


def parse_ncbi_taxonomy_phyla(xml_text: str) -> dict:
    root = ET.fromstring(xml_text)
    out = {}
    for taxon in root.findall(".//Taxon"):
        query_taxid = taxon.findtext("TaxId", default="").strip()
        if not query_taxid:
            continue
        phylum_name = None
        if taxon.findtext("Rank", default="").strip().lower() == "phylum":
            phylum_name = taxon.findtext("ScientificName", default="").strip()
        if not phylum_name:
            lineage_ex = taxon.find("LineageEx")
            if lineage_ex is not None:
                for lin_taxon in lineage_ex.findall("Taxon"):
                    if lin_taxon.findtext("Rank", default="").strip().lower() == "phylum":
                        phylum_name = lin_taxon.findtext("ScientificName", default="").strip()
                        break
        out[query_taxid] = phylum_name or "Unresolved phylum"
    return out


def get_taxid_to_phylum(taxids) -> dict:
    taxids = sorted({str(x) for x in taxids if pd.notna(x) and str(x).strip()})
    cached = load_cached_taxid_to_phylum()
    missing = [t for t in taxids if t not in cached]
    if missing:
        print(f"NCBI taxonomy phylum lookup: {len(missing):,} taxids not yet cached.")
    for i in range(0, len(missing), NCBI_BATCH_SIZE):
        batch = missing[i:i + NCBI_BATCH_SIZE]
        print(f"Fetching NCBI taxonomy {i + 1:,}-{min(i + len(batch), len(missing)):,} of {len(missing):,} missing taxids...")
        batch_mapping = parse_ncbi_taxonomy_phyla(fetch_ncbi_taxonomy_batch(batch))
        cached.update(batch_mapping)
        save_cached_taxid_to_phylum(cached)
        time.sleep(NCBI_DELAY_SECONDS)
    return {t: cached.get(t, "Unresolved phylum") for t in taxids}


def apply_phylum_name_style(name: str) -> str:
    if USE_MODERN_BACTERIAL_PHYLUM_NAMES:
        return MODERN_PHYLUM_RENAME.get(str(name), str(name))
    return str(name)


def build_presence_table(out: Path) -> pd.DataFrame:
    genome_text = get_kegg_text(KEGG_ENDPOINTS["genomes"])
    genomes = parse_bacterial_genome_list(genome_text)
    print(f"Bacterial KEGG genomes from list/genome/bacteria: {len(genomes):,}")

    taxonomy_text = get_kegg_text(KEGG_ENDPOINTS["taxonomy"])
    tax_links = parse_genome_taxonomy_links(taxonomy_text)
    print(f"Genome-to-taxonomy links parsed: {len(tax_links):,}")

    genomes = genomes.merge(tax_links, on="org_code", how="left")
    missing_tax = int(genomes["taxonomy_taxid"].isna().sum())
    if missing_tax:
        print(f"WARNING: {missing_tax:,} bacterial genomes lack taxonomy taxid links.")

    taxid_to_phylum = get_taxid_to_phylum(genomes["taxonomy_taxid"].dropna().astype(str).unique())
    genomes["phylum"] = genomes["taxonomy_taxid"].astype(str).map(taxid_to_phylum)
    genomes.loc[genomes["taxonomy_taxid"].isna(), "phylum"] = "Unresolved phylum"
    genomes["phylum"] = genomes["phylum"].fillna("Unresolved phylum").apply(apply_phylum_name_style)
    genomes["approx_species"] = genomes["organism_name"].apply(approx_species_name)

    all_bacterial_org_codes = set(genomes["org_code"])
    ko_link_tables = []
    for gene, ko in TARGET_KOS.items():
        link_text = get_kegg_text(f"link/genes/{ko}")
        links = parse_ko_gene_links(link_text)
        if links.empty:
            print(f"WARNING: no KEGG gene links parsed for {gene} ({ko}).")
            positive_org_codes = set()
        else:
            links["target_gene"] = gene
            links["target_ko"] = ko
            ko_link_tables.append(links)
            positive_org_codes = set(links["gene_prefix"]).intersection(all_bacterial_org_codes)
        genomes[f"{gene}_present"] = genomes["org_code"].isin(positive_org_codes)
        print(f"{gene} ({ko}): {int(genomes[f'{gene}_present'].sum()):,} positive bacterial KEGG genomes")

    if ko_link_tables:
        pd.concat(ko_link_tables, ignore_index=True).to_csv(out / "raw_kegg_ko_gene_links.tsv", sep="\t", index=False)
    copy_input_files_to_output(out)

    if COUNT_APPROXIMATE_SPECIES_INSTEAD_OF_KEGG_GENOMES:
        agg = {
            "phylum": "first",
            "taxonomy_taxid": "first",
            "organism_name": lambda x: "; ".join(sorted(set(map(str, x)))[:5]),
            "org_code": lambda x: ";".join(sorted(set(map(str, x)))),
            "t_number": lambda x: ";".join(sorted(set(map(str, x)))),
        }
        for gene in TARGET_KOS:
            agg[f"{gene}_present"] = "max"
        genomes = genomes.groupby("approx_species", as_index=False).agg(agg)
        genomes = genomes.rename(columns={"approx_species": "unit_id"})
        counting_unit = "approximate species"
    else:
        genomes["unit_id"] = genomes["t_number"]
        counting_unit = "KEGG bacterial genome"

    presence_cols = [f"{g}_present" for g in TARGET_KOS]
    genomes["All_three_present"] = genomes[presence_cols].all(axis=1)
    genomes["counting_unit"] = counting_unit

    first_cols = ["unit_id", "t_number", "org_code", "organism_name", "phylum", "taxonomy_taxid", "counting_unit"]
    return genomes[first_cols + [c for c in genomes.columns if c not in first_cols]]


def summarize_by_phylum(presence: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for phylum, sub in presence.groupby("phylum", dropna=False):
        row = {
            "phylum": phylum,
            "taxonomy_taxids": ";".join(sorted(set(map(str, sub["taxonomy_taxid"].dropna())))),
            "n_units": len(sub),
        }
        for gene in TARGET_KOS:
            n = int(sub[f"{gene}_present"].sum())
            row[f"{gene}_n_present"] = n
            row[f"{gene}_percent"] = 100 * n / len(sub)
        n_all = int(sub["All_three_present"].sum())
        row["All_three_n_present"] = n_all
        row["All_three_percent"] = 100 * n_all / len(sub)
        rows.append(row)
    return pd.DataFrame(rows).sort_values("n_units", ascending=False).reset_index(drop=True)


def export_excel(out: Path, presence: pd.DataFrame, summary: pd.DataFrame):
    excel_path = out / "kegg_tRNA_modifier_conservation_results.xlsx"
    settings = pd.DataFrame([
        ["data_source_mode", DATA_SOURCE_MODE],
        ["output_base_dir", str(output_base_folder())],
        ["cache_folder", str(cache_folder())],
        ["counting_unit", "approximate species" if COUNT_APPROXIMATE_SPECIES_INSTEAD_OF_KEGG_GENOMES else "KEGG bacterial genome representative"],
        ["genome_endpoint", endpoint_url(KEGG_ENDPOINTS["genomes"])],
        ["taxonomy_endpoint", endpoint_url(KEGG_ENDPOINTS["taxonomy"])],
        ["taxonomy_assignment", "NCBI Taxonomy phylum from KEGG organism-code taxonomy taxid"],
        ["MiaA", TARGET_KOS["MiaA"]],
        ["MnmE", TARGET_KOS["MnmE"]],
        ["MnmG", TARGET_KOS["MnmG"]],
        ["min_n_per_phylum_for_figure", MIN_GENOMES_PER_PHYLUM_FOR_FIGURE],
        ["max_phyla_in_figure", MAX_PHYLUMS_IN_FIGURE],
        ["include_all_three_in_figure", INCLUDE_ALL_THREE_IN_FIGURE],
        ["use_modern_bacterial_phylum_names", USE_MODERN_BACTERIAL_PHYLUM_NAMES],
        ["show_n_below_phylum_label", SHOW_N_BELOW_PHYLUM_LABEL],
        ["inter_phylum_spacing", INTER_PHYLUM_SPACING],
        ["bar_group_height", BAR_GROUP_HEIGHT],
        ["legend_center_mode", LEGEND_CENTER_MODE],
        ["y_tick_label_pad", Y_TICK_LABEL_PAD],
    ], columns=["setting", "value"])

    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="phylum_summary", index=False)
        presence.to_excel(writer, sheet_name="genome_presence", index=False)
        settings.to_excel(writer, sheet_name="settings", index=False)
    print(f"Saved: {excel_path}")


def plot_summary(out: Path, summary: pd.DataFrame):
    plot_df = summary[summary["n_units"] >= MIN_GENOMES_PER_PHYLUM_FOR_FIGURE].copy()
    if plot_df.empty:
        raise RuntimeError("No phyla pass MIN_GENOMES_PER_PHYLUM_FOR_FIGURE. Lower this setting.")
    if MAX_PHYLUMS_IN_FIGURE is not None:
        plot_df = plot_df.nlargest(int(MAX_PHYLUMS_IN_FIGURE), "n_units")
    plot_df = plot_df.sort_values("n_units", ascending=True).reset_index(drop=True)

    genes = list(TARGET_KOS.keys())
    if INCLUDE_ALL_THREE_IN_FIGURE:
        genes += ["All three"]

    y = np.arange(len(plot_df)) * INTER_PHYLUM_SPACING
    bar_h = BAR_GROUP_HEIGHT / len(genes)
    offsets = (np.arange(len(genes)) - (len(genes) - 1) / 2) * bar_h

    fig_h_cm = max(FIGURE_HEIGHT_CM, 0.62 * len(plot_df) * INTER_PHYLUM_SPACING + 4.5)
    fig, ax = plt.subplots(figsize=(cm_to_inch(FIGURE_WIDTH_CM), cm_to_inch(fig_h_cm)), dpi=DPI)

    for i, gene in enumerate(genes):
        col = "All_three_percent" if gene == "All three" else f"{gene}_percent"
        ax.barh(
            y + offsets[i],
            plot_df[col],
            height=bar_h * 0.92,
            label=gene,
            color=GENE_COLORS.get(gene),
            edgecolor="black",
            linewidth=BAR_EDGE_WIDTH,
        )

    if SHOW_N_BELOW_PHYLUM_LABEL:
        ylabels = [f"{p}\n({N_LABEL_PREFIX}{int(n)})" for p, n in zip(plot_df["phylum"], plot_df["n_units"])]
    else:
        ylabels = list(plot_df["phylum"])

    ax.set_yticks(y)
    ax.set_yticklabels(ylabels, fontsize=TICK_FONT_SIZE, linespacing=1.15)
    ax.tick_params(axis="x", labelsize=TICK_FONT_SIZE)
    ax.tick_params(axis="y", length=0, pad=Y_TICK_LABEL_PAD)
    ax.set_xlabel(X_AXIS_LABEL, fontsize=AXIS_LABEL_FONT_SIZE, labelpad=8)
    ax.set_ylabel(Y_AXIS_LABEL, fontsize=AXIS_LABEL_FONT_SIZE if Y_AXIS_LABEL else None)
    ax.set_xlim(0, 104)
    ax.grid(axis="x", linestyle=":", linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)

    if HIDE_TOP_AND_RIGHT_SPINES:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    handles, labels = ax.get_legend_handles_labels()
    fig.suptitle(TITLE, fontsize=TITLE_FONT_SIZE, y=TITLE_Y_POSITION)

    if LEGEND_CENTER_MODE.lower() == "axes":
        legend_x = (LEFT_MARGIN + RIGHT_MARGIN) / 2
    else:
        legend_x = 0.5

    fig.legend(
        handles,
        labels,
        frameon=False,
        fontsize=LEGEND_FONT_SIZE,
        loc="upper center",
        bbox_to_anchor=(legend_x, LEGEND_Y_POSITION),
        ncol=len(genes),
        handlelength=1.4,
        handletextpad=0.5,
        columnspacing=1.2,
    )

    fig.tight_layout(rect=[LEFT_MARGIN, BOTTOM_MARGIN, RIGHT_MARGIN, AXES_TOP_FOR_TITLE_AND_LEGEND])

    plotted_values_path = out / "figure_conservation_by_phylum_plotted_values.csv"
    plot_df.to_csv(plotted_values_path, index=False)
    print(f"Saved: {plotted_values_path}")

    for ext in EXPORT_FORMATS:
        path = out / f"figure_conservation_by_phylum.{ext}"
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
        print(f"Saved: {path}")

    if SHOW_FIGURE:
        plt.show()
    else:
        plt.close(fig)


def main():
    global OUTPUT_BASE_DIR, MANUAL_INPUT_DIR
    OUTPUT_BASE_DIR = select_directory("Select output directory")
    if DATA_SOURCE_MODE.strip().lower() == "manual":
        MANUAL_INPUT_DIR = select_directory("Select directory containing KEGG source tables")
    out = make_output_folder()
    print(f"Output folder: {out}")
    print(f"Cache folder: {cache_folder()}")
    print(f"Python: {sys.version.split()[0]}")
    print(f"Data source mode: {DATA_SOURCE_MODE}")

    presence = build_presence_table(out)
    summary = summarize_by_phylum(presence)
    export_excel(out, presence, summary)
    plot_summary(out, summary)

    print("\nDone.")
    print("Detailed values are available in the phylum_summary and genome_presence worksheets.")


if __name__ == "__main__":
    main()