## Conservation of tRNA-modification enzymes across bacterial phyla using KEGG and NCBI Taxonomy

### Overview

`tRNA_modifier_conservation_KEGG_NCBI.py` quantifies the conservation of selected tRNA-modification enzymes across bacterial phyla.

The script uses KEGG bacterial genome representatives and KEGG Orthology assignments to determine the fraction of genomes containing each target enzyme. KEGG organism codes are linked to NCBI Taxonomy identifiers, which are then used to assign each genome to a bacterial phylum.

The default analysis evaluates the conservation of MiaA, MnmE and MnmG, but the target KEGG Orthology identifiers can be modified in the configuration section.

### Data sources

The script retrieves data from:

* KEGG bacterial genome listings;
* KEGG genome-to-taxonomy links;
* KEGG Orthology-to-gene links; and
* NCBI Taxonomy E-utilities.

No genome FASTA files are required.

### Requirements

The script requires Python 3.9 or later and the following packages:

```bash
pip install pandas numpy matplotlib requests openpyxl
```

A graphical desktop environment is required for directory-selection dialogs.

### Usage

Run the script from a Python environment or command line:

```bash
python tRNA_modifier_conservation_KEGG_NCBI.py
```

At startup, the script prompts the user to select an output directory. A timestamped analysis folder is then created inside the selected directory.

The script supports two data-source modes:

```python
DATA_SOURCE_MODE = "online"
```

In online mode, KEGG and NCBI data are downloaded automatically and cached locally.

```python
DATA_SOURCE_MODE = "manual"
```

In manual mode, the user is prompted to select a folder containing previously downloaded KEGG tables.

### Required files for manual mode

The selected manual-input folder must contain:

```text
list_genome_bacteria.tsv
link_taxonomy_genome.tsv
link_genes_K00791.tsv
link_genes_K03650.tsv
link_genes_K03495.tsv
```

The required filenames change accordingly when different KEGG Orthology identifiers are specified.

### Default targets

The default target enzymes are:

* MiaA — KEGG Orthology `K00791`
* MnmE — KEGG Orthology `K03650`
* MnmG — KEGG Orthology `K03495`

The target dictionary can be extended to include additional enzymes.

### Counting units

By default, the analysis uses individual KEGG bacterial genome representatives as counting units.

Alternatively:

```python
COUNT_APPROXIMATE_SPECIES_INSTEAD_OF_KEGG_GENOMES = True
```

can be used to collapse genome representatives into approximate species-level units based on organism names. This option reduces overrepresentation of species with multiple KEGG genomes, but it is based on name parsing rather than formal taxonomic clustering.

### Phylum-level analysis

For each bacterial phylum, the script calculates:

* total number of analyzed genome representatives or approximate species;
* number and percentage containing each target KEGG Orthology identifier; and
* number and percentage containing all target enzymes.

Phyla represented by fewer than the configured minimum number of genomes are excluded from the figure but remain available in the exported tables.

### Main configuration options

The following parameters can be adjusted in the configuration section:

* target KEGG Orthology identifiers;
* minimum number of genomes required for figure inclusion;
* maximum number of phyla shown;
* inclusion of an “all targets present” category;
* use of modern or historical bacterial phylum names;
* counting by KEGG genome or approximate species;
* figure dimensions and font sizes;
* export formats;
* bar spacing, legend position and plot margins;
* use of cached downloads; and
* NCBI request batch size and delay.

Providing an email address in `NCBI_EUTILS_EMAIL` is recommended for reproducible use of NCBI E-utilities.

### Output files

The script generates:

* `kegg_tRNA_modifier_conservation_results.xlsx`
  Excel workbook containing the phylum-level summary, genome-level presence table and analysis settings.

* `raw_kegg_ko_gene_links.tsv`
  Combined KEGG Orthology-to-gene links used in the analysis.

* `figure_conservation_by_phylum_plotted_values.csv`
  Exact values included in the plotted figure.

* `figure_conservation_by_phylum.png`

* `figure_conservation_by_phylum.svg`

* `figure_conservation_by_phylum.pdf`

* copies of the KEGG input tables used for the analysis; and

* `_cache/` containing KEGG downloads and the cached NCBI taxid-to-phylum mapping.

### Interpretation and limitations

The results represent the proportion of KEGG genome representatives with an annotated gene assigned to the specified KEGG Orthology group. Absence of a KEGG Orthology assignment does not necessarily demonstrate biological absence of the corresponding enzyme.

Low apparent conservation may result from:

* incomplete genome annotation;
* incomplete or lineage-specific KEGG Orthology assignment;
* highly divergent homologues;
* alternative enzymes performing the same biochemical function; or
* unresolved taxonomy mappings.

Unexpected lineage-specific results should therefore be validated using complementary annotation resources or sequence-based homology searches.
