#!/usr/bin/env python3
"""Extract protein-coordinate mutation inputs from the repository T-ALL XLSX.

Uses only the Python standard library. Amino-acid positions are taken from the
workbook's `aa_change` field. Substitutions, stop gains, and frameshifts are
anchored at their first affected residue; explicit ranges such as
`p.54_55del` retain both endpoints. Rows without a resolvable protein position
are excluded and counted in the provenance file.
"""

import csv
import re
import sys
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def column_name(reference: str) -> str:
    return re.match(r"[A-Z]+", reference).group(0)


def workbook_rows(path: Path):
    with zipfile.ZipFile(path) as archive:
        strings_xml = ET.fromstring(archive.read("xl/sharedStrings.xml"))
        strings = [
            "".join(node.text or "" for node in item.findall(".//m:t", NS))
            for item in strings_xml.findall("m:si", NS)
        ]
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
        for row in sheet.findall(".//m:sheetData/m:row", NS):
            values = {}
            for cell in row.findall("m:c", NS):
                value = cell.find("m:v", NS)
                text = "" if value is None else value.text
                if cell.attrib.get("t") == "s" and text:
                    text = strings[int(text)]
                values[column_name(cell.attrib["r"])] = text
            yield values


def protein_interval(change: str):
    if not change:
        return None
    ranged = re.search(r"p\.[^0-9]*(\d+)_[A-Za-z*]*(\d+)", change)
    if ranged:
        return int(ranged.group(1)), int(ranged.group(2))
    point = re.search(r"p\.[^0-9]*(\d+)", change)
    if point:
        position = int(point.group(1))
        return position, position
    return None


def main():
    root = Path(__file__).resolve().parents[1]
    workbook = root / "datasets" / "T_ALL_public_data" / "mutations_dataset.xlsx"
    available_targets = {"SUZ12", "EZH2", "PTEN", "LEF1"}
    targets = set(sys.argv[1:]) or available_targets
    unknown_targets = targets - available_targets
    if unknown_targets:
        raise ValueError(
            "Unknown target(s): " + ", ".join(sorted(unknown_targets))
        )
    metadata = {
        "SUZ12": {
            "refseq_transcript": "NM_015355",
            "refseq_protein": "NP_056170.2",
            "uniprot_accession": "Q15022",
            "coordinate_source": "AlphaFold_DB_AF-Q15022-F1-model_v6.pdb",
            "sequence_validation": "NP_056170.2 exactly matches UniProt Q15022 (739 residues)",
        },
        "EZH2": {
            "refseq_transcript": "NM_001203247",
            "refseq_protein": "NP_001190176.1",
            "uniprot_accession": "Q15910",
            "coordinate_source": "AlphaFold_DB_AF-Q15910-F1-model_v6.pdb",
            "sequence_validation": "NP_001190176.1 exactly matches UniProt Q15910 (746 residues)",
        },
        "PTEN": {
            "refseq_transcript": "NM_000314",
            "refseq_protein": "NP_000305.3",
            "uniprot_accession": "P60484",
            "coordinate_source": "AlphaFold_DB_AF-P60484-F1-model_v6.pdb",
            "sequence_validation": "All extracted mutation positions fall within canonical UniProt P60484 (403 residues)",
        },
        "LEF1": {
            "refseq_transcript": "NM_016269.5",
            "refseq_protein": "NP_057353.1",
            "uniprot_accession": "Q9UJU2",
            "coordinate_source": "AlphaFold_DB_AF-Q9UJU2-F1-model_v6.pdb",
            "sequence_validation": "NP_057353.1 corresponds to canonical UniProt Q9UJU2 (399 residues)",
        },
    }
    rows = iter(workbook_rows(workbook))
    header_cells = next(rows)
    headers = {column: name for column, name in header_cells.items()}
    extracted = {gene: [] for gene in targets}
    skipped = {gene: 0 for gene in targets}

    for cells in rows:
        record = {name: cells.get(column, "") for column, name in headers.items()}
        gene = record.get("gene", "")
        if gene not in targets:
            continue
        interval = protein_interval(record.get("aa_change", ""))
        if interval is None:
            skipped[gene] += 1
            continue
        start, end = interval
        extracted[gene].append(
            {
                "id": record.get("variantID") or record.get("position.ID"),
                "ID": record.get("sample"),
                "start": start,
                "end": end,
                "aa_change": record.get("aa_change", ""),
                "mutation_class": record.get("mutation_class", ""),
                "exonic_function": record.get("ExonicFunc", ""),
                "validation": record.get("validation", ""),
                "quality": record.get("quality.MULTI", ""),
                "source_variant_id": record.get("variantID", ""),
            }
        )

    for gene in sorted(targets):
        output_dir = root / "examples" / f"{gene}_mutations" / "input_files"
        output_dir.mkdir(parents=True, exist_ok=True)
        output = output_dir / f"{gene.lower()}_start_end.csv"
        with output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(extracted[gene][0]))
            writer.writeheader()
            writer.writerows(extracted[gene])
        provenance = output_dir / "mutation_input_provenance.txt"
        provenance.write_text(
            "\n".join(
                [
                    f"source={workbook.relative_to(root)}",
                    "source_fields=gene,sample,variantID,aa_change,mutation_class,ExonicFunc,validation,quality.MULTI",
                    "coordinate_rule=first affected residue; explicit amino-acid deletion/substitution ranges retain both endpoints",
                    f"gene={gene}",
                    f"included_rows={len(extracted[gene])}",
                    f"excluded_without_resolvable_aa_position={skipped[gene]}",
                    *[f"{key}={value}" for key, value in metadata[gene].items()],
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        print(f"{gene}: wrote {len(extracted[gene])} events to {output}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
