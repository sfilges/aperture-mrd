"""CLI entry points for aperture-snv.

Each command corresponds to a step in the SNV integration pipeline and is
invoked by the Nextflow snv_integration.nf subworkflow.

Commands:
    aperture-snv extract  — Extract candidates from plasma CRAM at compendium sites
    aperture-snv filter   — Apply trained SVM to filter artifact reads
    aperture-snv score    — Compute tumor fraction and detection z-score
    aperture-snv noise    — Build noise profile from control plasma samples
"""

from __future__ import annotations

import csv
import json
from typing import TYPE_CHECKING, Annotated

import typer
from loguru import logger
from rich.console import Console
from rich.table import Table

if TYPE_CHECKING:
    from pathlib import Path

app = typer.Typer(
    name="aperture-snv",
    help="SNV candidate extraction, error suppression, and tumor fraction scoring.",
    no_args_is_help=True,
)
console = Console()


@app.command()
def extract(
    cram: Annotated[Path, typer.Option(help="Plasma CRAM/BAM file path")],
    compendium: Annotated[Path, typer.Option(help="Patient SNV compendium VCF")],
    reference: Annotated[Path, typer.Option(help="Reference FASTA path")],
    output: Annotated[Path, typer.Option(help="Output TSV file path")],
    min_mapq: Annotated[int, typer.Option(help="Min mapping quality")] = 0,
    min_bq: Annotated[int, typer.Option(help="Min base quality")] = 0,
) -> None:
    """Extract candidate reads from plasma CRAM at compendium loci."""
    from aperture_snv.extract import extract_all_candidates, resolve_paired_end_concordance

    logger.info("Extracting candidates from {} at {} compendium sites", cram, compendium)

    candidates = extract_all_candidates(
        cram_path=cram,
        compendium_path=compendium,
        reference_path=reference,
        min_mapq=min_mapq,
        min_bq=min_bq,
    )
    candidates = resolve_paired_end_concordance(candidates)

    with open(output, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "chrom", "pos", "ref", "alt", "read_name", "read_base",
            "vbq", "mrbq", "pir", "mapq", "is_read1",
            "is_concordant", "is_proper_pair", "insert_size",
        ])
        for _site, reads in candidates.items():
            for r in reads:
                writer.writerow([
                    r.chrom, r.pos, r.ref, r.alt, r.read_name, r.read_base,
                    r.vbq, f"{r.mrbq:.2f}", f"{r.pir:.4f}", r.mapq,
                    int(r.is_read1),
                    "" if r.is_concordant is None else int(r.is_concordant),
                    int(r.is_proper_pair), r.insert_size,
                ])

    n_sites = len(candidates)
    n_reads = sum(len(reads) for reads in candidates.values())
    logger.success("Extracted {} reads at {} compendium sites → {}", n_reads, n_sites, output)


@app.command()
def filter(
    input: Annotated[Path, typer.Option(help="Input candidates TSV")],
    model: Annotated[Path, typer.Option(help="Trained SVM model (pickle)")],
    output: Annotated[Path, typer.Option(help="Output filtered TSV")],
) -> None:
    """Apply trained SVM model to filter candidate reads."""
    import numpy as np

    from aperture_snv.svm import load_model

    logger.info("Filtering candidates from {} with model {}", input, model)
    svm_model = load_model(model)

    with open(input) as fin, open(output, "w", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        writer = csv.DictWriter(fout, fieldnames=reader.fieldnames, delimiter="\t")
        writer.writeheader()

        rows = list(reader)
        if not rows:
            logger.warning("No candidate reads in input file")
            return

        features = np.array([
            [
                float(r["vbq"]),
                float(r["mrbq"]),
                float(r["pir"]),
                0.5 if r["is_concordant"] == "" else float(r["is_concordant"]),
                float(r["mapq"]),
            ]
            for r in rows
        ])

        predictions = svm_model.predict(features)
        n_passed = 0
        for row, pred in zip(rows, predictions, strict=True):
            if pred == 1:
                writer.writerow(row)
                n_passed += 1

    logger.success("SVM filter: {}/{} reads passed", n_passed, len(rows))


@app.command()
def score(
    input: Annotated[Path, typer.Option(help="Filtered candidates TSV")],
    noise_profile: Annotated[Path, typer.Option(help="Noise profile JSON")],
    mean_coverage: Annotated[float, typer.Option(help="Mean coverage at compendium sites")],
    output: Annotated[Path, typer.Option(help="Output JSON with scores")],
    z_threshold: Annotated[float, typer.Option(help="z-score detection threshold")] = 1.2,
) -> None:
    """Compute tumor fraction and detection z-score from filtered candidates."""
    from aperture_snv.score import score_sample

    logger.info("Scoring candidates from {} against noise profile {}", input, noise_profile)

    with open(noise_profile) as f:
        noise = json.load(f)

    # Count detected variants (unique sites with alt-supporting reads)
    detected_sites: set[str] = set()
    total_reads = 0
    with open(input) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            key = f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}"
            if row["read_base"] == row["alt"]:
                detected_sites.add(key)
            total_reads += 1

    result = score_sample(
        detected_variants=len(detected_sites),
        total_sites=noise["compendium_size"],
        mean_coverage=mean_coverage,
        noise_rate=noise["noise_rate_per_read"],
        total_reads_evaluated=total_reads,
        noise_mean=noise["noise_mean"],
        noise_std=noise["noise_std"],
        z_threshold=z_threshold,
    )

    output_data = {
        "detected_variants": result.detected_variants,
        "total_sites": result.total_sites,
        "mean_coverage": result.mean_coverage,
        "tumor_fraction": result.tumor_fraction,
        "detection_rate": result.detection_rate,
        "z_score": result.z_score,
        "is_detected": result.is_detected,
        "z_threshold": z_threshold,
    }

    with open(output, "w") as f:
        json.dump(output_data, f, indent=2)

    # Rich output summary
    table = Table(title="SNV Detection Result")
    table.add_column("Metric", style="cyan")
    table.add_column("Value", style="bold")
    table.add_row("Detected variants", str(result.detected_variants))
    table.add_row("Compendium sites", str(result.total_sites))
    table.add_row("Tumor fraction", f"{result.tumor_fraction:.2e}")
    table.add_row("z-score", f"{result.z_score:.2f}")
    status = "[bold green]DETECTED[/]" if result.is_detected else "[bold red]NOT DETECTED[/]"
    table.add_row("Status", status)
    console.print(table)


@app.command()
def noise(
    compendium: Annotated[Path, typer.Option(help="Patient SNV compendium VCF")],
    controls: Annotated[list[Path], typer.Option(help="Control plasma CRAM/BAM files")],
    reference: Annotated[Path, typer.Option(help="Reference FASTA path")],
    model: Annotated[Path, typer.Option(help="Trained SVM model (pickle)")],
    output: Annotated[Path, typer.Option(help="Output noise profile JSON")],
) -> None:
    """Build noise profile from control plasma samples."""
    from aperture_snv.noise import build_noise_profile
    from aperture_snv.svm import load_model

    logger.info(
        "Building noise profile from {} control samples against {}",
        len(controls), compendium,
    )

    svm_model = load_model(model)
    profile = build_noise_profile(
        compendium_path=compendium,
        control_cram_paths=controls,
        reference_path=reference,
        svm_model=svm_model,
    )

    output_data = {
        "compendium_size": profile.compendium_size,
        "n_controls": profile.n_controls,
        "detection_rates": list(profile.detection_rates),
        "noise_mean": profile.noise_mean,
        "noise_std": profile.noise_std,
        "noise_rate_per_read": profile.noise_rate_per_read,
    }

    with open(output, "w") as f:
        json.dump(output_data, f, indent=2)

    logger.success(
        "Noise profile: μ={:.2e}, σ={:.2e}, n={} controls → {}",
        profile.noise_mean, profile.noise_std, profile.n_controls, output,
    )
