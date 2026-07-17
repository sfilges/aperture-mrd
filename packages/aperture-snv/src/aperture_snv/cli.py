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
from pathlib import Path  # noqa: TC003 — runtime import; typer resolves annotations at runtime
from typing import Annotated

import typer
from loguru import logger
from rich.console import Console
from rich.table import Table

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
    output: Annotated[Path, typer.Option(help="Output TSV file path")],
    reference: Annotated[
        Path | None, typer.Option(help="Reference FASTA path (required for CRAM)")
    ] = None,
    min_mapq: Annotated[int, typer.Option(help="Min mapping quality")] = 20,
    min_bq: Annotated[int, typer.Option(help="Min base quality at the variant site")] = 20,
    min_insert: Annotated[int, typer.Option(help="Min fragment insert size (bp)")] = 40,
    max_insert: Annotated[int, typer.Option(help="Max fragment insert size (bp)")] = 240,
) -> None:
    """Extract candidate fragments from a plasma CRAM/BAM at compendium loci."""
    from aperture_snv.extract import extract_all_fragments
    from aperture_snv.filters import FilterConfig

    logger.info("Extracting fragments from {} at {} compendium sites", cram, compendium)

    config = FilterConfig(
        min_mapq=min_mapq, min_bq=min_bq, min_insert=min_insert, max_insert=max_insert
    )
    candidates = extract_all_fragments(
        alignment_file_path=cram,
        compendium_vcf_path=compendium,
        reference_fasta_path=reference,
        config=config,
    )

    columns = [
        "chrom",
        "pos",
        "ref",
        "alt",
        "frag_id",
        "read_base",
        "supports_alt",
        "vbq",
        "mrbq",
        "n_low_bq",
        "mapq",
        "edit_distance",
        "pir",
        "dist_3prime",
        "concordance",
        "fragment_length",
        "softclip_len",
        "is_proper_pair",
        "n_mates_at_site",
    ]
    with open(output, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(columns)
        for _site, reads in candidates.items():
            for r in reads:
                writer.writerow(
                    [
                        r.chrom,
                        r.pos,
                        r.ref,
                        r.alt,
                        r.frag_id,
                        r.read_base,
                        int(r.supports_alt),
                        r.vbq,
                        f"{r.mrbq:.2f}",
                        r.n_low_bq,
                        r.mapq,
                        r.edit_distance,
                        f"{r.pir:.4f}",
                        r.dist_3prime,
                        r.concordance,
                        r.fragment_length,
                        r.softclip_len,
                        int(r.is_proper_pair),
                        r.n_mates_at_site,
                    ]
                )

    n_sites = len(candidates)
    n_reads = sum(len(reads) for reads in candidates.values())
    logger.success("Extracted {} fragments at {} compendium sites → {}", n_reads, n_sites, output)


@app.command()
def filter(
    input: Annotated[Path, typer.Option(help="Input candidates TSV")],
    model: Annotated[Path, typer.Option(help="Trained SVM model (pickle)")],
    output: Annotated[Path, typer.Option(help="Output filtered TSV")],
) -> None:
    """Apply trained SVM model to filter candidate reads."""
    import numpy as np

    from aperture_snv.models.svm import _CONCORDANCE_ENCODING, load_model

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

        features = np.array(
            [
                [
                    float(r["vbq"]),
                    float(r["mrbq"]),
                    float(r["pir"]),
                    _CONCORDANCE_ENCODING.get(r["concordance"], 0.5),
                    float(r["mapq"]),
                ]
                for r in rows
            ]
        )

        predictions = svm_model.predict(features)
        n_passed = 0
        for row, pred in zip(rows, predictions, strict=True):
            if pred == 1:
                writer.writerow(row)
                n_passed += 1

    logger.success("SVM filter: {}/{} reads passed", n_passed, len(rows))


@app.command()
def score(
    input: Annotated[Path, typer.Option(help="Extracted fragments TSV")],
    bam: Annotated[Path, typer.Option(help="Plasma CRAM/BAM (for the aligned-read count)")],
    noise_profile: Annotated[Path, typer.Option(help="Noise profile JSON")],
    output: Annotated[Path, typer.Option(help="Output JSON with scores")],
    reference: Annotated[
        Path | None, typer.Option(help="Reference FASTA path (required for CRAM)")
    ] = None,
    z_threshold: Annotated[float, typer.Option(help="z-score detection threshold")] = 1.2,
) -> None:
    """Compute the PPM mutational load and detection z-score from extracted fragments."""
    from aperture_snv.extract import count_mapped_reads
    from aperture_snv.score import score_sample

    logger.info("Scoring fragments from {} against noise profile {}", input, noise_profile)

    with open(noise_profile) as f:
        noise = json.load(f)

    # Numerator: alt-supporting fragments (weight 1.0 until a calibrated model exists).
    supporting_reads = 0.0
    with open(input) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["supports_alt"] == "1":
                supporting_reads += 1.0

    total_aligned_reads = count_mapped_reads(bam, reference)

    result = score_sample(
        supporting_reads=supporting_reads,
        total_aligned_reads=total_aligned_reads,
        noise_mean=noise["noise_mean"],
        noise_std=noise["noise_std"],
        z_threshold=z_threshold,
    )

    output_data = {
        "supporting_reads": result.supporting_reads,
        "total_aligned_reads": result.total_aligned_reads,
        "ppm": result.ppm,
        "z_score": result.z_score,
        "is_detected": result.is_detected,
        "z_threshold": z_threshold,
    }

    with open(output, "w") as f:
        json.dump(output_data, f, indent=2)

    table = Table(title="SNV Detection Result")
    table.add_column("Metric", style="cyan")
    table.add_column("Value", style="bold")
    table.add_row("Supporting reads", f"{result.supporting_reads:.0f}")
    table.add_row("Aligned reads", str(result.total_aligned_reads))
    table.add_row("PPM", f"{result.ppm:.3f}")
    table.add_row("z-score", f"{result.z_score:.2f}")
    status = "[bold green]DETECTED[/]" if result.is_detected else "[bold red]NOT DETECTED[/]"
    table.add_row("Status", status)
    console.print(table)


@app.command()
def noise(
    compendium: Annotated[Path, typer.Option(help="Patient SNV compendium VCF")],
    controls: Annotated[list[Path], typer.Option(help="Control plasma CRAM/BAM files")],
    output: Annotated[Path, typer.Option(help="Output noise profile JSON")],
    reference: Annotated[
        Path | None, typer.Option(help="Reference FASTA path (required for CRAM)")
    ] = None,
) -> None:
    """Build the background PPM distribution from control plasma samples."""
    from aperture_snv.noise import build_noise_profile

    logger.info(
        "Building noise profile from {} control samples against {}",
        len(controls),
        compendium,
    )

    profile = build_noise_profile(
        compendium_path=compendium,
        control_paths=controls,
        reference_path=reference,
    )

    output_data = {
        "n_controls": profile.n_controls,
        "control_ppms": list(profile.control_ppms),
        "noise_mean": profile.noise_mean,
        "noise_std": profile.noise_std,
    }

    with open(output, "w") as f:
        json.dump(output_data, f, indent=2)

    logger.success(
        "Noise profile: μ={:.3f} PPM, σ={:.3f}, n={} controls → {}",
        profile.noise_mean,
        profile.noise_std,
        profile.n_controls,
        output,
    )
