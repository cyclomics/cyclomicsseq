from pathlib import Path

from .algorithms.fraction import FractionCalibratedConsensus


def file_prep(
    hdf5_file,
    reference_file,
    output_dir,
    method,
    min_alignment_quality,
    vcf_file,
    skip_variants,
):
    if method == "fraction":
        calibrator = FractionCalibratedConsensus()
    else:
        raise ValueError(f"Unknown method: {method}")

    calibrator.prep_file(
        hdf5_file=Path(hdf5_file),
        reference_file=reference_file,
        output_dir=Path(output_dir),
        min_alignment_quality=min_alignment_quality,
        vcf_file=vcf_file,
        skip_variants=skip_variants,
    )


def model_fit(
    input_dir,
    output_dir,
    method,
):
    if method == "fraction":
        calibrator = FractionCalibratedConsensus()
    else:
        raise ValueError(f"Unknown method: {method}")

    calibrator.fit_model(
        input_dir=input_dir,
        output_dir=output_dir,
    )
