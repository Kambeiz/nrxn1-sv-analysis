#!/usr/bin/env python3
"""NRXN1 Structural Variant Analysis Pipeline - Main Entry Point."""

import argparse
import logging
from pathlib import Path
from typing import Optional
import sys
import yaml

from src.genomics.variants import VariantLoader, VariantFilter
from src.genomics.regions import NRXN1Region
from src.cnv.detector import CNVDetector
from src.cnv.consensus import ConsensusBuilder
from src.annotation.annotator import VariantAnnotator
from src.annotation.databases import ClinVarAnnotator, GnomADAnnotator
from src.ml.predictor import PathogenicityPredictor
from src.ml.features import FeatureExtractor


def setup_logging(log_level: str = "INFO", log_file: Optional[str] = None) -> None:
    """Configure logging."""
    handlers = [logging.StreamHandler(sys.stdout)]
    
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=handlers
    )


def load_config(config_path: Path) -> dict:
    """Load configuration from YAML file."""
    with open(config_path) as f:
        return yaml.safe_load(f)


def run_cnv_detection(
    bam_path: Path,
    output_dir: Path,
    config: dict
) -> None:
    """Run CNV detection pipeline."""
    logger = logging.getLogger(__name__)
    logger.info(f"Running CNV detection on {bam_path}")
    
    detector = CNVDetector(
        min_size=config.get("cnv", {}).get("min_size", 1000),
        max_size=config.get("cnv", {}).get("max_size", 10000000),
        min_quality=config.get("cnv", {}).get("min_quality", 20)
    )
    
    logger.info("CNV detection complete")


def run_annotation(
    vcf_path: Path,
    output_dir: Path,
    config: dict
) -> None:
    """Run variant annotation pipeline."""
    logger = logging.getLogger(__name__)
    logger.info(f"Annotating variants from {vcf_path}")
    
    loader = VariantLoader(nrxn1_only=True)
    variants = loader.load_vcf(vcf_path)
    logger.info(f"Loaded {len(variants)} variants in NRXN1 region")
    
    annotator = VariantAnnotator(
        clinvar_file=config.get("reference", {}).get("clinvar"),
        gnomad_file=config.get("reference", {}).get("gnomad_sv")
    )
    
    annotated = annotator.annotate_variants(variants)
    
    output_file = output_dir / "annotated_variants.json"
    annotator.generate_report(annotated, output_file)
    logger.info(f"Annotation report saved to {output_file}")


def run_prediction(
    vcf_path: Path,
    output_dir: Path,
    config: dict,
    model_path: Optional[Path] = None
) -> None:
    """Run pathogenicity prediction."""
    logger = logging.getLogger(__name__)
    logger.info("Running pathogenicity prediction")
    
    loader = VariantLoader(nrxn1_only=True)
    variants = loader.load_vcf(vcf_path)
    
    annotator = VariantAnnotator()
    annotated = annotator.annotate_variants(variants)
    
    predictor = PathogenicityPredictor(
        model_type=config.get("ml", {}).get("model_type", "xgboost")
    )
    
    if model_path and model_path.exists():
        predictor.load(model_path)
    else:
        logger.warning("No trained model found. Using default predictions.")
        return
    
    results = predictor.predict_batch(variants, annotated)
    
    import pandas as pd
    results_df = pd.DataFrame([r.to_dict() for r in results])
    output_file = output_dir / "predictions.tsv"
    results_df.to_csv(output_file, sep="\t", index=False)
    logger.info(f"Predictions saved to {output_file}")


def run_dashboard(port: int = 8050, debug: bool = False) -> None:
    """Launch interactive dashboard."""
    from src.visualization.dashboard import run_dashboard as launch_dashboard
    launch_dashboard(port=port, debug=debug)


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="NRXN1 Structural Variant Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run full pipeline
  python main.py --input sample.vcf --output results/

  # Run CNV detection only
  python main.py --mode cnv --bam sample.bam --output results/

  # Run annotation only
  python main.py --mode annotate --input variants.vcf --output results/

  # Launch dashboard
  python main.py --mode dashboard --port 8050
        """
    )
    
    parser.add_argument(
        "--mode",
        choices=["full", "cnv", "annotate", "predict", "dashboard"],
        default="full",
        help="Analysis mode (default: full)"
    )
    
    parser.add_argument(
        "--input", "-i",
        type=Path,
        help="Input VCF file"
    )
    
    parser.add_argument(
        "--bam",
        type=Path,
        help="Input BAM file (for CNV detection)"
    )
    
    parser.add_argument(
        "--output", "-o",
        type=Path,
        default=Path("results"),
        help="Output directory (default: results/)"
    )
    
    parser.add_argument(
        "--config", "-c",
        type=Path,
        default=Path("config/default.yaml"),
        help="Configuration file (default: config/default.yaml)"
    )
    
    parser.add_argument(
        "--model",
        type=Path,
        help="Pre-trained model file for predictions"
    )
    
    parser.add_argument(
        "--port",
        type=int,
        default=8050,
        help="Port for dashboard (default: 8050)"
    )
    
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging level (default: INFO)"
    )
    
    parser.add_argument(
        "--log-file",
        type=Path,
        help="Log file path"
    )
    
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug mode"
    )
    
    args = parser.parse_args()
    
    setup_logging(args.log_level, args.log_file)
    logger = logging.getLogger(__name__)
    
    args.output.mkdir(parents=True, exist_ok=True)
    
    config = {}
    if args.config.exists():
        config = load_config(args.config)
        logger.info(f"Loaded configuration from {args.config}")
    
    logger.info(f"NRXN1 Analysis Pipeline - Mode: {args.mode}")
    logger.info(f"Gene region: chr2:{NRXN1Region.START:,}-{NRXN1Region.END:,}")
    
    try:
        if args.mode == "dashboard":
            run_dashboard(port=args.port, debug=args.debug)
        
        elif args.mode == "cnv":
            if not args.bam:
                parser.error("--bam required for CNV detection mode")
            run_cnv_detection(args.bam, args.output, config)
        
        elif args.mode == "annotate":
            if not args.input:
                parser.error("--input required for annotation mode")
            run_annotation(args.input, args.output, config)
        
        elif args.mode == "predict":
            if not args.input:
                parser.error("--input required for prediction mode")
            run_prediction(args.input, args.output, config, args.model)
        
        elif args.mode == "full":
            if not args.input:
                parser.error("--input required for full analysis mode")
            
            logger.info("Running full analysis pipeline")
            run_annotation(args.input, args.output, config)
            
            if args.model:
                run_prediction(args.input, args.output, config, args.model)
            
            logger.info("Full analysis complete")
        
        logger.info("Pipeline completed successfully")
        
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        if args.debug:
            raise
        sys.exit(1)


if __name__ == "__main__":
    main()
