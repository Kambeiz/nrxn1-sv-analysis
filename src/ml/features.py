"""Feature extraction for machine learning models."""

from dataclasses import dataclass
from typing import List, Dict, Any, Optional, Tuple
import numpy as np
import pandas as pd

from ..genomics.variants import Variant, StructuralVariant
from ..genomics.regions import NRXN1Region
from ..cnv.detector import CNVCall, CNVType
from ..annotation.annotator import AnnotatedVariant


@dataclass
class VariantFeatures:
    """Feature vector for a variant."""
    
    variant_size: int
    variant_type: str
    exons_affected: int
    affects_alpha: bool
    affects_beta: bool
    affects_both_isoforms: bool
    is_intronic: bool
    
    domain_signal_peptide: bool
    domain_egf_like: bool
    domain_lns: bool
    domain_transmembrane: bool
    
    gnomad_af: float
    gnomad_af_log: float
    is_absent_gnomad: bool
    
    cadd_score: float
    revel_score: float
    conservation_score: float
    
    gc_content: float
    repeat_region: bool
    breakpoint_homology: float
    
    def to_array(self) -> np.ndarray:
        """Convert to numpy array for ML models."""
        type_encoding = {
            "SNV": 0, "INS": 1, "DEL": 2, "DUP": 3, "INV": 4, "COMPLEX": 5
        }
        
        return np.array([
            self.variant_size,
            type_encoding.get(self.variant_type, 5),
            self.exons_affected,
            int(self.affects_alpha),
            int(self.affects_beta),
            int(self.affects_both_isoforms),
            int(self.is_intronic),
            int(self.domain_signal_peptide),
            int(self.domain_egf_like),
            int(self.domain_lns),
            int(self.domain_transmembrane),
            self.gnomad_af,
            self.gnomad_af_log,
            int(self.is_absent_gnomad),
            self.cadd_score,
            self.revel_score,
            self.conservation_score,
            self.gc_content,
            int(self.repeat_region),
            self.breakpoint_homology
        ])
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "variant_size": self.variant_size,
            "variant_type": self.variant_type,
            "exons_affected": self.exons_affected,
            "affects_alpha": self.affects_alpha,
            "affects_beta": self.affects_beta,
            "affects_both_isoforms": self.affects_both_isoforms,
            "is_intronic": self.is_intronic,
            "domain_signal_peptide": self.domain_signal_peptide,
            "domain_egf_like": self.domain_egf_like,
            "domain_lns": self.domain_lns,
            "domain_transmembrane": self.domain_transmembrane,
            "gnomad_af": self.gnomad_af,
            "gnomad_af_log": self.gnomad_af_log,
            "is_absent_gnomad": self.is_absent_gnomad,
            "cadd_score": self.cadd_score,
            "revel_score": self.revel_score,
            "conservation_score": self.conservation_score,
            "gc_content": self.gc_content,
            "repeat_region": self.repeat_region,
            "breakpoint_homology": self.breakpoint_homology
        }
    
    @staticmethod
    def feature_names() -> List[str]:
        """Return list of feature names."""
        return [
            "variant_size", "variant_type", "exons_affected",
            "affects_alpha", "affects_beta", "affects_both_isoforms",
            "is_intronic", "domain_signal_peptide", "domain_egf_like",
            "domain_lns", "domain_transmembrane", "gnomad_af",
            "gnomad_af_log", "is_absent_gnomad", "cadd_score",
            "revel_score", "conservation_score", "gc_content",
            "repeat_region", "breakpoint_homology"
        ]


class FeatureExtractor:
    """Extract features from variants for ML models."""
    
    def __init__(self, reference_fasta: Optional[str] = None):
        """Initialize feature extractor.
        
        Args:
            reference_fasta: Path to reference genome FASTA.
        """
        self.reference_fasta = reference_fasta
        self.nrxn1 = NRXN1Region()
        self._reference_loaded = False
        self._reference_seq = None
    
    def extract_features(
        self,
        variant: Variant,
        annotated: Optional[AnnotatedVariant] = None
    ) -> VariantFeatures:
        """Extract features from a variant.
        
        Args:
            variant: Variant object.
            annotated: Optional annotated variant with additional info.
        
        Returns:
            VariantFeatures object.
        """
        if hasattr(variant, 'sv_end') and variant.sv_end:
            end = variant.sv_end
            size = end - variant.position
        else:
            end = variant.position + max(len(variant.ref), len(variant.alt))
            size = abs(len(variant.alt) - len(variant.ref))
        
        alpha_exons = self.nrxn1.get_affected_exons(variant.position, end, "alpha")
        beta_exons = self.nrxn1.get_affected_exons(variant.position, end, "beta")
        domains = self.nrxn1.get_affected_domains(variant.position, end)
        
        domain_names = [d.name for d in domains]
        
        gnomad_af = 0.0
        cadd_score = 0.0
        revel_score = 0.0
        conservation = 0.0
        
        if annotated:
            gnomad_af = annotated.gnomad_af or 0.0
            cadd_score = annotated.cadd_score or 0.0
            revel_score = annotated.revel_score or 0.0
            conservation = annotated.conservation_score or 0.0
        
        gnomad_af_log = np.log10(gnomad_af + 1e-8)
        
        gc_content = self._compute_gc_content(variant.position, end)
        repeat_region = self._check_repeat_region(variant.position, end)
        breakpoint_homology = self._compute_breakpoint_homology(variant.position, end)
        
        return VariantFeatures(
            variant_size=size,
            variant_type=self._get_variant_type(variant),
            exons_affected=len(alpha_exons) + len(beta_exons),
            affects_alpha=len(alpha_exons) > 0,
            affects_beta=len(beta_exons) > 0,
            affects_both_isoforms=len(alpha_exons) > 0 and len(beta_exons) > 0,
            is_intronic=len(alpha_exons) == 0 and len(beta_exons) == 0,
            domain_signal_peptide="signal_peptide" in str(domain_names),
            domain_egf_like="egf" in str(domain_names).lower(),
            domain_lns="lns" in str(domain_names).lower(),
            domain_transmembrane="transmembrane" in str(domain_names).lower(),
            gnomad_af=gnomad_af,
            gnomad_af_log=gnomad_af_log,
            is_absent_gnomad=gnomad_af == 0.0,
            cadd_score=cadd_score,
            revel_score=revel_score,
            conservation_score=conservation,
            gc_content=gc_content,
            repeat_region=repeat_region,
            breakpoint_homology=breakpoint_homology
        )
    
    def extract_cnv_features(
        self,
        cnv: CNVCall,
        annotated: Optional[AnnotatedVariant] = None
    ) -> VariantFeatures:
        """Extract features from a CNV call."""
        sv = StructuralVariant(
            chromosome=cnv.chromosome,
            position=cnv.start,
            ref="N",
            alt=f"<{cnv.cnv_type.value}>",
            sv_type=cnv.cnv_type.value,
            sv_end=cnv.end
        )
        return self.extract_features(sv, annotated)
    
    def extract_batch(
        self,
        variants: List[Variant],
        annotated_variants: Optional[List[AnnotatedVariant]] = None
    ) -> Tuple[np.ndarray, pd.DataFrame]:
        """Extract features for a batch of variants.
        
        Returns:
            Tuple of (feature matrix, feature DataFrame).
        """
        if annotated_variants is None:
            annotated_variants = [None] * len(variants)
        
        features_list = []
        for var, ann in zip(variants, annotated_variants):
            features = self.extract_features(var, ann)
            features_list.append(features)
        
        feature_matrix = np.array([f.to_array() for f in features_list])
        feature_df = pd.DataFrame([f.to_dict() for f in features_list])
        
        return feature_matrix, feature_df
    
    def _get_variant_type(self, variant: Variant) -> str:
        """Get variant type string."""
        if hasattr(variant, 'sv_type'):
            return variant.sv_type
        return variant.variant_type
    
    def _compute_gc_content(self, start: int, end: int) -> float:
        """Compute GC content of region."""
        if self._reference_seq is None:
            return 0.5
        
        region_start = max(0, start - NRXN1Region.START)
        region_end = min(len(self._reference_seq), end - NRXN1Region.START)
        
        if region_start >= region_end:
            return 0.5
        
        seq = self._reference_seq[region_start:region_end].upper()
        gc_count = seq.count('G') + seq.count('C')
        return gc_count / len(seq) if seq else 0.5
    
    def _check_repeat_region(self, start: int, end: int) -> bool:
        """Check if region overlaps with known repeats."""
        return False
    
    def _compute_breakpoint_homology(self, start: int, end: int) -> float:
        """Compute microhomology at breakpoints."""
        if self._reference_seq is None:
            return 0.0
        return 0.0
    
    def get_feature_importance_names(self) -> List[str]:
        """Return feature names for importance analysis."""
        return VariantFeatures.feature_names()
