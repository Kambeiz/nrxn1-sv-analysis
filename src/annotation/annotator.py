"""Variant annotation pipeline."""

from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional
from pathlib import Path
import pandas as pd
import subprocess
import json

from ..genomics.variants import Variant, StructuralVariant
from ..genomics.regions import NRXN1Region
from ..cnv.detector import CNVCall


@dataclass
class AnnotatedVariant:
    """Annotated variant with functional predictions."""
    
    variant: Variant
    gene: str = "NRXN1"
    transcript: Optional[str] = None
    consequence: Optional[str] = None
    impact: Optional[str] = None
    exon: Optional[str] = None
    intron: Optional[str] = None
    hgvsc: Optional[str] = None
    hgvsp: Optional[str] = None
    
    gnomad_af: Optional[float] = None
    gnomad_af_popmax: Optional[float] = None
    gnomad_homozygotes: Optional[int] = None
    
    clinvar_id: Optional[str] = None
    clinvar_significance: Optional[str] = None
    clinvar_review_status: Optional[str] = None
    
    cadd_score: Optional[float] = None
    revel_score: Optional[float] = None
    sift_prediction: Optional[str] = None
    polyphen_prediction: Optional[str] = None
    
    conservation_score: Optional[float] = None
    
    affected_isoforms: List[str] = field(default_factory=list)
    affected_domains: List[str] = field(default_factory=list)
    
    def __post_init__(self):
        """Compute isoform and domain impact."""
        nrxn1 = NRXN1Region()
        
        if hasattr(self.variant, 'end'):
            end = self.variant.end
        else:
            end = self.variant.position + 1
        
        alpha_exons = nrxn1.get_affected_exons(self.variant.position, end, "alpha")
        beta_exons = nrxn1.get_affected_exons(self.variant.position, end, "beta")
        
        if alpha_exons:
            self.affected_isoforms.append("alpha")
        if beta_exons:
            self.affected_isoforms.append("beta")
        
        domains = nrxn1.get_affected_domains(self.variant.position, end)
        self.affected_domains = [d.name for d in domains]
    
    @property
    def is_pathogenic(self) -> bool:
        """Check if variant is classified as pathogenic."""
        if self.clinvar_significance:
            return "pathogenic" in self.clinvar_significance.lower()
        return False
    
    @property
    def is_rare(self) -> bool:
        """Check if variant is rare (AF < 0.01)."""
        if self.gnomad_af is not None:
            return self.gnomad_af < 0.01
        return True
    
    @property
    def is_high_impact(self) -> bool:
        """Check if variant has high functional impact."""
        high_impact_consequences = [
            "frameshift", "stop_gained", "stop_lost",
            "start_lost", "splice_acceptor", "splice_donor"
        ]
        if self.consequence:
            return any(c in self.consequence.lower() for c in high_impact_consequences)
        return False
    
    @property
    def pathogenicity_score(self) -> float:
        """Compute aggregate pathogenicity score (0-1)."""
        score = 0.0
        weights_sum = 0.0
        
        if self.cadd_score is not None:
            cadd_norm = min(self.cadd_score / 40.0, 1.0)
            score += cadd_norm * 0.3
            weights_sum += 0.3
        
        if self.revel_score is not None:
            score += self.revel_score * 0.3
            weights_sum += 0.3
        
        if self.is_high_impact:
            score += 0.2
            weights_sum += 0.2
        
        if self.gnomad_af is not None:
            rarity_score = 1.0 if self.gnomad_af < 0.0001 else (
                0.8 if self.gnomad_af < 0.001 else (
                    0.5 if self.gnomad_af < 0.01 else 0.2
                )
            )
            score += rarity_score * 0.2
            weights_sum += 0.2
        
        if weights_sum > 0:
            return score / weights_sum
        return 0.5
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "chromosome": self.variant.chromosome,
            "position": self.variant.position,
            "ref": self.variant.ref,
            "alt": self.variant.alt,
            "gene": self.gene,
            "transcript": self.transcript,
            "consequence": self.consequence,
            "impact": self.impact,
            "exon": self.exon,
            "hgvsc": self.hgvsc,
            "hgvsp": self.hgvsp,
            "gnomad_af": self.gnomad_af,
            "gnomad_af_popmax": self.gnomad_af_popmax,
            "clinvar_id": self.clinvar_id,
            "clinvar_significance": self.clinvar_significance,
            "cadd_score": self.cadd_score,
            "revel_score": self.revel_score,
            "pathogenicity_score": self.pathogenicity_score,
            "affected_isoforms": ",".join(self.affected_isoforms),
            "affected_domains": ",".join(self.affected_domains),
            "is_pathogenic": self.is_pathogenic,
            "is_rare": self.is_rare,
            "is_high_impact": self.is_high_impact
        }


class VariantAnnotator:
    """Annotate variants with functional and clinical information."""
    
    def __init__(
        self,
        vep_cache_path: Optional[Path] = None,
        cadd_file: Optional[Path] = None,
        clinvar_file: Optional[Path] = None,
        gnomad_file: Optional[Path] = None
    ):
        """Initialize annotator with reference files.
        
        Args:
            vep_cache_path: Path to VEP cache directory.
            cadd_file: Path to CADD scores file.
            clinvar_file: Path to ClinVar VCF.
            gnomad_file: Path to gnomAD VCF.
        """
        self.vep_cache_path = vep_cache_path
        self.cadd_file = cadd_file
        self.clinvar_file = clinvar_file
        self.gnomad_file = gnomad_file
        
        self._clinvar_cache = {}
        self._gnomad_cache = {}
        self._cadd_cache = {}
    
    def annotate_variant(self, variant: Variant) -> AnnotatedVariant:
        """Annotate a single variant."""
        annotated = AnnotatedVariant(variant=variant)
        
        annotated = self._add_consequence_annotation(annotated)
        annotated = self._add_frequency_annotation(annotated)
        annotated = self._add_clinvar_annotation(annotated)
        annotated = self._add_pathogenicity_scores(annotated)
        
        return annotated
    
    def annotate_variants(self, variants: List[Variant]) -> List[AnnotatedVariant]:
        """Annotate a list of variants."""
        return [self.annotate_variant(v) for v in variants]
    
    def annotate_cnv(self, cnv: CNVCall) -> AnnotatedVariant:
        """Annotate a CNV call."""
        sv = StructuralVariant(
            chromosome=cnv.chromosome,
            position=cnv.start,
            ref="N",
            alt=f"<{cnv.cnv_type.value}>",
            sv_type=cnv.cnv_type.value,
            sv_end=cnv.end
        )
        
        annotated = AnnotatedVariant(variant=sv)
        
        nrxn1 = NRXN1Region()
        affected_exons = nrxn1.get_affected_exons(cnv.start, cnv.end, "alpha")
        
        if affected_exons:
            annotated.exon = f"{affected_exons[0].name}-{affected_exons[-1].name}"
            if len(affected_exons) >= 3:
                annotated.consequence = "exon_deletion"
                annotated.impact = "HIGH"
            else:
                annotated.consequence = "partial_exon_deletion"
                annotated.impact = "MODERATE"
        else:
            annotated.consequence = "intronic_deletion"
            annotated.impact = "LOW"
        
        return annotated
    
    def _add_consequence_annotation(
        self,
        annotated: AnnotatedVariant
    ) -> AnnotatedVariant:
        """Add consequence predictions using simple rules."""
        var = annotated.variant
        nrxn1 = NRXN1Region()
        
        if hasattr(var, 'sv_end') and var.sv_end:
            end = var.sv_end
        else:
            end = var.position + max(len(var.ref), len(var.alt))
        
        affected_exons = nrxn1.get_affected_exons(var.position, end, "alpha")
        
        if affected_exons:
            annotated.exon = affected_exons[0].name
            
            if var.is_indel:
                if var.length % 3 != 0:
                    annotated.consequence = "frameshift_variant"
                    annotated.impact = "HIGH"
                else:
                    annotated.consequence = "inframe_indel"
                    annotated.impact = "MODERATE"
            else:
                annotated.consequence = "missense_variant"
                annotated.impact = "MODERATE"
        else:
            annotated.consequence = "intron_variant"
            annotated.impact = "LOW"
        
        return annotated
    
    def _add_frequency_annotation(
        self,
        annotated: AnnotatedVariant
    ) -> AnnotatedVariant:
        """Add population frequency from gnomAD."""
        var = annotated.variant
        key = f"{var.chromosome}:{var.position}:{var.ref}:{var.alt}"
        
        if key in self._gnomad_cache:
            freq_data = self._gnomad_cache[key]
            annotated.gnomad_af = freq_data.get("af")
            annotated.gnomad_af_popmax = freq_data.get("af_popmax")
            annotated.gnomad_homozygotes = freq_data.get("homozygotes")
        
        return annotated
    
    def _add_clinvar_annotation(
        self,
        annotated: AnnotatedVariant
    ) -> AnnotatedVariant:
        """Add ClinVar clinical significance."""
        var = annotated.variant
        key = f"{var.chromosome}:{var.position}:{var.ref}:{var.alt}"
        
        if key in self._clinvar_cache:
            clinvar_data = self._clinvar_cache[key]
            annotated.clinvar_id = clinvar_data.get("id")
            annotated.clinvar_significance = clinvar_data.get("significance")
            annotated.clinvar_review_status = clinvar_data.get("review_status")
        
        return annotated
    
    def _add_pathogenicity_scores(
        self,
        annotated: AnnotatedVariant
    ) -> AnnotatedVariant:
        """Add pathogenicity prediction scores."""
        var = annotated.variant
        key = f"{var.chromosome}:{var.position}:{var.ref}:{var.alt}"
        
        if key in self._cadd_cache:
            annotated.cadd_score = self._cadd_cache[key]
        
        return annotated
    
    def load_clinvar_annotations(self, clinvar_path: Path) -> None:
        """Pre-load ClinVar annotations for NRXN1 region."""
        try:
            import cyvcf2
            vcf = cyvcf2.VCF(str(clinvar_path))
            
            region = f"{NRXN1Region.CHROMOSOME}:{NRXN1Region.START}-{NRXN1Region.END}"
            
            for record in vcf(region):
                for alt in record.ALT:
                    key = f"{record.CHROM}:{record.POS}:{record.REF}:{alt}"
                    self._clinvar_cache[key] = {
                        "id": record.ID,
                        "significance": record.INFO.get("CLNSIG"),
                        "review_status": record.INFO.get("CLNREVSTAT")
                    }
            
            vcf.close()
        except Exception:
            pass
    
    def to_dataframe(
        self,
        annotated_variants: List[AnnotatedVariant]
    ) -> pd.DataFrame:
        """Convert annotated variants to DataFrame."""
        return pd.DataFrame([av.to_dict() for av in annotated_variants])
    
    def generate_report(
        self,
        annotated_variants: List[AnnotatedVariant],
        output_path: Path
    ) -> None:
        """Generate annotation report."""
        df = self.to_dataframe(annotated_variants)
        
        pathogenic = df[df["is_pathogenic"] == True]
        high_impact = df[df["is_high_impact"] == True]
        rare = df[df["is_rare"] == True]
        
        report = {
            "summary": {
                "total_variants": len(df),
                "pathogenic_variants": len(pathogenic),
                "high_impact_variants": len(high_impact),
                "rare_variants": len(rare)
            },
            "variants": df.to_dict(orient="records")
        }
        
        with open(output_path, "w") as f:
            json.dump(report, f, indent=2, default=str)
