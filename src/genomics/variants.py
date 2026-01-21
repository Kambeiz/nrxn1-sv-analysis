"""Variant loading and filtering utilities."""

from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any, Iterator
from pathlib import Path
import pandas as pd
import numpy as np
from .regions import GenomicRegion, NRXN1Region


@dataclass
class Variant:
    """Represents a genomic variant."""
    
    chromosome: str
    position: int
    ref: str
    alt: str
    variant_id: Optional[str] = None
    quality: float = 0.0
    filter_status: str = "PASS"
    info: Dict[str, Any] = field(default_factory=dict)
    
    @property
    def end(self) -> int:
        """Calculate end position based on variant type."""
        if self.is_snv:
            return self.position + 1
        elif self.is_deletion:
            return self.position + len(self.ref)
        else:
            return self.position + max(len(self.ref), len(self.alt))
    
    @property
    def is_snv(self) -> bool:
        """Check if variant is a SNV."""
        return len(self.ref) == 1 and len(self.alt) == 1
    
    @property
    def is_insertion(self) -> bool:
        """Check if variant is an insertion."""
        return len(self.ref) < len(self.alt)
    
    @property
    def is_deletion(self) -> bool:
        """Check if variant is a deletion."""
        return len(self.ref) > len(self.alt)
    
    @property
    def is_indel(self) -> bool:
        """Check if variant is an indel."""
        return self.is_insertion or self.is_deletion
    
    @property
    def variant_type(self) -> str:
        """Determine variant type."""
        if self.is_snv:
            return "SNV"
        elif self.is_insertion:
            return "INS"
        elif self.is_deletion:
            return "DEL"
        else:
            return "COMPLEX"
    
    @property
    def length(self) -> int:
        """Calculate variant length."""
        return abs(len(self.alt) - len(self.ref))
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert variant to dictionary."""
        return {
            "chromosome": self.chromosome,
            "position": self.position,
            "end": self.end,
            "ref": self.ref,
            "alt": self.alt,
            "variant_id": self.variant_id,
            "quality": self.quality,
            "filter_status": self.filter_status,
            "variant_type": self.variant_type,
            "length": self.length,
            **self.info
        }


@dataclass
class StructuralVariant(Variant):
    """Represents a structural variant (CNV, inversion, etc.)."""
    
    sv_type: str = "DEL"
    sv_end: Optional[int] = None
    sv_length: Optional[int] = None
    confidence: float = 0.0
    caller: str = "unknown"
    
    @property
    def end(self) -> int:
        """Return SV end position."""
        if self.sv_end is not None:
            return self.sv_end
        if self.sv_length is not None:
            return self.position + abs(self.sv_length)
        return super().end
    
    @property
    def length(self) -> int:
        """Return SV length."""
        if self.sv_length is not None:
            return abs(self.sv_length)
        return self.end - self.position
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert SV to dictionary."""
        base = super().to_dict()
        base.update({
            "sv_type": self.sv_type,
            "sv_end": self.end,
            "sv_length": self.length,
            "confidence": self.confidence,
            "caller": self.caller
        })
        return base


class VariantLoader:
    """Load variants from various file formats."""
    
    def __init__(self, nrxn1_only: bool = True):
        """Initialize loader.
        
        Args:
            nrxn1_only: If True, filter to NRXN1 region only.
        """
        self.nrxn1_only = nrxn1_only
        self.nrxn1_region = NRXN1Region() if nrxn1_only else None
    
    def load_vcf(self, vcf_path: Path) -> List[Variant]:
        """Load variants from VCF file."""
        variants = []
        
        try:
            import cyvcf2
            vcf = cyvcf2.VCF(str(vcf_path))
            
            for record in vcf:
                if self.nrxn1_only:
                    chrom = record.CHROM.replace("chr", "")
                    if chrom != NRXN1Region.CHROMOSOME:
                        continue
                    if not (NRXN1Region.START <= record.POS <= NRXN1Region.END):
                        continue
                
                for alt in record.ALT:
                    info_dict = dict(record.INFO)
                    
                    if info_dict.get("SVTYPE"):
                        sv = StructuralVariant(
                            chromosome=record.CHROM,
                            position=record.POS,
                            ref=record.REF,
                            alt=alt or "<SV>",
                            variant_id=record.ID,
                            quality=record.QUAL or 0.0,
                            filter_status=";".join(record.FILTER) if record.FILTER else "PASS",
                            sv_type=info_dict.get("SVTYPE", "unknown"),
                            sv_end=info_dict.get("END"),
                            sv_length=info_dict.get("SVLEN"),
                            info=info_dict
                        )
                        variants.append(sv)
                    else:
                        var = Variant(
                            chromosome=record.CHROM,
                            position=record.POS,
                            ref=record.REF,
                            alt=alt,
                            variant_id=record.ID,
                            quality=record.QUAL or 0.0,
                            filter_status=";".join(record.FILTER) if record.FILTER else "PASS",
                            info=info_dict
                        )
                        variants.append(var)
            
            vcf.close()
            
        except ImportError:
            raise ImportError("cyvcf2 required for VCF parsing. Install with: pip install cyvcf2")
        
        return variants
    
    def load_bed(self, bed_path: Path) -> List[StructuralVariant]:
        """Load structural variants from BED file."""
        variants = []
        
        with open(bed_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                
                fields = line.strip().split("\t")
                if len(fields) < 3:
                    continue
                
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                name = fields[3] if len(fields) > 3 else None
                
                if self.nrxn1_only:
                    if chrom.replace("chr", "") != NRXN1Region.CHROMOSOME:
                        continue
                    if end < NRXN1Region.START or start > NRXN1Region.END:
                        continue
                
                sv = StructuralVariant(
                    chromosome=chrom,
                    position=start,
                    ref="N",
                    alt="<DEL>",
                    sv_type="DEL",
                    sv_end=end,
                    variant_id=name
                )
                variants.append(sv)
        
        return variants
    
    def load_dataframe(self, df: pd.DataFrame) -> List[Variant]:
        """Load variants from pandas DataFrame."""
        required_cols = ["chromosome", "position", "ref", "alt"]
        for col in required_cols:
            if col not in df.columns:
                raise ValueError(f"Missing required column: {col}")
        
        variants = []
        for _, row in df.iterrows():
            var = Variant(
                chromosome=str(row["chromosome"]),
                position=int(row["position"]),
                ref=str(row["ref"]),
                alt=str(row["alt"]),
                variant_id=row.get("variant_id"),
                quality=float(row.get("quality", 0.0)),
                filter_status=str(row.get("filter_status", "PASS"))
            )
            
            if self.nrxn1_only:
                if var.chromosome.replace("chr", "") != NRXN1Region.CHROMOSOME:
                    continue
                if not (NRXN1Region.START <= var.position <= NRXN1Region.END):
                    continue
            
            variants.append(var)
        
        return variants
    
    def to_dataframe(self, variants: List[Variant]) -> pd.DataFrame:
        """Convert list of variants to DataFrame."""
        return pd.DataFrame([v.to_dict() for v in variants])


class VariantFilter:
    """Filter variants based on various criteria."""
    
    def __init__(self):
        """Initialize filter."""
        self.nrxn1 = NRXN1Region()
    
    def filter_by_quality(
        self, 
        variants: List[Variant], 
        min_quality: float = 20.0
    ) -> List[Variant]:
        """Filter variants by quality score."""
        return [v for v in variants if v.quality >= min_quality]
    
    def filter_by_type(
        self, 
        variants: List[Variant], 
        variant_types: List[str]
    ) -> List[Variant]:
        """Filter variants by type (SNV, INS, DEL, etc.)."""
        return [v for v in variants if v.variant_type in variant_types]
    
    def filter_by_size(
        self, 
        variants: List[Variant], 
        min_size: int = 0, 
        max_size: int = float("inf")
    ) -> List[Variant]:
        """Filter variants by size."""
        return [v for v in variants if min_size <= v.length <= max_size]
    
    def filter_exonic(
        self, 
        variants: List[Variant], 
        isoform: str = "alpha"
    ) -> List[Variant]:
        """Filter to only exonic variants."""
        exonic = []
        for v in variants:
            affected = self.nrxn1.get_affected_exons(v.position, v.end, isoform)
            if affected:
                exonic.append(v)
        return exonic
    
    def filter_by_region(
        self, 
        variants: List[Variant], 
        region: GenomicRegion
    ) -> List[Variant]:
        """Filter variants overlapping a specific region."""
        filtered = []
        for v in variants:
            v_region = GenomicRegion(v.chromosome, v.position, v.end)
            if v_region.overlaps(region):
                filtered.append(v)
        return filtered
    
    def filter_pass_only(self, variants: List[Variant]) -> List[Variant]:
        """Keep only PASS variants."""
        return [v for v in variants if v.filter_status == "PASS"]
    
    def apply_filters(
        self,
        variants: List[Variant],
        min_quality: Optional[float] = None,
        variant_types: Optional[List[str]] = None,
        min_size: Optional[int] = None,
        max_size: Optional[int] = None,
        exonic_only: bool = False,
        pass_only: bool = True
    ) -> List[Variant]:
        """Apply multiple filters in sequence."""
        filtered = variants
        
        if pass_only:
            filtered = self.filter_pass_only(filtered)
        
        if min_quality is not None:
            filtered = self.filter_by_quality(filtered, min_quality)
        
        if variant_types is not None:
            filtered = self.filter_by_type(filtered, variant_types)
        
        if min_size is not None or max_size is not None:
            filtered = self.filter_by_size(
                filtered,
                min_size or 0,
                max_size or float("inf")
            )
        
        if exonic_only:
            filtered = self.filter_exonic(filtered)
        
        return filtered
