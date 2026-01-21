"""CNV detection module for NRXN1 structural variants."""

from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any, Tuple
from pathlib import Path
import numpy as np
import pandas as pd
from enum import Enum

from ..genomics.regions import NRXN1Region, GenomicRegion


class CNVType(Enum):
    """CNV type enumeration."""
    DELETION = "DEL"
    DUPLICATION = "DUP"
    INVERSION = "INV"
    TRANSLOCATION = "TRA"
    INSERTION = "INS"
    COMPLEX = "CPX"


@dataclass
class CNVCall:
    """Represents a CNV call."""
    
    chromosome: str
    start: int
    end: int
    cnv_type: CNVType
    caller: str
    quality: float = 0.0
    genotype: Optional[str] = None
    copy_number: Optional[int] = None
    log2_ratio: Optional[float] = None
    read_depth: Optional[float] = None
    split_reads: int = 0
    paired_end_reads: int = 0
    sample_id: Optional[str] = None
    info: Dict[str, Any] = field(default_factory=dict)
    
    @property
    def length(self) -> int:
        """Return CNV length."""
        return self.end - self.start
    
    @property
    def is_deletion(self) -> bool:
        """Check if CNV is a deletion."""
        return self.cnv_type == CNVType.DELETION
    
    @property
    def is_duplication(self) -> bool:
        """Check if CNV is a duplication."""
        return self.cnv_type == CNVType.DUPLICATION
    
    def overlaps(self, other: "CNVCall", min_reciprocal: float = 0.5) -> bool:
        """Check if two CNVs overlap with minimum reciprocal overlap."""
        if self.chromosome != other.chromosome:
            return False
        
        overlap_start = max(self.start, other.start)
        overlap_end = min(self.end, other.end)
        
        if overlap_start >= overlap_end:
            return False
        
        overlap_length = overlap_end - overlap_start
        overlap_self = overlap_length / self.length
        overlap_other = overlap_length / other.length
        
        return overlap_self >= min_reciprocal and overlap_other >= min_reciprocal
    
    def get_affected_exons(self, nrxn1: NRXN1Region) -> List[GenomicRegion]:
        """Get exons affected by this CNV."""
        return nrxn1.get_affected_exons(self.start, self.end, "alpha")
    
    def get_clinical_significance(self) -> str:
        """Preliminary clinical significance assessment."""
        nrxn1 = NRXN1Region()
        affected_exons = self.get_affected_exons(nrxn1)
        location = nrxn1.classify_variant_location(self.start, self.end)
        
        if not affected_exons:
            return "likely_benign"
        
        if self.is_deletion:
            if len(affected_exons) >= 3:
                return "likely_pathogenic"
            elif location == "alpha_specific":
                return "uncertain_significance"
            else:
                return "uncertain_significance"
        elif self.is_duplication:
            return "uncertain_significance"
        
        return "uncertain_significance"
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "chromosome": self.chromosome,
            "start": self.start,
            "end": self.end,
            "length": self.length,
            "cnv_type": self.cnv_type.value,
            "caller": self.caller,
            "quality": self.quality,
            "genotype": self.genotype,
            "copy_number": self.copy_number,
            "log2_ratio": self.log2_ratio,
            "split_reads": self.split_reads,
            "paired_end_reads": self.paired_end_reads,
            "sample_id": self.sample_id,
            **self.info
        }
    
    def to_bed(self) -> str:
        """Return BED format string."""
        return f"{self.chromosome}\t{self.start}\t{self.end}\t{self.cnv_type.value}\t{self.quality}"


class CNVDetector:
    """Multi-method CNV detector for NRXN1 region."""
    
    SUPPORTED_CALLERS = ["manta", "delly", "cnvnator", "lumpy", "gatk"]
    
    def __init__(
        self,
        min_size: int = 1000,
        max_size: int = 10_000_000,
        min_quality: float = 20.0,
        nrxn1_only: bool = True
    ):
        """Initialize CNV detector.
        
        Args:
            min_size: Minimum CNV size to report.
            max_size: Maximum CNV size to report.
            min_quality: Minimum quality score.
            nrxn1_only: Restrict to NRXN1 region.
        """
        self.min_size = min_size
        self.max_size = max_size
        self.min_quality = min_quality
        self.nrxn1_only = nrxn1_only
        self.nrxn1 = NRXN1Region() if nrxn1_only else None
    
    def detect_from_depth(
        self,
        depth_array: np.ndarray,
        positions: np.ndarray,
        sample_id: str = "sample",
        window_size: int = 1000,
        z_threshold: float = 3.0
    ) -> List[CNVCall]:
        """Detect CNVs from read depth data.
        
        Args:
            depth_array: Array of read depths.
            positions: Array of genomic positions.
            sample_id: Sample identifier.
            window_size: Sliding window size.
            z_threshold: Z-score threshold for calling.
        
        Returns:
            List of CNV calls.
        """
        calls = []
        
        mean_depth = np.mean(depth_array)
        std_depth = np.std(depth_array)
        
        if std_depth == 0:
            return calls
        
        z_scores = (depth_array - mean_depth) / std_depth
        log2_ratios = np.log2(depth_array / mean_depth + 0.001)
        
        in_cnv = False
        cnv_start = 0
        cnv_type = None
        cnv_depths = []
        
        for i, (pos, z, depth) in enumerate(zip(positions, z_scores, depth_array)):
            if not in_cnv:
                if z < -z_threshold:
                    in_cnv = True
                    cnv_start = pos
                    cnv_type = CNVType.DELETION
                    cnv_depths = [depth]
                elif z > z_threshold:
                    in_cnv = True
                    cnv_start = pos
                    cnv_type = CNVType.DUPLICATION
                    cnv_depths = [depth]
            else:
                is_same_type = (
                    (cnv_type == CNVType.DELETION and z < -z_threshold / 2) or
                    (cnv_type == CNVType.DUPLICATION and z > z_threshold / 2)
                )
                
                if is_same_type:
                    cnv_depths.append(depth)
                else:
                    cnv_end = positions[i - 1]
                    cnv_length = cnv_end - cnv_start
                    
                    if self.min_size <= cnv_length <= self.max_size:
                        avg_depth = np.mean(cnv_depths)
                        log2 = np.log2(avg_depth / mean_depth + 0.001)
                        
                        call = CNVCall(
                            chromosome=NRXN1Region.CHROMOSOME,
                            start=cnv_start,
                            end=cnv_end,
                            cnv_type=cnv_type,
                            caller="depth_analysis",
                            quality=abs(log2) * 10,
                            log2_ratio=log2,
                            read_depth=avg_depth,
                            sample_id=sample_id,
                            copy_number=self._estimate_copy_number(log2)
                        )
                        calls.append(call)
                    
                    in_cnv = False
                    cnv_depths = []
        
        return calls
    
    def _estimate_copy_number(self, log2_ratio: float) -> int:
        """Estimate copy number from log2 ratio."""
        ratio = 2 ** log2_ratio
        copy_number = round(ratio * 2)
        return max(0, copy_number)
    
    def parse_manta_vcf(self, vcf_path: Path, sample_id: str = None) -> List[CNVCall]:
        """Parse CNV calls from Manta VCF output."""
        calls = []
        
        try:
            import cyvcf2
            vcf = cyvcf2.VCF(str(vcf_path))
            
            for record in vcf:
                sv_type = record.INFO.get("SVTYPE")
                if sv_type not in ["DEL", "DUP", "INV", "INS"]:
                    continue
                
                chrom = record.CHROM.replace("chr", "")
                start = record.POS
                end = record.INFO.get("END", start)
                
                if self.nrxn1_only:
                    if chrom != NRXN1Region.CHROMOSOME:
                        continue
                    if end < NRXN1Region.START or start > NRXN1Region.END:
                        continue
                
                length = end - start
                if not (self.min_size <= length <= self.max_size):
                    continue
                
                cnv_type = CNVType(sv_type)
                
                call = CNVCall(
                    chromosome=chrom,
                    start=start,
                    end=end,
                    cnv_type=cnv_type,
                    caller="manta",
                    quality=record.QUAL or 0.0,
                    split_reads=record.INFO.get("SR", 0),
                    paired_end_reads=record.INFO.get("PR", 0),
                    sample_id=sample_id,
                    info=dict(record.INFO)
                )
                calls.append(call)
            
            vcf.close()
            
        except ImportError:
            raise ImportError("cyvcf2 required. Install with: pip install cyvcf2")
        
        return calls
    
    def parse_delly_vcf(self, vcf_path: Path, sample_id: str = None) -> List[CNVCall]:
        """Parse CNV calls from DELLY VCF output."""
        calls = []
        
        try:
            import cyvcf2
            vcf = cyvcf2.VCF(str(vcf_path))
            
            for record in vcf:
                sv_type = record.INFO.get("SVTYPE")
                if sv_type not in ["DEL", "DUP", "INV", "INS"]:
                    continue
                
                chrom = record.CHROM.replace("chr", "")
                start = record.POS
                end = record.INFO.get("END", start)
                
                if self.nrxn1_only:
                    if chrom != NRXN1Region.CHROMOSOME:
                        continue
                    if end < NRXN1Region.START or start > NRXN1Region.END:
                        continue
                
                length = end - start
                if not (self.min_size <= length <= self.max_size):
                    continue
                
                cnv_type = CNVType(sv_type)
                pe_support = record.INFO.get("PE", 0)
                sr_support = record.INFO.get("SR", 0)
                
                call = CNVCall(
                    chromosome=chrom,
                    start=start,
                    end=end,
                    cnv_type=cnv_type,
                    caller="delly",
                    quality=record.QUAL or 0.0,
                    split_reads=sr_support,
                    paired_end_reads=pe_support,
                    sample_id=sample_id,
                    info=dict(record.INFO)
                )
                calls.append(call)
            
            vcf.close()
            
        except ImportError:
            raise ImportError("cyvcf2 required. Install with: pip install cyvcf2")
        
        return calls
    
    def filter_calls(
        self,
        calls: List[CNVCall],
        min_quality: Optional[float] = None,
        cnv_types: Optional[List[CNVType]] = None,
        exonic_only: bool = False
    ) -> List[CNVCall]:
        """Filter CNV calls based on criteria."""
        filtered = calls
        
        if min_quality is not None:
            filtered = [c for c in filtered if c.quality >= min_quality]
        
        if cnv_types is not None:
            filtered = [c for c in filtered if c.cnv_type in cnv_types]
        
        if exonic_only and self.nrxn1:
            exonic = []
            for call in filtered:
                affected = call.get_affected_exons(self.nrxn1)
                if affected:
                    exonic.append(call)
            filtered = exonic
        
        return filtered
    
    def to_dataframe(self, calls: List[CNVCall]) -> pd.DataFrame:
        """Convert CNV calls to DataFrame."""
        if not calls:
            return pd.DataFrame()
        return pd.DataFrame([c.to_dict() for c in calls])
    
    def to_bed(self, calls: List[CNVCall], output_path: Path) -> None:
        """Write CNV calls to BED file."""
        with open(output_path, "w") as f:
            f.write("#chrom\tstart\tend\ttype\tquality\n")
            for call in calls:
                f.write(call.to_bed() + "\n")
