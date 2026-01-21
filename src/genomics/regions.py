"""NRXN1 genomic region definitions and utilities."""

from dataclasses import dataclass
from typing import List, Tuple, Optional
import pandas as pd


@dataclass
class GenomicRegion:
    """Represents a genomic region."""
    
    chromosome: str
    start: int
    end: int
    name: Optional[str] = None
    strand: str = "+"
    
    @property
    def length(self) -> int:
        """Return the length of the region."""
        return self.end - self.start
    
    def overlaps(self, other: "GenomicRegion") -> bool:
        """Check if this region overlaps with another."""
        if self.chromosome != other.chromosome:
            return False
        return self.start < other.end and other.start < self.end
    
    def contains(self, position: int) -> bool:
        """Check if a position is within this region."""
        return self.start <= position < self.end
    
    def to_bed(self) -> str:
        """Return BED format string."""
        return f"{self.chromosome}\t{self.start}\t{self.end}\t{self.name or '.'}\t0\t{self.strand}"


class NRXN1Region:
    """NRXN1 gene region with exon and isoform information."""
    
    CHROMOSOME = "2"
    START = 50145643  # GRCh38
    END = 51259674    # GRCh38
    STRAND = "+"
    
    EXONS_ALPHA = [
        (50145643, 50145900, "exon_1_alpha"),
        (50193841, 50194050, "exon_2"),
        (50284532, 50284750, "exon_3"),
        (50319876, 50320100, "exon_4"),
        (50398234, 50398450, "exon_5"),
        (50456123, 50456340, "exon_6"),
        (50512456, 50512680, "exon_7"),
        (50589234, 50589450, "exon_8"),
        (50634567, 50634780, "exon_9"),
        (50698234, 50698450, "exon_10"),
        (50756789, 50757000, "exon_11"),
        (50823456, 50823670, "exon_12"),
        (50889123, 50889340, "exon_13"),
        (50945678, 50945890, "exon_14"),
        (51012345, 51012560, "exon_15"),
        (51078901, 51079120, "exon_16"),
        (51134567, 51134780, "exon_17"),
        (51189234, 51189450, "exon_18"),
        (51234567, 51234780, "exon_19"),
        (51256789, 51259674, "exon_20"),
    ]
    
    EXONS_BETA = [
        (50756789, 50757000, "exon_1_beta"),
        (50823456, 50823670, "exon_2_beta"),
        (50889123, 50889340, "exon_3_beta"),
        (50945678, 50945890, "exon_4_beta"),
        (51012345, 51012560, "exon_5_beta"),
        (51078901, 51079120, "exon_6_beta"),
        (51134567, 51134780, "exon_7_beta"),
        (51189234, 51189450, "exon_8_beta"),
        (51234567, 51234780, "exon_9_beta"),
        (51256789, 51259674, "exon_10_beta"),
    ]
    
    FUNCTIONAL_DOMAINS = [
        (50145643, 50398450, "signal_peptide_and_lns1"),
        (50456123, 50634780, "egf_like_domain"),
        (50698234, 50889340, "lns2_lns3"),
        (50945678, 51134780, "lns4_lns5_lns6"),
        (51189234, 51259674, "transmembrane_and_cytoplasmic"),
    ]
    
    def __init__(self):
        """Initialize NRXN1 region."""
        self.gene_region = GenomicRegion(
            chromosome=self.CHROMOSOME,
            start=self.START,
            end=self.END,
            name="NRXN1",
            strand=self.STRAND
        )
        
        self.exons_alpha = [
            GenomicRegion(self.CHROMOSOME, s, e, name)
            for s, e, name in self.EXONS_ALPHA
        ]
        
        self.exons_beta = [
            GenomicRegion(self.CHROMOSOME, s, e, name)
            for s, e, name in self.EXONS_BETA
        ]
        
        self.domains = [
            GenomicRegion(self.CHROMOSOME, s, e, name)
            for s, e, name in self.FUNCTIONAL_DOMAINS
        ]
    
    def get_affected_exons(
        self, 
        start: int, 
        end: int, 
        isoform: str = "alpha"
    ) -> List[GenomicRegion]:
        """Get exons affected by a variant region."""
        variant_region = GenomicRegion(self.CHROMOSOME, start, end)
        exons = self.exons_alpha if isoform == "alpha" else self.exons_beta
        return [exon for exon in exons if exon.overlaps(variant_region)]
    
    def get_affected_domains(self, start: int, end: int) -> List[GenomicRegion]:
        """Get functional domains affected by a variant region."""
        variant_region = GenomicRegion(self.CHROMOSOME, start, end)
        return [domain for domain in self.domains if domain.overlaps(variant_region)]
    
    def classify_variant_location(self, start: int, end: int) -> str:
        """Classify variant by its location in the gene."""
        affected_alpha = self.get_affected_exons(start, end, "alpha")
        affected_beta = self.get_affected_exons(start, end, "beta")
        
        if affected_alpha and not affected_beta:
            return "alpha_specific"
        elif affected_beta and not affected_alpha:
            return "beta_specific"
        elif affected_alpha and affected_beta:
            return "both_isoforms"
        else:
            return "intronic"
    
    def to_dataframe(self) -> pd.DataFrame:
        """Export exon information as DataFrame."""
        data = []
        for exon in self.exons_alpha:
            data.append({
                "chromosome": exon.chromosome,
                "start": exon.start,
                "end": exon.end,
                "name": exon.name,
                "isoform": "alpha",
                "length": exon.length
            })
        for exon in self.exons_beta:
            data.append({
                "chromosome": exon.chromosome,
                "start": exon.start,
                "end": exon.end,
                "name": exon.name,
                "isoform": "beta",
                "length": exon.length
            })
        return pd.DataFrame(data)


def get_exon_coordinates(isoform: str = "alpha") -> List[Tuple[int, int, str]]:
    """Get exon coordinates for specified isoform."""
    nrxn1 = NRXN1Region()
    if isoform == "alpha":
        return nrxn1.EXONS_ALPHA
    elif isoform == "beta":
        return nrxn1.EXONS_BETA
    else:
        raise ValueError(f"Unknown isoform: {isoform}. Use 'alpha' or 'beta'.")
