"""Database connectors for variant annotation."""

from dataclasses import dataclass
from typing import Dict, Optional, List, Any
from pathlib import Path
import pandas as pd

from ..genomics.regions import NRXN1Region


@dataclass
class ClinVarEntry:
    """ClinVar database entry."""
    
    variation_id: str
    chromosome: str
    position: int
    ref: str
    alt: str
    clinical_significance: str
    review_status: str
    condition: Optional[str] = None
    gene: str = "NRXN1"
    molecular_consequence: Optional[str] = None
    last_evaluated: Optional[str] = None
    
    @property
    def stars(self) -> int:
        """Convert review status to star rating."""
        status_map = {
            "practice guideline": 4,
            "reviewed by expert panel": 3,
            "criteria provided, multiple submitters, no conflicts": 2,
            "criteria provided, single submitter": 1,
            "no assertion criteria provided": 0,
            "no assertion provided": 0
        }
        return status_map.get(self.review_status.lower(), 0)
    
    @property
    def is_pathogenic(self) -> bool:
        """Check if classified as pathogenic."""
        return "pathogenic" in self.clinical_significance.lower()
    
    @property
    def is_benign(self) -> bool:
        """Check if classified as benign."""
        return "benign" in self.clinical_significance.lower()
    
    @property
    def is_vus(self) -> bool:
        """Check if variant of uncertain significance."""
        return "uncertain" in self.clinical_significance.lower()


class ClinVarAnnotator:
    """Annotate variants with ClinVar clinical significance."""
    
    def __init__(self, clinvar_path: Optional[Path] = None):
        """Initialize ClinVar annotator.
        
        Args:
            clinvar_path: Path to ClinVar VCF or TSV file.
        """
        self.clinvar_path = clinvar_path
        self._cache: Dict[str, ClinVarEntry] = {}
        self._loaded = False
    
    def load(self) -> None:
        """Load ClinVar data for NRXN1 region."""
        if self.clinvar_path is None:
            return
        
        if str(self.clinvar_path).endswith(".vcf") or str(self.clinvar_path).endswith(".vcf.gz"):
            self._load_vcf()
        else:
            self._load_tsv()
        
        self._loaded = True
    
    def _load_vcf(self) -> None:
        """Load from VCF format."""
        try:
            import cyvcf2
            vcf = cyvcf2.VCF(str(self.clinvar_path))
            
            region = f"chr{NRXN1Region.CHROMOSOME}:{NRXN1Region.START}-{NRXN1Region.END}"
            
            try:
                records = vcf(region)
            except Exception:
                region = f"{NRXN1Region.CHROMOSOME}:{NRXN1Region.START}-{NRXN1Region.END}"
                records = vcf(region)
            
            for record in records:
                for alt in record.ALT:
                    key = self._make_key(record.CHROM, record.POS, record.REF, alt)
                    
                    entry = ClinVarEntry(
                        variation_id=record.ID or "",
                        chromosome=record.CHROM,
                        position=record.POS,
                        ref=record.REF,
                        alt=alt,
                        clinical_significance=record.INFO.get("CLNSIG", ""),
                        review_status=record.INFO.get("CLNREVSTAT", ""),
                        condition=record.INFO.get("CLNDN", ""),
                        molecular_consequence=record.INFO.get("MC", "")
                    )
                    self._cache[key] = entry
            
            vcf.close()
        except ImportError:
            pass
    
    def _load_tsv(self) -> None:
        """Load from TSV format."""
        df = pd.read_csv(self.clinvar_path, sep="\t")
        
        df = df[df["GeneSymbol"] == "NRXN1"]
        
        for _, row in df.iterrows():
            key = self._make_key(
                str(row.get("Chromosome", "")),
                int(row.get("Start", 0)),
                str(row.get("ReferenceAllele", "")),
                str(row.get("AlternateAllele", ""))
            )
            
            entry = ClinVarEntry(
                variation_id=str(row.get("VariationID", "")),
                chromosome=str(row.get("Chromosome", "")),
                position=int(row.get("Start", 0)),
                ref=str(row.get("ReferenceAllele", "")),
                alt=str(row.get("AlternateAllele", "")),
                clinical_significance=str(row.get("ClinicalSignificance", "")),
                review_status=str(row.get("ReviewStatus", "")),
                condition=str(row.get("PhenotypeList", ""))
            )
            self._cache[key] = entry
    
    def _make_key(self, chrom: str, pos: int, ref: str, alt: str) -> str:
        """Create lookup key."""
        chrom = str(chrom).replace("chr", "")
        return f"{chrom}:{pos}:{ref}:{alt}"
    
    def annotate(
        self,
        chromosome: str,
        position: int,
        ref: str,
        alt: str
    ) -> Optional[ClinVarEntry]:
        """Get ClinVar annotation for a variant."""
        if not self._loaded:
            self.load()
        
        key = self._make_key(chromosome, position, ref, alt)
        return self._cache.get(key)
    
    def get_pathogenic_variants(self) -> List[ClinVarEntry]:
        """Get all pathogenic variants in cache."""
        return [e for e in self._cache.values() if e.is_pathogenic]
    
    def get_statistics(self) -> Dict[str, int]:
        """Get statistics on cached variants."""
        pathogenic = sum(1 for e in self._cache.values() if e.is_pathogenic)
        benign = sum(1 for e in self._cache.values() if e.is_benign)
        vus = sum(1 for e in self._cache.values() if e.is_vus)
        
        return {
            "total": len(self._cache),
            "pathogenic": pathogenic,
            "benign": benign,
            "vus": vus
        }


@dataclass
class GnomADEntry:
    """gnomAD database entry."""
    
    chromosome: str
    position: int
    ref: str
    alt: str
    allele_frequency: float
    allele_count: int
    allele_number: int
    homozygote_count: int
    popmax_af: Optional[float] = None
    popmax_population: Optional[str] = None
    filters: str = "PASS"
    
    @property
    def is_rare(self) -> bool:
        """Check if variant is rare (AF < 0.01)."""
        return self.allele_frequency < 0.01
    
    @property
    def is_ultra_rare(self) -> bool:
        """Check if variant is ultra-rare (AF < 0.0001)."""
        return self.allele_frequency < 0.0001
    
    @property
    def is_singleton(self) -> bool:
        """Check if variant is singleton."""
        return self.allele_count == 1


class GnomADAnnotator:
    """Annotate variants with gnomAD population frequencies."""
    
    POPULATIONS = [
        "afr", "amr", "asj", "eas", "fin", "nfe", "oth", "sas"
    ]
    
    def __init__(self, gnomad_path: Optional[Path] = None):
        """Initialize gnomAD annotator.
        
        Args:
            gnomad_path: Path to gnomAD VCF or sites file.
        """
        self.gnomad_path = gnomad_path
        self._cache: Dict[str, GnomADEntry] = {}
        self._loaded = False
    
    def load(self) -> None:
        """Load gnomAD data for NRXN1 region."""
        if self.gnomad_path is None:
            return
        
        try:
            import cyvcf2
            vcf = cyvcf2.VCF(str(self.gnomad_path))
            
            region = f"chr{NRXN1Region.CHROMOSOME}:{NRXN1Region.START}-{NRXN1Region.END}"
            
            try:
                records = vcf(region)
            except Exception:
                region = f"{NRXN1Region.CHROMOSOME}:{NRXN1Region.START}-{NRXN1Region.END}"
                records = vcf(region)
            
            for record in records:
                for i, alt in enumerate(record.ALT):
                    key = self._make_key(record.CHROM, record.POS, record.REF, alt)
                    
                    af = record.INFO.get("AF")
                    if isinstance(af, (list, tuple)):
                        af = af[i] if i < len(af) else af[0]
                    
                    ac = record.INFO.get("AC")
                    if isinstance(ac, (list, tuple)):
                        ac = ac[i] if i < len(ac) else ac[0]
                    
                    entry = GnomADEntry(
                        chromosome=record.CHROM,
                        position=record.POS,
                        ref=record.REF,
                        alt=alt,
                        allele_frequency=float(af or 0),
                        allele_count=int(ac or 0),
                        allele_number=int(record.INFO.get("AN", 0)),
                        homozygote_count=int(record.INFO.get("nhomalt", 0)),
                        filters=";".join(record.FILTER) if record.FILTER else "PASS"
                    )
                    self._cache[key] = entry
            
            vcf.close()
            self._loaded = True
            
        except ImportError:
            pass
    
    def _make_key(self, chrom: str, pos: int, ref: str, alt: str) -> str:
        """Create lookup key."""
        chrom = str(chrom).replace("chr", "")
        return f"{chrom}:{pos}:{ref}:{alt}"
    
    def annotate(
        self,
        chromosome: str,
        position: int,
        ref: str,
        alt: str
    ) -> Optional[GnomADEntry]:
        """Get gnomAD annotation for a variant."""
        if not self._loaded:
            self.load()
        
        key = self._make_key(chromosome, position, ref, alt)
        return self._cache.get(key)
    
    def get_frequency(
        self,
        chromosome: str,
        position: int,
        ref: str,
        alt: str
    ) -> Optional[float]:
        """Get allele frequency for a variant."""
        entry = self.annotate(chromosome, position, ref, alt)
        return entry.allele_frequency if entry else None
    
    def get_rare_variants(self, threshold: float = 0.01) -> List[GnomADEntry]:
        """Get variants below frequency threshold."""
        return [e for e in self._cache.values() if e.allele_frequency < threshold]
    
    def get_statistics(self) -> Dict[str, Any]:
        """Get statistics on cached variants."""
        if not self._cache:
            return {"total": 0}
        
        frequencies = [e.allele_frequency for e in self._cache.values()]
        
        return {
            "total": len(self._cache),
            "rare": sum(1 for e in self._cache.values() if e.is_rare),
            "ultra_rare": sum(1 for e in self._cache.values() if e.is_ultra_rare),
            "singletons": sum(1 for e in self._cache.values() if e.is_singleton),
            "mean_af": sum(frequencies) / len(frequencies),
            "median_af": sorted(frequencies)[len(frequencies) // 2]
        }
