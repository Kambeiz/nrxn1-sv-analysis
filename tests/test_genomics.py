"""Unit tests for genomics modules."""

import pytest
import numpy as np
import pandas as pd

from src.genomics.regions import NRXN1Region, GenomicRegion, get_exon_coordinates
from src.genomics.variants import Variant, StructuralVariant, VariantLoader, VariantFilter


class TestGenomicRegion:
    """Tests for GenomicRegion class."""
    
    def test_region_creation(self):
        """Test basic region creation."""
        region = GenomicRegion(
            chromosome="2",
            start=50000000,
            end=50100000,
            name="test_region"
        )
        
        assert region.chromosome == "2"
        assert region.start == 50000000
        assert region.end == 50100000
        assert region.length == 100000
    
    def test_region_overlaps(self):
        """Test overlap detection."""
        region1 = GenomicRegion("2", 50000000, 50100000)
        region2 = GenomicRegion("2", 50050000, 50150000)
        region3 = GenomicRegion("2", 50200000, 50300000)
        region4 = GenomicRegion("3", 50000000, 50100000)
        
        assert region1.overlaps(region2) is True
        assert region1.overlaps(region3) is False
        assert region1.overlaps(region4) is False
    
    def test_region_contains(self):
        """Test position containment."""
        region = GenomicRegion("2", 50000000, 50100000)
        
        assert region.contains(50050000) is True
        assert region.contains(50000000) is True
        assert region.contains(49999999) is False
        assert region.contains(50100000) is False
    
    def test_to_bed(self):
        """Test BED format output."""
        region = GenomicRegion("2", 50000000, 50100000, "test")
        bed = region.to_bed()
        
        assert "2\t50000000\t50100000\ttest" in bed


class TestNRXN1Region:
    """Tests for NRXN1Region class."""
    
    def test_nrxn1_constants(self):
        """Test NRXN1 region constants."""
        assert NRXN1Region.CHROMOSOME == "2"
        assert NRXN1Region.START < NRXN1Region.END
        assert len(NRXN1Region.EXONS_ALPHA) == 20
        assert len(NRXN1Region.EXONS_BETA) == 10
    
    def test_nrxn1_initialization(self):
        """Test NRXN1Region initialization."""
        nrxn1 = NRXN1Region()
        
        assert len(nrxn1.exons_alpha) == 20
        assert len(nrxn1.exons_beta) == 10
        assert len(nrxn1.domains) == 5
    
    def test_get_affected_exons(self):
        """Test finding affected exons."""
        nrxn1 = NRXN1Region()
        
        first_exon = nrxn1.EXONS_ALPHA[0]
        affected = nrxn1.get_affected_exons(
            first_exon[0],
            first_exon[1],
            "alpha"
        )
        
        assert len(affected) >= 1
        assert affected[0].name == "exon_1_alpha"
    
    def test_classify_variant_location(self):
        """Test variant location classification."""
        nrxn1 = NRXN1Region()
        
        alpha_exon = nrxn1.EXONS_ALPHA[0]
        location = nrxn1.classify_variant_location(alpha_exon[0], alpha_exon[1])
        assert location in ["alpha_specific", "both_isoforms"]
        
        intronic_start = nrxn1.EXONS_ALPHA[0][1] + 1000
        intronic_end = nrxn1.EXONS_ALPHA[1][0] - 1000
        location = nrxn1.classify_variant_location(intronic_start, intronic_end)
        assert location == "intronic"
    
    def test_to_dataframe(self):
        """Test DataFrame export."""
        nrxn1 = NRXN1Region()
        df = nrxn1.to_dataframe()
        
        assert isinstance(df, pd.DataFrame)
        assert "chromosome" in df.columns
        assert "start" in df.columns
        assert "isoform" in df.columns
        assert len(df) == 30


class TestVariant:
    """Tests for Variant class."""
    
    def test_snv_creation(self):
        """Test SNV creation."""
        snv = Variant(
            chromosome="2",
            position=50200000,
            ref="A",
            alt="G"
        )
        
        assert snv.is_snv is True
        assert snv.is_indel is False
        assert snv.variant_type == "SNV"
        assert snv.length == 0
    
    def test_deletion_creation(self):
        """Test deletion creation."""
        deletion = Variant(
            chromosome="2",
            position=50200000,
            ref="ATCG",
            alt="A"
        )
        
        assert deletion.is_snv is False
        assert deletion.is_deletion is True
        assert deletion.variant_type == "DEL"
        assert deletion.length == 3
    
    def test_insertion_creation(self):
        """Test insertion creation."""
        insertion = Variant(
            chromosome="2",
            position=50200000,
            ref="A",
            alt="ATCG"
        )
        
        assert insertion.is_insertion is True
        assert insertion.variant_type == "INS"
        assert insertion.length == 3
    
    def test_to_dict(self):
        """Test dictionary conversion."""
        var = Variant("2", 50200000, "A", "G", variant_id="test_var")
        d = var.to_dict()
        
        assert d["chromosome"] == "2"
        assert d["position"] == 50200000
        assert d["variant_type"] == "SNV"


class TestStructuralVariant:
    """Tests for StructuralVariant class."""
    
    def test_sv_creation(self):
        """Test SV creation."""
        sv = StructuralVariant(
            chromosome="2",
            position=50200000,
            ref="N",
            alt="<DEL>",
            sv_type="DEL",
            sv_end=50300000
        )
        
        assert sv.end == 50300000
        assert sv.length == 100000
    
    def test_sv_with_length(self):
        """Test SV with explicit length."""
        sv = StructuralVariant(
            chromosome="2",
            position=50200000,
            ref="N",
            alt="<DEL>",
            sv_type="DEL",
            sv_length=-50000
        )
        
        assert sv.length == 50000


class TestVariantFilter:
    """Tests for VariantFilter class."""
    
    @pytest.fixture
    def sample_variants(self):
        """Create sample variants for testing."""
        return [
            Variant("2", 50200000, "A", "G", quality=30.0),
            Variant("2", 50200100, "AT", "A", quality=25.0),
            Variant("2", 50200200, "C", "T", quality=15.0),
            Variant("2", 50200300, "G", "GACGT", quality=40.0),
        ]
    
    def test_filter_by_quality(self, sample_variants):
        """Test quality filtering."""
        vf = VariantFilter()
        filtered = vf.filter_by_quality(sample_variants, min_quality=25.0)
        
        assert len(filtered) == 3
        assert all(v.quality >= 25.0 for v in filtered)
    
    def test_filter_by_type(self, sample_variants):
        """Test type filtering."""
        vf = VariantFilter()
        snvs = vf.filter_by_type(sample_variants, ["SNV"])
        
        assert len(snvs) == 2
        assert all(v.is_snv for v in snvs)
    
    def test_apply_multiple_filters(self, sample_variants):
        """Test multiple filters."""
        vf = VariantFilter()
        filtered = vf.apply_filters(
            sample_variants,
            min_quality=20.0,
            variant_types=["SNV"],
            pass_only=False
        )
        
        assert len(filtered) == 1
        assert filtered[0].position == 50200000


def test_get_exon_coordinates():
    """Test exon coordinate retrieval."""
    alpha_coords = get_exon_coordinates("alpha")
    beta_coords = get_exon_coordinates("beta")
    
    assert len(alpha_coords) == 20
    assert len(beta_coords) == 10
    
    with pytest.raises(ValueError):
        get_exon_coordinates("invalid")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
