"""Unit tests for CNV detection modules."""

import pytest
import numpy as np
import pandas as pd

from src.cnv.detector import CNVDetector, CNVCall, CNVType
from src.cnv.consensus import ConsensusBuilder, ConsensusCNV
from src.genomics.regions import NRXN1Region


class TestCNVType:
    """Tests for CNVType enum."""
    
    def test_cnv_types(self):
        """Test CNV type values."""
        assert CNVType.DELETION.value == "DEL"
        assert CNVType.DUPLICATION.value == "DUP"
        assert CNVType.INVERSION.value == "INV"


class TestCNVCall:
    """Tests for CNVCall class."""
    
    def test_cnv_call_creation(self):
        """Test CNV call creation."""
        cnv = CNVCall(
            chromosome="2",
            start=50200000,
            end=50300000,
            cnv_type=CNVType.DELETION,
            caller="manta",
            quality=50.0
        )
        
        assert cnv.chromosome == "2"
        assert cnv.length == 100000
        assert cnv.is_deletion is True
        assert cnv.is_duplication is False
    
    def test_cnv_overlap(self):
        """Test CNV overlap detection."""
        cnv1 = CNVCall("2", 50200000, 50300000, CNVType.DELETION, "manta")
        cnv2 = CNVCall("2", 50220000, 50320000, CNVType.DELETION, "delly")
        cnv3 = CNVCall("2", 50400000, 50500000, CNVType.DELETION, "manta")
        
        assert cnv1.overlaps(cnv2, min_reciprocal=0.5) is True
        assert cnv1.overlaps(cnv3, min_reciprocal=0.5) is False
    
    def test_get_affected_exons(self):
        """Test exon impact detection."""
        nrxn1 = NRXN1Region()
        first_exon = nrxn1.EXONS_ALPHA[0]
        
        cnv = CNVCall(
            chromosome="2",
            start=first_exon[0] - 100,
            end=first_exon[1] + 100,
            cnv_type=CNVType.DELETION,
            caller="test"
        )
        
        affected = cnv.get_affected_exons(nrxn1)
        assert len(affected) >= 1
    
    def test_clinical_significance(self):
        """Test clinical significance assessment."""
        nrxn1 = NRXN1Region()
        
        multi_exon_del = CNVCall(
            chromosome="2",
            start=nrxn1.EXONS_ALPHA[0][0],
            end=nrxn1.EXONS_ALPHA[4][1],
            cnv_type=CNVType.DELETION,
            caller="test"
        )
        
        sig = multi_exon_del.get_clinical_significance()
        assert sig in ["likely_pathogenic", "uncertain_significance"]
    
    def test_to_dict(self):
        """Test dictionary conversion."""
        cnv = CNVCall(
            chromosome="2",
            start=50200000,
            end=50300000,
            cnv_type=CNVType.DELETION,
            caller="manta",
            quality=50.0,
            sample_id="sample1"
        )
        
        d = cnv.to_dict()
        assert d["chromosome"] == "2"
        assert d["length"] == 100000
        assert d["cnv_type"] == "DEL"
    
    def test_to_bed(self):
        """Test BED format output."""
        cnv = CNVCall(
            chromosome="2",
            start=50200000,
            end=50300000,
            cnv_type=CNVType.DELETION,
            caller="manta",
            quality=50.0
        )
        
        bed = cnv.to_bed()
        assert "2\t50200000\t50300000\tDEL\t50.0" == bed


class TestCNVDetector:
    """Tests for CNVDetector class."""
    
    def test_detector_initialization(self):
        """Test detector initialization."""
        detector = CNVDetector(
            min_size=1000,
            max_size=1000000,
            min_quality=20.0
        )
        
        assert detector.min_size == 1000
        assert detector.max_size == 1000000
        assert detector.nrxn1_only is True
    
    def test_detect_from_depth(self):
        """Test depth-based CNV detection."""
        detector = CNVDetector(min_size=100, nrxn1_only=False)
        
        positions = np.arange(0, 10000, 10)
        depth = np.ones(len(positions)) * 30
        depth[300:350] = 5
        
        calls = detector.detect_from_depth(
            depth,
            positions,
            window_size=100,
            z_threshold=2.0
        )
        
        assert isinstance(calls, list)
    
    def test_estimate_copy_number(self):
        """Test copy number estimation."""
        detector = CNVDetector()
        
        assert detector._estimate_copy_number(-1.0) == 1
        assert detector._estimate_copy_number(0.0) == 2
        assert detector._estimate_copy_number(0.58) == 3
    
    def test_filter_calls(self):
        """Test CNV call filtering."""
        detector = CNVDetector()
        
        calls = [
            CNVCall("2", 50200000, 50300000, CNVType.DELETION, "manta", 50.0),
            CNVCall("2", 50400000, 50500000, CNVType.DUPLICATION, "manta", 30.0),
            CNVCall("2", 50600000, 50700000, CNVType.DELETION, "manta", 10.0),
        ]
        
        filtered = detector.filter_calls(calls, min_quality=25.0)
        assert len(filtered) == 2
        
        dels_only = detector.filter_calls(calls, cnv_types=[CNVType.DELETION])
        assert len(dels_only) == 2
    
    def test_to_dataframe(self):
        """Test DataFrame conversion."""
        detector = CNVDetector()
        
        calls = [
            CNVCall("2", 50200000, 50300000, CNVType.DELETION, "manta", 50.0),
            CNVCall("2", 50400000, 50500000, CNVType.DUPLICATION, "manta", 30.0),
        ]
        
        df = detector.to_dataframe(calls)
        
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 2
        assert "cnv_type" in df.columns


class TestConsensusBuilder:
    """Tests for ConsensusBuilder class."""
    
    @pytest.fixture
    def sample_calls(self):
        """Create sample calls from multiple callers."""
        return {
            "manta": [
                CNVCall("2", 50200000, 50300000, CNVType.DELETION, "manta", 50.0),
                CNVCall("2", 50500000, 50600000, CNVType.DUPLICATION, "manta", 40.0),
            ],
            "delly": [
                CNVCall("2", 50210000, 50290000, CNVType.DELETION, "delly", 45.0),
                CNVCall("2", 50800000, 50900000, CNVType.DELETION, "delly", 35.0),
            ],
            "cnvnator": [
                CNVCall("2", 50205000, 50295000, CNVType.DELETION, "cnvnator", 55.0),
            ]
        }
    
    def test_builder_initialization(self):
        """Test consensus builder initialization."""
        builder = ConsensusBuilder(
            min_reciprocal_overlap=0.5,
            min_callers=2
        )
        
        assert builder.min_reciprocal_overlap == 0.5
        assert builder.min_callers == 2
    
    def test_build_consensus(self, sample_calls):
        """Test consensus building."""
        builder = ConsensusBuilder(min_callers=2)
        consensus = builder.build_consensus(sample_calls)
        
        assert isinstance(consensus, list)
        assert len(consensus) >= 1
        
        for c in consensus:
            assert c.num_callers >= 2
    
    def test_calls_overlap(self):
        """Test overlap checking."""
        builder = ConsensusBuilder()
        
        call1 = CNVCall("2", 50200000, 50300000, CNVType.DELETION, "manta")
        call2 = CNVCall("2", 50220000, 50320000, CNVType.DELETION, "delly")
        call3 = CNVCall("2", 50400000, 50500000, CNVType.DELETION, "manta")
        
        assert builder._calls_overlap(call1, call2) is True
        assert builder._calls_overlap(call1, call3) is False
    
    def test_summarize(self, sample_calls):
        """Test summary statistics."""
        builder = ConsensusBuilder(min_callers=2)
        consensus = builder.build_consensus(sample_calls)
        summary = builder.summarize(consensus)
        
        assert "total_calls" in summary
        assert "deletions" in summary
        assert "duplications" in summary
        assert "mean_size" in summary
    
    def test_to_dataframe(self, sample_calls):
        """Test DataFrame export."""
        builder = ConsensusBuilder(min_callers=2)
        consensus = builder.build_consensus(sample_calls)
        df = builder.to_dataframe(consensus)
        
        assert isinstance(df, pd.DataFrame)
        if len(consensus) > 0:
            assert "num_callers" in df.columns


class TestConsensusCNV:
    """Tests for ConsensusCNV class."""
    
    def test_consensus_cnv_creation(self):
        """Test consensus CNV creation."""
        merged_calls = [
            CNVCall("2", 50200000, 50300000, CNVType.DELETION, "manta"),
            CNVCall("2", 50210000, 50290000, CNVType.DELETION, "delly"),
        ]
        
        consensus = ConsensusCNV(
            chromosome="2",
            start=50205000,
            end=50295000,
            cnv_type=CNVType.DELETION,
            supporting_callers=["manta", "delly"],
            merged_calls=merged_calls,
            consensus_quality=50.0
        )
        
        assert consensus.num_callers == 2
        assert consensus.length == 90000
        assert "delly" in consensus.caller_string
        assert "manta" in consensus.caller_string


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
