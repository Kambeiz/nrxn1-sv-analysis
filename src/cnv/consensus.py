"""Consensus CNV calling from multiple callers."""

from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
import numpy as np
import pandas as pd

from .detector import CNVCall, CNVType


@dataclass
class ConsensusCNV:
    """Consensus CNV call from multiple callers."""
    
    chromosome: str
    start: int
    end: int
    cnv_type: CNVType
    supporting_callers: List[str]
    merged_calls: List[CNVCall]
    consensus_quality: float
    
    @property
    def length(self) -> int:
        """Return consensus CNV length."""
        return self.end - self.start
    
    @property
    def num_callers(self) -> int:
        """Return number of supporting callers."""
        return len(self.supporting_callers)
    
    @property
    def caller_string(self) -> str:
        """Return comma-separated list of callers."""
        return ",".join(sorted(self.supporting_callers))
    
    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return {
            "chromosome": self.chromosome,
            "start": self.start,
            "end": self.end,
            "length": self.length,
            "cnv_type": self.cnv_type.value,
            "num_callers": self.num_callers,
            "supporting_callers": self.caller_string,
            "consensus_quality": self.consensus_quality
        }


class ConsensusBuilder:
    """Build consensus CNV calls from multiple callers."""
    
    def __init__(
        self,
        min_reciprocal_overlap: float = 0.5,
        min_callers: int = 2,
        same_type_required: bool = True
    ):
        """Initialize consensus builder.
        
        Args:
            min_reciprocal_overlap: Minimum reciprocal overlap to merge calls.
            min_callers: Minimum number of callers required for consensus.
            same_type_required: Require same CNV type to merge.
        """
        self.min_reciprocal_overlap = min_reciprocal_overlap
        self.min_callers = min_callers
        self.same_type_required = same_type_required
    
    def build_consensus(
        self,
        caller_results: Dict[str, List[CNVCall]]
    ) -> List[ConsensusCNV]:
        """Build consensus calls from multiple caller results.
        
        Args:
            caller_results: Dictionary mapping caller name to list of calls.
        
        Returns:
            List of consensus CNV calls.
        """
        all_calls = []
        for caller, calls in caller_results.items():
            for call in calls:
                call.caller = caller
                all_calls.append(call)
        
        if not all_calls:
            return []
        
        all_calls.sort(key=lambda c: (c.chromosome, c.start))
        
        clusters = self._cluster_overlapping_calls(all_calls)
        
        consensus_calls = []
        for cluster in clusters:
            if len(set(c.caller for c in cluster)) >= self.min_callers:
                consensus = self._merge_cluster(cluster)
                if consensus:
                    consensus_calls.append(consensus)
        
        return consensus_calls
    
    def _cluster_overlapping_calls(
        self,
        calls: List[CNVCall]
    ) -> List[List[CNVCall]]:
        """Cluster overlapping CNV calls."""
        if not calls:
            return []
        
        clusters = []
        current_cluster = [calls[0]]
        
        for call in calls[1:]:
            overlaps_cluster = False
            
            for cluster_call in current_cluster:
                if self._calls_overlap(call, cluster_call):
                    if not self.same_type_required or call.cnv_type == cluster_call.cnv_type:
                        overlaps_cluster = True
                        break
            
            if overlaps_cluster:
                current_cluster.append(call)
            else:
                if current_cluster:
                    clusters.append(current_cluster)
                current_cluster = [call]
        
        if current_cluster:
            clusters.append(current_cluster)
        
        return clusters
    
    def _calls_overlap(self, call1: CNVCall, call2: CNVCall) -> bool:
        """Check if two calls overlap with minimum reciprocal overlap."""
        if call1.chromosome != call2.chromosome:
            return False
        
        overlap_start = max(call1.start, call2.start)
        overlap_end = min(call1.end, call2.end)
        
        if overlap_start >= overlap_end:
            return False
        
        overlap_length = overlap_end - overlap_start
        overlap1 = overlap_length / call1.length if call1.length > 0 else 0
        overlap2 = overlap_length / call2.length if call2.length > 0 else 0
        
        return (overlap1 >= self.min_reciprocal_overlap and 
                overlap2 >= self.min_reciprocal_overlap)
    
    def _merge_cluster(self, cluster: List[CNVCall]) -> Optional[ConsensusCNV]:
        """Merge a cluster of calls into a consensus call."""
        if not cluster:
            return None
        
        callers = list(set(c.caller for c in cluster))
        
        starts = [c.start for c in cluster]
        ends = [c.end for c in cluster]
        consensus_start = int(np.median(starts))
        consensus_end = int(np.median(ends))
        
        type_counts = {}
        for c in cluster:
            type_counts[c.cnv_type] = type_counts.get(c.cnv_type, 0) + 1
        consensus_type = max(type_counts, key=type_counts.get)
        
        qualities = [c.quality for c in cluster if c.quality > 0]
        consensus_quality = np.mean(qualities) if qualities else 0.0
        
        consensus_quality += len(callers) * 10
        
        return ConsensusCNV(
            chromosome=cluster[0].chromosome,
            start=consensus_start,
            end=consensus_end,
            cnv_type=consensus_type,
            supporting_callers=callers,
            merged_calls=cluster,
            consensus_quality=consensus_quality
        )
    
    def to_dataframe(self, consensus_calls: List[ConsensusCNV]) -> pd.DataFrame:
        """Convert consensus calls to DataFrame."""
        if not consensus_calls:
            return pd.DataFrame()
        return pd.DataFrame([c.to_dict() for c in consensus_calls])
    
    def summarize(self, consensus_calls: List[ConsensusCNV]) -> Dict:
        """Generate summary statistics for consensus calls."""
        if not consensus_calls:
            return {
                "total_calls": 0,
                "deletions": 0,
                "duplications": 0,
                "mean_size": 0,
                "mean_callers": 0
            }
        
        deletions = sum(1 for c in consensus_calls if c.cnv_type == CNVType.DELETION)
        duplications = sum(1 for c in consensus_calls if c.cnv_type == CNVType.DUPLICATION)
        
        return {
            "total_calls": len(consensus_calls),
            "deletions": deletions,
            "duplications": duplications,
            "mean_size": np.mean([c.length for c in consensus_calls]),
            "median_size": np.median([c.length for c in consensus_calls]),
            "mean_callers": np.mean([c.num_callers for c in consensus_calls]),
            "size_range": (
                min(c.length for c in consensus_calls),
                max(c.length for c in consensus_calls)
            )
        }
