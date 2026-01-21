"""Machine learning modules for pathogenicity prediction."""

from .features import FeatureExtractor
from .predictor import PathogenicityPredictor

__all__ = ["FeatureExtractor", "PathogenicityPredictor"]
