"""Variant annotation modules."""

from .annotator import VariantAnnotator, AnnotatedVariant
from .databases import ClinVarAnnotator, GnomADAnnotator

__all__ = ["VariantAnnotator", "AnnotatedVariant", "ClinVarAnnotator", "GnomADAnnotator"]
