"""Pathogenicity prediction models."""

from dataclasses import dataclass
from typing import List, Dict, Any, Optional, Tuple, Union
from pathlib import Path
import numpy as np
import pandas as pd
import json
import pickle

from .features import FeatureExtractor, VariantFeatures
from ..genomics.variants import Variant
from ..annotation.annotator import AnnotatedVariant


@dataclass
class PredictionResult:
    """Prediction result for a variant."""
    
    pathogenicity_score: float
    classification: str
    confidence: float
    feature_contributions: Dict[str, float]
    
    @property
    def is_pathogenic(self) -> bool:
        """Check if predicted pathogenic."""
        return self.classification in ["Pathogenic", "Likely_Pathogenic"]
    
    @property
    def is_benign(self) -> bool:
        """Check if predicted benign."""
        return self.classification in ["Benign", "Likely_Benign"]
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "pathogenicity_score": self.pathogenicity_score,
            "classification": self.classification,
            "confidence": self.confidence,
            "is_pathogenic": self.is_pathogenic,
            "top_features": dict(sorted(
                self.feature_contributions.items(),
                key=lambda x: abs(x[1]),
                reverse=True
            )[:5])
        }


class PathogenicityPredictor:
    """Machine learning model for variant pathogenicity prediction."""
    
    CLASSIFICATION_THRESHOLDS = {
        "Pathogenic": 0.9,
        "Likely_Pathogenic": 0.7,
        "Uncertain_Significance": 0.3,
        "Likely_Benign": 0.1,
        "Benign": 0.0
    }
    
    def __init__(
        self,
        model_type: str = "xgboost",
        model_path: Optional[Path] = None
    ):
        """Initialize predictor.
        
        Args:
            model_type: Type of model ('xgboost', 'lightgbm', 'random_forest').
            model_path: Path to pre-trained model file.
        """
        self.model_type = model_type
        self.model_path = model_path
        self.model = None
        self.feature_extractor = FeatureExtractor()
        self.is_trained = False
        self._feature_names = VariantFeatures.feature_names()
    
    def train(
        self,
        X: np.ndarray,
        y: np.ndarray,
        validation_split: float = 0.2,
        **kwargs
    ) -> Dict[str, float]:
        """Train the pathogenicity prediction model.
        
        Args:
            X: Feature matrix (n_samples, n_features).
            y: Labels (0=benign, 1=pathogenic).
            validation_split: Fraction for validation.
            **kwargs: Additional model parameters.
        
        Returns:
            Dictionary with training metrics.
        """
        from sklearn.model_selection import train_test_split
        from sklearn.metrics import (
            accuracy_score, precision_score, recall_score,
            f1_score, roc_auc_score
        )
        
        X_train, X_val, y_train, y_val = train_test_split(
            X, y, test_size=validation_split, random_state=42, stratify=y
        )
        
        if self.model_type == "xgboost":
            self.model = self._create_xgboost_model(**kwargs)
        elif self.model_type == "lightgbm":
            self.model = self._create_lightgbm_model(**kwargs)
        elif self.model_type == "random_forest":
            self.model = self._create_random_forest_model(**kwargs)
        else:
            raise ValueError(f"Unknown model type: {self.model_type}")
        
        self.model.fit(X_train, y_train)
        
        y_pred = self.model.predict(X_val)
        y_proba = self.model.predict_proba(X_val)[:, 1]
        
        metrics = {
            "accuracy": accuracy_score(y_val, y_pred),
            "precision": precision_score(y_val, y_pred, zero_division=0),
            "recall": recall_score(y_val, y_pred, zero_division=0),
            "f1_score": f1_score(y_val, y_pred, zero_division=0),
            "roc_auc": roc_auc_score(y_val, y_proba) if len(np.unique(y_val)) > 1 else 0.0,
            "n_train": len(y_train),
            "n_val": len(y_val)
        }
        
        self.is_trained = True
        return metrics
    
    def _create_xgboost_model(self, **kwargs):
        """Create XGBoost classifier."""
        try:
            from xgboost import XGBClassifier
            
            default_params = {
                "n_estimators": 100,
                "max_depth": 6,
                "learning_rate": 0.1,
                "subsample": 0.8,
                "colsample_bytree": 0.8,
                "random_state": 42,
                "eval_metric": "logloss",
                "use_label_encoder": False
            }
            default_params.update(kwargs)
            
            return XGBClassifier(**default_params)
        except ImportError:
            raise ImportError("XGBoost required. Install with: pip install xgboost")
    
    def _create_lightgbm_model(self, **kwargs):
        """Create LightGBM classifier."""
        try:
            from lightgbm import LGBMClassifier
            
            default_params = {
                "n_estimators": 100,
                "max_depth": 6,
                "learning_rate": 0.1,
                "subsample": 0.8,
                "colsample_bytree": 0.8,
                "random_state": 42,
                "verbose": -1
            }
            default_params.update(kwargs)
            
            return LGBMClassifier(**default_params)
        except ImportError:
            raise ImportError("LightGBM required. Install with: pip install lightgbm")
    
    def _create_random_forest_model(self, **kwargs):
        """Create Random Forest classifier."""
        from sklearn.ensemble import RandomForestClassifier
        
        default_params = {
            "n_estimators": 100,
            "max_depth": 10,
            "min_samples_split": 5,
            "min_samples_leaf": 2,
            "random_state": 42,
            "n_jobs": -1
        }
        default_params.update(kwargs)
        
        return RandomForestClassifier(**default_params)
    
    def predict(
        self,
        variant: Variant,
        annotated: Optional[AnnotatedVariant] = None
    ) -> PredictionResult:
        """Predict pathogenicity for a single variant.
        
        Args:
            variant: Variant to predict.
            annotated: Optional annotated variant.
        
        Returns:
            PredictionResult object.
        """
        if not self.is_trained and self.model is None:
            raise RuntimeError("Model not trained. Call train() or load() first.")
        
        features = self.feature_extractor.extract_features(variant, annotated)
        X = features.to_array().reshape(1, -1)
        
        proba = self.model.predict_proba(X)[0, 1]
        
        classification = self._score_to_classification(proba)
        
        confidence = abs(proba - 0.5) * 2
        
        feature_contributions = self._get_feature_contributions(X[0])
        
        return PredictionResult(
            pathogenicity_score=proba,
            classification=classification,
            confidence=confidence,
            feature_contributions=feature_contributions
        )
    
    def predict_batch(
        self,
        variants: List[Variant],
        annotated_variants: Optional[List[AnnotatedVariant]] = None
    ) -> List[PredictionResult]:
        """Predict pathogenicity for multiple variants."""
        if annotated_variants is None:
            annotated_variants = [None] * len(variants)
        
        return [
            self.predict(var, ann)
            for var, ann in zip(variants, annotated_variants)
        ]
    
    def _score_to_classification(self, score: float) -> str:
        """Convert probability score to ACMG-like classification."""
        if score >= self.CLASSIFICATION_THRESHOLDS["Pathogenic"]:
            return "Pathogenic"
        elif score >= self.CLASSIFICATION_THRESHOLDS["Likely_Pathogenic"]:
            return "Likely_Pathogenic"
        elif score >= self.CLASSIFICATION_THRESHOLDS["Uncertain_Significance"]:
            return "Uncertain_Significance"
        elif score >= self.CLASSIFICATION_THRESHOLDS["Likely_Benign"]:
            return "Likely_Benign"
        else:
            return "Benign"
    
    def _get_feature_contributions(self, X: np.ndarray) -> Dict[str, float]:
        """Get feature contributions using SHAP or feature importance."""
        contributions = {}
        
        if hasattr(self.model, 'feature_importances_'):
            importances = self.model.feature_importances_
            for name, importance, value in zip(
                self._feature_names, importances, X
            ):
                contributions[name] = float(importance * value)
        
        return contributions
    
    def cross_validate(
        self,
        X: np.ndarray,
        y: np.ndarray,
        n_folds: int = 5
    ) -> Dict[str, Any]:
        """Perform cross-validation.
        
        Returns:
            Dictionary with CV metrics.
        """
        from sklearn.model_selection import StratifiedKFold, cross_val_score
        
        if self.model is None:
            if self.model_type == "xgboost":
                self.model = self._create_xgboost_model()
            elif self.model_type == "lightgbm":
                self.model = self._create_lightgbm_model()
            else:
                self.model = self._create_random_forest_model()
        
        cv = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=42)
        
        scores = {
            "accuracy": cross_val_score(self.model, X, y, cv=cv, scoring="accuracy"),
            "precision": cross_val_score(self.model, X, y, cv=cv, scoring="precision"),
            "recall": cross_val_score(self.model, X, y, cv=cv, scoring="recall"),
            "f1": cross_val_score(self.model, X, y, cv=cv, scoring="f1"),
            "roc_auc": cross_val_score(self.model, X, y, cv=cv, scoring="roc_auc")
        }
        
        return {
            metric: {
                "mean": float(np.mean(values)),
                "std": float(np.std(values)),
                "values": values.tolist()
            }
            for metric, values in scores.items()
        }
    
    def get_feature_importance(self) -> pd.DataFrame:
        """Get feature importance ranking."""
        if self.model is None or not hasattr(self.model, 'feature_importances_'):
            return pd.DataFrame()
        
        importance_df = pd.DataFrame({
            "feature": self._feature_names,
            "importance": self.model.feature_importances_
        })
        
        return importance_df.sort_values("importance", ascending=False)
    
    def save(self, path: Path) -> None:
        """Save model to file."""
        model_data = {
            "model_type": self.model_type,
            "feature_names": self._feature_names,
            "thresholds": self.CLASSIFICATION_THRESHOLDS
        }
        
        with open(path.with_suffix(".json"), "w") as f:
            json.dump(model_data, f, indent=2)
        
        with open(path.with_suffix(".pkl"), "wb") as f:
            pickle.dump(self.model, f)
    
    def load(self, path: Path) -> None:
        """Load model from file."""
        with open(path.with_suffix(".json"), "r") as f:
            model_data = json.load(f)
        
        self.model_type = model_data["model_type"]
        self._feature_names = model_data["feature_names"]
        
        with open(path.with_suffix(".pkl"), "rb") as f:
            self.model = pickle.load(f)
        
        self.is_trained = True
    
    def explain_prediction(
        self,
        variant: Variant,
        annotated: Optional[AnnotatedVariant] = None
    ) -> Dict[str, Any]:
        """Generate detailed explanation for a prediction."""
        result = self.predict(variant, annotated)
        features = self.feature_extractor.extract_features(variant, annotated)
        
        explanation = {
            "prediction": result.to_dict(),
            "features": features.to_dict(),
            "reasoning": []
        }
        
        if features.exons_affected > 0:
            explanation["reasoning"].append(
                f"Variant affects {features.exons_affected} exon(s), "
                "increasing pathogenicity likelihood."
            )
        
        if features.is_absent_gnomad:
            explanation["reasoning"].append(
                "Variant absent from gnomAD population database, "
                "suggesting it may be deleterious."
            )
        
        if features.affects_both_isoforms:
            explanation["reasoning"].append(
                "Variant affects both alpha and beta isoforms, "
                "potentially disrupting multiple protein functions."
            )
        
        if features.cadd_score > 20:
            explanation["reasoning"].append(
                f"High CADD score ({features.cadd_score:.1f}) suggests "
                "strong deleteriousness."
            )
        
        return explanation
