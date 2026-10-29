"""Find metabolic pathways affected in public gene expression microarray studies."""

from .geo import ExpressionDataset, load_series

__all__ = ["ExpressionDataset", "load_series"]
