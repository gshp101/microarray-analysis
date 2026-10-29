"""Find metabolic pathways affected in public gene expression microarray studies."""

__version__ = "0.2.0"

from .geo import ExpressionDataset, load_series

__all__ = ["ExpressionDataset", "load_series"]
