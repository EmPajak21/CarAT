"""Data processing components for CarAT.

The processing package provides tools for handling and transforming data:

- `PostProcessConfig` and `PreProcessConfig`: Dataclasses for pre- and post-processing.
- `DataPreprocessor`: Class takes value chain data and prepares it for LPFormulator.
"""

from .datatypes import PostProcessConfig, PreProcessConfig
from .preprocessor import DataPreprocessor
