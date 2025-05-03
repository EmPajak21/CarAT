"""Visualization components for CarAT.

The vis package provides tools for creating visual representations CarAT outputs:

- `SankeyDiagramGenerator`: Generates Sankey diagrams to visualize carbon flows.
- `mermaid_plot`: Creates value chain graphs using Mermaid.js syntax.
"""

from .sankey import SankeyDiagramGenerator
from .vc_graph import mermaid_plot
