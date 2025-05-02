import networkx as nx
import base64
import io
import re
import requests
from typing import Tuple, Optional
from IPython.display import Image, display
from PIL import Image as PILImage
import matplotlib.pyplot as plt


def nx_to_mermaid(G: nx.Graph) -> str:
    """
    Convert a NetworkX bipartite graph to Mermaid graph syntax.

    Parameters:
    G : NetworkX Graph
        A bipartite graph where nodes have a 'bipartite' attribute (0 or 1)
        indicating which set they belong to.

    Returns:
    str : Mermaid graph representation
    """

    # Get the two sets of nodes
    duplets = {n for n, d in G.nodes(data=True) if n.startswith("d")}

    # Start the Mermaid graph definition
    mermaid = "graph LR\n"

    # Add nodes to the graph definition
    for node, data in G.nodes(data=True):
        node_text = data.get(
            "label", "tr"
        )  # Default to node ID if 'txt' is not available
        if node in duplets:
            mermaid += f"    {node}([{node_text}\n{node}])\n"
        else:
            mermaid += f"    {node}({node_text}\n{node})\n"

    # Add edges
    for u, v in G.edges():
        mermaid += f"    {u} --> {v}\n"

    return mermaid


def mm(
    graph: str,
    output_filename: str = "mermaid_diagram.png",
    figsize: Tuple[int, int] = (12, 10),
    dpi: int = 300,
) -> bool:
    """
    Renders a Mermaid diagram using the mermaid.ink service and saves it as an image.

    Args:
        graph: The Mermaid diagram code as a string
        output_filename: Filename to save the rendered diagram (default: 'mermaid_diagram.png')
        figsize: Figure size as (width, height) in inches (default: (12, 10))
        dpi: Resolution of the output image (default: 300)

    Returns:
        bool: True if the diagram was successfully rendered and saved, False otherwise
    """
    # Clean up the input graph - strip any leading/trailing whitespace
    graph = graph.strip()

    # URL-safe base64 encoding
    graphbytes = graph.encode("utf8")
    base64_bytes = base64.b64encode(graphbytes)
    base64_string = base64_bytes.decode("ascii")

    # Make the base64 string URL-safe
    base64_string = base64_string.replace("+", "-").replace("/", "_").rstrip("=")

    # Fetch the image from mermaid.ink
    url = "https://mermaid.ink/img/" + base64_string
    response = requests.get(url)

    if response.status_code == 200:
        img = PILImage.open(io.BytesIO(response.content))

        # Solution: Use only matplotlib for rendering and display
        plt.figure(figsize=figsize)
        plt.imshow(img)
        plt.axis("off")
        plt.savefig(output_filename, dpi=1000, bbox_inches="tight")
        plt.close()  # Close the figure to prevent double display

        # Use IPython's display for showing the saved image
        display(Image(output_filename))
        return True
    else:
        print(f"Error fetching diagram: {response.status_code}")
        print(f"URL attempted: {url}")
        return False


def mermaid_plot(
    graph,
    output_filename: Optional[str] = None,
    figsize: Tuple[int, int] = (12, 10),
    dpi: int = 300,
) -> bool:
    """
    Converts a networkx graph to mermaid code and renders it.

    Args:
        graph: The networkx graph to convert and render
        output_filename: Optional custom filename for the output image
                        (default: 'mermaid_diagram.png')
        figsize: Figure size as (width, height) in inches (default: (12, 10))
        dpi: Resolution of the output image (default: 300)

    Returns:
        bool: True if the diagram was successfully rendered and saved, False otherwise
    """
    # Convert networkx graph to mermaid code
    mermaid_code = nx_to_mermaid(graph)

    # If no output filename specified, use default
    if output_filename is None:
        output_filename = "mermaid_diagram.png"

    # Generate the diagram
    mm(mermaid_code, output_filename, figsize, dpi)
