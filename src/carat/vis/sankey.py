"""Sankey diagram generation for CarAT.

Main class:
- SankeyDiagramGenerator: create and render Sankey diagrams from graph data.
"""

import os
import webbrowser

import numpy as np
import pandas as pd
import plotly.graph_objects as go


class SankeyDiagramGenerator:
    """Generate Sankey diagrams from graph data with biogenic content information."""

    def __init__(self, graph, beta_d_res, beta_t_res):
        """
        Initialize the Sankey diagram generator.

        Args:
            graph_path (str): Path to the pickled graph file
            beta_d_res (pd.DataFrame): DataFrame with duplet data
            beta_t_res (pd.DataFrame): DataFrame with triplet data
        """
        self.graph = graph
        self.beta_d_res = beta_d_res
        self.beta_t_res = beta_t_res
        self.sub_node_info = {}
        self.node_labels = {}

    def prepare_data(self):
        """
        Prepare the data for the Sankey diagram.

        This includes creating biogenic content dictionaries and preparing
        the data structure for the visualization.
        """
        # Create dictionaries for biogenic content by duplet and triplet
        self.beta_d_res["DUPLET"] = (
            "d:" + self.beta_d_res["COCD"] + "," + self.beta_d_res["PBG"]
        )
        self.beta_t_res["TRIPLET"] = (
            "t:"
            + self.beta_t_res["COCD"]
            + ","
            + self.beta_t_res["BP"]
            + ","
            + self.beta_t_res["PBG_GR"]
        )

        # Get biogenic content for duplets
        self._process_biogenic_content(
            self.beta_d_res[self.beta_d_res["ATTRIBUTE"] == "biogenic"][
                ["SMILES", "OPTIMAL", "DUPLET"]
            ],
            "DUPLET",
        )

        # Get biogenic content for triplets
        self._process_biogenic_content(
            self.beta_t_res[self.beta_t_res["ATTRIBUTE"] == "biogenic"][
                ["SMILES", "OPTIMAL", "TRIPLET"]
            ],
            "TRIPLET",
        )

        # Set node labels
        self.node_labels = {
            node: self.graph.nodes[node].get("label", node)
            for node in self.graph.nodes()
        }

    def _process_biogenic_content(self, df, key_column):
        """
        Process biogenic content data from a DataFrame.

        Args:
            df (pd.DataFrame): DataFrame containing biogenic content data
            key_column (str): Column name for the node identifier
        """
        for _, row in df.iterrows():
            node_key = row[key_column]
            smiles = row["SMILES"]
            optimal = row["OPTIMAL"]

            # Create the outer dictionary key if it doesn't exist
            if node_key not in self.sub_node_info:
                self.sub_node_info[node_key] = {}

            # Add the inner dictionary key-value pair
            self.sub_node_info[node_key][smiles] = optimal

    def create_edge_dataframe(self):
        """
        Create a DataFrame from graph edges with additional biogenic content.

        Returns:
            pd.DataFrame: Processed DataFrame with expanded substream data
        """
        # Create DataFrame from edges
        edges = list(self.graph.edges(data=True))
        df_trial = pd.DataFrame(edges, columns=["source", "target", "data"])

        # Use random values for ratio (replace with actual ratios if available)
        df_trial["ratio"] = np.random.rand(len(df_trial))

        # Drop the original 'data' column if no longer needed
        df_trial = df_trial.drop("data", axis=1)

        # Create 'substreams' column with expanded data
        df_trial["substreams"] = df_trial["source"].apply(self._expand_substreams)

        # Explode 'substreams' column to create separate rows
        df_trial = df_trial.explode("substreams").reset_index(drop=True)

        # Extract data from 'substreams' column into separate columns
        df_trial = pd.concat(
            [
                df_trial.drop(["substreams"], axis=1),
                df_trial["substreams"].apply(pd.Series),
            ],
            axis=1,
        )

        # Aggregate ratios for each source-target-substream combination
        df_agg = df_trial.groupby(
            ["source", "target", "substream"], as_index=False
        ).agg(
            {
                "ratio": "sum",
                "biogenic content": "first",
            }
        )

        return df_agg

    def _expand_substreams(self, node):
        """
        Expand substreams and values for a node.
        Returns a default "No carbon data" entry if no substreams are found.

        Args:
            node (str): Node identifier

        Returns:
            list: List of dictionaries with substream data
        """
        if node in self.sub_node_info and self.sub_node_info[node]:
            substreams_data = self.sub_node_info[node]
            return [
                {"substream": substream, "biogenic content": value}
                for substream, value in substreams_data.items()
            ]
        else:
            # Return a default "No carbon data" entry to ensure the node stays connected
            return [{"substream": "No carbon data", "biogenic content": 0.0}]

    def generate_sankey_diagram(
        self, title="Biogenic Carbon Content in TDI Value Chain", width=1000, height=600
    ):
        """
        Generate a Sankey diagram visualization.

        Args:
            title (str): Title for the diagram
            width (int): Width of the diagram in pixels
            height (int): Height of the diagram in pixels

        Returns:
            plotly.graph_objects.Figure: The generated Sankey diagram figure
        """
        # Create edge dataframe
        df_agg = self.create_edge_dataframe()

        # Initialize nodes
        nodes = list(set(df_agg["source"]).union(set(df_agg["target"])))

        # Create indices for the nodes
        node_indices = {node: i for i, node in enumerate(nodes)}

        # Prepare the Sankey diagram data
        links = self._prepare_links(df_agg, node_indices)

        # Update node labels list
        node_labels_list = [self.node_labels.get(node, node) for node in nodes]

        # Create colors for nodes based on starting character
        node_colors = [
            "lightblue" if node.startswith("t") else "darkblue" for node in nodes
        ]

        # Create the Sankey diagram with conditional hover templates
        fig = go.Figure(
            data=[
                go.Sankey(
                    node=dict(
                        pad=10,
                        thickness=10,
                        line=dict(color="black", width=0.5),
                        label=node_labels_list,
                        color=node_colors,
                    ),
                    link=dict(
                        source=[link["source"] for link in links],
                        target=[link["target"] for link in links],
                        value=[link["value"] for link in links],
                        label=[link["label"] for link in links],
                        hoverinfo="all",
                        hoverlabel=dict(bgcolor="white"),
                        hovertemplate=(
                            "Substream: %{label}<br>"
                            "Biogenic Content: %{customdata:.2f}<extra></extra>"
                        ),
                        color=[link["color"] for link in links],
                        customdata=[link["biogenic content"] for link in links],
                    ),
                    visible=True,
                )
            ]
        )

        # Update layout
        fig.update_layout(
            title=title,
            font=dict(size=14, color="black"),
            width=width,
            height=height,
        )

        return fig

    def _prepare_links(self, df_agg, node_indices):
        """
        Prepare links data for the Sankey diagram.

        Args:
            df_agg (pd.DataFrame): Aggregated DataFrame with edge data
            node_indices (dict): Dictionary mapping node names to indices

        Returns:
            list: List of link dictionaries for the Sankey diagram
        """
        links = []
        for _, row in df_agg.iterrows():
            source = row["source"]
            target = row["target"]
            substream = row["substream"]
            ratio = row["ratio"]
            biogenic_content = row["biogenic content"]
            source_idx = node_indices[source]
            target_idx = node_indices[target]

            # Set color based on the type of stream
            if substream == "No carbon data":
                # Pale yellow for non-carbon streams
                color = "rgba(255, 255, 200, 0.8)"  # Pale yellow with fixed opacity
            elif biogenic_content == 0:
                # Grey for zero biogenic content streams that aren't "No carbon data"
                color = "rgba(128, 128, 128, 0.8)"  # Grey with fixed opacity
            else:
                # Light green with calculated opacity for streams with biogenic content
                color = f"rgba(144, 238, 144, {biogenic_content**3 + 0.1:.2f})"

            # Skip links involving nodes with 'N/A' text
            if (
                self.node_labels.get(source) == "N/A"
                or self.node_labels.get(target) == "N/A"
            ):
                continue

            # Add link for each substream
            links.append(
                {
                    "source": source_idx,
                    "target": target_idx,
                    "value": ratio,
                    "label": substream,
                    "color": color,
                    "biogenic content": biogenic_content,
                }
            )

        return links

    def save_diagram(self, fig, filename="sankey_diagram.html", scale=2):
        """
        Save the Sankey diagram to a file.

        Args:
            fig (plotly.graph_objects.Figure): The Sankey diagram figure
            filename (str): Output filename
            scale (int): Resolution scale factor
        """
        fig.write_html(filename)
        print(f"Diagram saved as {filename}")

    def display_diagram(self, fig):
        """
        Display the Sankey diagram.

        Args:
            fig (plotly.graph_objects.Figure): The Sankey diagram figure
        """
        fig.show()

    @staticmethod
    def view_saved_diagram(html_file_path):
        """
        Open a previously saved Sankey diagram HTML file in the default web browser.

        Args:
            html_file_path (str): Path to the saved HTML file

        Returns:
            bool: True if successful, False otherwise
        """
        if not os.path.exists(html_file_path):
            print(f"Error: File not found at {html_file_path}")
            return False

        try:
            print(f"Opening {html_file_path} in default web browser...")
            webbrowser.open(f"file://{os.path.abspath(html_file_path)}")
            return True
        except Exception as e:
            print(f"Error opening diagram in browser: {str(e)}")
            return False
