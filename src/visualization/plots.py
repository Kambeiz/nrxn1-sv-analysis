"""Plotting utilities for NRXN1 variant visualization."""

from typing import List, Optional, Dict, Any
import numpy as np
import pandas as pd

from ..genomics.regions import NRXN1Region, GenomicRegion
from ..cnv.detector import CNVCall, CNVType
from ..annotation.annotator import AnnotatedVariant


class CNVPlotter:
    """Create visualizations for CNV analysis."""
    
    def __init__(self):
        """Initialize plotter."""
        self.nrxn1 = NRXN1Region()
    
    def plot_cnv_landscape(
        self,
        cnv_calls: List[CNVCall],
        title: str = "NRXN1 CNV Landscape"
    ) -> Dict[str, Any]:
        """Create CNV landscape plot data.
        
        Args:
            cnv_calls: List of CNV calls.
            title: Plot title.
        
        Returns:
            Plotly figure dictionary.
        """
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        
        fig = make_subplots(
            rows=2, cols=1,
            row_heights=[0.7, 0.3],
            shared_xaxes=True,
            vertical_spacing=0.05,
            subplot_titles=("CNV Calls", "Gene Structure")
        )
        
        deletions = [c for c in cnv_calls if c.cnv_type == CNVType.DELETION]
        duplications = [c for c in cnv_calls if c.cnv_type == CNVType.DUPLICATION]
        
        for i, cnv in enumerate(deletions):
            fig.add_trace(
                go.Scatter(
                    x=[cnv.start, cnv.end, cnv.end, cnv.start, cnv.start],
                    y=[i, i, i + 0.8, i + 0.8, i],
                    fill="toself",
                    fillcolor="rgba(255, 0, 0, 0.3)",
                    line=dict(color="red", width=1),
                    name=f"DEL {cnv.length // 1000}kb",
                    hovertemplate=(
                        f"Type: Deletion<br>"
                        f"Start: {cnv.start:,}<br>"
                        f"End: {cnv.end:,}<br>"
                        f"Size: {cnv.length:,} bp<br>"
                        f"Quality: {cnv.quality:.1f}"
                    )
                ),
                row=1, col=1
            )
        
        offset = len(deletions)
        for i, cnv in enumerate(duplications):
            fig.add_trace(
                go.Scatter(
                    x=[cnv.start, cnv.end, cnv.end, cnv.start, cnv.start],
                    y=[offset + i, offset + i, offset + i + 0.8, offset + i + 0.8, offset + i],
                    fill="toself",
                    fillcolor="rgba(0, 0, 255, 0.3)",
                    line=dict(color="blue", width=1),
                    name=f"DUP {cnv.length // 1000}kb",
                    hovertemplate=(
                        f"Type: Duplication<br>"
                        f"Start: {cnv.start:,}<br>"
                        f"End: {cnv.end:,}<br>"
                        f"Size: {cnv.length:,} bp<br>"
                        f"Quality: {cnv.quality:.1f}"
                    )
                ),
                row=1, col=1
            )
        
        for exon in self.nrxn1.exons_alpha:
            fig.add_trace(
                go.Scatter(
                    x=[exon.start, exon.end, exon.end, exon.start, exon.start],
                    y=[0, 0, 1, 1, 0],
                    fill="toself",
                    fillcolor="rgba(0, 128, 0, 0.5)",
                    line=dict(color="green", width=1),
                    name=exon.name,
                    showlegend=False
                ),
                row=2, col=1
            )
        
        fig.add_trace(
            go.Scatter(
                x=[self.nrxn1.START, self.nrxn1.END],
                y=[0.5, 0.5],
                mode="lines",
                line=dict(color="black", width=2),
                showlegend=False
            ),
            row=2, col=1
        )
        
        fig.update_layout(
            title=title,
            xaxis2_title="Genomic Position (chr2)",
            yaxis_title="CNV Calls",
            height=600,
            showlegend=True
        )
        
        fig.update_xaxes(range=[self.nrxn1.START - 10000, self.nrxn1.END + 10000])
        
        return fig.to_dict()
    
    def plot_cnv_size_distribution(
        self,
        cnv_calls: List[CNVCall]
    ) -> Dict[str, Any]:
        """Plot CNV size distribution."""
        import plotly.graph_objects as go
        
        del_sizes = [c.length for c in cnv_calls if c.cnv_type == CNVType.DELETION]
        dup_sizes = [c.length for c in cnv_calls if c.cnv_type == CNVType.DUPLICATION]
        
        fig = go.Figure()
        
        if del_sizes:
            fig.add_trace(go.Histogram(
                x=del_sizes,
                name="Deletions",
                marker_color="red",
                opacity=0.7
            ))
        
        if dup_sizes:
            fig.add_trace(go.Histogram(
                x=dup_sizes,
                name="Duplications",
                marker_color="blue",
                opacity=0.7
            ))
        
        fig.update_layout(
            title="CNV Size Distribution",
            xaxis_title="CNV Size (bp)",
            yaxis_title="Count",
            barmode="overlay"
        )
        
        return fig.to_dict()
    
    def plot_pathogenicity_scores(
        self,
        annotated_variants: List[AnnotatedVariant]
    ) -> Dict[str, Any]:
        """Plot pathogenicity score distribution."""
        import plotly.graph_objects as go
        
        scores = [av.pathogenicity_score for av in annotated_variants]
        classifications = []
        for s in scores:
            if s >= 0.9:
                classifications.append("Pathogenic")
            elif s >= 0.7:
                classifications.append("Likely Pathogenic")
            elif s >= 0.3:
                classifications.append("VUS")
            elif s >= 0.1:
                classifications.append("Likely Benign")
            else:
                classifications.append("Benign")
        
        colors = {
            "Pathogenic": "red",
            "Likely Pathogenic": "orange",
            "VUS": "gray",
            "Likely Benign": "lightblue",
            "Benign": "green"
        }
        
        fig = go.Figure()
        
        fig.add_trace(go.Histogram(
            x=scores,
            marker_color=[colors.get(c, "gray") for c in classifications],
            name="Pathogenicity Scores"
        ))
        
        for threshold, label in [(0.9, "P"), (0.7, "LP"), (0.3, "VUS"), (0.1, "LB")]:
            fig.add_vline(x=threshold, line_dash="dash", line_color="gray")
            fig.add_annotation(x=threshold, y=1, text=label, showarrow=False)
        
        fig.update_layout(
            title="Pathogenicity Score Distribution",
            xaxis_title="Pathogenicity Score",
            yaxis_title="Count"
        )
        
        return fig.to_dict()


class GeneViewer:
    """Interactive gene structure viewer."""
    
    def __init__(self):
        """Initialize gene viewer."""
        self.nrxn1 = NRXN1Region()
    
    def create_gene_track(
        self,
        variants: Optional[List[AnnotatedVariant]] = None,
        cnvs: Optional[List[CNVCall]] = None,
        show_domains: bool = True
    ) -> Dict[str, Any]:
        """Create interactive gene track visualization.
        
        Args:
            variants: Optional list of variants to display.
            cnvs: Optional list of CNVs to display.
            show_domains: Whether to show functional domains.
        
        Returns:
            Plotly figure dictionary.
        """
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        
        n_rows = 3 if show_domains else 2
        if cnvs:
            n_rows += 1
        if variants:
            n_rows += 1
        
        fig = make_subplots(
            rows=n_rows, cols=1,
            shared_xaxes=True,
            vertical_spacing=0.02,
            row_heights=[1] * n_rows
        )
        
        current_row = 1
        
        for exon in self.nrxn1.exons_alpha:
            fig.add_trace(
                go.Scatter(
                    x=[exon.start, exon.end, exon.end, exon.start, exon.start],
                    y=[0.2, 0.2, 0.8, 0.8, 0.2],
                    fill="toself",
                    fillcolor="rgba(0, 100, 0, 0.7)",
                    line=dict(color="darkgreen", width=1),
                    name=exon.name,
                    hovertemplate=f"{exon.name}<br>Size: {exon.length} bp"
                ),
                row=current_row, col=1
            )
        
        fig.add_trace(
            go.Scatter(
                x=[self.nrxn1.START, self.nrxn1.END],
                y=[0.5, 0.5],
                mode="lines",
                line=dict(color="black", width=3),
                name="NRXN1 Alpha",
                showlegend=True
            ),
            row=current_row, col=1
        )
        current_row += 1
        
        for exon in self.nrxn1.exons_beta:
            fig.add_trace(
                go.Scatter(
                    x=[exon.start, exon.end, exon.end, exon.start, exon.start],
                    y=[0.2, 0.2, 0.8, 0.8, 0.2],
                    fill="toself",
                    fillcolor="rgba(0, 0, 100, 0.7)",
                    line=dict(color="darkblue", width=1),
                    name=exon.name,
                    showlegend=False
                ),
                row=current_row, col=1
            )
        
        fig.add_trace(
            go.Scatter(
                x=[self.nrxn1.exons_beta[0].start, self.nrxn1.END],
                y=[0.5, 0.5],
                mode="lines",
                line=dict(color="navy", width=3),
                name="NRXN1 Beta",
                showlegend=True
            ),
            row=current_row, col=1
        )
        current_row += 1
        
        if show_domains:
            domain_colors = {
                "signal_peptide_and_lns1": "rgba(255, 165, 0, 0.5)",
                "egf_like_domain": "rgba(255, 0, 255, 0.5)",
                "lns2_lns3": "rgba(0, 255, 255, 0.5)",
                "lns4_lns5_lns6": "rgba(255, 255, 0, 0.5)",
                "transmembrane_and_cytoplasmic": "rgba(128, 0, 128, 0.5)"
            }
            
            for domain in self.nrxn1.domains:
                fig.add_trace(
                    go.Scatter(
                        x=[domain.start, domain.end, domain.end, domain.start, domain.start],
                        y=[0.1, 0.1, 0.9, 0.9, 0.1],
                        fill="toself",
                        fillcolor=domain_colors.get(domain.name, "rgba(128,128,128,0.5)"),
                        line=dict(width=1),
                        name=domain.name.replace("_", " ").title(),
                        hovertemplate=f"{domain.name}<br>Size: {domain.length:,} bp"
                    ),
                    row=current_row, col=1
                )
            current_row += 1
        
        if cnvs:
            for i, cnv in enumerate(cnvs):
                color = "red" if cnv.cnv_type == CNVType.DELETION else "blue"
                fig.add_trace(
                    go.Scatter(
                        x=[cnv.start, cnv.end, cnv.end, cnv.start, cnv.start],
                        y=[i % 5, i % 5, (i % 5) + 0.8, (i % 5) + 0.8, i % 5],
                        fill="toself",
                        fillcolor=f"rgba({color}, 0.3)",
                        line=dict(color=color, width=1),
                        name=f"{cnv.cnv_type.value}",
                        showlegend=False
                    ),
                    row=current_row, col=1
                )
            current_row += 1
        
        if variants:
            positions = [v.variant.position for v in variants]
            scores = [v.pathogenicity_score for v in variants]
            
            fig.add_trace(
                go.Scatter(
                    x=positions,
                    y=scores,
                    mode="markers",
                    marker=dict(
                        size=8,
                        color=scores,
                        colorscale="RdYlGn_r",
                        showscale=True,
                        colorbar=dict(title="Pathogenicity")
                    ),
                    name="Variants",
                    hovertemplate="Position: %{x}<br>Score: %{y:.2f}"
                ),
                row=current_row, col=1
            )
        
        fig.update_layout(
            title="NRXN1 Gene Structure and Variants",
            height=150 * n_rows,
            showlegend=True,
            xaxis=dict(range=[self.nrxn1.START - 50000, self.nrxn1.END + 50000])
        )
        
        fig.update_xaxes(title_text="Genomic Position (chr2)", row=n_rows, col=1)
        
        return fig.to_dict()
    
    def export_svg(self, figure_dict: Dict[str, Any], output_path: str) -> None:
        """Export figure to SVG."""
        import plotly.graph_objects as go
        
        fig = go.Figure(figure_dict)
        fig.write_image(output_path, format="svg")
    
    def export_html(self, figure_dict: Dict[str, Any], output_path: str) -> None:
        """Export figure to interactive HTML."""
        import plotly.graph_objects as go
        
        fig = go.Figure(figure_dict)
        fig.write_html(output_path)
