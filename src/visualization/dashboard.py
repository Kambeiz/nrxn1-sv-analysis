"""Interactive dashboard for NRXN1 analysis."""

from typing import List, Optional, Dict, Any
from pathlib import Path
import json


def create_dashboard(
    port: int = 8050,
    debug: bool = False
) -> Any:
    """Create and return a Dash dashboard application.
    
    Args:
        port: Port to run the dashboard on.
        debug: Enable debug mode.
    
    Returns:
        Dash application instance.
    """
    try:
        import dash
        from dash import dcc, html, dash_table
        from dash.dependencies import Input, Output, State
        import plotly.graph_objects as go
    except ImportError:
        raise ImportError(
            "Dash required for dashboard. Install with: pip install dash"
        )
    
    from .plots import CNVPlotter, GeneViewer
    from ..genomics.regions import NRXN1Region
    
    app = dash.Dash(
        __name__,
        title="NRXN1 Structural Variant Analysis",
        suppress_callback_exceptions=True
    )
    
    app.layout = html.Div([
        html.Div([
            html.H1("NRXN1 Structural Variant Analysis Dashboard"),
            html.P("Interactive visualization of CNVs and variants in NRXN1")
        ], className="header"),
        
        html.Div([
            dcc.Tabs([
                dcc.Tab(label="Gene Overview", children=[
                    html.Div([
                        html.H3("NRXN1 Gene Structure"),
                        html.P(
                            "NRXN1 is located on chromosome 2 and spans approximately "
                            "1.1 Mb. It encodes neurexin 1, a synaptic adhesion molecule "
                            "critical for synapse formation and function."
                        ),
                        html.Div([
                            html.Div([
                                html.H4("Gene Statistics"),
                                html.Ul([
                                    html.Li(f"Chromosome: 2"),
                                    html.Li(f"Start: {NRXN1Region.START:,}"),
                                    html.Li(f"End: {NRXN1Region.END:,}"),
                                    html.Li(f"Size: {(NRXN1Region.END - NRXN1Region.START):,} bp"),
                                    html.Li(f"Alpha exons: {len(NRXN1Region.EXONS_ALPHA)}"),
                                    html.Li(f"Beta exons: {len(NRXN1Region.EXONS_BETA)}"),
                                ])
                            ], style={"width": "30%", "display": "inline-block"}),
                            html.Div([
                                dcc.Graph(id="gene-structure-plot")
                            ], style={"width": "70%", "display": "inline-block"})
                        ])
                    ])
                ]),
                
                dcc.Tab(label="CNV Analysis", children=[
                    html.Div([
                        html.H3("CNV Detection Results"),
                        html.Div([
                            html.Label("Upload VCF file:"),
                            dcc.Upload(
                                id="upload-vcf",
                                children=html.Div([
                                    "Drag and Drop or ",
                                    html.A("Select VCF File")
                                ]),
                                style={
                                    "width": "100%",
                                    "height": "60px",
                                    "lineHeight": "60px",
                                    "borderWidth": "1px",
                                    "borderStyle": "dashed",
                                    "borderRadius": "5px",
                                    "textAlign": "center"
                                }
                            )
                        ]),
                        html.Div([
                            html.Div([
                                html.H4("Filter Options"),
                                html.Label("Minimum Size (bp):"),
                                dcc.Input(
                                    id="min-size",
                                    type="number",
                                    value=1000
                                ),
                                html.Label("Maximum Size (bp):"),
                                dcc.Input(
                                    id="max-size",
                                    type="number",
                                    value=10000000
                                ),
                                html.Label("CNV Type:"),
                                dcc.Checklist(
                                    id="cnv-types",
                                    options=[
                                        {"label": "Deletions", "value": "DEL"},
                                        {"label": "Duplications", "value": "DUP"}
                                    ],
                                    value=["DEL", "DUP"]
                                ),
                                html.Button("Apply Filters", id="apply-filters")
                            ], style={"width": "25%", "display": "inline-block", "verticalAlign": "top"}),
                            html.Div([
                                dcc.Graph(id="cnv-landscape-plot"),
                                dcc.Graph(id="cnv-size-distribution")
                            ], style={"width": "75%", "display": "inline-block"})
                        ])
                    ])
                ]),
                
                dcc.Tab(label="Variant Annotation", children=[
                    html.Div([
                        html.H3("Variant Annotation Results"),
                        html.Div([
                            html.Label("Search variant:"),
                            dcc.Input(
                                id="variant-search",
                                type="text",
                                placeholder="chr2:50200000:A:G"
                            ),
                            html.Button("Search", id="search-btn")
                        ]),
                        html.Div(id="variant-details"),
                        dash_table.DataTable(
                            id="variant-table",
                            columns=[
                                {"name": "Position", "id": "position"},
                                {"name": "Ref", "id": "ref"},
                                {"name": "Alt", "id": "alt"},
                                {"name": "Consequence", "id": "consequence"},
                                {"name": "gnomAD AF", "id": "gnomad_af"},
                                {"name": "ClinVar", "id": "clinvar_significance"},
                                {"name": "Pathogenicity", "id": "pathogenicity_score"}
                            ],
                            page_size=20,
                            sort_action="native",
                            filter_action="native",
                            style_table={"overflowX": "auto"}
                        )
                    ])
                ]),
                
                dcc.Tab(label="ML Predictions", children=[
                    html.Div([
                        html.H3("Pathogenicity Predictions"),
                        html.Div([
                            dcc.Graph(id="pathogenicity-distribution"),
                            dcc.Graph(id="feature-importance")
                        ]),
                        html.Div([
                            html.H4("Prediction Details"),
                            html.Div(id="prediction-explanation")
                        ])
                    ])
                ]),
                
                dcc.Tab(label="Reports", children=[
                    html.Div([
                        html.H3("Generate Reports"),
                        html.Div([
                            html.Label("Report Type:"),
                            dcc.Dropdown(
                                id="report-type",
                                options=[
                                    {"label": "Summary Report", "value": "summary"},
                                    {"label": "Clinical Report", "value": "clinical"},
                                    {"label": "Full Analysis", "value": "full"}
                                ],
                                value="summary"
                            ),
                            html.Button("Generate Report", id="generate-report"),
                            html.Div(id="report-output")
                        ])
                    ])
                ])
            ])
        ], className="main-content"),
        
        dcc.Store(id="cnv-data-store"),
        dcc.Store(id="variant-data-store")
    ])
    
    @app.callback(
        Output("gene-structure-plot", "figure"),
        Input("gene-structure-plot", "id")
    )
    def update_gene_structure(_):
        viewer = GeneViewer()
        fig_dict = viewer.create_gene_track()
        return go.Figure(fig_dict)
    
    return app


def run_dashboard(
    port: int = 8050,
    host: str = "127.0.0.1",
    debug: bool = False
) -> None:
    """Run the dashboard server.
    
    Args:
        port: Port number.
        host: Host address.
        debug: Enable debug mode.
    """
    app = create_dashboard(port=port, debug=debug)
    app.run_server(host=host, port=port, debug=debug)


def create_streamlit_app() -> None:
    """Create a Streamlit-based dashboard (alternative to Dash)."""
    try:
        import streamlit as st
        import plotly.graph_objects as go
    except ImportError:
        raise ImportError(
            "Streamlit required. Install with: pip install streamlit"
        )
    
    from .plots import CNVPlotter, GeneViewer
    from ..genomics.regions import NRXN1Region
    
    st.set_page_config(
        page_title="NRXN1 SV Analysis",
        page_icon="🧬",
        layout="wide"
    )
    
    st.title("NRXN1 Structural Variant Analysis")
    st.markdown(
        "Interactive analysis of structural variants in the NRXN1 gene"
    )
    
    tab1, tab2, tab3, tab4 = st.tabs([
        "Gene Overview", "CNV Analysis", "Variant Annotation", "ML Predictions"
    ])
    
    with tab1:
        st.header("NRXN1 Gene Structure")
        
        col1, col2 = st.columns([1, 3])
        
        with col1:
            st.subheader("Gene Statistics")
            st.write(f"**Chromosome:** 2")
            st.write(f"**Start:** {NRXN1Region.START:,}")
            st.write(f"**End:** {NRXN1Region.END:,}")
            st.write(f"**Size:** {(NRXN1Region.END - NRXN1Region.START):,} bp")
            st.write(f"**Alpha exons:** {len(NRXN1Region.EXONS_ALPHA)}")
            st.write(f"**Beta exons:** {len(NRXN1Region.EXONS_BETA)}")
        
        with col2:
            viewer = GeneViewer()
            fig_dict = viewer.create_gene_track()
            st.plotly_chart(go.Figure(fig_dict), use_container_width=True)
    
    with tab2:
        st.header("CNV Detection Results")
        
        uploaded_file = st.file_uploader("Upload VCF file", type=["vcf", "vcf.gz"])
        
        col1, col2 = st.columns([1, 3])
        
        with col1:
            st.subheader("Filter Options")
            min_size = st.number_input("Minimum Size (bp)", value=1000)
            max_size = st.number_input("Maximum Size (bp)", value=10000000)
            cnv_types = st.multiselect(
                "CNV Types",
                ["DEL", "DUP"],
                default=["DEL", "DUP"]
            )
        
        with col2:
            st.info("Upload a VCF file to visualize CNV calls")
    
    with tab3:
        st.header("Variant Annotation")
        
        variant_input = st.text_input(
            "Enter variant (chr:pos:ref:alt)",
            placeholder="2:50200000:A:G"
        )
        
        if st.button("Annotate Variant"):
            st.info("Variant annotation will be displayed here")
    
    with tab4:
        st.header("ML Pathogenicity Predictions")
        st.info("Upload variants to generate pathogenicity predictions")


if __name__ == "__main__":
    run_dashboard(debug=True)
