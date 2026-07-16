# Storage-Water Quality Roadmap

Progress updates: [PROGRESS_UPDATES.md](PROGRESS_UPDATES.md)

## Goal

Use the finalized HJA storage-paper framework to explain long-term stream
chemistry patterns across watersheds.

## Current Work

Move from exploratory storage-chemistry checks toward a first paper structure:
questions, prelim figures, and a results outline. The paper outline itself
lives in the Google Doc.

Keep this paper focused on long-term chemistry fingerprints, storage metrics,
annual chemistry structure, and pairwise chemical synchrony. Fire, flood,
landslide, and post-2020 disturbance checks belong in the separate disturbance
project.

## What Exists Now

- The main analysis script is `run_active_workflow.R`.
- Storage-chemistry synthesis tables are written to
  `/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality/outputs/05_synthesis/storage_chemistry_links`.
- Prelim main-paper figures are written to
  `/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality/exploratory_plots/06_figures/main_paper_candidates`.
- The prelim figure set was regenerated on 2026-07-16 using
  `HJA_RUN_STEPS=links,figures Rscript run_active_workflow.R`.
- Current strongest synthesis signal: storage metrics relate clearly to
  chemistry profiles and annual chemistry structure, while pairwise synchrony
  signals are present but more variable.
- Fire, flood, landslide, and post-2020 disturbance checks belong in the
  disturbance project unless they are clearly framed as future work.

## Prelim Figures

- `Fig1_storage_chemistry_response_matrix.png`: strongest overview figure for
  showing which storage metrics align with each chemistry response family.
- `Fig2_top_storage_links_by_response.png`: useful for narrowing the story to
  the top storage-chemistry links.
- `Fig3_pairwise_storage_similarity_scatter.png`: supports the pairwise
  synchrony/agreement argument; needs final axis-label cleanup before manuscript
  use.
- `Fig4_annual_stream_chemistry_storage_pca.png`: prelim annual chemistry
  ordination figure.
- `FigS1`-`FigS5`: supplemental ordination and pairwise heatmap options.

## Current Next Steps

1. Review `Fig1` and `Fig2` as the likely main synthesis figures.
2. Decide whether `Fig3` belongs in the main paper or supplement.
3. Use the storage-chemistry response matrix to write the first results outline.
4. Clean figure labels and captions after the figure shortlist is chosen.
5. Keep dated notes in [Progress updates](PROGRESS_UPDATES.md).

## Useful Files

- `run_active_workflow.R`: main analysis script.
- `05_synthesis/5b_synthesize_storage_chemistry_links.R`: synthesis tables and
  response matrix.
- `06_figures/6a_candidate_main_figures.R`: prelim main-paper figures.
- `PROGRESS_UPDATES.md`: dated decisions and next actions.
