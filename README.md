# Storage-Water Quality Relationship

Project planning notes are maintained separately.

## Main Idea

Use the finalized storage-paper results to ask how watershed storage is related
to stream chemistry at HJA.

The storage paper showed that HJA watersheds do not fall along one simple
high-storage to low-storage ranking. Different storage measures describe
different parts of watershed behavior. This paper asks whether those differences
also show up in stream chemistry.

The paper outline and story notes live in the Google Doc, not in this project
folder.

## Questions

- Do HJA catchments have different chemistry fingerprints?
- Do storage metrics align with chemistry behavior?
- Are catchments with more similar storage metrics more chemically synchronized?

## What To Run

From this project folder:

```r
Rscript run_active_workflow.R
```

This imports the finalized storage-paper results, runs the active chemistry and
storage checks, and writes outputs here:

`outputs/05_synthesis/storage_chemistry_links`

Prelim figures are written here:

`exploratory_plots/06_figures/prelim_main_figures`

## Files To Edit

- `05_synthesis/5b_synthesize_storage_chemistry_links.R`: storage-chemistry
  summary tables and results outline.
- `06_figures/6a_prelim_main_figures.R`: prelim figure labels, figure list, and
  plotting choices.

Generated tables and figures in Box should be recreated by running the scripts,
not edited by hand.

## What Belongs In This Paper

- Long-term storage and chemistry patterns.
- Site differences in chemistry.
- Whether sites with similar storage behave similarly in chemistry.
- Results based on the finalized storage-paper data.

## What Does Not Belong Here

- Fire, flood, and landslide response.
- Post-2020 disturbance checks.
- Older storage proxy work that has been replaced by the finalized storage-paper
  results.

## Archive Notes

Older exploratory and deprecated storage-proxy work is kept in
[`_archive/`](_archive/) for reference. The scripts listed above are the active
paper analysis.

## Next Steps

Current next steps are maintained separately.
