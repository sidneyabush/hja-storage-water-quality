# Storage-Water Quality Relationship

Start here:

- [Roadmap](ROADMAP.md)
- [Progress updates](PROGRESS_UPDATES.md)

## Main Idea

Use the finalized storage-paper results to ask how watershed storage is related
to stream chemistry at HJA.

The storage paper showed that HJA watersheds do not fall along one simple
high-storage to low-storage ranking. Different storage measures describe
different parts of watershed behavior. This paper asks whether those differences
also show up in stream chemistry.

The paper outline and story notes live in the Google Doc, not in this GitHub
repo.

## Questions

- Do HJA catchments have different chemistry fingerprints?
- Do storage metrics align with chemistry behavior?
- Are catchments with more similar storage metrics more chemically synchronized?

## What To Run

From the repo root:

```r
Rscript run_active_workflow.R
```

This imports the finalized storage-paper results, runs the active chemistry and
storage checks, and writes outputs here:

`outputs/05_synthesis/storage_chemistry_links`

Candidate figures are written here:

`exploratory_plots/06_figures/main_paper_candidates`

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

Those pieces belong in the separate HJA disturbance-water quality repo.

## Archive Notes

Older exploratory and deprecated storage-proxy work is kept in
[`_archive/`](_archive/) for reference. The active paper workflow is still the
main set of scripts listed above.

## Next Steps

Current next steps are tracked in [Roadmap](ROADMAP.md).
