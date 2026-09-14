# Publication sources

`weekly_correlation_report.tex` and `weekly_correlation_slides.tex` are the editable report and 15-slide Beamer sources. They share the bibliography, generated tables, and vector data charts. The report's appendix covers all 91 pairs in the same ordering as MATLAB's column-major lower triangle.

Run from the repository root after MATLAB has written its numerical results:

```bash
python3 scripts/build_publication.py --data outputs/weekly_research/data --pilot PATH/TO/pilot_summary.json
python3 scripts/check_publication.py
```

The builder does not launch estimation. It validates dimensions, Friday spacing, the input fingerprint when the source is available, return reconstruction, and rolling-output coverage. It writes the PDFs to `output/pdf/` and supporting material under `research/generated/`. The baseline publication expects the `asof` closure policy. Masked-model results need separate sensitivity interpretation.

Dependencies: Python 3 with NumPy, pandas, Matplotlib, and Pillow; a TeX installation with `latexmk`, `pdflatex`, Beamer, and the standard packages listed in the sources; Poppler's `pdfinfo`, `pdftotext`, and `pdftoppm` for verification. On the supplied Mac, system `python3` has the chart libraries and TeX is under `/Library/TeX/texbin`. No PowerPoint or LibreOffice conversion is required.

`--no-compile` regenerates evidence only. Omitting `--pilot` creates a deliberately labeled draft with no timing or posterior claims. A zero-completed-iteration pilot is handled explicitly, with no fabricated runtime projection. Final documents require a saved pilot summary and visual inspection of the regenerated PDF pages.

The publication manifest records the numerical CSV fingerprints and pilot path. Prior sensitivity comes from `data/prior_summary.json`. Pilot timings, iteration counts, proposal counts, and analytical memory/storage comparisons come from the saved sampler summary. Those analytical sizes are not measured peak resident memory.

When a sibling `resource_observations.json` accompanies the pilot summary, the builder also reports its sampled MATLAB resident-set memory and fingerprints that file. This optional record contains positive `rss_kib` measurements with observation timestamps and process elapsed times. The builder converts KiB to GiB and labels every sample as point-in-time memory, never as a measured peak. Such a sample includes the MATLAB runtime and excludes helper processes.
