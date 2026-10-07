#!/usr/bin/env python3
"""Check compiled PDF structure/text and render every page for visual review."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def tool(name):
    path = shutil.which(name)
    if path:
        return path
    path = Path("/opt/homebrew/bin") / name
    if path.exists():
        return str(path)
    raise RuntimeError(f"Install Poppler: required tool {name} is unavailable")


def check_publication(run_id, no_render=False):
    with tempfile.TemporaryDirectory(prefix="dsc-publication-check-") as temporary:
        return check_in_workspace(run_id, no_render, Path(temporary))


def check_in_workspace(run_id, no_render, qa):
    output = ROOT / "report"
    checks = {}
    for label in ("report", "figure"):
        stem = f"{label}_{run_id}"
        pdf = output / (stem + ".pdf")
        info = subprocess.check_output([tool("pdfinfo"), str(pdf)], text=True)
        count = int(re.search(r"^Pages:\s+(\d+)", info, re.M).group(1))
        content = subprocess.check_output([tool("pdftotext"), "-layout", str(pdf), "-"], text=True)
        if "??" in content:
            raise AssertionError(f"Unresolved reference in {pdf.name}")
        if label == "report":
            forbidden = ["Appendix", "Selected Bayesian correlation paths",
                         "Selected smoothed conditional innovation-correlation paths",
                         "All 91 pairwise correlation paths", "Complete stationarity diagnostics",
                         "Empirical-Bayes prior calibration"]
            present = [phrase for phrase in forbidden if phrase in content]
            if present:
                raise AssertionError(f"Unexpected appendix or sample-path content in report: {present}")
        if label == "figure":
            forbidden = ["All Bayesian correlation paths", "Conditional on posterior-mean",
                         "Bands exclude parameter uncertainty", "Convergence not established"]
            present = [phrase for phrase in forbidden if phrase in content]
            if present:
                raise AssertionError(f"Unexpected provenance header text in companion figures: {present}")
        # The builder checks LaTeX logs before removing its temporary workspace.
        checks[stem] = {"pages": count, "text_characters": len(content)}
        if not no_render:
            page_dir = qa / stem
            page_dir.mkdir(exist_ok=True)
            subprocess.run([tool("pdftoppm"), "-scale-to", "1400", "-png", str(pdf), str(page_dir / "page")], check=True)
            # A contact sheet is only a navigation aid; inspect full page PNGs as needed.
            from PIL import Image, ImageDraw
            paths = sorted((p for p in page_dir.glob("page-*.png") if int(p.stem.split("-")[-1]) <= count),
                           key=lambda p: int(p.stem.split("-")[-1]))
            if len(paths) != count or [int(p.stem.split("-")[-1]) for p in paths] != list(range(1, count + 1)):
                raise AssertionError("Rendered page inventory disagrees with the current PDF")
            checks[stem]["rendered_pages"] = [p.name for p in paths]
            checks[stem]["contact_sheet_count"] = (count + 11) // 12
            for block in range(0, len(paths), 12):
                group = paths[block:block+12]
                sheet = Image.new("RGB", (1200, ((len(group)+2)//3)*340), "#DDE2E7")
                draw = ImageDraw.Draw(sheet)
                for k, path in enumerate(group):
                    im = Image.open(path).convert("RGB")
                    im.thumbnail((382, 312))
                    x, y = (k % 3)*400 + (400-im.width)//2, (k//3)*340 + 18
                    sheet.paste(im, (x, y))
                    draw.text(((k%3)*400+10, (k//3)*340+3), path.stem, fill="black")
                sheet.save(qa / f"{stem}_contact_{block//12+1:02d}.png")
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-render", action="store_true")
    parser.add_argument("--run-id", required=True, help="Run-folder ID used in the two PDF filenames")
    args = parser.parse_args()
    print(json.dumps(check_publication(args.run_id, args.no_render), indent=2))


if __name__ == "__main__":
    main()
