#!/usr/bin/env python3
"""Check compiled PDF structure/text and render every page for visual review."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def tool(name):
    path = shutil.which(name)
    if path:
        return path
    path = Path("/opt/homebrew/bin") / name
    if path.exists():
        return str(path)
    raise RuntimeError(f"Install Poppler: required tool {name} is unavailable")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-render", action="store_true")
    parser.add_argument("--report-only", action="store_true", help="Check only the report PDF")
    args = parser.parse_args()
    output = ROOT / "output/pdf"
    qa = ROOT / "research/generated/qa"
    qa.mkdir(parents=True, exist_ok=True)
    checks = {}
    stems = ["weekly_correlation_report"] if args.report_only else ["weekly_correlation_report", "weekly_correlation_slides"]
    for stem in stems:
        pdf = output / (stem + ".pdf")
        info = subprocess.check_output([tool("pdfinfo"), str(pdf)], text=True)
        count = int(re.search(r"^Pages:\s+(\d+)", info, re.M).group(1))
        content = subprocess.check_output([tool("pdftotext"), "-layout", str(pdf), "-"], text=True)
        if "??" in content:
            raise AssertionError(f"Unresolved reference in {pdf.name}")
        if "slides" in stem and not 12 <= count <= 15:
            raise AssertionError(f"Expected 12--15 slides, got {count}")
        if "report" in stem and "All 91 pairwise correlation paths" not in content:
            raise AssertionError("Missing complete pair appendix")
        publication_manifest = ROOT / "research/generated/publication_manifest.json"
        if "report" in stem and publication_manifest.exists():
            manifest = json.loads(publication_manifest.read_text())
            if manifest.get("posterior_run") and "All 91 Bayesian correlation paths" not in content:
                raise AssertionError("Missing Bayesian correlation-path appendix from the selected sampler run")
        log = (output / (stem + ".log")).read_text(errors="replace")
        issues = [line for line in log.splitlines() if re.search(r"Overfull \\[hv]box|undefined|Citation .* undefined", line)]
        if issues:
            raise AssertionError(f"LaTeX layout/reference issues in {stem}: {issues}")
        checks[stem] = {"pages": count, "text_characters": len(content), "latex_issues": issues}
        if not args.no_render:
            page_dir = qa / stem
            page_dir.mkdir(exist_ok=True)
            # These are this checker's reproducible QA images, not source artifacts.
            # Remove only surplus pages/contact sheets from an older PDF build.
            for stale in page_dir.glob("page-*.png"):
                match = re.fullmatch(r"page-(\d+)\.png", stale.name)
                if match and int(match.group(1)) > count:
                    stale.unlink()
            for stale in qa.glob(stem + "_contact_*.png"):
                if re.fullmatch(re.escape(stem) + r"_contact_\d+\.png", stale.name):
                    stale.unlink()
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
    (qa / "checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    main()
