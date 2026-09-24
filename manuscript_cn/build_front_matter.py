#!/usr/bin/env python3
"""Build the Chinese abstract and two conceptual figures, without experiments.

From the repository root: python manuscript_cn/build_front_matter.py
Requires Pandoc, XeTeX, the existing single-column template's packages/fonts,
and PyMuPDF. Intermediate files stay in ignored build/manuscript_cn/.
"""
from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path
from urllib.parse import quote, urlsplit

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
BUILD = REPO / "build" / "manuscript_cn"


def run(args: list[str], log: str) -> None:
    result = subprocess.run(args, cwd=BUILD, capture_output=True, text=True)
    output = result.stdout + result.stderr
    (BUILD / log).write_text(output)
    if result.returncode:
        raise RuntimeError(f"Build failed; see {BUILD / log}\n{output[-4000:]}")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    BUILD.mkdir(parents=True, exist_ok=True)
    abstract = (ROOT / "abstract.md").read_text().strip()
    title, subtitle, date, body = abstract.split("\n\n", 3)
    date = "中文摘要与总览图 · " + date.rsplit(" · ", 1)[-1]
    abstract = "\n\n".join([title, subtitle, date, body])
    title = title.removeprefix("# ")
    subtitle = subtitle.strip("*")
    captions = (ROOT / "figures" / "captions.md").read_text().strip()
    blocks = re.findall(r"^## 图 [AB]：[^\n]+\n.*?(?=^## |\Z)", captions, flags=re.M | re.S)
    if len(blocks) != 2:
        raise ValueError("captions.md must contain level-2 figure sections named 图 A and 图 B")

    def rebase(match: re.Match) -> str:
        label, target = match.groups()
        if urlsplit(target).scheme or target.startswith("#"):
            return match.group(0)
        return f"{label}(figures/{target})"

    blocks = [re.sub(r"(!?\[[^\]]*\])\(([^)]+)\)", rebase, b.strip()) for b in blocks]
    reading = abstract + "\n\n" + "\n\n".join(blocks) + "\n"
    (ROOT / "front_matter.md").write_text(reading)
    pdf_source = body + "\n\n\\clearpage\n\n" + "\n\n\\clearpage\n\n".join(blocks)
    inputs = [ROOT / "abstract.md", ROOT / "figures" / "captions.md"]

    def pdf_target(match: re.Match) -> str:
        label, target = match.groups()
        if urlsplit(target).scheme or target.startswith("#"):
            return match.group(0)
        path = (ROOT / target).resolve()
        if label.startswith("!") and path.suffix == ".png":
            path = path.with_suffix(".pdf")  # Embed the vector source in the PDF.
        if not path.is_file():
            raise FileNotFoundError(path)
        inputs.append(path)
        if label.startswith("!"):
            return f"{label}({path})"
        url = "https://github.com/maris205/rsm/blob/main/" + quote(path.relative_to(REPO).as_posix())
        return f"{label}({url})"

    pdf_source = re.sub(r"(!?\[[^\]]*\])\(([^)]+)\)", pdf_target, pdf_source)
    source = BUILD / "front_matter.md"
    source.write_text(pdf_source)
    template_source = REPO / "manuscript_v01" / "typeset" / "paper.tex"
    template = BUILD / "template.tex"
    template.write_text(template_source.read_text().replace("中文完整初稿", "中文摘要与总览图"))
    run(["pandoc", str(source), "--from=markdown-implicit_figures", "--to=latex-smart",
         "--standalone", "--template", str(template), "--metadata", f"title={title}",
         "--metadata", f"subtitle={subtitle}", "--metadata", f"date={date}",
         "--output", "front_matter.tex"], "pandoc.log")
    if executable := shutil.which("xelatex"):
        compiler = [executable]
    else:
        executable = shutil.which("xetex")
        if not executable:
            raise RuntimeError("XeTeX is required; no system packages are installed by this script")
        fmt = BUILD / "xelatex.fmt"
        if not fmt.exists():
            run([executable, "-ini", "-etex", "-interaction=nonstopmode", "-halt-on-error",
                 "-jobname=xelatex", "-progname=xelatex", "xelatex.ini"], "format_build.log")
        compiler = [executable, "-progname=xelatex", f"-fmt={fmt}"]
    for index in (1, 2):
        run(compiler + ["-interaction=nonstopmode", "-halt-on-error", "-no-shell-escape",
                        "front_matter.tex"], f"compile{index}.log")
    log = (BUILD / "front_matter.log").read_text(errors="replace")
    issues = [line for line in log.splitlines()
              if line.startswith(("Missing character:", "Overfull \\"))]
    if issues:
        raise RuntimeError("Typesetting needs correction: " + "\n".join(issues))

    import fitz
    pdf = BUILD / "front_matter.pdf"
    with fitz.open(pdf) as document:
        extracted = "\n".join(page.get_text() for page in document)
        compact = re.sub(r"\s+", "", extracted)
        for phrase in ("黎曼标准模型", "摘要", "关键词", "中心95%", "图A", "图B"):
            if phrase not in compact:
                raise RuntimeError(f"Missing expected PDF content: {phrase}")
        pages = len(document)
        (BUILD / "pdf_text.txt").write_text(extracted)
    shutil.copyfile(pdf, ROOT / "front_matter.pdf")
    report = {"pages": pages, "typesetting_issues": issues,
              "inputs_sha256": {p.relative_to(REPO).as_posix(): sha(p) for p in inputs},
              "template_sha256": sha(template_source), "renderer_sha256": sha(Path(__file__)),
              "pdf_sha256": sha(pdf), "new_scientific_computation": False}
    (BUILD / "build_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"pdf": str(ROOT / "front_matter.pdf"), "pages": pages,
                      "typesetting_issues": issues}, ensure_ascii=False))


if __name__ == "__main__":
    main()
