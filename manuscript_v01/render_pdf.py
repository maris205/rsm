#!/usr/bin/env python3
"""Render paper.md as a single-column Chinese PDF without system installs.

Run from any directory: python manuscript_v01/render_pdf.py
Requires Pandoc, XeTeX, fontspec, installed AR PL SungtiL GB / DejaVu Sans
Mono / Tinos fonts, and PyMuPDF for post-build validation. If xelatex is absent, a
local xelatex format is built and cached under build/manuscript_typeset/.
This script never assembles or edits scientific source Markdown.
"""
from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path
from urllib.parse import quote, unquote, urlsplit

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
BUILD = REPO / "build" / "manuscript_typeset"
SOURCE = ROOT / "paper.md"
OUTPUT = ROOT / "paper.pdf"
TEMPLATE = ROOT / "typeset" / "paper.tex"
PUBLIC_BASE = "https://github.com/maris205/rsm/blob/main/"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(command: list[str], *, log: str | None = None, stdin: str | None = None) -> str:
    result = subprocess.run(command, cwd=BUILD, input=stdin, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if log:
        (BUILD / log).write_text(result.stdout)
    if result.returncode:
        raise RuntimeError(f"Command failed ({result.returncode}): {command}\n{result.stdout[-6000:]}")
    return result.stdout


def inline_text(value) -> str:
    if isinstance(value, list):
        return "".join(inline_text(v) for v in value)
    if not isinstance(value, dict):
        return ""
    kind, contents = value.get("t"), value.get("c")
    if kind == "Str":
        return contents
    if kind in ("Space", "SoftBreak", "LineBreak"):
        return " "
    if kind in ("Code", "Math"):
        return contents[-1]
    return inline_text(contents)


def prepare_ast(source_text: str):
    ast = json.loads(run(["pandoc", "--from=markdown", "--to=json"], stdin=source_text))
    blocks = ast["blocks"]
    if blocks[0]["t"] != "Header" or blocks[1]["t"] != "Para" or blocks[2]["t"] != "Para":
        raise ValueError("Expected title, English subtitle and date as the first three blocks")
    title = blocks[0]["c"][2]
    subtitle = blocks[1]["c"]
    # The source English subtitle is emphasized; the template owns its style.
    if len(subtitle) == 1 and subtitle[0]["t"] == "Emph":
        subtitle = subtitle[0]["c"]
    ast["meta"].update({"title": {"t": "MetaInlines", "c": title},
                        "subtitle": {"t": "MetaInlines", "c": subtitle},
                        "date": {"t": "MetaInlines", "c": blocks[2]["c"]}})
    ast["blocks"] = blocks[3:]
    headers = {inline_text(b["c"][2]): b["c"][1][0]
               for b in ast["blocks"] if b["t"] == "Header"}
    image_sources = []
    remapped_links = []

    def rewrite(node):
        if isinstance(node, list):
            return [rewrite(item) for item in node]
        if not isinstance(node, dict):
            return node
        node = {key: rewrite(value) for key, value in node.items()}
        if node.get("t") in ("RawBlock", "RawInline") and node["c"][0] == "html":
            anchor = re.fullmatch(r'<a\s+id="([^"]+)"\s*>(?:\s*</a>)?', node["c"][1].strip())
            if anchor:
                node["c"] = ["latex", r"\hypertarget{" + anchor.group(1) + "}{}"]
            elif node["c"][1].strip() == "</a>":
                node["c"] = ["latex", ""]
            else:
                raise ValueError(f"Unsupported raw HTML must not be silently lost: {node['c'][1]}")
        if node.get("t") == "Image":
            target = node["c"][2][0]
            if urlsplit(target).scheme:
                raise ValueError("Images must be local reproducible artifacts")
            path = (ROOT / unquote(target)).resolve()
            if not path.is_file():
                raise FileNotFoundError(path)
            image_sources.append(dict(path=str(path.relative_to(REPO)), sha256=digest(path)))
            node["c"][2][0] = str(path)
        if node.get("t") == "Link":
            target = node["c"][2][0]
            split = urlsplit(target)
            if target.startswith("#") or split.scheme:
                return node
            path = (ROOT / unquote(split.path)).resolve()
            if not path.is_file():
                raise FileNotFoundError(f"Broken manuscript link: {target}")
            new_target = None
            if path.parent == ROOT / "chapters" and path.suffix == ".md":
                chapter_title = next((line[2:].strip() for line in path.read_text().splitlines()
                                      if line.startswith("# ")), None)
                if chapter_title in headers:
                    new_target = "#" + headers[chapter_title]
            if new_target is None:
                relative = path.relative_to(REPO).as_posix()
                new_target = PUBLIC_BASE + quote(relative, safe="/")
                if split.fragment:
                    new_target += "#" + quote(split.fragment)
            node["c"][2][0] = new_target
            remapped_links.append(dict(source=target, rendered=new_target))
        return node

    ast = rewrite(ast)
    first_chapter = next(index for index, block in enumerate(ast["blocks"])
                         if block["t"] == "Header" and block["c"][0] == 1)
    ast["blocks"].insert(first_chapter, {"t": "RawBlock", "c": ["latex",
                                r"\clearpage\tableofcontents\clearpage"]})
    if len(image_sources) != 3:
        raise ValueError(f"Expected the three original figures, found {len(image_sources)}")
    return ast, image_sources, remapped_links


def compiler() -> list[str]:
    if executable := shutil.which("xelatex"):
        return [executable]
    executable = shutil.which("xetex")
    if not executable:
        raise RuntimeError("XeTeX is required; no system packages will be installed")
    fmt = BUILD / "xelatex.fmt"
    if not fmt.exists():
        run([executable, "-ini", "-etex", "-interaction=nonstopmode", "-halt-on-error",
             "-jobname=xelatex", "-progname=xelatex", "xelatex.ini"], log="format_build.log")
    return [executable, "-progname=xelatex", f"-fmt={fmt}"]


def validate_pdf(path: Path, source_text: str) -> dict:
    import fitz
    with fitz.open(path) as document:
        text = "\n".join(page.get_text() for page in document)
        destinations = document.get_toc()
        links = [link for page in document for link in page.get_links()]
        bad_links = [link for link in links if link.get("kind") in (fitz.LINK_LAUNCH, fitz.LINK_GOTOR)
                     or link.get("uri", "").startswith(("file:", "/", "../", "chapters/"))]
        tags = re.findall(r"\\tag\{([^}]+)\}", source_text)
        compact_text = re.sub(r"\s+", "", text)
        missing_tags = [tag for tag in tags if f"({tag})" not in compact_text]
        checks = dict(single_column_template=True, original_figures=len({image[0] for page in document
                                                                        for image in page.get_images()}),
                      toc_entries=len(destinations), clickable_links=len(links), bad_links=bad_links,
                      equation_tags=len(tags), missing_equation_tags=missing_tags,
                      title_present="黎曼结构启发的逆对数平方动力学" in compact_text,
                      no_invented_author=not bool(document.metadata.get("author")))
        if bad_links or missing_tags or checks["original_figures"] != 3 or not checks["title_present"]:
            raise RuntimeError(f"PDF content validation failed: {checks}")
        (BUILD / "paper_text.txt").write_text(text)
        return dict(pages=len(document), checks=checks)


def main():
    BUILD.mkdir(parents=True, exist_ok=True)
    source_hash = digest(SOURCE)
    source_text = SOURCE.read_text()
    ast, images, links = prepare_ast(source_text)
    ast_path = BUILD / "paper_ast.json"
    ast_path.write_text(json.dumps(ast, ensure_ascii=False))
    tex_path = BUILD / "paper.tex"
    run(["pandoc", str(ast_path), "--from=json", "--to=latex", "--standalone",
         "--template", str(TEMPLATE), "--output", str(tex_path)], log="pandoc.log")
    tex = tex_path.read_text()
    # paper.md supplies complete, explicitly numbered captions immediately
    # after each image. Keep each image with that original caption, instead
    # of generating a second figure number from the Markdown alt text.
    figure_pattern = (r"\\begin\{figure\}\n\\centering\n(\\includegraphics[^\n]+)\n"
                      r"\\caption\{([^\n]+)\}\n\\end\{figure\}\n\n\\emph\{([^\n]+)\}")
    def figure_group(match):
        return ("\\begin{figure}[H]\n\\centering\n" + match.group(1)
                + "\n\\caption*{" + match.group(3) + "}\n\\end{figure}")
    tex, caption_count = re.subn(figure_pattern, figure_group, tex)
    if caption_count != 3:
        raise ValueError("Expected each of the three figures to have its original full caption")
    # The installed Chinese font lacks some Latin diacritics. Retain each
    # character with a local Latin font, without altering scientific text.
    tex = re.sub(r"[\u00c0-\u024f]", lambda match: r"{\latinfont " + match.group(0) + "}", tex)
    # Inline math is indivisible; legal breakpoints outside it prevent long
    # numeric constants from forcing adjacent Chinese text past the margin.
    tex = tex.replace(r"\(", r"\allowbreak{}\(").replace(r"\)", r"\)\allowbreak{}")
    # Reflow the two paired canonical equations into one equation per line.
    # This is a line-breaking change only; the single original tag is kept.
    tex = tex.replace(r"\dot q_\ell&=p_\ell/a,&", r"\dot q_\ell&=p_\ell/a,\\")
    tex = tex.replace(r"\dot R&=P_R/I,&", r"\dot R&=P_R/I,\\")
    tex_path.write_text(tex)
    command = compiler()
    for pass_number in range(1, 4):
        run(command + ["-interaction=nonstopmode", "-halt-on-error", "-no-shell-escape",
                       "paper.tex"], log=f"compile_pass{pass_number}.log")
    final_log = (BUILD / "paper.log").read_text(errors="replace")
    missing_glyphs = [line for line in final_log.splitlines() if line.startswith("Missing character:")]
    if missing_glyphs:
        raise RuntimeError("Missing glyphs in output: " + "\n".join(missing_glyphs))
    report = validate_pdf(BUILD / "paper.pdf", source_text)
    if digest(SOURCE) != source_hash:
        raise RuntimeError("paper.md changed during rendering; run again after editing finishes")
    shutil.copyfile(BUILD / "paper.pdf", OUTPUT)
    report.update(source_sha256=source_hash, renderer_sha256=digest(Path(__file__)),
                  template_sha256=digest(TEMPLATE), output_sha256=digest(OUTPUT),
                  compiler=command, figures=images, rewritten_links=links,
                  overfull_boxes=re.findall(r"Overfull \\[hv]box[^\n]*", final_log),
                  local_format_only=True)
    (BUILD / "build_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(dict(pdf=str(OUTPUT), pages=report["pages"], checks=report["checks"],
                          overfull_boxes=report["overfull_boxes"], report=str(BUILD / "build_report.json")),
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
