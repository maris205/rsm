#!/usr/bin/env python3
"""Assemble completed Chinese chapters and render a cumulative single-column PDF.

From the repository root: python manuscript_cn/build_manuscript.py
Requires Pandoc, XeTeX, the fonts/packages used by manuscript_v01/typeset/paper.tex,
and PyMuPDF. No system installation, network access, or scientific computation is
performed. The original three-page front matter and scientific sources are left
untouched; temporary files stay in ignored build/manuscript_cn_full/.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
from pathlib import Path
from urllib.parse import quote, unquote, urlsplit

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
BUILD = REPO / "build" / "manuscript_cn_full"
TEMPLATE_SOURCE = REPO / "manuscript_v01" / "typeset" / "paper.tex"
PUBLIC_BASE = "https://github.com/maris205/rsm/blob/main/"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(command: list[str], *, log: str, stdin: str | None = None) -> str:
    result = subprocess.run(command, cwd=BUILD, input=stdin, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (BUILD / log).write_text(result.stdout)
    if result.returncode:
        raise RuntimeError(f"Build failed; see {BUILD / log}\n{result.stdout[-5000:]}")
    return result.stdout


def rebase(text: str, source: Path) -> str:
    """Resolve relative Markdown links at their source location before merging."""
    def replace(match: re.Match) -> str:
        label, target = match.groups()
        split = urlsplit(target)
        if split.scheme or target.startswith("#"):
            return match.group(0)
        path = (source.parent / unquote(split.path)).resolve()
        if path == ROOT / "references.md" and split.fragment:
            return f"{label}(#{split.fragment})"
        relative = Path(os.path.relpath(path, ROOT)).as_posix()
        suffix = ("?" + split.query if split.query else "") + ("#" + split.fragment if split.fragment else "")
        return f"{label}({quote(relative, safe='/._-')}{suffix})"
    return re.sub(r"(!?\[[^\]]*\])\(([^)]+)\)", replace, text)


def assemble() -> tuple[str, list[Path], list[Path]]:
    abstract_path = ROOT / "abstract.md"
    captions_path = ROOT / "figures" / "captions.md"
    references_path = ROOT / "references.md"
    chapters = sorted((ROOT / "chapters").glob("[0-9][0-9]_*.md"))
    if not chapters or not references_path.is_file():
        raise FileNotFoundError("At least one completed chapter and references.md are required")
    abstract = abstract_path.read_text().strip()
    title, subtitle, date, body = abstract.split("\n\n", 3)
    chapter_numbers = [str(int(path.name[:2])) for path in chapters]
    coverage = "第 " + "、".join(chapter_numbers) + " 章"
    date = re.sub(r" · 摘要(?= ·|$)", f" · 累计稿（已完成{coverage}）", date)
    captions = captions_path.read_text()
    figures = re.findall(r"^## 图 [AB]：[^\n]+\n.*?(?=^## |\Z)", captions, flags=re.M | re.S)
    if len(figures) != 2:
        raise ValueError("Expected exactly the two front-matter figures A and B")
    pieces = ["\n\n".join([title, subtitle, date, body])]
    pieces.extend(rebase(block.strip(), captions_path) for block in figures)
    pieces.extend(rebase(path.read_text().strip(), path) for path in chapters)
    pieces.append(rebase(references_path.read_text().strip(), references_path))
    text = "\n\n".join(pieces) + "\n"
    (ROOT / "paper.md").write_text(text)
    return text, [abstract_path, captions_path, *chapters, references_path], chapters


def inline_text(value) -> str:
    if isinstance(value, list):
        return "".join(inline_text(item) for item in value)
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


def prepare_ast(text: str, chapters: list[Path]):
    ast = json.loads(run(["pandoc", "--from=markdown-implicit_figures", "--to=json"],
                         stdin=text, log="parse.log"))
    blocks = ast["blocks"]
    if [block["t"] for block in blocks[:3]] != ["Header", "Para", "Para"]:
        raise ValueError("Expected the original title, subtitle, and date metadata")
    subtitle = blocks[1]["c"]
    if len(subtitle) == 1 and subtitle[0]["t"] == "Emph":
        subtitle = subtitle[0]["c"]
    ast["meta"].update(title={"t": "MetaInlines", "c": blocks[0]["c"][2]},
                       subtitle={"t": "MetaInlines", "c": subtitle},
                       date={"t": "MetaInlines", "c": blocks[2]["c"]})
    ast["blocks"] = blocks[3:]
    headers = {inline_text(b["c"][2]): b["c"][1][0]
               for b in ast["blocks"] if b["t"] == "Header"}
    images, links, anchors = [], [], []

    def rewrite(node):
        if isinstance(node, list):
            return [rewrite(item) for item in node]
        if not isinstance(node, dict):
            return node
        node = {key: rewrite(value) for key, value in node.items()}
        if node.get("t") == "Table":
            # Chinese pipe tables have no word spaces for Pandoc to infer
            # wrapping widths. Supply paragraph-column widths, not text edits.
            specifications = node["c"][2]
            count = len(specifications)
            widths = {2: [0.27, 0.73], 3: [0.22, 0.39, 0.39]}.get(count, [1 / count] * count)
            if all(isinstance(value, (int, float)) for value in specifications):
                node["c"][2] = widths  # Pandoc API before 1.21.
            else:
                for specification, width in zip(specifications, widths):
                    specification[1] = {"t": "ColWidth", "c": width}
        if node.get("t") in ("RawBlock", "RawInline") and node["c"][0] == "html":
            anchor = re.fullmatch(r'<a\s+id="([^"]+)"\s*>(?:\s*</a>)?', node["c"][1].strip())
            if anchor:
                anchors.append(anchor.group(1))
                node["c"] = ["latex", r"\hypertarget{" + anchor.group(1) + "}{}"]
            elif node["c"][1].strip() == "</a>":
                node["c"] = ["latex", ""]
            else:
                raise ValueError(f"Unsupported HTML must not silently disappear: {node['c'][1]}")
        if node.get("t") == "Image":
            target = node["c"][2][0]
            if urlsplit(target).scheme:
                raise ValueError("Figures must be local reproducible artifacts")
            path = (ROOT / unquote(target)).resolve()
            if path.suffix == ".png" and path.with_suffix(".pdf").is_file():
                path = path.with_suffix(".pdf")
            if not path.is_file():
                raise FileNotFoundError(path)
            images.append({"path": path.relative_to(REPO).as_posix(), "sha256": digest(path)})
            node["c"][2][0] = str(path)
        if node.get("t") == "Link":
            target = node["c"][2][0]
            split = urlsplit(target)
            if split.scheme or target.startswith("#"):
                return node
            path = (ROOT / unquote(split.path)).resolve()
            if not path.is_file():
                raise FileNotFoundError(f"Broken manuscript link: {target}")
            replacement = None
            if path in chapters:
                title = next(line[2:].strip() for line in path.read_text().splitlines()
                             if line.startswith("# "))
                if title in headers:
                    replacement = "#" + headers[title]
            if replacement is None:
                replacement = PUBLIC_BASE + quote(path.relative_to(REPO).as_posix())
                if split.fragment:
                    replacement += "#" + quote(split.fragment)
            links.append({"source": target, "rendered": replacement})
            node["c"][2][0] = replacement
        return node

    ast = rewrite(ast)
    laid_out = []
    inserted_toc = False
    for block in ast["blocks"]:
        if block["t"] == "Header":
            heading = inline_text(block["c"][2])
            if heading.startswith(("图 A：", "图 B：")):
                laid_out.append({"t": "RawBlock", "c": ["latex", r"\clearpage"]})
            elif block["c"][0] == 1:
                if not inserted_toc:
                    laid_out.append({"t": "RawBlock", "c": ["latex", r"\clearpage\tableofcontents\clearpage"]})
                    inserted_toc = True
                else:
                    laid_out.append({"t": "RawBlock", "c": ["latex", r"\clearpage"]})
        laid_out.append(block)
    ast["blocks"] = laid_out
    if len(images) != 2 or not inserted_toc:
        raise ValueError("Expected two front-matter figures and at least one chapter")
    return ast, images, links, anchors


def compiler() -> list[str]:
    if executable := shutil.which("xelatex"):
        return [executable]
    executable = shutil.which("xetex")
    if not executable:
        raise RuntimeError("XeTeX is required; this script does not install system packages")
    fmt = BUILD / "xelatex.fmt"
    if not fmt.exists():
        run([executable, "-ini", "-etex", "-interaction=nonstopmode", "-halt-on-error",
             "-jobname=xelatex", "-progname=xelatex", "xelatex.ini"], log="format_build.log")
    return [executable, "-progname=xelatex", f"-fmt={fmt}"]


def validate_pdf(path: Path, text: str, chapters: list[Path], anchors: list[str]) -> dict:
    import fitz
    with fitz.open(path) as document:
        extracted = "\n".join(page.get_text() for page in document)
        compact = re.sub(r"\s+", "", extracted)
        links = [link for page in document for link in page.get_links()]
        bad_links = [link for link in links if link.get("kind") in (fitz.LINK_LAUNCH, fitz.LINK_GOTOR)
                     or link.get("uri", "").startswith(("file:", "/", "../", "chapters/"))]
        tags = re.findall(r"\\tag\{([^}]+)\}", text)
        missing_tags = [tag for tag in tags if f"({tag})" not in compact]
        expected = ["黎曼标准模型", "摘要", "关键词", "图A", "图B", "目录", "参考文献"]
        for chapter in chapters:
            expected.extend(re.sub(r"^#+\s+", "", line).replace(" ", "")
                            for line in chapter.read_text().splitlines()
                            if re.match(r"^#{1,2} ", line))
        missing_headings = [heading for heading in expected if heading not in compact]
        # XeTeX embeds each PDF figure as a reusable form XObject, preserving vectors.
        xobjects = {entry[0] for page in document for entry in page.get_xobjects()}
        tex = (BUILD / "paper.tex").read_text()
        missing_anchors = [anchor for anchor in anchors
                           if r"\hypertarget{" + anchor + "}" not in tex]
        unresolved_anchors = [anchor for anchor in anchors
                              if document.resolve_link("#" + anchor)[0] < 0]
        figure_pages = [number + 1 for number, page in enumerate(document)
                        if page.get_xobjects()]
        checks = dict(completed_chapters=[path.stem for path in chapters],
                      equation_tags=len(tags), missing_equation_tags=missing_tags,
                      missing_headings=missing_headings, reference_anchors=len(anchors),
                      missing_reference_anchors=missing_anchors,
                      unresolved_reference_anchors=unresolved_anchors,
                      toc_entries=len(document.get_toc()), clickable_links=len(links),
                      bad_links=bad_links, vector_figure_xobjects=len(xobjects),
                      front_matter_figure_pages=figure_pages,
                      no_invented_author=not bool(document.metadata.get("author")))
        if (bad_links or missing_tags or missing_headings or missing_anchors or unresolved_anchors
                or len(xobjects) < 2 or figure_pages != [2, 3]):
            raise RuntimeError(f"PDF validation failed: {checks}")
        (BUILD / "paper_text.txt").write_text(extracted)
        return {"pages": len(document), "checks": checks}


def main() -> None:
    BUILD.mkdir(parents=True, exist_ok=True)
    text, inputs, chapters = assemble()
    input_hashes = {path: digest(path) for path in inputs}
    ast, images, links, anchors = prepare_ast(text, chapters)
    ast_path = BUILD / "paper_ast.json"
    ast_path.write_text(json.dumps(ast, ensure_ascii=False))
    template = BUILD / "template.tex"
    template_text = TEMPLATE_SOURCE.read_text().replace("中文完整初稿", "中文逐章修订累计稿")
    # The legacy \hbar macro takes its overbar from the text Roman family,
    # whose Chinese font lacks U+00AF. Use the existing AMS mathematical
    # glyph instead; this preserves \hbar in all manuscript source files.
    template_text = template_text.replace(
        r"\usepackage{amsmath,amssymb,mathtools,mathrsfs}",
        r"\usepackage{amsmath,amssymb,mathtools,mathrsfs}" + "\n"
        + r'\DeclareMathSymbol{\rsmhbar}{\mathord}{AMSb}{"7E}' + "\n"
        + r"\AtBeginDocument{\renewcommand{\hbar}{\rsmhbar}}")
    template.write_text(template_text)
    run(["pandoc", str(ast_path), "--from=json", "--to=latex", "--standalone",
         "--template", str(template), "--output", "paper.tex"], log="pandoc.log")
    tex_path = BUILD / "paper.tex"
    tex = tex_path.read_text()
    # Preserve Latin diacritics absent from the Chinese body font.
    tex = re.sub(r"[\u00c0-\u024f]", lambda match: r"{\latinfont " + match.group(0) + "}", tex)
    tex = tex.replace(r"\(", r"\allowbreak{}\(").replace(r"\)", r"\)\allowbreak{}")
    tex_path.write_text(tex)
    command = compiler()
    for pass_number in range(1, 4):
        run(command + ["-interaction=nonstopmode", "-halt-on-error", "-no-shell-escape", "paper.tex"],
            log=f"compile_pass{pass_number}.log")
    log = (BUILD / "paper.log").read_text(errors="replace")
    issues = [line for line in log.splitlines()
              if line.startswith(("Missing character:", "Overfull \\"))]
    if issues:
        raise RuntimeError("Typesetting requires correction: " + "\n".join(issues))
    report = validate_pdf(BUILD / "paper.pdf", text, chapters, anchors)
    if any(digest(path) != value for path, value in input_hashes.items()):
        raise RuntimeError("A manuscript source changed during rendering; rerun on finished sources")
    shutil.copyfile(BUILD / "paper.pdf", ROOT / "paper.pdf")
    report.update(inputs_sha256={p.relative_to(REPO).as_posix(): v for p, v in input_hashes.items()},
                  source_sha256=digest(ROOT / "paper.md"), renderer_sha256=digest(Path(__file__)),
                  template_sha256=digest(TEMPLATE_SOURCE), output_sha256=digest(ROOT / "paper.pdf"),
                  figures=images, rewritten_links=links, typesetting_issues=issues,
                  new_scientific_computation=False)
    (BUILD / "build_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"pdf": str(ROOT / "paper.pdf"), "pages": report["pages"],
                      "checks": report["checks"], "typesetting_issues": issues}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
