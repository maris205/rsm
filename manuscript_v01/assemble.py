#!/usr/bin/env python3
"""Join the editable manuscript sections, preserving local Markdown links."""
from pathlib import Path
import os
import re
from urllib.parse import urlsplit

ROOT = Path(__file__).resolve().parent
PARTS = [ROOT / "abstract.md"] + sorted((ROOT / "chapters").glob("0[1-7]_*.md"))
PARTS += sorted((ROOT / "chapters").glob("appendix_*.md")) + [ROOT / "references.md"]
LINK = re.compile(r"(!?\[[^\]\n]*\]\()([^\s)]+)(\))")


def rebase_links(source: str, folder: Path) -> str:
    def replace(match: re.Match) -> str:
        target = match.group(2)
        if target.startswith("#") or urlsplit(target).scheme:
            return match.group(0)
        name, marker, fragment = target.partition("#")
        path = (folder / name).resolve()
        if path == ROOT / "references.md" and marker:
            revised = "#" + fragment
        else:
            revised = Path(os.path.relpath(path, ROOT)).as_posix()
            if marker:
                revised += "#" + fragment
        return match.group(1) + revised + match.group(3)
    return LINK.sub(replace, source)


def main() -> None:
    chapters = [p for p in PARTS if p.parent.name == "chapters" and p.name[0].isdigit()]
    if len(chapters) != 7 or [p.name[:2] for p in chapters] != [f"{n:02d}" for n in range(1, 8)]:
        raise ValueError("Expected one complete source file for each of chapters 1–7")
    sections = [rebase_links(p.read_text(encoding="utf-8"), p.parent).strip() for p in PARTS]
    output = ROOT / "paper.md"
    output.write_text("\n\n".join(sections) + "\n", encoding="utf-8")
    print(f"Assembled {len(PARTS)} source files into {output.name}")


if __name__ == "__main__":
    main()
