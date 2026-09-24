"""Shared Helvetica typography for all manuscript figures."""

from functools import lru_cache
from pathlib import Path
from matplotlib import font_manager


@lru_cache(maxsize=1)
def figure_font():
    """Require genuine Helvetica; never substitute Nimbus Sans or DejaVu."""
    font_dir = Path(__file__).resolve().parents[2] / "fonts" / "Helvetica"
    if font_dir.is_dir():
        for path in sorted(font_dir.iterdir()):
            if path.suffix.lower() in {".ttf", ".otf", ".ttc"}:
                font_manager.fontManager.addfont(str(path))
    for family in ("Helvetica", "Helvetica LT Std"):
        try:
            path = font_manager.findfont(
                font_manager.FontProperties(family=[family]),
                fallback_to_default=False,
            )
        except ValueError:
            continue
        # Inspect the font itself rather than trusting an OS font alias.
        if font_manager.FontProperties(fname=path).get_name() == family:
            return family
    raise RuntimeError(
        "Helvetica is required for figure text. Install Helvetica or place "
        "your Helvetica TTF/OTF files in fonts/Helvetica/ (see fonts/README.md), "
        "then restart Python. Font substitution is disabled."
    )


def font_rcparams():
    """Style ordinary text while preserving each figure's math renderer."""
    family = figure_font()
    return {
        "font.family": "sans-serif",
        "font.sans-serif": [family] if family == "Helvetica" else [family, "Helvetica"],
    }
