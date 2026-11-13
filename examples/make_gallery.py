"""Regenerate the README gallery (docs/images/*.png) with the real ChemDraw application.

Every image is the PNG that ChemDraw itself exports for a document the MCP tools drew, trimmed and scaled for the web.
The CDXML sources are written next to them (docs/images/cdxml/) so each figure can be opened and edited in ChemDraw.

Run (macOS or Windows with ChemDraw; ChemDraw stays in the background):

    .venv/bin/python examples/make_gallery.py            # everything
    .venv/bin/python examples/make_gallery.py sn2 suzuki  # only these entries
"""
from __future__ import annotations

import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
OUT = ROOT / "docs" / "images"
WORK = tempfile.mkdtemp(prefix="chemdrawmcp-gallery-")
os.environ["CHEMDRAWMCP_OUTPUT_DIR"] = WORK
os.environ.setdefault("CHEMDRAWMCP_ONLINE", "0")

from PIL import Image, ImageChops  # noqa: E402

from chemdrawmcp import server  # noqa: E402

MAX_WIDTH = 1500
PAD = 36


def entries():
    """name -> (tool, kwargs).  The prompt that would produce each one is in README.md."""
    return {
        "sn2": ("draw_mechanism", dict(
            title="SN2: hydroxide + chloromethane",
            stages=[
                {"structures": ["[OH-]", "CCl"], "arrows": ["lp:0:0->0:0-1:0", "1:0-1:1->1:1"], "label": "backside attack, inversion"},
                {"structures": ["CO", "[Cl-]"], "label": "methanol + chloride"},
            ])),
        "cyanohydrin": ("draw_mechanism", dict(
            title="Cyanohydrin formation",
            stages=[
                {"structures": ["[C-]#N", "CC(C)=O"], "arrows": ["lp:0:0->0:0-1:1", "1:1=1:3->1:3"], "label": "nucleophilic addition"},
                {"structures": ["N#CC(C)(C)[O-]", "[H][OH2+]"], "arrows": ["lp:0:5->0:5-1:0", "1:0-1:1->1:1"], "label": "proton transfer"},
                {"structures": ["N#CC(C)(C)O", "O"], "label": "acetone cyanohydrin"},
            ])),
        "esterification": ("draw_reaction", dict(
            reaction="acetic acid + ethanol --[H2SO4 (cat.)]--> ethyl acetate + water",
            conditions_below=["reflux, 4 h"], product_labels=["85%", None], arrow="equilibrium")),
        "ibuprofen": ("draw_scheme", dict(
            title="BHC synthesis of ibuprofen", max_width=1000,
            steps=[
                {"reactants": ["CC(C)Cc1ccccc1"], "products": ["CC(=O)c1ccc(CC(C)C)cc1"], "above": ["Ac2O"], "below": ["HF"]},
                {"reactants": ["CC(=O)c1ccc(CC(C)C)cc1"], "products": ["CC(O)c1ccc(CC(C)C)cc1"], "above": ["H2"], "below": ["Raney Ni"]},
                {"reactants": ["CC(O)c1ccc(CC(C)C)cc1"], "products": ["CC(C)Cc1ccc(C(C)C(=O)O)cc1"], "above": ["CO"], "below": ["Pd cat."]},
            ])),
        "retro": ("draw_reaction", dict(
            reactants=["CC(=O)Nc1ccc(O)cc1"], products=["Nc1ccc(O)cc1", "CC(=O)OC(C)=O"], arrow="retrosynthetic",
            reagents_above=[], conditions_below=[])),
        "suzuki": ("apply_named_reaction", dict(reaction="suzuki", reactants=["4-bromoanisole", "phenylboronic acid"])),
        "drugs": ("draw_structures", dict(
            structures=["aspirin", "ibuprofen", "testosterone", "caffeine", "naproxen", "sildenafil"],
            labels=["aspirin", "ibuprofen", "testosterone", "caffeine", "naproxen", "sildenafil"], columns=3)),
        "nmr": ("predict_nmr", dict(structure="ethyl acetate", nucleus="1H", keep_open=False)),
        "styles_nature": ("draw_reaction", dict(
            reaction="benzaldehyde + nitromethane --[NaOH]--> (E)-beta-nitrostyrene + water",
            conditions_below=["MeOH, 0 °C"], style="nature")),
    }


CROP_BOTTOM = {"nmr": 0.53}      # keep only the structure and the spectrum of the NMR window


def trim(path: Path, keep: float = 1.0) -> Image.Image:
    im = Image.open(path)
    if im.mode in ("RGBA", "LA", "P"):
        im = im.convert("RGBA")
        bg = Image.new("RGBA", im.size, (255, 255, 255, 255))
        im = Image.alpha_composite(bg, im)
    im = im.convert("RGB")
    if keep < 1.0:
        im = im.crop((0, 0, im.width, round(im.height * keep)))
    bbox = ImageChops.difference(im, Image.new("RGB", im.size, (255, 255, 255))).getbbox()
    if bbox:
        im = im.crop(bbox)
    canvas = Image.new("RGB", (im.width + 2 * PAD, im.height + 2 * PAD), (255, 255, 255))
    canvas.paste(im, (PAD, PAD))
    if canvas.width > MAX_WIDTH:
        canvas = canvas.resize((MAX_WIDTH, round(canvas.height * MAX_WIDTH / canvas.width)), Image.LANCZOS)
    return canvas


def main(argv) -> int:
    wanted = argv or list(entries())
    (OUT / "cdxml").mkdir(parents=True, exist_ok=True)
    failures = 0
    for name in wanted:
        tool, kwargs = entries()[name]
        raw = getattr(server, tool)(return_image=False, **kwargs)
        data = json.loads(raw) if isinstance(raw, str) else raw
        data = data.get("drawn", data)          # apply_named_reaction nests the drawing under "drawn"
        render = data.get("render") or {}
        if data.get("error") or render.get("rendered_by") != "ChemDraw":
            print("FAILED %-16s %s" % (name, data.get("error") or render.get("rendered_by")))
            failures += 1
            continue
        trim(Path(render["path"]), CROP_BOTTOM.get(name, 1.0)).save(OUT / (name + ".png"), optimize=True)
        if data.get("cdxml_path"):
            shutil.copy(data["cdxml_path"], OUT / "cdxml" / (name + ".cdxml"))
        print("ok     %-16s %s" % (name, render.get("size_px")))
        try:
            if data.get("window_id"):
                server.close_document(str(data["window_id"]))
        except Exception:
            pass
    shutil.rmtree(WORK, ignore_errors=True)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
