#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

WHITE = (255, 255, 255)
DARKRED = (139, 0, 0)
DARKBLUE = (31, 72, 115)
DIAG = (228, 228, 225)

CELL, PAD, S = 132, 26, 2
FONT = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
FONT_B = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"


def ramp(v, dark, vmax):
    t = max(0.0, min(vmax, v)) / vmax
    return tuple(round(a + (b - a) * t) for a, b in zip(WHITE, dark))


def _ink(rgb):
    """Whichever of near-black / white has more contrast on this fill."""
    r, g, b = [c / 255 for c in rgb]
    lum = .2126 * r + .7152 * g + .0722 * b
    return (20, 20, 20) if (lum + .05) / .05 > 1.05 / (lum + .05) else WHITE


def draw(panels, labels, out: Path, fmt=lambda v: f"{v:.2f}",
         caption="1/sqrt(MSD)   higher = assays agree", values=False):
    n = len(labels)
    cell, pad = CELL * S, PAD * S
    rowlab, head, title = int(cell * 1.02), int(cell * .58), int(cell * .40)
    keyh = int(cell * .96)
    panelw = rowlab + n * cell + int(cell * .34)
    W = pad * 2 + panelw * len(panels)
    H = pad * 2 + title + head + n * cell + keyh

    img = Image.new("RGB", (W, H), WHITE)
    d = ImageDraw.Draw(img)
    f_val = ImageFont.truetype(FONT, int(cell * .205))
    f_lab = ImageFont.truetype(FONT, int(cell * .155))
    f_ttl = ImageFont.truetype(FONT_B, int(cell * .235))
    grey, rule = (108, 112, 118), (206, 206, 202)
    lh = lambda f: (d.textbbox((0, 0), "Hg", font=f)[3]
                    - d.textbbox((0, 0), "Hg", font=f)[1])

    def block(txt, box, font, fill, anchor="c"):
        x0, y0, x1, y1 = box
        lines = txt.split("\n")
        step = lh(font) * 1.42
        oy = y0 + (y1 - y0 - step * len(lines)) / 2
        for k, line in enumerate(lines):
            l, t, r, b = d.textbbox((0, 0), line, font=font)
            x = (x0 + (x1 - x0 - (r - l)) / 2 - l if anchor == "c"
                 else x1 - (r - l) - l)
            d.text((x, oy + k * step - t), line, font=font, fill=fill)

    for pi, (name, dark, vmax, ktick, vals) in enumerate(panels):
        px = pad + pi * panelw
        d.text((px, pad), name, font=f_ttl, fill=(24, 24, 28))
        hy = pad + title
        for j in range(n):
            block(labels[j].replace(" ", "\n", 1),
                  (px + rowlab + j * cell, hy,
                   px + rowlab + (j + 1) * cell, hy + head), f_lab, grey)
        gy = hy + head
        for i in range(n):
            block(labels[i].replace(" ", "\n", 1),
                  (px, gy + i * cell, px + rowlab - int(cell * .12),
                   gy + (i + 1) * cell), f_lab, grey, anchor="r")
            for j in range(n):
                x0, y0 = px + rowlab + j * cell, gy + i * cell
                box = (x0, y0, x0 + cell, y0 + cell)
                if i == j:
                    d.rectangle(box, fill=DIAG, outline=rule, width=max(1, S))
                    continue
                if j < i:
                    continue
                col = ramp(vals[(i, j)], dark, vmax)
                d.rectangle(box, fill=col, outline=rule, width=max(1, S))
                if values:
                    txt = fmt(vals[(i, j)])
                    l, t, r, b = d.textbbox((0, 0), txt, font=f_val)
                    d.text((x0 + (cell - (r - l)) / 2 - l,
                            y0 + (cell - (b - t)) / 2 - t),
                           txt, font=f_val, fill=_ink(col))

        ky = gy + n * cell + int(cell * .22)
        kw, kh = int(cell * 1.9), int(cell * .17)
        kx = px + rowlab
        d.text((kx, ky - lh(f_lab) - int(cell * .07)), caption, font=f_lab, fill=grey)
        for p in range(kw):
            d.line([(kx + p, ky), (kx + p, ky + kh)],
                   fill=ramp(vmax * p / (kw - 1), dark, vmax))
        d.rectangle((kx, ky, kx + kw, ky + kh), outline=rule, width=max(1, S))
        d.text((kx, ky + kh + int(cell * .06)), "0", font=f_lab, fill=grey)
        l, t, r, b = d.textbbox((0, 0), ktick, font=f_lab)
        d.text((kx + kw - (r - l), ky + kh + int(cell * .06)), ktick,
               font=f_lab, fill=grey)

    out.parent.mkdir(parents=True, exist_ok=True)
    img.resize((W // S, H // S), Image.LANCZOS).save(out, dpi=(300, 300))
    print(f"wrote {out}  ({W//S}x{H//S} px)")
