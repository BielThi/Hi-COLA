"""Slide versions of the notebook figures, for a LaTeX beamer talk (4:3 frame, 128 x 96 mm).

    import slides
    slides.show(fig, 'compare_Abacus_Small_pk', slide=True)

slides.show always shows the figure as it is. With slide=True it also saves it, unchanged, as
SLIDE_DIR/<name>.pdf (vector) and SLIDE_DIR/<name>.png (PNG_DPI), scaled as a whole (text, lines,
panels together, aspect ratio kept) to fit in SLIDE_SIZE. The slide file is therefore the notebook
figure, only magnified or shrunk: a figure made at 6.6 in wide with 10.5 pt labels and saved at
4.7 in wide has 7.5 pt labels on the slide. An existing file of the same name is overwritten.

Include it at its natural size, without width= or height=:

    \\includegraphics{slides/<name>.pdf}
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_mixed import MixedModeRenderer
from matplotlib.backends.backend_pdf import Op, PdfFile, RendererPdf

import hicola_paths as hp

SLIDE_SIZE = (4.7, 3.1)    # inches (119 x 79 mm): the 128 x 96 mm frame minus side margins and the frame title
PNG_DPI = 300              # pixels per inch of the PNG at its slide size
SLIDE_DIR = hp.NOTEBOOK_DIR / "slides"


def _save_pdf(fig, path, scale):
    # Matplotlib's own PDF writing (FigureCanvasPdf.print_pdf) with one change: the page is
    # `scale` times the figure size and everything is drawn under a uniform scale transform.
    fig.canvas.draw()   # finish the layout as on screen
    saved_dpi = fig.dpi
    fig.dpi = 72        # 72 pdf points per inch, as in print_pdf
    width, height = fig.get_size_inches()
    # Fonts embedded as TrueType (42), not matplotlib's default Type 3, which many PDF viewers
    # draw without hinting (soft, blurry text on screen and projectors)
    try:
        with mpl.rc_context({"pdf.fonttype": 42}):
            _write_pdf(fig, path, scale, width, height)
    finally:
        fig.dpi = saved_dpi


def _write_pdf(fig, path, scale, width, height):
    file = PdfFile(str(path))
    try:
        file.newPage(width * scale, height * scale)   # writes the page size (MediaBox)
        file.width, file.height = width, height        # drawing (clipping) coordinates stay those of the figure
        file.output(scale, 0, 0, scale, 0, 0, Op.concat_matrix)
        renderer = MixedModeRenderer(fig, width, height, PNG_DPI / scale, RendererPdf(file, PNG_DPI / scale, height, width))
        fig.draw(renderer)
        renderer.finalize()
        file.finalize()
    finally:
        file.close()


def save(fig, name):
    """Save the figure as SLIDE_DIR/<name>.pdf and .png, scaled as a whole to fit in SLIDE_SIZE."""
    SLIDE_DIR.mkdir(parents=True, exist_ok=True)
    width, height = fig.get_size_inches()
    scale = min(SLIDE_SIZE[0] / width, SLIDE_SIZE[1] / height)
    pdf, png = SLIDE_DIR / f"{name}.pdf", SLIDE_DIR / f"{name}.png"
    _save_pdf(fig, pdf, scale)
    # PNG: PNG_DPI pixels per inch of the slide size, and the file says so (natural size = slide size)
    fig.savefig(png, dpi=PNG_DPI * scale, pil_kwargs={"dpi": (PNG_DPI, PNG_DPI)})
    print(f"Saved {pdf}, {png}: {width:.2f} x {height:.2f} in scaled by {scale:.3f} "
          f"to {width * scale:.2f} x {height * scale:.2f} in (a 10 pt label becomes {10 * scale:.1f} pt)")
    return pdf, png


def show(fig, name, slide):
    """Show the figure; with slide=True, also save it for the slides as <name>."""
    if slide:
        save(fig, name)
    plt.show()
