#!/usr/bin/env python3
"""
Create a PowerPoint slide with benchmark comparison tables for Supplementary Figure 4.
"""

from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from lxml import etree
from pptx.oxml.ns import qn


def set_cell_format(cell, text, font_size=10, bold=False, font_color=None,
                    fill_color=None, alignment=PP_ALIGN.CENTER):
    """Format a table cell with specified properties."""
    cell.text = ""
    p = cell.text_frame.paragraphs[0]
    p.alignment = alignment
    run = p.add_run()
    run.text = text
    run.font.size = Pt(font_size)
    run.font.name = "Arial"
    if bold:
        run.font.bold = True
    if font_color:
        run.font.color.rgb = font_color
    if fill_color:
        cell.fill.solid()
        cell.fill.fore_color.rgb = fill_color
    # Vertical centering
    cell.vertical_anchor = MSO_ANCHOR.MIDDLE
    # Reduce margins for compact layout
    cell.margin_left = Emu(45720)   # ~0.05 inch
    cell.margin_right = Emu(45720)
    cell.margin_top = Emu(27432)    # ~0.03 inch
    cell.margin_bottom = Emu(27432)


def set_cell_border(cell, color_str="BFBFBF", width=Pt(0.5)):
    """Set thin borders on a cell using low-level XML manipulation."""
    tc = cell._tc
    tcPr = tc.get_or_add_tcPr()

    for border_name in ['lnL', 'lnR', 'lnT', 'lnB']:
        ln = tcPr.find(qn(f'a:{border_name}'))
        if ln is not None:
            tcPr.remove(ln)
        ln = etree.SubElement(tcPr, qn(f'a:{border_name}'))
        ln.set('w', str(int(width)))
        solidFill = etree.SubElement(ln, qn('a:solidFill'))
        srgbClr = etree.SubElement(solidFill, qn('a:srgbClr'))
        srgbClr.set('val', color_str)


def create_benchmark_slide():
    prs = Presentation()
    # Widescreen dimensions
    prs.slide_width = Inches(13.33)
    prs.slide_height = Inches(7.5)

    slide_layout = prs.slide_layouts[6]  # Blank layout
    slide = prs.slides.add_slide(slide_layout)

    # Colors
    DARK_BLUE = RGBColor(0x1F, 0x3A, 0x5F)
    WHITE = RGBColor(0xFF, 0xFF, 0xFF)
    LIGHT_GRAY = RGBColor(0xF2, 0xF2, 0xF2)
    ROW_WHITE = RGBColor(0xFF, 0xFF, 0xFF)
    BORDER_STR = "BFBFBF"
    BLACK = RGBColor(0x00, 0x00, 0x00)
    DARK_GRAY_TEXT = RGBColor(0x33, 0x33, 0x33)

    # ========== TITLE ==========
    left = Inches(0.5)
    top = Inches(0.25)
    width = Inches(12.33)
    height = Inches(0.6)
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.LEFT
    run = p.add_run()
    run.text = "Cross-Dataset Benchmarking of Pretrained Nuclei Segmentation Models"
    run.font.size = Pt(18)
    run.font.bold = True
    run.font.name = "Arial"
    run.font.color.rgb = DARK_BLUE

    # ========== TABLE 1: Cross-dataset benchmark ==========
    # Subtitle for table 1
    sub1_top = Inches(0.9)
    txBox1 = slide.shapes.add_textbox(Inches(0.5), sub1_top, Inches(10), Inches(0.35))
    tf1 = txBox1.text_frame
    p1 = tf1.paragraphs[0]
    p1.alignment = PP_ALIGN.LEFT
    r1 = p1.add_run()
    r1.text = "A. Cross-dataset evaluation (pretrained models, no fine-tuning)"
    r1.font.size = Pt(12)
    r1.font.bold = True
    r1.font.name = "Arial"
    r1.font.color.rgb = DARK_GRAY_TEXT

    rows, cols = 4, 5  # header + 3 data rows, 5 columns
    table1_left = Inches(0.5)
    table1_top = Inches(1.3)
    table1_width = Inches(12.33)
    table1_height = Inches(1.6)

    table1_shape = slide.shapes.add_table(rows, cols, table1_left, table1_top,
                                          table1_width, table1_height)
    table1 = table1_shape.table

    # Set column widths
    col_widths = [Inches(1.8), Inches(1.8), Inches(2.91), Inches(2.91), Inches(2.91)]
    for i, w in enumerate(col_widths):
        table1.columns[i].width = w

    # Header row
    headers1 = ["Model", "Training Data",
                "NuInsSeg (n=665)\nDice / AJI / PQ",
                "MoNuSeg (n=32)\nDice / AJI / PQ",
                "CryoNuSeg (n=30)\nDice / AJI / PQ"]

    for j, hdr in enumerate(headers1):
        cell = table1.cell(0, j)
        set_cell_format(cell, hdr, font_size=10, bold=True,
                        font_color=WHITE, fill_color=DARK_BLUE)
        set_cell_border(cell, BORDER_STR)

    # Data rows
    data1 = [
        ["HoVerNet", "PanNuke", "0.497 / 0.313 / 0.275", "0.790 / 0.442 / 0.410", "0.777 / 0.524 / 0.424"],
        ["StarDist", "MoNuSeg+TNBC", "0.452 / 0.284 / 0.278", "0.753 / 0.425 / 0.411", "0.733 / 0.501 / 0.416"],
        ["CellViT-SAM-H", "PanNuke", "0.659 / 0.458 / 0.403", "0.803 / 0.448 / 0.426", "0.792 / 0.546 / 0.453"],
    ]

    # CellViT is row index 2 (0-based in data) - best in all datasets
    bold_row_idx = 2

    for i, row_data in enumerate(data1):
        row_fill = LIGHT_GRAY if i % 2 == 0 else ROW_WHITE
        is_best = (i == bold_row_idx)
        for j, val in enumerate(row_data):
            cell = table1.cell(i + 1, j)
            # Bold the metric columns for CellViT (best model) and model name
            cell_bold = is_best and j >= 2
            name_bold = is_best and j == 0
            set_cell_format(cell, val, font_size=10, bold=(cell_bold or name_bold),
                            font_color=BLACK, fill_color=row_fill,
                            alignment=PP_ALIGN.CENTER if j >= 2 else PP_ALIGN.LEFT)
            set_cell_border(cell, BORDER_STR)

    # ========== TABLE 2: Skin-specific comparison ==========
    sub2_top = Inches(3.15)
    txBox2 = slide.shapes.add_textbox(Inches(0.5), sub2_top, Inches(10), Inches(0.35))
    tf2 = txBox2.text_frame
    p2 = tf2.paragraphs[0]
    p2.alignment = PP_ALIGN.LEFT
    r2 = p2.add_run()
    r2.text = "B. Skin tissue evaluation (qMAP dataset)"
    r2.font.size = Pt(12)
    r2.font.bold = True
    r2.font.name = "Arial"
    r2.font.color.rgb = DARK_GRAY_TEXT

    rows2, cols2 = 5, 6  # header + 4 data rows, 6 columns
    table2_left = Inches(0.5)
    table2_top = Inches(3.55)
    table2_width = Inches(9.0)
    table2_height = Inches(2.0)

    table2_shape = slide.shapes.add_table(rows2, cols2, table2_left, table2_top,
                                          table2_width, table2_height)
    table2 = table2_shape.table

    # Set column widths for table 2
    col_widths2 = [Inches(2.8), Inches(1.24), Inches(1.24), Inches(1.24), Inches(1.24), Inches(1.24)]
    for i, w in enumerate(col_widths2):
        table2.columns[i].width = w

    # Header row
    headers2 = ["Model", "Dice", "AJI", "DQ", "SQ", "PQ"]
    for j, hdr in enumerate(headers2):
        cell = table2.cell(0, j)
        set_cell_format(cell, hdr, font_size=10, bold=True,
                        font_color=WHITE, fill_color=DARK_BLUE)
        set_cell_border(cell, BORDER_STR)

    # Data rows
    data2 = [
        ["StarDist (pretrained)", "0.71", "0.66", "0.57", "0.62", "0.59"],
        ["CellViT (pretrained)", "0.72", "0.88", "0.70", "0.77", "0.88"],
        ["HoVerNet (pretrained)", "0.66", "0.40", "0.53", "0.76", "0.41"],
        ["HoVerNet-skin (retrained)", "0.73", "0.90", "0.72", "0.80", "0.90"],
    ]

    # Find best values per column
    best_per_col = {}
    for j in range(1, 6):
        max_val = max(float(data2[i][j]) for i in range(4))
        best_per_col[j] = max_val

    for i, row_data in enumerate(data2):
        row_fill = LIGHT_GRAY if i % 2 == 0 else ROW_WHITE
        for j, val in enumerate(row_data):
            cell = table2.cell(i + 1, j)
            # Bold the best value in each metric column
            cell_bold = False
            if j >= 1:
                if abs(float(val) - best_per_col[j]) < 0.001:
                    cell_bold = True
            # Also bold the model name if it's the retrained HoVerNet (overall best)
            if j == 0 and i == 3:
                cell_bold = True
            set_cell_format(cell, val, font_size=10, bold=cell_bold,
                            font_color=BLACK, fill_color=row_fill,
                            alignment=PP_ALIGN.CENTER if j >= 1 else PP_ALIGN.LEFT)
            set_cell_border(cell, BORDER_STR)

    # ========== CAPTION / FOOTNOTE ==========
    caption_top = Inches(5.8)
    caption_left = Inches(0.5)
    caption_width = Inches(12.33)
    caption_height = Inches(1.4)
    txBox3 = slide.shapes.add_textbox(caption_left, caption_top, caption_width, caption_height)
    tf3 = txBox3.text_frame
    tf3.word_wrap = True
    p3 = tf3.paragraphs[0]
    p3.alignment = PP_ALIGN.LEFT
    p3.space_before = Pt(0)
    p3.space_after = Pt(0)
    r3 = p3.add_run()
    r3.text = (
        "All pretrained models evaluated using publicly released weights without fine-tuning "
        "on evaluation datasets. MoNuSeg restricted to 32 TCGA images (1000\u00d71000); 50 additional "
        "256\u00d7256 patches excluded due to asymmetric padding requirements across models. "
        "Metrics: Dice = binary overlap; AJI = Aggregated Jaccard Index; PQ = Panoptic Quality "
        "(IoU > 0.5). Bold = best performance per dataset."
    )
    r3.font.size = Pt(8)
    r3.font.name = "Arial"
    r3.font.color.rgb = RGBColor(0x55, 0x55, 0x55)
    r3.font.italic = True

    # Save
    output_path = "/home/kyu_insilica_co/qMAP/previously_revised_manuscript/supp_fig4_benchmark_table.pptx"
    prs.save(output_path)
    print(f"Saved to: {output_path}")
    print("Slide contents:")
    print(f"  - Title: Cross-Dataset Benchmarking of Pretrained Nuclei Segmentation Models")
    print(f"  - Table A: 3 models x 3 datasets (Dice/AJI/PQ)")
    print(f"  - Table B: 4 models x 5 metrics (skin-specific)")
    print(f"  - Caption footnote with methodology notes")
    print(f"  - Slide dimensions: 13.33 x 7.5 inches (widescreen)")


if __name__ == "__main__":
    create_benchmark_slide()
