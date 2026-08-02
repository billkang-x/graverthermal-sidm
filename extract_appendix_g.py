"""Extract Appendix G from the PDF (rescaling symmetry)."""
import fitz

pdf_path = r"D:\graverthermal-sidm\2606.19428v1.pdf"
doc = fitz.open(pdf_path)

# Search for "Appendix G" in the last pages
for i in range(15, len(doc)):
    page = doc[i]
    text = page.get_text()
    if 'Appendix G' in text or 'Rescaling' in text or 'rescaling symmetry' in text.lower():
        print(f"\n===== PAGE {i+1} =====")
        print(text)
