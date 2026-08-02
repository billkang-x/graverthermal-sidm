"""Extract detailed setup sections from Schmidt et al. 2026."""
import fitz
import re

pdf_path = r"D:\graverthermal-sidm\2606.19428v1.pdf"
doc = fitz.open(pdf_path)

full_text = ""
for i in range(len(doc)):
    page = doc[i]
    text = page.get_text()
    full_text += f"\n\n===== PAGE {i+1} =====\n\n" + text

# Save full text for searching
with open(r"D:\graverthermal-sidm\notes\schmidt2026_fulltext.txt", "w", encoding="utf-8") as f:
    f.write(full_text)

# Find the key equation 19 context
print("="*80)
print("COOLING RATE EQUATION 19 CONTEXT:")
print("="*80)
idx = full_text.find("Eq. (19)")
if idx < 0:
    idx = full_text.find("8 √π")
if idx >= 0:
    print(full_text[max(0,idx-500):idx+1500])

print("\n\n" + "="*80)
print("APPENDIX D CONTEXT:")
print("="*80)
idx = full_text.find("Appendix D")
while idx >= 0:
    print("\n--- match at", idx, "---")
    print(full_text[max(0,idx-200):idx+800])
    idx = full_text.find("Appendix D", idx+1)
    if idx > 30000:
        break

print("\n\n" + "="*80)
print("SIGMA_T AND R_DISS VALUES:")
print("="*80)
for kw in ['sigma_T/m', 'r_diss =', 'r_diss=', 'rdiss', '3.6 kpc', '7.09', 't_core', 't_coll']:
    for m in re.finditer(re.escape(kw), full_text):
        s = max(0, m.start()-100)
        e = min(len(full_text), m.end()+200)
        print(f"\n[{kw}] ...{' '.join(full_text[s:e].split())}...")
