"""Extract key simulation setup sections from Schmidt et al. 2026."""
import sys
import fitz  # PyMuPDF

pdf_path = r"D:\graverthermal-sidm\2606.19428v1.pdf"

doc = fitz.open(pdf_path)
print(f"Total pages: {len(doc)}")

# Extract text from pages 1-15 (intro + setup)
full_text = ""
for i in range(min(len(doc), 25)):
    page = doc[i]
    text = page.get_text()
    full_text += f"\n\n===== PAGE {i+1} =====\n\n" + text

# Search for key terms
import re
keywords = ['NFW', 'scale radius', 'r_s', 'rho_s', 'sigma_T', 'r_diss',
            'initial condition', 'Appendix D', 'cooling rate',
            'M_sun', 'kpc', 'km/s', 't_coll', 't_core',
            'n_shells', 'fluid', 'gravothermal']

print("\n" + "="*80)
print("Searching for key setup parameters...")
print("="*80)

for kw in keywords:
    # find all occurrences with surrounding context
    pattern = re.compile(r'.{0,100}' + re.escape(kw) + r'.{0,100}', re.IGNORECASE | re.DOTALL)
    matches = pattern.findall(full_text)
    if matches:
        print(f"\n--- {kw} ({len(matches)} matches) ---")
        for m in matches[:3]:  # show first 3
            # clean up whitespace
            clean = ' '.join(m.split())
            print(f"  ...{clean}...")
