from Bio import SeqIO
import glob, os

# GBK_DIR = r"C:\Users\anacl\Desktop\Programa_Corredor_Donas\POST_NATURE_DOCUMENTOS\Z_Arbol\MIBIG"
GBK_DIR = r"C:\Users\a.."
OUTPUT  = r"C:\Users\a.."

# awk '
#   /^>/{
#     h=$0
#     seen[h]++
#     if(seen[h]>1) print h "_dup" seen[h]
#     else print h
#     next
#   }
#   {print}
# ' terpene_core.aln.trim.faa > terpene_core.aln.trim.uniq.faa
#finally run the iq tree command
# iqtree2 \
#   -s terpene_core.aln.trim.uniq.faa \
#   -m MFP \
#   -bb 1000 \
#   -alrt 1000 \
#   -nt AUTO

# things antiSMASH uses inside gene_functions / sec_met_domain
KEYWORDS = [
    "phytoene_synt",
    "Lycopene_cycl",
    "terpene_cyclase",
    "Terpene_synth",
    "Terpene_synth_C",
    "trichodiene_synth",
    "NapT7",
    "TRI5",
    "fung_ggpps",
    "fung_ggpps2"
]

def is_terpene_core(feature) -> bool:
    gf = " ".join(feature.qualifiers.get("gene_functions", [])).lower()
    sd = " ".join(feature.qualifiers.get("sec_met_domain", [])).lower()
    text = gf + " " + sd
    return any(k in text for k in KEYWORDS)

n_files = 0
n_hits = 0

with open(OUTPUT, "w") as out:
    for gbk in glob.iglob(os.path.join(GBK_DIR, "**", "*.gbk"), recursive=True):
        n_files += 1
        try:
            for record in SeqIO.parse(gbk, "genbank"):
                for feat in record.features:
                    if feat.type != "CDS":
                        continue

                    if not is_terpene_core(feat):
                        continue

                    translation = feat.qualifiers.get("translation")
                    if not translation:
                        continue

                    locus = feat.qualifiers.get("locus_tag", ["unknown"])[0]
                    gf = " ".join(feat.qualifiers.get("gene_functions", [])).replace(" ", "_")
                    sd = " ".join(feat.qualifiers.get("sec_met_domain", [])).split(" ")[0].replace(" ", "_")

                    header = f">PREDICTED|{os.path.basename(gbk)}|{record.id}|{locus}|{sd}|{gf}"
                    out.write(header + "\n")
                    out.write(translation[0] + "\n")
                    n_hits += 1

        except Exception as e:
            print(f"Skipping {gbk}: {e}")

print(f"Processed {n_files} GBK files")
print(f"Extracted {n_hits} terpene core protein sequences")
print(f"Wrote: {OUTPUT}")
