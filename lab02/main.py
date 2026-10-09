import random
import tkinter as tk
from tkinter import ttk
from tkinter.scrolledtext import ScrolledText

# Standard RNA codon table: RNA codon -> amino acid (3-letter name).
CODON_TABLE = {
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "UAU": "Tyr", "UAC": "Tyr", "UAA": "STOP", "UAG": "STOP",
    "UGU": "Cys", "UGC": "Cys", "UGA": "STOP", "UGG": "Trp",
    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}
BASES = "UCAG"
STOP_CODONS = {"UAA", "UAG", "UGA"}
START_CODON = "AUG"


def random_non_stop_codon():
    """Choose a codon that won't accidentally end/start an ORF."""
    choices = [
        codon for codon, aa in CODON_TABLE.items()
        if codon not in STOP_CODONS and codon != START_CODON
    ]
    return random.choice(choices)



def generate_sequence(length=100):
    """Generate a truly random RNA sequence of the requested length."""
    return "".join(random.choices(BASES, k=length))

def translate_orfs(sequence):
    """
    Find complete proteins in reading frame 0.
    An ORF begins at AUG and ends at the first in-frame stop codon.
    After a stop codon, scanning resumes for the next AUG.
    """
    proteins = []
    codons = [sequence[i:i + 3] for i in range(0, len(sequence) - 2, 3)]
    i = 0

    while i < len(codons):
        if codons[i] != START_CODON:
            i += 1
            continue

        start_index = i
        amino_acids = ["Met"]
        i += 1
        while i < len(codons):
            codon = codons[i]
            amino_acid = CODON_TABLE[codon]
            if amino_acid == "STOP":
                proteins.append({
                    "start": start_index * 3 + 1,
                    "end": i * 3 + 3,
                    "amino_acids": amino_acids,
                    "stop": codon,
                    "codons": codons[start_index:i + 1],
                })
                i += 1
                break
            amino_acids.append(amino_acid)
            i += 1
        else:
            # Ignore an AUG without a stop codon: it is not a complete protein.
            pass

    return proteins


class RNAProteinApp:
    def __init__(self, root):
        self.root = root
        self.root.title("RNA → Protein Translator")
        self.root.geometry("850x700")
        self.root.minsize(680, 560)

        outer = ttk.Frame(root, padding=18)
        outer.pack(fill="both", expand=True)

        ttk.Label(
            outer, text="RNA → Protein Translator", font=("TkDefaultFont", 18, "bold")
        ).pack(anchor="w")
        ttk.Label(
            outer,
            text=("Generate a random RNA sequence (A, U, C, G) and find complete "
                  "proteins that start with AUG and end with UAA, UAG, or UGA."),
            wraplength=790,
        ).pack(anchor="w", pady=(5, 14))

        controls = ttk.Frame(outer)
        controls.pack(fill="x", pady=(0, 8))
        ttk.Button(controls, text="Regenerate sequence", command=self.regenerate).pack(side="left")
        self.length_label = ttk.Label(controls, text="")
        self.length_label.pack(side="right")

        ttk.Label(outer, text="RNA sequence", font=("TkDefaultFont", 11, "bold")).pack(anchor="w")
        self.sequence_box = ScrolledText(outer, height=5, wrap="word", font=("Consolas", 12))
        self.sequence_box.pack(fill="x", pady=(4, 12))
        self.sequence_box.configure(state="disabled")

        ttk.Label(outer, text="Translated proteins", font=("TkDefaultFont", 11, "bold")).pack(anchor="w")
        self.results_box = ScrolledText(outer, height=14, wrap="word", font=("Consolas", 11))
        self.results_box.pack(fill="both", expand=True, pady=(4, 10))
        self.results_box.configure(state="disabled")

        ttk.Label(
            outer,
            text=("Note: translation is read in frame 1 (groups of three bases from the "
                  "first character). A final incomplete codon is ignored. Amino-acid "
                  "chains are shown using three-letter abbreviations; stop codons mark "
                  "the end and are not part of the chain."),
            wraplength=790,
        ).pack(anchor="w")

        self.regenerate()

    def regenerate(self):
        sequence = generate_sequence(100)
        proteins = translate_orfs(sequence)

        self.sequence_box.configure(state="normal")
        self.sequence_box.delete("1.0", "end")

        # Keep the sequence unspaced so nucleotide positions match text positions.
        self.sequence_box.insert("1.0", sequence)

        # Give each recognized complete protein its own highlight color.
        protein_colors = [
            ("#ffe08a", "#332600"),
            ("#a8e6a3", "#123b16"),
            ("#a8d8ff", "#102f4a"),
            ("#f6b3d2", "#4a1730"),
        ]
        for index, protein in enumerate(proteins):
            background, foreground = protein_colors[index % len(protein_colors)]
            tag = f"protein_{index + 1}"
            self.sequence_box.tag_configure(
                tag, background=background, foreground=foreground
            )
            # Protein positions are 1-based and include the stop codon.
            start_offset = protein["start"] - 1
            end_offset = protein["end"]
            start_text_index = f"1.0 + {start_offset} chars"
            end_text_index = f"1.0 + {end_offset} chars"
            self.sequence_box.tag_add(tag, start_text_index, end_text_index)

        self.sequence_box.configure(state="disabled")
        self.length_label.configure(text=f"{len(sequence)} nucleotides")

        self.results_box.configure(state="normal")
        self.results_box.delete("1.0", "end")
        if not proteins:
            self.results_box.insert("end", "No complete proteins found. Regenerate to try again.\n")
        else:
            self.results_box.insert("end", f"Found {len(proteins)} complete protein(s):\n\n")
            for index, protein in enumerate(proteins, start=1):
                aa = " – ".join(protein["amino_acids"])
                codons = " ".join(protein["codons"])
                self.results_box.insert(
                    "end",
                    f"Protein {index} (nucleotides {protein['start']}–{protein['end']})\n"
                    f"  Amino acids: {aa}\n"
                    f"  Codons:      {codons}\n"
                    f"  Stop codon:  {protein['stop']}\n\n"
                )
        self.results_box.configure(state="disabled")


if __name__ == "__main__":
    root = tk.Tk()
    app = RNAProteinApp(root)
    root.mainloop()
