
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from collections import Counter


# Standard RNA codon table reused from the previous lab.
# RNA codon -> amino acid (three-letter abbreviation).
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
}

# Convert the RNA codon table into a DNA codon table.
# Example: AUG becomes ATG, while its amino acid remains Met.
DNA_CODON_TABLE = {
    codon.replace("U", "T"): amino_acid
    for codon, amino_acid in CODON_TABLE.items()
}


def read_fasta(filepath):
    """Read a FASTA file and return its sequences."""
    sequences = []
    current_sequence = []

    with open(filepath, "r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_sequence:
                    sequences.append("".join(current_sequence))
                    current_sequence = []
            else:
                # Remove whitespace and normalize lowercase bases.
                current_sequence.append("".join(line.split()).upper())

    # Save the final FASTA record.
    if current_sequence:
        sequences.append("".join(current_sequence))

    if not sequences:
        raise ValueError("The selected file contains no FASTA sequences.")

    return sequences


def count_codons(sequences):
    """
    Count every overlapping DNA triplet.

    The window advances by one nucleotide.
    Windows containing ambiguous or invalid bases are skipped.
    """
    counts = Counter()
    valid_windows = 0
    total_nucleotides = sum(len(seq) for seq in sequences)

    for sequence in sequences:
        for i in range(len(sequence) - 2):
            codon = sequence[i:i + 3]

            if codon in DNA_CODON_TABLE:
                counts[codon] += 1
                valid_windows += 1

    return counts, valid_windows, total_nucleotides


def get_top_codons(counts, total_windows, limit=10):
    """Return the most frequent codons and their amino acids."""
    results = []

    # Sort by count descending, then alphabetically for ties.
    for codon, count in sorted(
        counts.items(),
        key=lambda item: (-item[1], item[0])
    )[:limit]:
        rna_codon = codon.replace("T", "U")
        amino_acid = CODON_TABLE[rna_codon]
        percentage = count / total_windows * 100

        results.append(
            (codon, rna_codon, amino_acid, count, percentage)
        )

    return results


class GenomeCodonApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Genome Codon Frequency Analyzer")
        self.root.geometry("850x560")
        self.root.minsize(700, 450)

        main = ttk.Frame(root, padding=18)
        main.pack(fill="both", expand=True)

        ttk.Label(
            main,
            text="Genome Codon Frequency Analyzer",
            font=("TkDefaultFont", 17, "bold")
        ).pack(anchor="w")

        ttk.Label(
            main,
            text=(
                "Load a DNA genome in FASTA format to find the 10 most "
                "frequent overlapping codons and their amino acids."
            ),
            wraplength=800
        ).pack(anchor="w", pady=(6, 14))

        controls = ttk.Frame(main)
        controls.pack(fill="x", pady=(0, 12))

        ttk.Button(
            controls,
            text="Open FASTA file",
            command=self.open_fasta
        ).pack(side="left")

        self.file_label = ttk.Label(
            controls, text="No file selected"
        )
        self.file_label.pack(side="left", padx=12)

        self.summary = ttk.Label(
            main,
            text="Load a FASTA file to begin."
        )
        self.summary.pack(anchor="w", pady=(0, 12))

        columns = (
            "rank", "dna", "rna", "amino_acid",
            "count", "frequency"
        )

        self.table = ttk.Treeview(
            main,
            columns=columns,
            show="headings",
            height=10
        )

        headings = {
            "rank": "Rank",
            "dna": "DNA codon",
            "rna": "RNA codon",
            "amino_acid": "Amino acid",
            "count": "Occurrences",
            "frequency": "Frequency (%)"
        }

        widths = {
            "rank": 60,
            "dna": 100,
            "rna": 100,
            "amino_acid": 140,
            "count": 130,
            "frequency": 140
        }

        for column in columns:
            self.table.heading(column, text=headings[column])
            self.table.column(
                column,
                width=widths[column],
                anchor="center"
            )

        scrollbar = ttk.Scrollbar(
            main,
            orient="vertical",
            command=self.table.yview
        )
        self.table.configure(yscrollcommand=scrollbar.set)

        self.table.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        # Keep the footer below the results table.
        footer = ttk.Frame(main)
        footer.pack(side="bottom", fill="x", pady=(12, 0))

        ttk.Label(
            footer,
            text=(
                "Method: overlapping windows of 3 nucleotides, "
                "sliding by 1 position. Incomplete triplets and "
                "windows containing ambiguous bases are ignored."
            ),
            wraplength=800
        ).pack(anchor="w")

    def open_fasta(self):
        filepath = filedialog.askopenfilename(
            title="Select a genome FASTA file",
            filetypes=[
                ("FASTA files", "*.fasta *.fa *.fna *.fas"),
                ("Text files", "*.txt"),
                ("All files", "*.*")
            ]
        )

        if not filepath:
            return

        try:
            sequences = read_fasta(filepath)
            counts, total_windows, total_bases = count_codons(
                sequences
            )

            if total_windows == 0:
                raise ValueError(
                    "No valid DNA triplets were found in the file."
                )

            results = get_top_codons(counts, total_windows)

            # Clear the previous results.
            for item in self.table.get_children():
                self.table.delete(item)

            for rank, result in enumerate(results, start=1):
                dna, rna, amino_acid, count, percentage = result

                self.table.insert(
                    "",
                    "end",
                    values=(
                        rank,
                        dna,
                        rna,
                        amino_acid,
                        count,
                        f"{percentage:.4f}%"
                    )
                )

            self.file_label.configure(
                text=filepath.split("/")[-1].split("\\")[-1]
            )

            self.summary.configure(
                text=(
                    f"FASTA records: {len(sequences)} | "
                    f"Nucleotides: {total_bases:,} | "
                    f"Valid triplets counted: {total_windows:,}"
                )
            )

        except (OSError, UnicodeError, ValueError) as error:
            messagebox.showerror(
                "Error reading genome",
                str(error)
            )


if __name__ == "__main__":
    root = tk.Tk()
    app = GenomeCodonApp(root)
    root.mainloop()
