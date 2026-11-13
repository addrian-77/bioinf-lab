import math
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox


# name -> (recognized site 5'->3', cut offset on the top strand)
# All five sites are palindromes, so scanning the top strand is enough.
ENZYMES = {
    "EcoRI":   ("GAATTC", 1),   # G^AATTC
    "BamHI":   ("GGATCC", 1),   # G^GATCC
    "HindIII": ("AAGCTT", 1),   # A^AGCTT
    "TaqI":    ("TCGA", 1),     # T^CGA
    "HaeIII":  ("GGCC", 2),     # GG^CC (blunt)
}
ALL = "All 5"
DEFAULT_FASTA = Path(__file__).parent.parent / "lab07" / "ecoli_k12_mg1655.fasta"

# agarose gel resolves roughly this range; anything outside piles at the edges
GEL_MIN_BP, GEL_MAX_BP = 100, 20000
LADDER = [250, 500, 750, 1000, 1500, 2000, 2500, 3000, 4000, 5000,
          6000, 8000, 10000, 20000]


def read_fasta(filepath):
    """Read a FASTA file and return its first sequence (reused from lab05)."""
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
                current_sequence.append("".join(line.split()).upper())

    if current_sequence:
        sequences.append("".join(current_sequence))

    if not sequences:
        raise ValueError("The selected file contains no FASTA sequences.")

    return sequences[0]


# ---------------------------------------------------------------- Digestion

def find_cuts(sequence, site, offset):
    """0-based positions where the top strand is cut (cut is before index)."""
    cuts = []
    i = sequence.find(site)

    while i != -1:
        cuts.append(i + offset)
        i = sequence.find(site, i + 1)

    return cuts


def fragments(length, cuts):
    """Fragment lengths of a linear molecule cut at the given positions."""
    bounds = [0] + sorted(set(cuts)) + [length]
    return [b - a for a, b in zip(bounds, bounds[1:])]


def digest(sequence):
    """Return {lane name: cut positions} for every enzyme plus all of them."""
    result = {name: find_cuts(sequence, site, offset)
              for name, (site, offset) in ENZYMES.items()}
    result[ALL] = sorted(set().union(*result.values()))
    return result


# ---------------------------------------------------------------- Interface

class Lab07App:
    def __init__(self, root):
        self.root = root
        self.root.title("Restriction Digest")
        self.root.geometry("1150x720")
        self.root.minsize(900, 600)

        self.sequence = ""
        self.cuts = {}

        main = ttk.Frame(root, padding=18)
        main.pack(fill="both", expand=True)

        ttk.Label(
            main,
            text="Restriction Digest & Gel Electrophoresis",
            font=("TkDefaultFont", 17, "bold")
        ).pack(anchor="w")

        ttk.Label(
            main,
            text=(
                "Digests a DNA sequence with EcoRI, BamHI, HindIII, TaqI and "
                "HaeIII. The sequence is treated as linear. Positions are "
                "1-based: a cut at position p falls between bases p and p+1."
            ),
            wraplength=1100
        ).pack(anchor="w", pady=(6, 10))

        controls = ttk.Frame(main)
        controls.pack(fill="x", pady=(0, 8))
        ttk.Button(controls, text="Open FASTA",
                   command=self.open_fasta).pack(side="left")
        self.info = ttk.Label(controls, text="")
        self.info.pack(side="left", padx=12)

        body = ttk.Frame(main)
        body.pack(fill="both", expand=True)

        left = ttk.Frame(body)
        left.pack(side="left", fill="both", expand=True)

        self.summary = self.make_table(
            left, "Enzymes (select one to list its fragments)",
            ("enzyme", "site", "cleavage", "cuts", "fragments"), height=6,
        )
        self.summary.bind("<<TreeviewSelect>>", self.show_fragments)

        self.detail = self.make_table(
            left, "Fragments", ("#", "start", "end", "length"), height=12,
        )

        gel_frame = ttk.LabelFrame(body, text="Simulated agarose gel",
                                   padding=6)
        gel_frame.pack(side="left", fill="both", padx=(8, 0))
        self.gel = tk.Canvas(gel_frame, width=520, height=560,
                             background="#111", highlightthickness=0)
        self.gel.pack(fill="both", expand=True)

        if DEFAULT_FASTA.exists():
            self.load(DEFAULT_FASTA)

    def make_table(self, parent, title, columns, height):
        frame = ttk.LabelFrame(parent, text=title, padding=6)
        frame.pack(fill="both", expand=True, pady=(0, 8))
        table = ttk.Treeview(frame, columns=columns, show="headings",
                             height=height)
        scroll = ttk.Scrollbar(frame, command=table.yview)
        table.configure(yscrollcommand=scroll.set)

        for column in columns:
            table.heading(column, text=column.capitalize())
            table.column(column, width=90, anchor="center")

        scroll.pack(side="right", fill="y")
        table.pack(fill="both", expand=True)
        return table

    def open_fasta(self):
        filepath = filedialog.askopenfilename(
            title="Select a DNA FASTA file",
            filetypes=[
                ("FASTA files", "*.fasta *.fa *.fna *.fas"),
                ("Text files", "*.txt"),
                ("All files", "*.*")
            ]
        )

        if filepath:
            self.load(filepath)

    def load(self, filepath):
        try:
            self.sequence = read_fasta(filepath)
        except (OSError, UnicodeError, ValueError) as error:
            messagebox.showerror("Error reading sequence", str(error))
            return

        self.cuts = digest(self.sequence)
        self.info.configure(
            text=f"{Path(filepath).name}   |   length: "
                 f"{len(self.sequence):,} b"
        )

        self.summary.delete(*self.summary.get_children())
        for name, cuts in self.cuts.items():
            if name == ALL:
                site = cleavage = "-"
            else:
                site, offset = ENZYMES[name]
                cleavage = site[:offset] + "^" + site[offset:]
            self.summary.insert("", "end", iid=name, values=(
                name, site, cleavage, len(cuts), len(cuts) + 1,
            ))

        self.summary.selection_set("EcoRI")
        self.draw_gel()

    def show_fragments(self, _event=None):
        selected = self.summary.selection()
        if not selected:
            return

        self.detail.delete(*self.detail.get_children())
        bounds = [0] + self.cuts[selected[0]] + [len(self.sequence)]

        for number, (a, b) in enumerate(zip(bounds, bounds[1:]), 1):
            self.detail.insert("", "end", values=(number, a + 1, b, b - a))

    def draw_gel(self):
        canvas = self.gel
        canvas.delete("all")
        width, height = int(canvas["width"]), int(canvas["height"])
        top, bottom = 50, height - 20
        lanes = ["Ladder"] + list(self.cuts)
        margin = 40   # room for ladder size labels
        lane_w = (width - margin) / len(lanes)

        def y_of(size):
            # migration distance is ~ linear in log(size)
            size = min(max(size, GEL_MIN_BP), GEL_MAX_BP)
            span = math.log10(GEL_MAX_BP) - math.log10(GEL_MIN_BP)
            return top + (math.log10(GEL_MAX_BP) - math.log10(size)) \
                / span * (bottom - top)

        for index, name in enumerate(lanes):
            x0 = margin + index * lane_w + 8
            x1 = margin + (index + 1) * lane_w - 8
            canvas.create_rectangle(x0, top - 14, x1, top - 8,
                                    outline="#555")   # well
            canvas.create_text((x0 + x1) / 2, top - 28, text=name,
                               fill="white", font=("TkDefaultFont", 9))

            if name == "Ladder":
                sizes = LADDER
                for size in sizes:
                    canvas.create_text(x0 - 2, y_of(size), text=str(size),
                                       fill="#888", anchor="e",
                                       font=("TkDefaultFont", 7))
            else:
                sizes = fragments(len(self.sequence), self.cuts[name])

            # brightness ~ total DNA mass landing on each pixel row
            rows = {}
            for size in sizes:
                row = round(y_of(size))
                rows[row] = rows.get(row, 0) + size
            heaviest = max(rows.values())

            for row, mass in rows.items():
                level = 60 + int(195 * math.sqrt(mass / heaviest))
                color = f"#{level // 3:02x}{level:02x}{level // 2:02x}"
                canvas.create_line(x0, row, x1, row, fill=color, width=2)


if __name__ == "__main__":
    root = tk.Tk()
    app = Lab07App(root)
    root.mainloop()
