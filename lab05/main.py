
import math
import random
import tkinter as tk
from tkinter import ttk, filedialog, messagebox


MIN_SAMPLE, MAX_SAMPLE = 100, 3000
LADDER = [100, 200, 300, 500, 700, 1000, 1500, 2000, 3000]


def read_fasta(filepath):
    """Read a FASTA file and return its first sequence (reused from lab03)."""
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


# ---------------------------------------------------------------- Exercise 1

def take_samples(sequence, count=10):
    """
    Take random fragments of the sequence, between 100 and 3000 bases.

    A fragment cannot be longer than the sequence itself, so the upper
    limit is the sequence length.
    """
    upper = min(MAX_SAMPLE, len(sequence))

    if upper < MIN_SAMPLE:
        raise ValueError(
            f"The sequence must have at least {MIN_SAMPLE} nucleotides."
        )

    samples = []

    for _ in range(count):
        length = random.randint(MIN_SAMPLE, upper)
        start = random.randint(0, len(sequence) - length)
        samples.append(sequence[start:start + length])

    return samples


def migration(length):
    """
    Relative distance travelled through the gel, from 0 (well) to 1 (bottom).

    Distance is roughly linear in log(size): short fragments meet less
    friction and travel further than long ones.
    """
    span = math.log10(MAX_SAMPLE) - math.log10(MIN_SAMPLE)
    return (math.log10(MAX_SAMPLE) - math.log10(length)) / span


# ---------------------------------------------------------------- Exercise 2

def smallest_period(unit):
    """Smallest p such that the unit is made of a repeated p-letter block."""
    for p in range(1, len(unit) + 1):
        if len(unit) % p == 0 and unit[:p] * (len(unit) // p) == unit:
            return p


def find_repeats(sequence, min_unit=3, max_unit=5, min_len=6, max_len=10):
    """
    Find tandem repeats: a unit of 3-5 bases repeated back to back,
    covering 6-10 bases in total.

    Only maximal runs are reported (the run cannot be extended left or
    right), and units that are themselves repetitive (e.g. AAA, ATAT)
    are skipped because they belong to a shorter unit.
    """
    repeats = []

    for k in range(min_unit, max_unit + 1):
        for i in range(len(sequence) - 2 * k + 1):
            unit = sequence[i:i + k]

            if smallest_period(unit) != k or set(unit) - set("ACGT"):
                continue

            # The periodic region must start here; otherwise this is the
            # same repeat read from a shifted frame (CAGCAG vs AGCAG...).
            if i > 0 and sequence[i - 1] == sequence[i - 1 + k]:
                continue

            copies = 1
            while sequence[i + copies * k:i + (copies + 1) * k] == unit:
                copies += 1

            length = copies * k

            if copies >= 2 and min_len <= length <= max_len:
                repeats.append((i + 1, unit, copies, length))

    repeats.sort()
    return repeats


# ---------------------------------------------------------------- Interface

class Lab05App:
    def __init__(self, root):
        self.root = root
        self.root.title("Gel Electrophoresis & Repeat Finder")
        self.root.geometry("950x640")
        self.root.minsize(800, 520)

        self.sequence = ""

        main = ttk.Frame(root, padding=18)
        main.pack(fill="both", expand=True)

        ttk.Label(
            main,
            text="Gel Electrophoresis & Repeat Finder",
            font=("TkDefaultFont", 17, "bold")
        ).pack(anchor="w")

        ttk.Label(
            main,
            text=(
                "Load a DNA sequence from NCBI in FASTA format "
                "(1000-3000 nucleotides)."
            ),
            wraplength=900
        ).pack(anchor="w", pady=(6, 14))

        controls = ttk.Frame(main)
        controls.pack(fill="x", pady=(0, 12))

        ttk.Button(
            controls, text="Open FASTA file", command=self.open_fasta
        ).pack(side="left")

        self.file_label = ttk.Label(controls, text="No file selected")
        self.file_label.pack(side="left", padx=12)

        notebook = ttk.Notebook(main)
        notebook.pack(fill="both", expand=True)

        # Exercise 1 tab
        gel_tab = ttk.Frame(notebook, padding=10)
        notebook.add(gel_tab, text="Ex 1: Electrophoresis")

        gel_controls = ttk.Frame(gel_tab)
        gel_controls.pack(fill="x", pady=(0, 8))

        ttk.Button(
            gel_controls, text="New random samples",
            command=self.run_gel
        ).pack(side="left")

        self.samples_label = ttk.Label(gel_controls, text="")
        self.samples_label.pack(side="left", padx=12)

        self.canvas = tk.Canvas(gel_tab, background="#1b1b2f")
        self.canvas.pack(fill="both", expand=True)
        self.canvas.bind("<Configure>", lambda _event: self.draw_gel())
        self.samples = []

        # Exercise 2 tab
        repeat_tab = ttk.Frame(notebook, padding=10)
        notebook.add(repeat_tab, text="Ex 2: Repeats")

        self.repeat_summary = ttk.Label(
            repeat_tab, text="Load a FASTA file to begin."
        )
        self.repeat_summary.pack(anchor="w", pady=(0, 8))

        columns = ("position", "unit", "unit_len", "copies", "length",
                   "region")
        headings = {
            "position": "Position", "unit": "Repeated unit",
            "unit_len": "Unit length", "copies": "Copies",
            "length": "Total length", "region": "Repeated region",
        }

        table_frame = ttk.Frame(repeat_tab)
        table_frame.pack(fill="both", expand=True)

        self.table = ttk.Treeview(
            table_frame, columns=columns, show="headings"
        )

        for column in columns:
            self.table.heading(column, text=headings[column])
            self.table.column(column, width=120, anchor="center")

        scrollbar = ttk.Scrollbar(
            table_frame, orient="vertical", command=self.table.yview
        )
        self.table.configure(yscrollcommand=scrollbar.set)

        self.table.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

    def open_fasta(self):
        filepath = filedialog.askopenfilename(
            title="Select a DNA FASTA file",
            filetypes=[
                ("FASTA files", "*.fasta *.fa *.fna *.fas"),
                ("Text files", "*.txt"),
                ("All files", "*.*")
            ]
        )

        if not filepath:
            return

        try:
            sequence = read_fasta(filepath)

            if not 1000 <= len(sequence) <= 3000:
                messagebox.showwarning(
                    "Sequence length",
                    f"The sequence has {len(sequence)} nucleotides; the "
                    "exercise asks for 1000-3000. Continuing anyway."
                )

            self.sequence = sequence
            self.file_label.configure(
                text=f"{filepath.split('/')[-1].split(chr(92))[-1]} "
                     f"({len(sequence):,} nt)"
            )
            self.run_gel()
            self.run_repeats()

        except (OSError, UnicodeError, ValueError) as error:
            messagebox.showerror("Error reading sequence", str(error))

    def run_gel(self):
        if not self.sequence:
            messagebox.showinfo("No sequence", "Load a FASTA file first.")
            return

        self.samples = take_samples(self.sequence)
        self.samples_label.configure(
            text="Sample sizes (bp): "
                 + ", ".join(str(len(s)) for s in self.samples)
        )
        self.draw_gel()

    def draw_gel(self):
        c = self.canvas
        c.delete("all")

        width, height = c.winfo_width(), c.winfo_height()
        lanes = ["Ladder"] + [f"S{i}" for i in range(1, 11)]
        lane_w = width / (len(lanes) + 1)
        top, bottom = 50, height - 30

        for index, name in enumerate(lanes):
            x = lane_w * (index + 1)
            c.create_text(x, 15, text=name, fill="white")
            # Well
            c.create_rectangle(x - lane_w * 0.35, 28, x + lane_w * 0.35,
                               36, outline="#888")

            if index == 0:
                sizes = LADDER
            elif self.samples:
                sizes = [len(self.samples[index - 1])]
            else:
                continue

            for size in sizes:
                y = top + migration(size) * (bottom - top)
                c.create_rectangle(x - lane_w * 0.35, y - 3,
                                   x + lane_w * 0.35, y + 3,
                                   fill="#7fffd4", outline="")

                if index == 0:
                    c.create_text(x - lane_w * 0.45, y, text=str(size),
                                  anchor="e", fill="white",
                                  font=("TkDefaultFont", 8))

    def run_repeats(self):
        repeats = find_repeats(self.sequence)

        for item in self.table.get_children():
            self.table.delete(item)

        for position, unit, copies, length in repeats:
            region = self.sequence[position - 1:position - 1 + length]
            self.table.insert("", "end", values=(
                position, unit, len(unit), copies, length, region
            ))

        self.repeat_summary.configure(
            text=f"Found {len(repeats)} tandem repeats "
                 "(unit 3-5 bases, total 6-10 bases). Positions are 1-based."
        )


if __name__ == "__main__":
    root = tk.Tk()
    app = Lab05App(root)
    root.mainloop()
