
import random
import tkinter as tk
from tkinter import ttk, filedialog, messagebox


IR_LEN = 10        # length of each terminal inverted repeat
TSD_LEN = 5        # target site duplication (direct repeat) length
CORE_LEN = (20, 40)
HOST_LEN = (100, 140)   # keeps the final sequence within 200-400 b
COMPLEMENT = str.maketrans("ACGT", "TGCA")
MODES = ("Separate", "Nested", "Intersecting")


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def random_dna(length):
    return "".join(random.choice("ACGT") for _ in range(length))


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


# ---------------------------------------------------------------- Exercise 1

def insert_transposon(sequence, elements, position, name):
    """
    Insert a new transposon at `position`, like in the figure:
    the TSD_LEN bases at the target site are duplicated on both sides,
    and the element is flanked by an inverted repeat pair.

    `elements` holds dicts with 0-based [start, end) of every previously
    inserted element; they are shifted / stretched in place.
    """
    ir = random_dna(IR_LEN)
    transposon = ir + random_dna(random.randint(*CORE_LEN)) \
        + reverse_complement(ir)
    added = len(transposon) + TSD_LEN
    cut = position + TSD_LEN

    # Bases before `cut` keep their index; bases from `position` on move
    # right by `added` (the TSD bases exist in both copies).
    for element in elements:
        start, end = element["start"], element["end"]

        if position <= start:
            element["start"], element["end"] = start + added, end + added
        elif position < end:
            element["end"] = end + added
            # Landing on one of the inverted repeats destroys the element.
            if cut < start + IR_LEN or position > end - IR_LEN:
                element["intact"] = False

    elements.append({
        "name": name, "start": cut, "end": cut + len(transposon),
        "intact": True,
    })

    return sequence[:cut] + transposon + sequence[position:]


def build_sequence(count=3, mode="Separate"):
    """
    Build an artificial 200-400 b sequence with `count` transposons.

    Separate:     target sites only in host DNA, between elements.
    Nested:       each new element lands inside the core of the previous.
    Intersecting: each new element lands on the right inverted repeat of
                  the previous one, cutting 4-5 of its bases away.
    """
    sequence = random_dna(random.randint(*HOST_LEN))
    elements = []

    for number in range(1, count + 1):
        previous = elements[-1] if elements else None

        if mode == "Nested" and previous:
            position = random.randint(previous["start"] + IR_LEN,
                                      previous["end"] - IR_LEN - TSD_LEN)
        elif mode == "Intersecting" and previous:
            position = random.randint(previous["end"] - IR_LEN + 4,
                                      previous["end"] - TSD_LEN)
        else:
            free = [p for p in range(len(sequence) - TSD_LEN + 1)
                    if all(p + TSD_LEN <= e["start"] or p >= e["end"]
                           for e in elements)]
            position = random.choice(free)

        sequence = insert_transposon(sequence, elements, position,
                                     f"TE{number}")

    return sequence, elements


# ---------------------------------------------------------------- Exercise 2

def find_transposons(sequence, ir_len=IR_LEN, tsd_len=TSD_LEN):
    """
    Detect transposons without knowing their sequences:
    a k-mer whose reverse complement appears further downstream
    (inverted repeats), with identical tsd_len bases just outside both
    ends (direct repeats from the target site duplication).

    Returns 0-based (start, end) pairs, end exclusive, IRs included.
    """
    found = []

    for i in range(tsd_len, len(sequence) - ir_len + 1):
        left_ir = sequence[i:i + ir_len]
        right_ir = reverse_complement(left_ir)
        tsd = sequence[i - tsd_len:i]
        j = sequence.find(right_ir, i + ir_len)

        while j != -1:
            end = j + ir_len
            if sequence[end:end + tsd_len] == tsd:
                found.append((i, end))
            j = sequence.find(right_ir, j + 1)

    return found


# ---------------------------------------------------------------- Interface

class Lab06App:
    def __init__(self, root):
        self.root = root
        self.root.title("Transposable Element Finder")
        self.root.geometry("1000x680")
        self.root.minsize(820, 560)

        self.sequence = ""
        self.elements = []

        main = ttk.Frame(root, padding=18)
        main.pack(fill="both", expand=True)

        ttk.Label(
            main,
            text="Transposable Element Finder",
            font=("TkDefaultFont", 17, "bold")
        ).pack(anchor="w")

        ttk.Label(
            main,
            text=(
                f"Generates a 200-400 b DNA sequence with 3-4 transposons "
                f"({IR_LEN} b inverted repeats, {TSD_LEN} b target site "
                "duplications) and detects them back from the sequence "
                "alone. Positions are 1-based, inclusive."
            ),
            wraplength=950
        ).pack(anchor="w", pady=(6, 14))

        controls = ttk.Frame(main)
        controls.pack(fill="x", pady=(0, 10))

        ttk.Label(controls, text="Transposons:").pack(side="left")
        self.count = tk.IntVar(value=3)
        ttk.Spinbox(controls, from_=3, to=4, width=4,
                    textvariable=self.count, state="readonly"
                    ).pack(side="left", padx=(4, 12))

        ttk.Label(controls, text="Placement:").pack(side="left")
        self.mode = tk.StringVar(value=MODES[0])
        ttk.Combobox(controls, values=MODES, width=12,
                     textvariable=self.mode, state="readonly"
                     ).pack(side="left", padx=(4, 12))

        ttk.Button(controls, text="Generate",
                   command=self.generate).pack(side="left")
        ttk.Button(controls, text="Save FASTA",
                   command=self.save_fasta).pack(side="left", padx=6)
        ttk.Button(controls, text="Open FASTA",
                   command=self.open_fasta).pack(side="left")

        self.info = ttk.Label(main, text="")
        self.info.pack(anchor="w", pady=(0, 6))

        self.text = tk.Text(main, height=8, wrap="char",
                            font=("TkFixedFont", 11))
        self.text.pack(fill="x")
        self.text.tag_configure("te", background="#d9c8f0")
        self.text.tag_configure("ir", background="#9b6fd6",
                                foreground="white")
        self.text.tag_configure("tsd", foreground="#d62828",
                                underline=True)

        ttk.Label(
            main,
            text="Legend: purple = detected transposon, dark = inverted "
                 "repeats, red underlined = direct repeats (TSD).",
        ).pack(anchor="w", pady=(4, 10))

        tables = ttk.Frame(main)
        tables.pack(fill="both", expand=True)

        self.truth = self.make_table(
            tables, "Inserted (ground truth)",
            ("name", "start", "end", "length", "status"),
        )
        self.found = self.make_table(
            tables, "Detected",
            ("start", "end", "length", "match"),
        )

    def make_table(self, parent, title, columns):
        frame = ttk.LabelFrame(parent, text=title, padding=6)
        frame.pack(side="left", fill="both", expand=True, padx=4)
        table = ttk.Treeview(frame, columns=columns, show="headings")

        for column in columns:
            table.heading(column, text=column.capitalize())
            table.column(column, width=80, anchor="center")

        table.pack(fill="both", expand=True)
        return table

    def generate(self):
        self.sequence, self.elements = build_sequence(
            self.count.get(), self.mode.get()
        )
        self.show()

    def save_fasta(self):
        if not self.sequence:
            messagebox.showinfo("No sequence", "Generate a sequence first.")
            return

        filepath = filedialog.asksaveasfilename(
            defaultextension=".fasta",
            filetypes=[("FASTA files", "*.fasta")]
        )

        if filepath:
            with open(filepath, "w", encoding="utf-8") as file:
                file.write(">artificial_transposons\n")
                for i in range(0, len(self.sequence), 60):
                    file.write(self.sequence[i:i + 60] + "\n")

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
            self.sequence = read_fasta(filepath)
            self.elements = []   # no ground truth for loaded files
            self.show()
        except (OSError, UnicodeError, ValueError) as error:
            messagebox.showerror("Error reading sequence", str(error))

    def show(self):
        detected = find_transposons(self.sequence)
        truth = {(e["start"], e["end"]): e["name"] for e in self.elements}

        self.text.delete("1.0", "end")
        self.text.insert("1.0", self.sequence)

        for start, end in detected:
            spans = {
                "te": (start, end),
                "ir": [(start, start + IR_LEN), (end - IR_LEN, end)],
                "tsd": [(start - TSD_LEN, start), (end, end + TSD_LEN)],
            }
            for tag, ranges in spans.items():
                for a, b in ([ranges] if tag == "te" else ranges):
                    self.text.tag_add(tag, f"1.0+{a}c", f"1.0+{b}c")

        for table in (self.truth, self.found):
            table.delete(*table.get_children())

        for e in self.elements:
            hit = (e["start"], e["end"]) in detected
            status = ("detected" if hit else "missed") if e["intact"] \
                else ("broken, detected" if hit else "broken (IR cut)")
            self.truth.insert("", "end", values=(
                e["name"], e["start"] + 1, e["end"], e["end"] - e["start"],
                status,
            ))

        for start, end in detected:
            self.found.insert("", "end", values=(
                start + 1, end, end - start,
                truth.get((start, end), "false positive"
                          if self.elements else "-"),
            ))

        self.info.configure(
            text=f"Sequence length: {len(self.sequence)} b   |   "
                 f"inserted: {len(self.elements)}   |   "
                 f"detected: {len(detected)}"
        )


if __name__ == "__main__":
    root = tk.Tk()
    app = Lab06App(root)
    root.mainloop()
