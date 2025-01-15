import tkinter as tk
from tkinter import messagebox
import graphviz
import os
import re

class FiniteAutomata:
    """Finite Automata for detecting genetic mutation patterns in DNA sequences."""

    def __init__(self, mutation_patterns):
        self.patterns = mutation_patterns
        self.automata = {pattern: self._build_states(pattern) for pattern in mutation_patterns}
        self.history = {}  # To store state transitions for each pattern

    def _build_states(self, pattern):
        """Builds a state machine for a specific pattern."""
        states = [{} for _ in range(len(pattern) + 1)]
        for i in range(len(pattern)):
            for nucleotide in 'ATGC':
                if nucleotide == pattern[i]:
                    states[i][nucleotide] = i + 1
                else:
                    states[i][nucleotide] = states[0].get(nucleotide, 0)
            states[i + 1] = states[i].copy()
        return states

    def recognize(self, sequence):
        """Checks for each mutation pattern in the sequence."""
        results = {}
        self.history = {}
        for pattern, states in self.automata.items():
            current_state = 0
            self.history[pattern] = [current_state]
            for char in sequence:
                current_state = states[current_state].get(char, 0)
                self.history[pattern].append(current_state)
                if current_state == len(pattern):
                    results[pattern] = True
                    break
            else:
                results[pattern] = False
        return results

    def visualize(self):
        """Visualizes the finite automaton using Graphviz with separate diagrams for each mutation pattern."""
        output_paths = []  # List to store file paths

        for pattern in self.patterns:
            # Create a new graph for the current pattern
            dot = graphviz.Digraph(comment=f'Finite Automaton for Pattern: {pattern}')
            
            # Get states for the current pattern
            states = self.automata[pattern]

            # Add nodes and transitions for each state
            for i in range(len(states)):
                # Add the start state and accepting state
                if i == 0:
                    dot.node(str(i), "Start State", shape="ellipse")
                elif i == len(states) - 1:
                    dot.node(str(i), "Accepting State", shape="doublecircle")
                else:
                    dot.node(str(i), f"State {i}", shape="circle")
                
                # Add transitions for each nucleotide
                for nucleotide, next_state in states[i].items():
                    dot.edge(str(i), str(next_state), label=nucleotide)
            
            # Save the diagram for the current pattern
            output_path = f'finite_automaton_{pattern}'
            dot.render(output_path, format='pdf', cleanup=True)
            output_paths.append(f'{output_path}.pdf')  # Append the generated file path

        return output_paths

def validate_sequence(sequence):
    """Validates that the sequence contains only valid DNA nucleotides (A, T, G, C)."""
    valid_nucleotides = {'A', 'T', 'G', 'C'}
    return all(nucleotide in valid_nucleotides for nucleotide in sequence)

# Global limits
MAX_SEQUENCE_LENGTH = 30000  # Maximum DNA sequence length
MAX_PATTERN_COUNT = 15        # Maximum number of patterns
MAX_PATTERN_LENGTH = 100      # Maximum length of a single pattern

mutation_disease_database = {
    "GAGGAGGAGGAG": "Sickle Cell Anemia",  # Mutation in the HBB gene causing sickle cell anemia
    "CGGCGGCAT": "Fragile X Syndrome",    # Mutation in the FMR1 gene causing Fragile X syndrome
    "TGGAGGAG": "Cystic Fibrosis",        # Mutation in the CFTR gene causing cystic fibrosis
    "AAAGGAGGAAG": "Huntington's Disease", # Mutation in the HTT gene causing Huntington's disease
    "CTGCTGCTG": "Myotonic Dystrophy",    # Mutation in the DMPK gene causing myotonic dystrophy
    "AGCAGCAGC": "Alpha-Thalassemia",     # Mutation in the HBA gene causing alpha-thalassemia
    "CAGGAGAGG": "Beta-Thalassemia",      # Mutation in the HBB gene causing beta-thalassemia
    "GTTGTTGTTG": "Tay-Sachs Disease",    # Mutation in the HEXA gene causing Tay-Sachs disease
    "GAAAGGAG": "Marfan Syndrome",        # Mutation in the FBN1 gene causing Marfan syndrome
    "CCAGGAG": "Retinitis Pigmentosa",    # Mutation in the RPGR gene causing retinitis pigmentosa
    "GGAGAGGA": "Phenylketonuria",        # Mutation in the PAH gene causing phenylketonuria
    "AGGAGGAAG": "Wilson's Disease",      # Mutation in the ATP7B gene causing Wilson's disease
    "TGTGTGAG": "Charcot-Marie-Tooth Disease", # Mutation in the PMP22 gene causing CMT disease
    "CGTGCAT": "Alzheimer's Disease",     # Mutation in the APP gene associated with Alzheimer's disease
    "AGAGAGAGAG": "Hemophilia A",         # Mutation in the F8 gene causing Hemophilia A
    "GCGCGCAG": "Prader-Willi Syndrome",  # Mutation in the 15q11-q13 region causing Prader-Willi syndrome
}

disease_info = {
    "Sickle Cell Anemia": "A blood disorder affecting hemoglobin, causing red blood cells to become misshapen.",
    "Fragile X Syndrome": "A genetic condition causing intellectual disability and behavioral challenges.",
    "Cystic Fibrosis": "A genetic disorder affecting the lungs and digestive system due to mucus buildup.",
    "Huntington's Disease": "A progressive brain disorder that causes uncontrolled movements and cognitive decline.",
    "Myotonic Dystrophy": "A muscle disorder causing prolonged muscle contractions and weakness.",
    "Alpha-Thalassemia": "A blood disorder reducing hemoglobin production, leading to anemia.",
    "Beta-Thalassemia": "A genetic blood disorder causing reduced hemoglobin and red blood cells.",
    "Tay-Sachs Disease": "A fatal disorder affecting the central nervous system, leading to motor and cognitive decline.",
    "Marfan Syndrome": "A connective tissue disorder affecting the heart, eyes, blood vessels, and skeleton.",
    "Retinitis Pigmentosa": "A group of eye diseases causing retinal degeneration and vision loss.",
    "Phenylketonuria": "A metabolic disorder where the body cannot break down the amino acid phenylalanine.",
    "Wilson's Disease": "A condition causing copper buildup in the body, affecting the liver and brain.",
    "Charcot-Marie-Tooth Disease": "A neurological disorder causing peripheral nerve damage and muscle weakness.",
    "Alzheimer's Disease": "A progressive brain disorder causing memory loss and cognitive decline.",
    "Hemophilia A": "A bleeding disorder caused by a deficiency in clotting factor VIII.",
    "Prader-Willi Syndrome": "A genetic disorder causing obesity, intellectual disability, and short stature."
}


def detect_mutations():
    """Detects mutations, displays results, and outputs analysis to a file."""
    dna_sequence = dna_entry.get().strip().upper()
    mutation_patterns_input = mutation_patterns_entry.get().strip().upper()

    # Split patterns using commas or newline as delimiters
    mutation_patterns = re.split(r"[,\n]+", mutation_patterns_input)
    mutation_patterns = [pattern.strip() for pattern in mutation_patterns if pattern.strip()]

    errors = []  # Collect error messages here

    # Check if both inputs are entered
    if not dna_sequence:
        errors.append("DNA sequence is required.")
    if not mutation_patterns:
        errors.append("At least one mutation pattern is required.")

    # Validate DNA sequence
    if dna_sequence and len(dna_sequence) > MAX_SEQUENCE_LENGTH:
        errors.append(f"DNA sequence is too long. Max length is {MAX_SEQUENCE_LENGTH} nucleotides.")
    elif dna_sequence and not validate_sequence(dna_sequence):
        errors.append("Invalid DNA sequence. Use only A, T, G, and C.")

    # Validate mutation patterns
    if mutation_patterns and len(mutation_patterns) > MAX_PATTERN_COUNT:
        errors.append(f"Too many mutation patterns. Max count is {MAX_PATTERN_COUNT}.")
    elif mutation_patterns and any(len(pattern) > MAX_PATTERN_LENGTH for pattern in mutation_patterns):
        errors.append(f"One or more mutation patterns are too long. Max length is {MAX_PATTERN_LENGTH} nucleotides.")
    elif mutation_patterns and not all(validate_sequence(pattern) for pattern in mutation_patterns):
        errors.append("Invalid mutation pattern(s). Use only A, T, G, and C.")

    # Display errors if any
    if errors:
        result_text.set("Error(s):\n" + "\n".join(errors))
        return

    # If all validations pass, proceed with detection
    fa = FiniteAutomata(mutation_patterns)
    results = fa.recognize(dna_sequence)

    # Check for overall mutation detection
    overall_mutation_detected = any(results.values())

    # Cross-reference the results with the disease database
    detected_diseases = []
    disease_descriptions = []
    for pattern in results:
        if results[pattern]:  # If the pattern is detected
            disease = mutation_disease_database.get(pattern, None)
            if disease:
                detected_diseases.append(disease)
                disease_descriptions.append(f"{disease}: {disease_info[disease]}")

    # Display the results in the GUI
    result_text.set(
        f"DNA Sequence: {dna_sequence}\n\n"
        f"Mutation Analysis:\n" +
        "\n".join(
            f"Pattern '{pattern}': {'Detected' if detected else 'Not Detected'}\n"
            f"State Transitions: {' -> '.join(map(str, fa.history[pattern]))}"
            for pattern, detected in results.items()
        ) +
        f"\n\nOverall Mutation Detected: {'Yes' if overall_mutation_detected else 'No'}\n\n"
        f"Detected Diseases/Syndromes:\n" +
        ("\n".join(disease_descriptions) if detected_diseases else "None")
    )

    # Generate FA diagrams for each pattern
    output_paths = fa.visualize()

    # Open all generated PDF files
    for path in output_paths:
        if os.name == 'posix':  # For Linux and macOS
            os.system(f'open {path}')  # macOS
        elif os.name == 'nt':  # For Windows
            os.startfile(path)  # Windows

    # Output results to a .txt file
    output_filename = 'mutation_analysis_results.txt'
    with open(output_filename, 'w') as file:
        file.write(f"DNA Sequence: {dna_sequence}\n\n")
        file.write("Mutation Analysis:\n")
        for pattern, detected in results.items():
            file.write(f"Pattern '{pattern}': {'Detected' if detected else 'Not Detected'}\n")
            file.write(f"State Transitions: {' -> '.join(map(str, fa.history[pattern]))}\n")
        file.write(f"\nOverall Mutation Detected: {'Yes' if overall_mutation_detected else 'No'}\n\n")
        file.write("Detected Diseases/Syndromes:\n")
        if detected_diseases:
            for disease, description in zip(detected_diseases, disease_descriptions):
                file.write(f"{disease}: {description}\n")
        else:
            file.write("None\n")

    # Notify the user that the results were saved
    messagebox.showinfo("Mutation Analysis", f"Results have been saved to {output_filename}")

# GUI setup
root = tk.Tk()
root.title("Genetic Mutation Detection")

frame = tk.Frame(root, padx=10, pady=10)
frame.pack()

# Input limitations description
limitations_label = tk.Label(
    frame,
    text=( 
        "Input Limitations:\n"
        "- Max DNA Sequence Length: 30,000 nucleotides\n"
        "- Max Mutation Pattern Length: 100 nucleotides\n"
        "- Max Number of Mutation Patterns: 15\n"
        "Note: Use only A, T, G, and C for DNA sequences and patterns."
    ),
    justify="left",
    fg="blue"
)
limitations_label.grid(row=0, columnspan=2, sticky="w", pady=5)

# DNA sequence input with nucleotide counter
tk.Label(frame, text="DNA Sequence:").grid(row=1, column=0, sticky="w", pady=5)

# Create a frame for the DNA sequence entry and counter
dna_frame = tk.Frame(frame)
dna_frame.grid(row=1, column=1, pady=5, sticky="w")

# DNA sequence entry
dna_entry = tk.Entry(dna_frame, width=50)
dna_entry.grid(row=0, column=0, pady=5)

# Nucleotide counter label
nucleotide_count = tk.StringVar(value="Nucleotides: 0")
nucleotide_count_label = tk.Label(dna_frame, textvariable=nucleotide_count, fg="green")
nucleotide_count_label.grid(row=0, column=1, padx=10)

# Function to update nucleotide count
def update_nucleotide_count(*args):
    dna_sequence = dna_entry.get().strip().upper()
    count = len(dna_sequence)
    if count > MAX_SEQUENCE_LENGTH:
        nucleotide_count.set(f"Nucleotides: {count} (Exceeds max length!)")
        nucleotide_count_label.config(fg="red")
    else:
        nucleotide_count.set(f"Nucleotides: {count}")
        nucleotide_count_label.config(fg="green")

# Bind the update function to the DNA entry field
dna_entry.bind("<KeyRelease>", update_nucleotide_count)

# Mutation patterns input
tk.Label(frame, text="Mutation Patterns (comma-separated):").grid(row=2, column=0, sticky="w", pady=5)
mutation_patterns_entry = tk.Entry(frame, width=50)
mutation_patterns_entry.grid(row=2, column=1, pady=5)

# Detect button
detect_button = tk.Button(frame, text="Detect Mutations", command=detect_mutations)
detect_button.grid(row=3, columnspan=2, pady=10)

# Result display
result_text = tk.StringVar()  # Initialize the StringVar to hold the result text
result_label = tk.Label(frame, textvariable=result_text, justify="left", anchor="w", wraplength=400)
result_label.grid(row=4, columnspan=2, sticky="w")

# Run the GUI
root.mainloop()
