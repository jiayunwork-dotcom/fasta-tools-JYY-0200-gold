// codon_table.go extends the translate package with alternative genetic code
// tables and codon optimization utilities.
package translate

// GeneticCode represents an alternative genetic code (e.g., mitochondrial).
type GeneticCode struct {
	Name  string
	ID    int
	Table map[string]byte
}

// MitochondrialCode returns the vertebrate mitochondrial genetic code.
// Differences from standard: UGA=W, AGA=*, AGG=*, AUA=M.
func MitochondrialCode() *GeneticCode {
	table := make(map[string]byte)
	// Copy standard table.
	for k, v := range standardCodonTable {
		table[k] = v
	}
	// Apply mitochondrial differences.
	table["UGA"] = 'W' // Trp instead of Stop
	table["AGA"] = '*' // Stop instead of Arg
	table["AGG"] = '*' // Stop instead of Arg
	table["AUA"] = 'M' // Met instead of Ile
	return &GeneticCode{Name: "Vertebrate Mitochondrial", ID: 2, Table: table}
}

// YeastMitoCode returns the yeast mitochondrial genetic code.
func YeastMitoCode() *GeneticCode {
	table := make(map[string]byte)
	for k, v := range standardCodonTable {
		table[k] = v
	}
	table["CUU"] = 'T'
	table["CUC"] = 'T'
	table["CUA"] = 'T'
	table["CUG"] = 'T'
	table["UGA"] = 'W'
	return &GeneticCode{Name: "Yeast Mitochondrial", ID: 3, Table: table}
}

// TranslateWithCode translates an RNA sequence using a custom genetic code.
func TranslateWithCode(rna string, code *GeneticCode) (string, error) {
	if len(rna) == 0 {
		return "", ErrEmptyInput
	}
	upper := toUpper(rna)
	var protein []byte
	for i := 0; i+3 <= len(upper); i += 3 {
		codon := upper[i : i+3]
		aa, ok := code.Table[codon]
		if !ok {
			return "", ErrInvalidBase
		}
		protein = append(protein, aa)
	}
	return string(protein), nil
}

// ReverseTranslate returns all possible RNA codons for a given amino acid using
// the standard genetic code.
func ReverseTranslate(aminoAcid byte) []string {
	var codons []string
	for codon, aa := range standardCodonTable {
		if aa == aminoAcid {
			codons = append(codons, codon)
		}
	}
	return codons
}

// IsStartCodon checks if the given DNA codon is a standard start codon.
func IsStartCodon(codon string) bool {
	return toUpper(codon) == "ATG"
}

// IsStopCodon checks if the given DNA codon is a standard stop codon.
func IsStopCodon(codon string) bool {
	upper := toUpper(codon)
	rna := DNAToRNA(upper)
	return rna == "UAA" || rna == "UAG" || rna == "UGA"
}

// CodonCount returns the total number of codons (complete triplets) in a seq.
func CodonCount(seq string) int {
	return len(seq) / 3
}

// Degenerate returns the number of degenerate (synonymous) codons for an amino
// acid using the standard genetic code.
func Degenerate(aminoAcid byte) int {
	return len(ReverseTranslate(aminoAcid))
}
