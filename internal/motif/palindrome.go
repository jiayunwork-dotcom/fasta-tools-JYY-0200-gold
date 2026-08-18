// palindrome.go extends the motif package with palindrome detection. DNA
// palindromes are sequences that read the same on both strands (the sequence
// equals its own reverse complement). These are important for identifying
// restriction enzyme sites and hairpin structures.
package motif

import "strings"

// Palindrome represents a palindromic site found in a sequence.
type Palindrome struct {
	Start  int
	End    int
	Length int
	Seq    string
}

// FindPalindromes finds all palindromic sequences of length between minLen and
// maxLen in the given DNA sequence. Only even-length palindromes are considered
// (as DNA palindromes are by definition even-length).
func FindPalindromes(seq string, minLen, maxLen int) []Palindrome {
	upper := strings.ToUpper(seq)
	var results []Palindrome
	for length := minLen; length <= maxLen; length += 2 {
		for i := 0; i+length <= len(upper); i++ {
			sub := upper[i : i+length]
			if isPalindrome(sub) {
				results = append(results, Palindrome{
					Start:  i,
					End:    i + length,
					Length: length,
					Seq:    seq[i : i+length],
				})
			}
		}
	}
	return results
}

// isPalindrome checks if a DNA sequence equals its reverse complement.
func isPalindrome(s string) bool {
	n := len(s)
	for i := 0; i < n/2; i++ {
		if !isComplement(s[i], s[n-1-i]) {
			return false
		}
	}
	return true
}

// isComplement returns true if a and b are Watson-Crick complements.
func isComplement(a, b byte) bool {
	switch a {
	case 'A':
		return b == 'T'
	case 'T':
		return b == 'A'
	case 'C':
		return b == 'G'
	case 'G':
		return b == 'C'
	default:
		return false
	}
}

// LongestPalindrome returns the longest palindromic site, or nil if none found.
func LongestPalindrome(palindromes []Palindrome) *Palindrome {
	if len(palindromes) == 0 {
		return nil
	}
	best := &palindromes[0]
	for i := 1; i < len(palindromes); i++ {
		if palindromes[i].Length > best.Length {
			best = &palindromes[i]
		}
	}
	return best
}

// CountByLength groups palindromes by their length and returns counts.
func CountByLength(palindromes []Palindrome) map[int]int {
	counts := make(map[int]int)
	for _, p := range palindromes {
		counts[p.Length]++
	}
	return counts
}

// HasHairpotential checks if a palindromic site could form a hairpin loop
// structure. A minimum stem length of 4 is required.
func HasHairpotential(p Palindrome) bool {
	return p.Length >= 8 // At least 4bp stem on each side
}

// InvertedRepeats finds inverted repeat pairs within a sequence. An inverted
// repeat is a palindromic sequence where the two halves are separated by a
// spacer region.
func InvertedRepeats(seq string, armLen, minSpacer, maxSpacer int) []InvertedRepeat {
	upper := strings.ToUpper(seq)
	var results []InvertedRepeat
	for i := 0; i+armLen <= len(upper); i++ {
		arm1 := upper[i : i+armLen]
		arm1RC := revCompStr(arm1)
		for spacer := minSpacer; spacer <= maxSpacer; spacer++ {
			j := i + armLen + spacer
			if j+armLen > len(upper) {
				break
			}
			arm2 := upper[j : j+armLen]
			if arm2 == arm1RC {
				results = append(results, InvertedRepeat{
					Start1:    i,
					End1:      i + armLen,
					Start2:    j,
					End2:      j + armLen,
					ArmLength: armLen,
					Spacer:    spacer,
				})
			}
		}
	}
	return results
}

// InvertedRepeat represents a pair of inverted repeats.
type InvertedRepeat struct {
	Start1    int
	End1      int
	Start2    int
	End2      int
	ArmLength int
	Spacer    int
}

// revCompStr returns the reverse complement of a DNA string.
func revCompStr(s string) string {
	out := make([]byte, len(s))
	for i := range s {
		switch s[len(s)-1-i] {
		case 'A':
			out[i] = 'T'
		case 'T':
			out[i] = 'A'
		case 'C':
			out[i] = 'G'
		case 'G':
			out[i] = 'C'
		default:
			out[i] = 'N'
		}
	}
	return string(out)
}
