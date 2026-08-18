// Package digest provides restriction enzyme digestion simulation. Given a DNA
// sequence and a set of restriction enzyme cut sites, it finds all recognition
// sites and reports the resulting fragments.
package digest

import (
	"errors"
	"sort"
	"strings"
)

// ErrNoEnzyme is returned when no enzyme is specified.
var ErrNoEnzyme = errors.New("digest: no enzyme specified")

// ErrEmptySeq is returned when the sequence is empty.
var ErrEmptySeq = errors.New("digest: empty sequence")

// Enzyme represents a restriction enzyme with its recognition site and cut
// positions (offset from the 5' end of the recognition site).
type Enzyme struct {
	Name    string
	Site    string // recognition sequence (uppercase, supports IUPAC)
	CutPos  int   // cut position on the forward strand (0-based from site start)
}

// Fragment is a piece of DNA resulting from enzyme digestion.
type Fragment struct {
	Start  int    // 0-based start in the original sequence
	End    int    // 0-based end (exclusive)
	Length int
	Seq    string
}

// CommonEnzymes returns a map of well-known restriction enzymes.
func CommonEnzymes() map[string]Enzyme {
	return map[string]Enzyme{
		"EcoRI":   {Name: "EcoRI", Site: "GAATTC", CutPos: 1},
		"BamHI":   {Name: "BamHI", Site: "GGATCC", CutPos: 1},
		"HindIII": {Name: "HindIII", Site: "AAGCTT", CutPos: 1},
		"NotI":    {Name: "NotI", Site: "GCGGCCGC", CutPos: 2},
		"XhoI":    {Name: "XhoI", Site: "CTCGAG", CutPos: 1},
		"SalI":    {Name: "SalI", Site: "GTCGAC", CutPos: 1},
		"PstI":    {Name: "PstI", Site: "CTGCAG", CutPos: 5},
		"SmaI":    {Name: "SmaI", Site: "CCCGGG", CutPos: 3},
		"KpnI":    {Name: "KpnI", Site: "GGTACC", CutPos: 5},
		"SacI":    {Name: "SacI", Site: "GAGCTC", CutPos: 5},
	}
}

// CutSite represents a position where an enzyme cuts the sequence.
type CutSite struct {
	Enzyme   string
	Position int // position in the sequence where the cut occurs
	SiteStart int // start of the recognition site
}

// FindSites locates all recognition sites for the given enzyme in the sequence.
func FindSites(seq string, enz Enzyme) ([]CutSite, error) {
	if len(seq) == 0 {
		return nil, ErrEmptySeq
	}
	upper := strings.ToUpper(seq)
	site := strings.ToUpper(enz.Site)
	var sites []CutSite
	for i := 0; i+len(site) <= len(upper); i++ {
		if upper[i:i+len(site)] == site {
			sites = append(sites, CutSite{
				Enzyme:    enz.Name,
				Position:  i + enz.CutPos,
				SiteStart: i,
			})
		}
	}
	return sites, nil
}

// Digest cuts the sequence with the given enzyme and returns the resulting
// fragments in order.
func Digest(seq string, enz Enzyme) ([]Fragment, error) {
	if len(seq) == 0 {
		return nil, ErrEmptySeq
	}
	sites, err := FindSites(seq, enz)
	if err != nil {
		return nil, err
	}
	if len(sites) == 0 {
		return []Fragment{{Start: 0, End: len(seq), Length: len(seq), Seq: seq}}, nil
	}
	// Collect cut positions and sort.
	cuts := make([]int, 0, len(sites))
	for _, s := range sites {
		cuts = append(cuts, s.Position)
	}
	sort.Ints(cuts)
	// Build fragments.
	var frags []Fragment
	prev := 0
	for _, c := range cuts {
		if c > prev && c <= len(seq) {
			frags = append(frags, Fragment{
				Start:  prev,
				End:    c,
				Length: c - prev,
				Seq:    seq[prev:c],
			})
			prev = c
		}
	}
	if prev < len(seq) {
		frags = append(frags, Fragment{
			Start:  prev,
			End:    len(seq),
			Length: len(seq) - prev,
			Seq:    seq[prev:],
		})
	}
	return frags, nil
}

// DigestMulti cuts the sequence with multiple enzymes simultaneously.
func DigestMulti(seq string, enzymes []Enzyme) ([]Fragment, error) {
	if len(seq) == 0 {
		return nil, ErrEmptySeq
	}
	if len(enzymes) == 0 {
		return nil, ErrNoEnzyme
	}
	var allCuts []int
	for _, enz := range enzymes {
		sites, err := FindSites(seq, enz)
		if err != nil {
			return nil, err
		}
		for _, s := range sites {
			allCuts = append(allCuts, s.Position)
		}
	}
	sort.Ints(allCuts)
	// Deduplicate.
	if len(allCuts) > 0 {
		j := 0
		for i := 1; i < len(allCuts); i++ {
			if allCuts[i] != allCuts[j] {
				j++
				allCuts[j] = allCuts[i]
			}
		}
		allCuts = allCuts[:j+1]
	}
	var frags []Fragment
	prev := 0
	for _, c := range allCuts {
		if c > prev && c <= len(seq) {
			frags = append(frags, Fragment{
				Start:  prev,
				End:    c,
				Length: c - prev,
				Seq:    seq[prev:c],
			})
			prev = c
		}
	}
	if prev < len(seq) {
		frags = append(frags, Fragment{
			Start:  prev,
			End:    len(seq),
			Length: len(seq) - prev,
			Seq:    seq[prev:],
		})
	}
	return frags, nil
}

// FragmentLengths returns just the lengths of fragments, useful for gel
// electrophoresis simulation.
func FragmentLengths(frags []Fragment) []int {
	lengths := make([]int, len(frags))
	for i, f := range frags {
		lengths[i] = f.Length
	}
	sort.Sort(sort.Reverse(sort.IntSlice(lengths)))
	return lengths
}
