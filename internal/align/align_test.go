package align

import (
	"math"
	"testing"
)

func TestGlobalIdentical(t *testing.T) {
	res, err := Global("ACGT", "ACGT", DefaultScoring())
	if err != nil {
		t.Fatalf("global: %v", err)
	}
	if res.AlignedA != "ACGT" || res.AlignedB != "ACGT" {
		t.Fatalf("unexpected alignment: %q vs %q", res.AlignedA, res.AlignedB)
	}
	if res.Identity != 1.0 {
		t.Fatalf("identity = %f, want 1.0", res.Identity)
	}
}

func TestGlobalWithGaps(t *testing.T) {
	res, err := Global("AGTACG", "ACGT", DefaultScoring())
	if err != nil {
		t.Fatalf("global: %v", err)
	}
	// With mismatches and gaps the score can be negative for dissimilar seqs.
	// Just verify we get a valid alignment.
	if res.Length < 4 {
		t.Fatalf("alignment length too short: %d", res.Length)
	}
	if len(res.AlignedA) != len(res.AlignedB) {
		t.Fatal("aligned sequences differ in length")
	}
}

func TestGlobalEmptySequence(t *testing.T) {
	_, err := Global("", "ACGT", DefaultScoring())
	if err != ErrEmptySequence {
		t.Fatalf("expected ErrEmptySequence, got %v", err)
	}
}

func TestLocalBasic(t *testing.T) {
	res, err := Local("GGACGTGG", "ACGT", DefaultScoring())
	if err != nil {
		t.Fatalf("local: %v", err)
	}
	// The local alignment should find the matching ACGT region.
	if res.Score < 8 {
		t.Fatalf("local score too low: %d", res.Score)
	}
	if res.AlignedA != "ACGT" || res.AlignedB != "ACGT" {
		t.Fatalf("unexpected local alignment: %q vs %q", res.AlignedA, res.AlignedB)
	}
}

func TestLocalNoMatch(t *testing.T) {
	sc := Scoring{Match: 1, Mismatch: -3, GapOpen: -5, GapExt: -2}
	res, err := Local("AAAA", "CCCC", sc)
	if err != nil {
		t.Fatalf("local: %v", err)
	}
	if res.Score != 0 {
		t.Fatalf("score should be 0, got %d", res.Score)
	}
}

func TestHammingDistance(t *testing.T) {
	d, err := HammingDistance("ACGT", "ACGA")
	if err != nil {
		t.Fatalf("hamming: %v", err)
	}
	if d != 1 {
		t.Fatalf("hamming = %d, want 1", d)
	}
}

func TestHammingDifferentLengths(t *testing.T) {
	_, err := HammingDistance("AC", "ACGT")
	if err == nil {
		t.Fatal("expected error for different lengths")
	}
}

func TestEditDistance(t *testing.T) {
	cases := []struct {
		a, b string
		want int
	}{
		{"", "ABC", 3},
		{"ABC", "", 3},
		{"ACGT", "ACGT", 0},
		{"ACGT", "ACG", 1},
		{"kitten", "sitting", 3},
	}
	for _, c := range cases {
		got := EditDistance(c.a, c.b)
		if got != c.want {
			t.Fatalf("EditDistance(%q, %q) = %d, want %d", c.a, c.b, got, c.want)
		}
	}
}

func TestGlobalIdentity(t *testing.T) {
	res, _ := Global("ACGT", "AGGT", DefaultScoring())
	// One mismatch out of 4 => identity 0.75.
	if math.Abs(res.Identity-0.75) > 0.01 {
		t.Fatalf("identity = %f, want 0.75", res.Identity)
	}
}
