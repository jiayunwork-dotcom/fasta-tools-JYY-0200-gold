package seq

import (
	"math"
	"testing"
)

func TestReverseComplement(t *testing.T) {
	got := ReverseComplement("GATTACA")
	want := "TGTAATC"
	if got != want {
		t.Fatalf("ReverseComplement(GATTACA)=%q, want %q", got, want)
	}
}

func TestReverseComplementCase(t *testing.T) {
	got := ReverseComplement("aCGT")
	want := "ACGt"
	if got != want {
		t.Fatalf("ReverseComplement(aCGT)=%q, want %q", got, want)
	}
}

func TestReverseComplementEmpty(t *testing.T) {
	if got := ReverseComplement(""); got != "" {
		t.Fatalf("want empty, got %q", got)
	}
}

func TestGCContent(t *testing.T) {
	cases := []struct {
		in   string
		want float64
	}{
		{"ACGT", 50.0},
		{"GG", 100.0},
		{"AT", 0.0},
		{"ACGTacgt", 50.0},
	}
	for _, c := range cases {
		if got := GCContent(c.in); got != c.want {
			t.Fatalf("GCContent(%q)=%v, want %v", c.in, got, c.want)
		}
	}
}

func TestGCContentEmpty(t *testing.T) {
	if got := GCContent(""); got != 0 {
		t.Fatalf("GCContent(\"\")=%v, want 0", got)
	}
}

func TestGCContentFraction(t *testing.T) {
	got := GCContent("acgtACGTN")
	want := 4.0 / 9.0 * 100.0
	if math.Abs(got-want) > 1e-9 {
		t.Fatalf("GCContent(acgtACGTN)=%v, want %v", got, want)
	}
}
