package stats

import (
	"math"
	"testing"

	"fasta-tools/internal/fasta"
)

func makeRecords(seqs ...string) []fasta.Record {
	recs := make([]fasta.Record, len(seqs))
	for i, s := range seqs {
		recs[i] = fasta.Record{Header: "seq" + string(rune('0'+i)), Sequence: s}
	}
	return recs
}

func TestComputeBasic(t *testing.T) {
	recs := makeRecords("ACGT", "AAGGCC", "TTT")
	s, err := Compute(recs)
	if err != nil {
		t.Fatalf("compute: %v", err)
	}
	if s.NumRecords != 3 {
		t.Fatalf("num = %d", s.NumRecords)
	}
	if s.TotalBases != 13 {
		t.Fatalf("total = %d", s.TotalBases)
	}
	if s.MinLength != 3 {
		t.Fatalf("min = %d", s.MinLength)
	}
	if s.MaxLength != 6 {
		t.Fatalf("max = %d", s.MaxLength)
	}
}

func TestComputeEmpty(t *testing.T) {
	_, err := Compute(nil)
	if err != ErrNoRecords {
		t.Fatalf("expected ErrNoRecords, got %v", err)
	}
}

func TestN50(t *testing.T) {
	// Sequences of lengths: 2, 3, 4, 5, 6 = total 20, half = 10.
	// From longest: 6 (cum=6), 5 (cum=11) >= 10 → N50=5.
	recs := makeRecords("AA", "CCC", "GGGG", "TTTTT", "AAAAAA")
	s, _ := Compute(recs)
	if s.N50 != 5 {
		t.Fatalf("N50 = %d, want 5", s.N50)
	}
}

func TestGCPercent(t *testing.T) {
	recs := makeRecords("GCGCGC")
	s, _ := Compute(recs)
	if math.Abs(s.GCPercent-100) > 0.01 {
		t.Fatalf("GC = %f, want 100", s.GCPercent)
	}
}

func TestLengthHistogram(t *testing.T) {
	recs := makeRecords("A", "AA", "AAA", "AAAA", "AAAAA")
	bins := LengthHistogram(recs, 2)
	if len(bins) == 0 {
		t.Fatal("empty histogram")
	}
	total := 0
	for _, b := range bins {
		total += b.Count
	}
	if total != 5 {
		t.Fatalf("total count = %d, want 5", total)
	}
}

func TestPercentile(t *testing.T) {
	recs := makeRecords("A", "AA", "AAA", "AAAA", "AAAAA")
	p50 := Percentile(recs, 50)
	if p50 != 3.0 {
		t.Fatalf("p50 = %f, want 3.0", p50)
	}
}

func TestFilterByLength(t *testing.T) {
	recs := makeRecords("A", "AA", "AAA", "AAAA")
	filtered := FilterByLength(recs, 2, 3)
	if len(filtered) != 2 {
		t.Fatalf("filtered = %d, want 2", len(filtered))
	}
}
