package orf

import "testing"

func TestFindBasic(t *testing.T) {
	// ATG + 30 codons + TAA = 93nt (< 100 default).
	// Make a longer one: ATG + 34 codons + TAA = 105nt.
	seq := "ATG" + "ACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACG" + "TAA"
	opts := DefaultOptions()
	opts.MinLength = 50
	orfs, err := Find(seq, opts)
	if err != nil {
		t.Fatalf("find: %v", err)
	}
	if len(orfs) == 0 {
		t.Fatal("expected at least one ORF")
	}
	if orfs[0].Start != 0 {
		t.Fatalf("ORF start = %d, want 0", orfs[0].Start)
	}
}

func TestFindNoORF(t *testing.T) {
	// No ATG at all.
	seq := "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
	opts := DefaultOptions()
	opts.MinLength = 10
	orfs, err := Find(seq, opts)
	if err != nil {
		t.Fatalf("find: %v", err)
	}
	if len(orfs) != 0 {
		t.Fatalf("expected 0 ORFs, got %d", len(orfs))
	}
}

func TestFindBothStrands(t *testing.T) {
	// Forward ORF: ATG...TAA
	seq := "NNATGAAAAAAAAATAAN"
	opts := Options{MinLength: 10, BothStrands: true, AllowNested: true}
	orfs, err := Find(seq, opts)
	if err != nil {
		t.Fatalf("find: %v", err)
	}
	if len(orfs) == 0 {
		t.Fatal("expected at least one ORF with both strands")
	}
}

func TestFindEmpty(t *testing.T) {
	_, err := Find("", DefaultOptions())
	if err != ErrEmptySequence {
		t.Fatalf("expected ErrEmptySequence, got %v", err)
	}
}

func TestLongest(t *testing.T) {
	orfs := []ORF{
		{Length: 50},
		{Length: 200},
		{Length: 100},
	}
	best := Longest(orfs)
	if best.Length != 200 {
		t.Fatalf("longest = %d, want 200", best.Length)
	}
}

func TestLongestEmpty(t *testing.T) {
	if Longest(nil) != nil {
		t.Fatal("longest(nil) should be nil")
	}
}

func TestFilterByLength(t *testing.T) {
	orfs := []ORF{
		{Length: 50},
		{Length: 200},
		{Length: 100},
	}
	filtered := FilterByLength(orfs, 100)
	if len(filtered) != 2 {
		t.Fatalf("expected 2, got %d", len(filtered))
	}
}

func TestCountByFrame(t *testing.T) {
	orfs := []ORF{
		{Frame: 0}, {Frame: 0}, {Frame: 1}, {Frame: 2},
	}
	counts := CountByFrame(orfs)
	if counts[0] != 2 || counts[1] != 1 || counts[2] != 1 {
		t.Fatalf("counts = %v", counts)
	}
}
