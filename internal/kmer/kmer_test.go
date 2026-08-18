package kmer

import (
	"testing"
)

func TestCountBasic(t *testing.T) {
	counts, err := Count("ACG", 2)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if len(counts) != 2 || counts["AC"] != 1 || counts["CG"] != 1 {
		t.Fatalf("bad counts: %v", counts)
	}
}

func TestCountOverlap(t *testing.T) {
	counts, err := Count("AAAA", 2)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if counts["AA"] != 3 {
		t.Fatalf("want AA=3, got %v", counts)
	}
}

func TestCountKOne(t *testing.T) {
	counts, err := Count("A", 1)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if counts["A"] != 1 {
		t.Fatalf("want A=1, got %v", counts)
	}
}

func TestCountKTooSmall(t *testing.T) {
	_, err := Count("AC", 0)
	if err == nil {
		t.Fatal("expected error for k < 1")
	}
}

func TestCountShortSeq(t *testing.T) {
	counts, err := Count("AC", 3)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if counts == nil || len(counts) != 0 {
		t.Fatalf("want non-nil empty map, got %v", counts)
	}
}

func TestCountEmptySeq(t *testing.T) {
	counts, err := Count("", 1)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if counts == nil || len(counts) != 0 {
		t.Fatalf("want non-nil empty map, got %v", counts)
	}
}
