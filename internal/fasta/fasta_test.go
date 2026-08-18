package fasta

import (
	"strings"
	"testing"
)

func TestParseValid(t *testing.T) {
	in := ">seq1\nACGT\n>seq2\nacgtN\n"
	recs, err := Parse(strings.NewReader(in))
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if len(recs) != 2 {
		t.Fatalf("want 2 records, got %d", len(recs))
	}
	if recs[0].Header != "seq1" || recs[0].Sequence != "ACGT" {
		t.Fatalf("bad record 0: %+v", recs[0])
	}
	if recs[1].Header != "seq2" || recs[1].Sequence != "acgtN" {
		t.Fatalf("bad record 1: %+v", recs[1])
	}
}

func TestParseMixedLines(t *testing.T) {
	in := ">r\nACGT\nACGU\nN\n"
	recs, err := Parse(strings.NewReader(in))
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if len(recs) != 1 || recs[0].Sequence != "ACGTACGUN" {
		t.Fatalf("bad merge: %+v", recs)
	}
}

func TestParseBadChar(t *testing.T) {
	in := ">r\nACGTZ\n"
	_, err := Parse(strings.NewReader(in))
	if err == nil {
		t.Fatal("expected error for invalid character")
	}
}

func TestParseMissingHeader(t *testing.T) {
	in := "ACGT\n"
	_, err := Parse(strings.NewReader(in))
	if err == nil {
		t.Fatal("expected error for sequence before header")
	}
}

func TestParseEmptyFile(t *testing.T) {
	_, err := Parse(strings.NewReader(""))
	if err == nil {
		t.Fatal("expected error for empty input")
	}
}

func TestParseEmptySequence(t *testing.T) {
	in := ">\n"
	recs, err := Parse(strings.NewReader(in))
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if len(recs) != 1 || recs[0].Sequence != "" {
		t.Fatalf("bad empty-seq record: %+v", recs)
	}
}
