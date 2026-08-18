package convert

import (
	"bytes"
	"strings"
	"testing"

	"fasta-tools/internal/fasta"
)

func testRecords() []fasta.Record {
	return []fasta.Record{
		{Header: "seq1", Sequence: "ACGTACGT"},
		{Header: "seq2", Sequence: "NNNN"},
	}
}

func TestWriteFASTA(t *testing.T) {
	var buf bytes.Buffer
	err := WriteFASTA(&buf, testRecords(), 4)
	if err != nil {
		t.Fatalf("write: %v", err)
	}
	out := buf.String()
	if !strings.Contains(out, ">seq1\n") {
		t.Fatal("missing header")
	}
	// With width 4, "ACGTACGT" should be split into "ACGT\nACGT\n".
	if !strings.Contains(out, "ACGT\nACGT\n") {
		t.Fatalf("wrapping: %q", out)
	}
}

func TestWriteFASTAInvalidWidth(t *testing.T) {
	var buf bytes.Buffer
	err := WriteFASTA(&buf, testRecords(), 0)
	if err != ErrInvalidWidth {
		t.Fatalf("expected ErrInvalidWidth, got %v", err)
	}
}

func TestWriteTab(t *testing.T) {
	var buf bytes.Buffer
	err := WriteTab(&buf, testRecords())
	if err != nil {
		t.Fatalf("write: %v", err)
	}
	lines := strings.Split(strings.TrimSpace(buf.String()), "\n")
	if len(lines) != 2 {
		t.Fatalf("expected 2 lines, got %d", len(lines))
	}
	if !strings.Contains(lines[0], "seq1\tACGTACGT") {
		t.Fatalf("line0: %q", lines[0])
	}
}

func TestWriteJSON(t *testing.T) {
	var buf bytes.Buffer
	err := WriteJSON(&buf, testRecords())
	if err != nil {
		t.Fatalf("write: %v", err)
	}
	if !strings.Contains(buf.String(), `"header":"seq1"`) {
		t.Fatalf("missing JSON: %q", buf.String())
	}
}

func TestWriteSingle(t *testing.T) {
	var buf bytes.Buffer
	err := WriteSingle(&buf, testRecords())
	if err != nil {
		t.Fatalf("write: %v", err)
	}
	if !strings.Contains(buf.String(), ">seq1\nACGTACGT\n") {
		t.Fatalf("wrong format: %q", buf.String())
	}
}

func TestToUpperCase(t *testing.T) {
	recs := []fasta.Record{{Header: "a", Sequence: "acgt"}}
	upper := ToUpperCase(recs)
	if upper[0].Sequence != "ACGT" {
		t.Fatalf("upper: %q", upper[0].Sequence)
	}
}

func TestFilterByHeader(t *testing.T) {
	recs := testRecords()
	filtered := FilterByHeader(recs, "seq1")
	if len(filtered) != 1 || filtered[0].Header != "seq1" {
		t.Fatalf("filter: %v", filtered)
	}
}

func TestRemoveDuplicates(t *testing.T) {
	recs := []fasta.Record{
		{Header: "a", Sequence: "ACGT"},
		{Header: "b", Sequence: "ACGT"},
		{Header: "c", Sequence: "TTTT"},
	}
	deduped := RemoveDuplicates(recs)
	if len(deduped) != 2 {
		t.Fatalf("dedup: %d", len(deduped))
	}
}

func TestSplitByCount(t *testing.T) {
	recs := make([]fasta.Record, 7)
	chunks := SplitByCount(recs, 3)
	if len(chunks) != 3 {
		t.Fatalf("chunks = %d, want 3", len(chunks))
	}
	if len(chunks[2]) != 1 {
		t.Fatalf("last chunk = %d, want 1", len(chunks[2]))
	}
}
