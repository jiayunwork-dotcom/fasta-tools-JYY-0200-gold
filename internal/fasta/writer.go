package fasta

import (
	"fmt"
	"io"
	"strings"
)

// DefaultLineWidth is the default number of characters per sequence line when
// writing FASTA output.
const DefaultLineWidth = 80

// Writer outputs FASTA records with configurable formatting.
type Writer struct {
	w         io.Writer
	lineWidth int
	written   int
}

// NewWriter creates a FASTA writer with the default line width.
func NewWriter(w io.Writer) *Writer {
	return &Writer{w: w, lineWidth: DefaultLineWidth}
}

// NewWriterWidth creates a FASTA writer with a custom line width.
func NewWriterWidth(w io.Writer, width int) *Writer {
	if width <= 0 {
		width = DefaultLineWidth
	}
	return &Writer{w: w, lineWidth: width}
}

// Write outputs a single record. The header line begins with '>' followed by
// the header text; the sequence is wrapped at the configured line width.
func (fw *Writer) Write(r Record) error {
	if _, err := fmt.Fprintf(fw.w, ">%s\n", r.Header); err != nil {
		return err
	}
	seq := r.Sequence
	for i := 0; i < len(seq); i += fw.lineWidth {
		end := i + fw.lineWidth
		if end > len(seq) {
			end = len(seq)
		}
		if _, err := fmt.Fprintf(fw.w, "%s\n", seq[i:end]); err != nil {
			return err
		}
	}
	if len(seq) == 0 {
		if _, err := fmt.Fprint(fw.w, "\n"); err != nil {
			return err
		}
	}
	fw.written++
	return nil
}

// WriteAll writes multiple records sequentially.
func (fw *Writer) WriteAll(records []Record) error {
	for _, r := range records {
		if err := fw.Write(r); err != nil {
			return err
		}
	}
	return nil
}

// Written returns the number of records written so far.
func (fw *Writer) Written() int { return fw.written }

// LineWidth returns the configured line width.
func (fw *Writer) LineWidth() int { return fw.lineWidth }

// Merge concatenates the sequences of all records into a single record with
// the given header.
func Merge(records []Record, header string) Record {
	var sb strings.Builder
	for _, r := range records {
		sb.WriteString(r.Sequence)
	}
	return Record{Header: header, Sequence: sb.String()}
}

// Split breaks a single record into multiple records of at most chunkSize
// nucleotides each, naming them header_1, header_2, etc.
func Split(r Record, chunkSize int) []Record {
	if chunkSize <= 0 {
		chunkSize = 1
	}
	var out []Record
	for i := 0; i < len(r.Sequence); i += chunkSize {
		end := i + chunkSize
		if end > len(r.Sequence) {
			end = len(r.Sequence)
		}
		hdr := fmt.Sprintf("%s_%d", r.Header, len(out)+1)
		out = append(out, Record{Header: hdr, Sequence: r.Sequence[i:end]})
	}
	return out
}

// Subset returns records at the given 0-based indices.
func Subset(records []Record, indices []int) []Record {
	var out []Record
	for _, idx := range indices {
		if idx >= 0 && idx < len(records) {
			out = append(out, records[idx])
		}
	}
	return out
}
