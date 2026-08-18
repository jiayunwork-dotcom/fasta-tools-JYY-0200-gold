package main

import (
	"flag"
	"fmt"
	"os"
	"strconv"
	"strings"

	"fasta-tools/internal/fasta"
	"fasta-tools/internal/kmer"
	"fasta-tools/internal/seq"
)

func main() {
	if len(os.Args) < 2 {
		usage()
		os.Exit(2)
	}

	cmd, rest := os.Args[1], os.Args[2:]
	var err error
	switch cmd {
	case "gc":
		err = runGC(rest)
	case "rc":
		err = runRC(rest)
	case "kmer":
		err = runKmer(rest)
	case "-h", "--help", "help":
		usage()
		os.Exit(0)
	default:
		usage()
		os.Exit(2)
	}

	if err != nil {
		fmt.Fprintln(os.Stderr, "error:", err)
		os.Exit(1)
	}
}

// reorder moves flags in front of positional arguments so that
// flag.Parse can handle flags placed after positionals.
func reorder(args []string) []string {
	var flags, pos []string
	i := 0
	for i < len(args) {
		a := args[i]
		if strings.HasPrefix(a, "-") && a != "-" && a != "--" {
			flags = append(flags, a)
			if !strings.Contains(a, "=") && i+1 < len(args) && !strings.HasPrefix(args[i+1], "-") {
				flags = append(flags, args[i+1])
				i++
			}
		} else {
			pos = append(pos, a)
		}
		i++
	}
	return append(flags, pos...)
}

func readFASTA(path string) ([]fasta.Record, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	return fasta.Parse(f)
}

func runGC(args []string) error {
	fs := flag.NewFlagSet("gc", flag.ContinueOnError)
	if err := fs.Parse(reorder(args)); err != nil {
		return err
	}
	paths := fs.Args()
	if len(paths) != 1 {
		return fmt.Errorf("gc requires exactly one FASTA file")
	}
	records, err := readFASTA(paths[0])
	if err != nil {
		return err
	}
	for _, r := range records {
		fmt.Printf("%s\t%.2f\n", r.Header, seq.GCContent(r.Sequence))
	}
	return nil
}

func runRC(args []string) error {
	fs := flag.NewFlagSet("rc", flag.ContinueOnError)
	if err := fs.Parse(reorder(args)); err != nil {
		return err
	}
	paths := fs.Args()
	if len(paths) != 1 {
		return fmt.Errorf("rc requires exactly one FASTA file")
	}
	records, err := readFASTA(paths[0])
	if err != nil {
		return err
	}
	for _, r := range records {
		fmt.Printf(">%s\n%s\n", r.Header, seq.ReverseComplement(r.Sequence))
	}
	return nil
}

func runKmer(args []string) error {
	fs := flag.NewFlagSet("kmer", flag.ContinueOnError)
	k := fs.Int("k", 0, "k-mer size")
	if err := fs.Parse(reorder(args)); err != nil {
		return err
	}
	paths := fs.Args()
	if len(paths) < 1 {
		return fmt.Errorf("kmer requires a FASTA file")
	}
	kk := *k
	if kk == 0 && len(paths) >= 2 {
		n, err := strconv.Atoi(paths[1])
		if err != nil {
			return fmt.Errorf("invalid k %q: %v", paths[1], err)
		}
		kk = n
	}
	if kk < 1 {
		return fmt.Errorf("k must be >= 1")
	}
	records, err := readFASTA(paths[0])
	if err != nil {
		return err
	}
	for _, r := range records {
		counts, err := kmer.Count(r.Sequence, kk)
		if err != nil {
			return err
		}
		for km, c := range counts {
			fmt.Printf("%s\t%s\t%d\n", r.Header, km, c)
		}
	}
	return nil
}

func usage() {
	fmt.Fprintln(os.Stderr, "usage: fasta-tools <command> [flags] <file>")
	fmt.Fprintln(os.Stderr, "commands:")
	fmt.Fprintln(os.Stderr, "  gc    <file>           GC content per record")
	fmt.Fprintln(os.Stderr, "  rc    <file>           reverse complement per record")
	fmt.Fprintln(os.Stderr, "  kmer  <file> [k|-k k]  k-mer counts per record")
}
