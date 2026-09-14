// Command bamroundtrip performs a no-op bamdriver decode/encode round-trip
// and compares BAMs by canonical header and record content.
package main

import (
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/rainoffallingstar/bamdriver/pkg/bamnative"
)

const usageText = `Usage:
  bamroundtrip roundtrip --input INPUT.bam --output OUTPUT.bam --report REPORT.json [--index]
  bamroundtrip compare --left LEFT.bam --right RIGHT.bam --report REPORT.json

The comparator intentionally evaluates decoded BAM structure, not compressed
bytes. It requires exact header and record ordering for no-op round-trips.
`

type comparisonReport struct {
	SchemaVersion     string `json:"schema_version"`
	Mode              string `json:"mode"`
	LeftPath          string `json:"left_path"`
	RightPath         string `json:"right_path"`
	CreatedAtUTC      string `json:"created_at_utc"`
	HeadersEqual      bool   `json:"headers_equal"`
	RecordsCompared   int64  `json:"records_compared"`
	RecordsEqual      bool   `json:"records_equal"`
	Equal             bool   `json:"equal"`
	FirstDifference   string `json:"first_difference,omitempty"`
	LeftRecordDigest  string `json:"left_record_digest"`
	RightRecordDigest string `json:"right_record_digest"`
}

type roundTripReport struct {
	SchemaVersion   string `json:"schema_version"`
	InputPath       string `json:"input_path"`
	OutputPath      string `json:"output_path"`
	CreatedAtUTC    string `json:"created_at_utc"`
	RecordsWritten  int64  `json:"records_written"`
	InputHeaderHash string `json:"input_header_hash"`
	IndexCreated    bool   `json:"index_created"`
}

func main() {
	if len(os.Args) < 2 {
		failUsage("missing command")
	}

	var err error
	switch os.Args[1] {
	case "roundtrip":
		err = runRoundTrip(os.Args[2:])
	case "compare":
		err = runComparison(os.Args[2:])
	case "--help", "-h", "help":
		fmt.Print(usageText)
		return
	default:
		failUsage(fmt.Sprintf("unknown command %q", os.Args[1]))
	}
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

func failUsage(message string) {
	fmt.Fprintln(os.Stderr, message)
	fmt.Fprint(os.Stderr, usageText)
	os.Exit(2)
}

func runRoundTrip(arguments []string) error {
	flagSet := flag.NewFlagSet("roundtrip", flag.ContinueOnError)
	flagSet.SetOutput(io.Discard)
	inputPath := flagSet.String("input", "", "input BAM path")
	outputPath := flagSet.String("output", "", "output BAM path")
	reportPath := flagSet.String("report", "", "round-trip report JSON path")
	buildIndex := flagSet.Bool("index", false, "create a BAI index for the round-trip BAM")
	if err := flagSet.Parse(arguments); err != nil {
		return fmt.Errorf("parse roundtrip flags: %w", err)
	}
	if *inputPath == "" || *outputPath == "" || *reportPath == "" {
		return errors.New("roundtrip requires --input, --output, and --report")
	}
	if filepath.Clean(*inputPath) == filepath.Clean(*outputPath) {
		return errors.New("input and output BAM paths must differ")
	}

	inputFile, err := os.Open(*inputPath)
	if err != nil {
		return fmt.Errorf("open input BAM: %w", err)
	}
	defer inputFile.Close()

	inputReader, err := bamnative.NewReader(inputFile)
	if err != nil {
		return fmt.Errorf("open BAM reader: %w", err)
	}
	outputWriter, err := bamnative.NewWriter(*outputPath, inputReader.Header())
	if err != nil {
		return fmt.Errorf("open BAM writer: %w", err)
	}

	var recordsWritten int64
	for {
		record, readErr := inputReader.Read()
		if errors.Is(readErr, io.EOF) {
			break
		}
		if readErr != nil {
			_ = outputWriter.Close()
			return fmt.Errorf("read input record %d: %w", recordsWritten, readErr)
		}
		if err := outputWriter.Write(record); err != nil {
			_ = outputWriter.Close()
			return fmt.Errorf("write output record %d: %w", recordsWritten, err)
		}
		recordsWritten++
	}
	if err := outputWriter.Close(); err != nil {
		return fmt.Errorf("close output BAM: %w", err)
	}

	indexCreated := false
	if *buildIndex {
		if err := bamnative.BuildIndex(*outputPath); err != nil {
			return fmt.Errorf("build output BAM index: %w", err)
		}
		indexCreated = true
	}

	report := roundTripReport{
		SchemaVersion:   "gate6.bamdriver-roundtrip/v1",
		InputPath:       *inputPath,
		OutputPath:      *outputPath,
		CreatedAtUTC:    time.Now().UTC().Format(time.RFC3339),
		RecordsWritten:  recordsWritten,
		InputHeaderHash: hashHeader(inputReader.Header()),
		IndexCreated:    indexCreated,
	}
	if err := writeJSONReport(*reportPath, report); err != nil {
		return fmt.Errorf("write round-trip report: %w", err)
	}
	return nil
}

func runComparison(arguments []string) error {
	flagSet := flag.NewFlagSet("compare", flag.ContinueOnError)
	flagSet.SetOutput(io.Discard)
	leftPath := flagSet.String("left", "", "left BAM path")
	rightPath := flagSet.String("right", "", "right BAM path")
	reportPath := flagSet.String("report", "", "comparison report JSON path")
	if err := flagSet.Parse(arguments); err != nil {
		return fmt.Errorf("parse compare flags: %w", err)
	}
	if *leftPath == "" || *rightPath == "" || *reportPath == "" {
		return errors.New("compare requires --left, --right, and --report")
	}

	report, err := compareBAMFiles(*leftPath, *rightPath)
	if err != nil {
		return err
	}
	report.CreatedAtUTC = time.Now().UTC().Format(time.RFC3339)
	if err := writeJSONReport(*reportPath, report); err != nil {
		return fmt.Errorf("write comparison report: %w", err)
	}
	if !report.Equal {
		return fmt.Errorf("canonical BAM comparison failed: %s", report.FirstDifference)
	}
	return nil
}

func compareBAMFiles(leftPath string, rightPath string) (comparisonReport, error) {
	report := comparisonReport{
		SchemaVersion: "gate6.bamdriver-compare/v1",
		Mode:          "ordered-canonical-record-comparison",
		LeftPath:      leftPath,
		RightPath:     rightPath,
		RecordsEqual:  true,
	}

	leftFile, err := os.Open(leftPath)
	if err != nil {
		return report, fmt.Errorf("open left BAM: %w", err)
	}
	defer leftFile.Close()
	rightFile, err := os.Open(rightPath)
	if err != nil {
		return report, fmt.Errorf("open right BAM: %w", err)
	}
	defer rightFile.Close()

	leftReader, err := bamnative.NewReader(leftFile)
	if err != nil {
		return report, fmt.Errorf("open left BAM reader: %w", err)
	}
	rightReader, err := bamnative.NewReader(rightFile)
	if err != nil {
		return report, fmt.Errorf("open right BAM reader: %w", err)
	}

	report.HeadersEqual = headersEqual(leftReader.Header(), rightReader.Header())
	if !report.HeadersEqual {
		report.RecordsEqual = false
		report.FirstDifference = "decoded BAM headers differ"
		report.Equal = false
		return report, nil
	}

	leftDigest := sha256.New()
	rightDigest := sha256.New()
	for recordOrdinal := int64(0); ; recordOrdinal++ {
		leftRecord, leftReadError := leftReader.Read()
		rightRecord, rightReadError := rightReader.Read()
		if errors.Is(leftReadError, io.EOF) || errors.Is(rightReadError, io.EOF) {
			if errors.Is(leftReadError, io.EOF) && errors.Is(rightReadError, io.EOF) {
				break
			}
			report.RecordsEqual = false
			report.FirstDifference = fmt.Sprintf("record count differs at ordinal %d", recordOrdinal)
			break
		}
		if leftReadError != nil {
			return report, fmt.Errorf("read left record %d: %w", recordOrdinal, leftReadError)
		}
		if rightReadError != nil {
			return report, fmt.Errorf("read right record %d: %w", recordOrdinal, rightReadError)
		}

		leftCanonical := canonicalRecord(leftRecord)
		rightCanonical := canonicalRecord(rightRecord)
		_, _ = leftDigest.Write([]byte(leftCanonical))
		_, _ = leftDigest.Write([]byte{'\n'})
		_, _ = rightDigest.Write([]byte(rightCanonical))
		_, _ = rightDigest.Write([]byte{'\n'})
		report.RecordsCompared++

		if leftCanonical != rightCanonical && report.FirstDifference == "" {
			report.RecordsEqual = false
			report.FirstDifference = fmt.Sprintf("record %d differs: left=%s right=%s", recordOrdinal, recordIdentity(leftRecord), recordIdentity(rightRecord))
		}
	}

	report.LeftRecordDigest = hex.EncodeToString(leftDigest.Sum(nil))
	report.RightRecordDigest = hex.EncodeToString(rightDigest.Sum(nil))
	report.Equal = report.HeadersEqual && report.RecordsEqual && report.LeftRecordDigest == report.RightRecordDigest
	if !report.Equal && report.FirstDifference == "" {
		report.FirstDifference = "canonical record stream digest differs"
	}
	return report, nil
}

func headersEqual(leftHeader *bamnative.Header, rightHeader *bamnative.Header) bool {
	return hashHeader(leftHeader) == hashHeader(rightHeader)
}

func hashHeader(header *bamnative.Header) string {
	headerJSON, err := json.Marshal(header)
	if err != nil {
		panic(fmt.Sprintf("marshal BAM header: %v", err))
	}
	digest := sha256.Sum256(headerJSON)
	return hex.EncodeToString(digest[:])
}

func canonicalRecord(record *bamnative.Record) string {
	var builder strings.Builder
	fmt.Fprintf(&builder, "%q|%d|%d|%d|%d|", record.Name, record.Flags, record.RefID, record.Pos, record.MapQ)
	for _, operation := range record.Cigar {
		fmt.Fprintf(&builder, "%d%c", operation.Len, operation.Op)
	}
	fmt.Fprintf(&builder, "|%d|%d|%d|%q|%x|", record.MateRefID, record.MatePos, record.TLen, record.Seq, record.Qual)
	for _, auxiliaryField := range record.Aux {
		fmt.Fprintf(&builder, "%s:%c:%c:%s;", auxiliaryField.Tag, auxiliaryField.Type, auxiliaryField.ArrayType, canonicalAuxiliaryValue(auxiliaryField.Value))
	}
	return builder.String()
}

func canonicalAuxiliaryValue(value interface{}) string {
	switch typedValue := value.(type) {
	case float32:
		return fmt.Sprintf("float32:%08x", math.Float32bits(typedValue))
	case float64:
		return fmt.Sprintf("float64:%016x", math.Float64bits(typedValue))
	case []byte:
		return "bytes:" + hex.EncodeToString(typedValue)
	case string:
		return fmt.Sprintf("string:%q", typedValue)
	default:
		return fmt.Sprintf("%T:%v", value, value)
	}
}

func recordIdentity(record *bamnative.Record) string {
	return fmt.Sprintf("qname=%q flag=%d ref_id=%d pos=%d cigar=%q", record.Name, record.Flags, record.RefID, record.Pos, cigarString(record.Cigar))
}

func cigarString(cigar []bamnative.CigarOp) string {
	var builder strings.Builder
	for _, operation := range cigar {
		fmt.Fprintf(&builder, "%d%c", operation.Len, operation.Op)
	}
	return builder.String()
}

func writeJSONReport(path string, report interface{}) error {
	reportBytes, err := json.MarshalIndent(report, "", "  ")
	if err != nil {
		return err
	}
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	return os.WriteFile(path, append(reportBytes, '\n'), 0o644)
}
