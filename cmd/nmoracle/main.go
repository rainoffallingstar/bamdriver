// Command nmoracle independently computes conventional and strand-aware
// bisulfite edit-distance values from BAM records and a reference FASTA.
package main

import (
	"bufio"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/rainoffallingstar/bamdriver/pkg/bamnative"
)

type scoreComponents struct {
	Mismatches                  int64 `json:"mismatches"`
	BisulfiteConversionsIgnored int64 `json:"bisulfite_conversions_ignored"`
	InsertionBases              int64 `json:"insertion_bases"`
	DeletionBases               int64 `json:"deletion_bases"`
	SoftClipBases               int64 `json:"soft_clip_bases"`
}

type oracleSummary struct {
	SchemaVersion               string          `json:"schema_version"`
	InputPath                   string          `json:"input_path"`
	ReferencePath               string          `json:"reference_path"`
	CreatedAtUTC                string          `json:"created_at_utc"`
	RecordsRead                 int64           `json:"records_read"`
	MappedRecordsProcessed      int64           `json:"mapped_records_processed"`
	UnmappedRecordsSkipped      int64           `json:"unmapped_records_skipped"`
	ReferenceLookupFailures     int64           `json:"reference_lookup_failures"`
	MalformedAlignmentFailures  int64           `json:"malformed_alignment_failures"`
	ConventionalNMTotal         int64           `json:"conventional_nm_total"`
	XenofilxBisulfiteNMTotal    int64           `json:"xenofilx_bisulfite_nm_total"`
	XenofilxClassificationTotal int64           `json:"xenofilx_classification_score_total"`
	Components                  scoreComponents `json:"components"`
}

type alignmentScore struct {
	ConventionalNM      int64
	XenofilxBisulfiteNM int64
	XenofilxScore       int64
	Components          scoreComponents
}

func main() {
	flagSet := flag.NewFlagSet("nmoracle", flag.ExitOnError)
	inputPath := flagSet.String("input", "", "input BAM path")
	referencePath := flagSet.String("reference", "", "reference FASTA path")
	reportPath := flagSet.String("report", "", "summary JSON path")
	recordsPath := flagSet.String("records", "", "optional per-record TSV path")
	limit := flagSet.Int64("limit", 0, "maximum records to read; zero scans all records")
	flagSet.Parse(os.Args[1:])
	if *inputPath == "" || *referencePath == "" || *reportPath == "" {
		fmt.Fprintln(os.Stderr, "nmoracle requires --input, --reference, and --report")
		os.Exit(2)
	}
	if *limit < 0 {
		fmt.Fprintln(os.Stderr, "--limit must be zero or positive")
		os.Exit(2)
	}

	summary, err := scanBAM(*inputPath, *referencePath, *recordsPath, *limit)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	if err := writeSummary(*reportPath, summary); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

func scanBAM(inputPath string, referencePath string, recordsPath string, limit int64) (oracleSummary, error) {
	summary := oracleSummary{
		SchemaVersion: "gate6.nm-oracle/v1",
		InputPath:     inputPath,
		ReferencePath: referencePath,
		CreatedAtUTC:  time.Now().UTC().Format(time.RFC3339),
	}
	inputFile, err := os.Open(inputPath)
	if err != nil {
		return summary, fmt.Errorf("open BAM: %w", err)
	}
	defer inputFile.Close()
	reader, err := bamnative.NewReader(inputFile)
	if err != nil {
		return summary, fmt.Errorf("create BAM reader: %w", err)
	}
	referenceReader, err := bamnative.NewFastaReader(referencePath)
	if err != nil {
		return summary, fmt.Errorf("open reference FASTA: %w", err)
	}
	defer referenceReader.Close()

	var recordsFile *os.File
	var recordsWriter *bufio.Writer
	if recordsPath != "" {
		if err := os.MkdirAll(filepath.Dir(recordsPath), 0o755); err != nil {
			return summary, fmt.Errorf("create records report directory: %w", err)
		}
		recordsFile, err = os.Create(recordsPath)
		if err != nil {
			return summary, fmt.Errorf("create records report: %w", err)
		}
		defer recordsFile.Close()
		recordsWriter = bufio.NewWriterSize(recordsFile, 1<<20)
		defer recordsWriter.Flush()
		if _, err := fmt.Fprintln(recordsWriter, "record_ordinal\tqname\tflag\treference\tposition_0_based\tcigar\tconventional_nm\txenofilx_bisulfite_nm\tinsertions\tdeletions\tsoft_clips\tbisulfite_conversions_ignored\tclassification_score"); err != nil {
			return summary, fmt.Errorf("write records report header: %w", err)
		}
	}

	referenceNames := make(map[int32]string, len(reader.Header().References))
	for _, reference := range reader.Header().References {
		referenceNames[reference.ID] = reference.Name
	}

	for limit == 0 || summary.RecordsRead < limit {
		record, readErr := reader.Read()
		if errors.Is(readErr, io.EOF) {
			break
		}
		if readErr != nil {
			return summary, fmt.Errorf("read record %d: %w", summary.RecordsRead, readErr)
		}
		summary.RecordsRead++
		if record.IsUnmapped() || record.RefID < 0 {
			summary.UnmappedRecordsSkipped++
			continue
		}

		referenceName, found := referenceNames[record.RefID]
		if !found {
			summary.ReferenceLookupFailures++
			continue
		}
		referenceSpan, spanErr := calculateReferenceSpan(record.Cigar)
		if spanErr != nil {
			summary.MalformedAlignmentFailures++
			continue
		}
		referenceSequence, found := referenceReader.GetRegion(referenceName, int64(record.Pos), int64(record.Pos)+int64(referenceSpan))
		if !found {
			summary.ReferenceLookupFailures++
			continue
		}
		score, scoreErr := scoreAlignment(record, referenceSequence)
		if scoreErr != nil {
			summary.MalformedAlignmentFailures++
			continue
		}
		summary.MappedRecordsProcessed++
		summary.ConventionalNMTotal += score.ConventionalNM
		summary.XenofilxBisulfiteNMTotal += score.XenofilxBisulfiteNM
		summary.XenofilxClassificationTotal += score.XenofilxScore
		addComponents(&summary.Components, score.Components)
		if recordsWriter != nil {
			if _, err := fmt.Fprintf(
				recordsWriter,
				"%d\t%s\t%d\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
				summary.RecordsRead-1,
				escapeTSV(record.Name),
				record.Flags,
				escapeTSV(referenceName),
				record.Pos,
				cigarString(record.Cigar),
				score.ConventionalNM,
				score.XenofilxBisulfiteNM,
				score.Components.InsertionBases,
				score.Components.DeletionBases,
				score.Components.SoftClipBases,
				score.Components.BisulfiteConversionsIgnored,
				score.XenofilxScore,
			); err != nil {
				return summary, fmt.Errorf("write records report for record %d: %w", summary.RecordsRead-1, err)
			}
		}
	}
	return summary, nil
}

func scoreAlignment(record *bamnative.Record, referenceSequence []byte) (alignmentScore, error) {
	if record == nil {
		return alignmentScore{}, errors.New("record is nil")
	}
	readSequence := []byte(strings.ToUpper(record.Seq))
	readPosition := 0
	referencePosition := 0
	components := scoreComponents{}

	for operationIndex, operation := range record.Cigar {
		if operation.Len <= 0 {
			return alignmentScore{}, fmt.Errorf("CIGAR operation %d has invalid length %d", operationIndex, operation.Len)
		}
		switch operation.Op {
		case bamnative.CigarMatch:
			if err := validateSpan(readPosition, operation.Len, len(readSequence), "read", operation); err != nil {
				return alignmentScore{}, err
			}
			if err := validateSpan(referencePosition, operation.Len, len(referenceSequence), "reference", operation); err != nil {
				return alignmentScore{}, err
			}
			for baseOffset := 0; baseOffset < operation.Len; baseOffset++ {
				referenceBase := toUpper(referenceSequence[referencePosition+baseOffset])
				readBase := readSequence[readPosition+baseOffset]
				if referenceBase == readBase {
					continue
				}
				components.Mismatches++
				if isStrandAwareBisulfiteConversion(record, referenceBase, readBase) {
					components.BisulfiteConversionsIgnored++
				}
			}
			readPosition += operation.Len
			referencePosition += operation.Len
		case bamnative.CigarMismatch:
			if err := validateSpan(readPosition, operation.Len, len(readSequence), "read", operation); err != nil {
				return alignmentScore{}, err
			}
			if err := validateSpan(referencePosition, operation.Len, len(referenceSequence), "reference", operation); err != nil {
				return alignmentScore{}, err
			}
			components.Mismatches += int64(operation.Len)
			for baseOffset := 0; baseOffset < operation.Len; baseOffset++ {
				referenceBase := toUpper(referenceSequence[referencePosition+baseOffset])
				readBase := readSequence[readPosition+baseOffset]
				if isStrandAwareBisulfiteConversion(record, referenceBase, readBase) {
					components.BisulfiteConversionsIgnored++
				}
			}
			readPosition += operation.Len
			referencePosition += operation.Len
		case bamnative.CigarEqual:
			if err := validateSpan(readPosition, operation.Len, len(readSequence), "read", operation); err != nil {
				return alignmentScore{}, err
			}
			if err := validateSpan(referencePosition, operation.Len, len(referenceSequence), "reference", operation); err != nil {
				return alignmentScore{}, err
			}
			readPosition += operation.Len
			referencePosition += operation.Len
		case bamnative.CigarInsertion:
			if err := validateSpan(readPosition, operation.Len, len(readSequence), "read", operation); err != nil {
				return alignmentScore{}, err
			}
			components.InsertionBases += int64(operation.Len)
			readPosition += operation.Len
		case bamnative.CigarDeletion:
			if err := validateSpan(referencePosition, operation.Len, len(referenceSequence), "reference", operation); err != nil {
				return alignmentScore{}, err
			}
			components.DeletionBases += int64(operation.Len)
			referencePosition += operation.Len
		case bamnative.CigarSkip:
			if err := validateSpan(referencePosition, operation.Len, len(referenceSequence), "reference", operation); err != nil {
				return alignmentScore{}, err
			}
			referencePosition += operation.Len
		case bamnative.CigarSoftClip:
			if err := validateSpan(readPosition, operation.Len, len(readSequence), "read", operation); err != nil {
				return alignmentScore{}, err
			}
			components.SoftClipBases += int64(operation.Len)
			readPosition += operation.Len
		case bamnative.CigarHardClip, bamnative.CigarPadding:
		default:
			return alignmentScore{}, fmt.Errorf("unsupported CIGAR operation %q", operation.Op)
		}
	}
	if readPosition != len(readSequence) {
		return alignmentScore{}, fmt.Errorf("CIGAR consumes %d read bases; sequence has %d", readPosition, len(readSequence))
	}

	conventionalNM := components.Mismatches + components.InsertionBases + components.DeletionBases
	xenofilxBisulfiteNM := conventionalNM - components.BisulfiteConversionsIgnored
	xenofilxScore := xenofilxBisulfiteNM + components.InsertionBases + components.SoftClipBases
	return alignmentScore{
		ConventionalNM:      conventionalNM,
		XenofilxBisulfiteNM: xenofilxBisulfiteNM,
		XenofilxScore:       xenofilxScore,
		Components:          components,
	}, nil
}

func calculateReferenceSpan(cigar []bamnative.CigarOp) (int, error) {
	span := 0
	for _, operation := range cigar {
		if operation.Len <= 0 {
			return 0, fmt.Errorf("CIGAR operation %q has invalid length %d", operation.Op, operation.Len)
		}
		switch operation.Op {
		case bamnative.CigarMatch, bamnative.CigarEqual, bamnative.CigarMismatch, bamnative.CigarDeletion, bamnative.CigarSkip:
			span += operation.Len
		case bamnative.CigarInsertion, bamnative.CigarSoftClip, bamnative.CigarHardClip, bamnative.CigarPadding:
		default:
			return 0, fmt.Errorf("unsupported CIGAR operation %q", operation.Op)
		}
	}
	if span == 0 {
		return 0, errors.New("CIGAR consumes no reference bases")
	}
	return span, nil
}

func validateSpan(start int, length int, available int, source string, operation bamnative.CigarOp) error {
	if start < 0 || length < 0 || start > available-length {
		return fmt.Errorf("CIGAR %d%c exceeds %s sequence", operation.Len, operation.Op, source)
	}
	return nil
}

func isStrandAwareBisulfiteConversion(record *bamnative.Record, referenceBase byte, readBase byte) bool {
	conversionContext := bismarkGenomeConversionContext(record)
	if conversionContext == "GA" {
		return referenceBase == 'G' && readBase == 'A'
	}
	if conversionContext == "CT" {
		return referenceBase == 'C' && readBase == 'T'
	}

	if record.IsReverse() {
		return referenceBase == 'G' && readBase == 'A'
	}
	return referenceBase == 'C' && readBase == 'T'
}

func bismarkGenomeConversionContext(record *bamnative.Record) string {
	if record == nil {
		return ""
	}
	auxiliaryField := record.GetAuxField("XG")
	if auxiliaryField == nil || auxiliaryField.Type != bamnative.AuxTypeString {
		return ""
	}
	conversionContext, ok := auxiliaryField.Value.(string)
	if !ok || (conversionContext != "CT" && conversionContext != "GA") {
		return ""
	}
	return conversionContext
}

func toUpper(base byte) byte {
	if base >= 'a' && base <= 'z' {
		return base - ('a' - 'A')
	}
	return base
}

func addComponents(target *scoreComponents, value scoreComponents) {
	target.Mismatches += value.Mismatches
	target.BisulfiteConversionsIgnored += value.BisulfiteConversionsIgnored
	target.InsertionBases += value.InsertionBases
	target.DeletionBases += value.DeletionBases
	target.SoftClipBases += value.SoftClipBases
}

func escapeTSV(value string) string {
	return strings.NewReplacer("\t", "\\t", "\n", "\\n", "\r", "\\r").Replace(value)
}

func cigarString(cigar []bamnative.CigarOp) string {
	var builder strings.Builder
	for _, operation := range cigar {
		fmt.Fprintf(&builder, "%d%c", operation.Len, operation.Op)
	}
	return builder.String()
}

func writeSummary(reportPath string, summary oracleSummary) error {
	reportBytes, err := json.MarshalIndent(summary, "", "  ")
	if err != nil {
		return err
	}
	if err := os.MkdirAll(filepath.Dir(reportPath), 0o755); err != nil {
		return err
	}
	if err := os.WriteFile(reportPath, append(reportBytes, '\n'), 0o644); err != nil {
		return fmt.Errorf("write report: %w", err)
	}
	return nil
}
