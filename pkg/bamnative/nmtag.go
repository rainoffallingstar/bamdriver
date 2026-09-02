package bamnative

import (
	"bytes"
	"fmt"
)

// CalculateNM calculates NM while preserving the historical no-error API.
func CalculateNM(record *Record, reference []byte, isBisulfite bool) int {
	nm, err := CalculateNMChecked(record, reference, isBisulfite)
	if err != nil {
		return 0
	}
	return nm
}

// CalculateNMChecked calculates NM and rejects malformed CIGAR or truncated input.
func CalculateNMChecked(record *Record, reference []byte, isBisulfite bool) (int, error) {
	return CalculateNMCheckedWindow(record, reference, 0, isBisulfite)
}

// CalculateNMCheckedWindow calculates NM against a reference slice that begins
// at referenceStart in the record's reference coordinate system. It preserves
// CalculateNMChecked semantics while allowing callers to fetch only the span
// consumed by a record's CIGAR.
func CalculateNMCheckedWindow(record *Record, reference []byte, referenceStart int64, isBisulfite bool) (int, error) {
	if record == nil {
		return 0, fmt.Errorf("record is nil")
	}
	if referenceStart < 0 {
		return 0, fmt.Errorf("reference window has negative start %d", referenceStart)
	}
	if record.RefID < 0 || record.Flags&FlagUnmapped != 0 {
		return 0, nil
	}
	if record.Pos < 0 {
		return 0, fmt.Errorf("record has negative reference position %d", record.Pos)
	}
	if int64(record.Pos) < referenceStart {
		return 0, fmt.Errorf("record position %d precedes reference window start %d", record.Pos, referenceStart)
	}

	readSequence := []byte(record.Seq)
	readIndex := 0
	referenceIndex := int(int64(record.Pos) - referenceStart)
	nm := 0

	ensureReadAvailable := func(length int, operation byte) error {
		if length < 0 || readIndex > len(readSequence)-length {
			return fmt.Errorf("CIGAR %d%c exceeds read length at offset %d", length, operation, readIndex)
		}
		return nil
	}
	ensureReferenceAvailable := func(length int, operation byte) error {
		if length < 0 || referenceIndex > len(reference)-length {
			return fmt.Errorf("CIGAR %d%c exceeds reference length at offset %d", length, operation, referenceIndex)
		}
		return nil
	}

	for operationIndex, operation := range record.Cigar {
		if operation.Len <= 0 {
			return 0, fmt.Errorf("CIGAR operation %d has invalid length %d", operationIndex, operation.Len)
		}
		switch operation.Op {
		case CigarMatch:
			if err := ensureReadAvailable(operation.Len, operation.Op); err != nil {
				return 0, err
			}
			if err := ensureReferenceAvailable(operation.Len, operation.Op); err != nil {
				return 0, err
			}
			for baseOffset := 0; baseOffset < operation.Len; baseOffset++ {
				referenceBase := toUpper(reference[referenceIndex+baseOffset])
				readBase := toUpper(readSequence[readIndex+baseOffset])
				if referenceBase != readBase && !isBisulfiteConversion(record, referenceBase, readBase, isBisulfite) {
					nm++
				}
			}
			readIndex += operation.Len
			referenceIndex += operation.Len
		case CigarMismatch:
			if err := ensureReadAvailable(operation.Len, operation.Op); err != nil {
				return 0, err
			}
			if err := ensureReferenceAvailable(operation.Len, operation.Op); err != nil {
				return 0, err
			}
			if !isBisulfite {
				nm += operation.Len
				readIndex += operation.Len
				referenceIndex += operation.Len
				break
			}
			nm += operation.Len
			for baseOffset := 0; baseOffset < operation.Len; baseOffset++ {
				referenceBase := toUpper(reference[referenceIndex+baseOffset])
				readBase := toUpper(readSequence[readIndex+baseOffset])
				if isBisulfiteConversion(record, referenceBase, readBase, true) {
					nm--
				}
			}
			readIndex += operation.Len
			referenceIndex += operation.Len
		case CigarEqual:
			if err := ensureReadAvailable(operation.Len, operation.Op); err != nil {
				return 0, err
			}
			if err := ensureReferenceAvailable(operation.Len, operation.Op); err != nil {
				return 0, err
			}
			readIndex += operation.Len
			referenceIndex += operation.Len
		case CigarInsertion:
			if err := ensureReadAvailable(operation.Len, operation.Op); err != nil {
				return 0, err
			}
			nm += operation.Len
			readIndex += operation.Len
		case CigarDeletion:
			if err := ensureReferenceAvailable(operation.Len, operation.Op); err != nil {
				return 0, err
			}
			nm += operation.Len
			referenceIndex += operation.Len
		case CigarSkip:
			if err := ensureReferenceAvailable(operation.Len, operation.Op); err != nil {
				return 0, err
			}
			referenceIndex += operation.Len
		case CigarSoftClip:
			if err := ensureReadAvailable(operation.Len, operation.Op); err != nil {
				return 0, err
			}
			readIndex += operation.Len
		case CigarHardClip, CigarPadding:
		default:
			return 0, fmt.Errorf("unknown CIGAR operation %q", operation.Op)
		}
	}

	if readIndex != len(readSequence) {
		return 0, fmt.Errorf("CIGAR consumes %d read bases, sequence has %d", readIndex, len(readSequence))
	}
	return nm, nil
}

func isBisulfiteConversion(record *Record, referenceBase byte, readBase byte, isBisulfite bool) bool {
	if !isBisulfite {
		return false
	}

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

func bismarkGenomeConversionContext(record *Record) string {
	if record == nil {
		return ""
	}
	auxiliaryField := record.GetAuxField("XG")
	if auxiliaryField == nil || auxiliaryField.Type != AuxTypeString {
		return ""
	}
	conversionContext, ok := auxiliaryField.Value.(string)
	if !ok || (conversionContext != "CT" && conversionContext != "GA") {
		return ""
	}
	return conversionContext
}

// toUpper converts byte to uppercase
func toUpper(b byte) byte {
	if b >= 'a' && b <= 'z' {
		return b - 32
	}
	return b
}

// HasNM checks if a record has the NM tag
func HasNM(record *Record, tagName string) bool {
	aux := record.GetAuxField(tagName)
	return aux != nil
}

// FastqToSeq converts a FASTQ-style sequence string to byte slice (for testing)
func FastqToSeq(seq string) []byte {
	return bytes.ToUpper([]byte(seq))
}
