from pathlib import Path
import dnaio
import csv

## VIASH START
par = {

}
## VIASH END

def resolve_header_name_to_index(header_entries, column_name):
    try:
        return header_entries.index(column_name)
    except ValueError as e:
        raise ValueError(f"Column name '{column_name}' could not "
                         "be found in the header of the CSV file.") from e


def csv_records(csv_file, delimiter, quote_character, 
                header, sequence_column, name_column,
                sequence_column_index, name_column_index):
    with open(csv_file, newline='') as csvfile:
        csv_reader = csv.reader(csvfile,
                                delimiter=delimiter,
                                quotechar=quote_character)
        for linenum, line in enumerate(csv_reader):
            if not linenum: # First row
                num_columns = len(line)
                if header:
                    if sequence_column:
                        sequence_column_index = resolve_header_name_to_index(line, sequence_column)
                    if name_column:
                        name_column_index = resolve_header_name_to_index(line, name_column)
                    continue
            if not (linenum - header): # First 'data' line
                if (not sequence_column_index and not name_column_index and len(line) == 2):
                    name_column_index, sequence_column_index = 0, 1
                if sequence_column_index == name_column_index:
                    raise ValueError("The same columns were selected for both the FASTQ sequences and "
                                     "headers.")
                if sequence_column_index is None:
                    raise ValueError("Either 'sequence_column_index' or 'sequence_column' needs "
                                     "to be specified.")
                if name_column_index is None:
                    raise ValueError("Either 'name_column' or 'name_column_index' needs to "
                                     "be specified.")
                if name_column_index >= num_columns:
                    raise ValueError(f"Requested to use column number {name_column_index} "
                                     f"(0 based) for the FASTA headers, but only {num_columns} "
                                     "were found on the first line.")
                if sequence_column_index >= num_columns:
                    raise ValueError(f"Requested to use column number {sequence_column_index} "
                                     f"(0 based) for the FASTA sequences, but only {num_columns} "
                                     "were found on the first line.") 
            if len(line) != num_columns:
                raise ValueError(f"Number of columns ({len(line)}) found on line {linenum+1} "
                                 "is different compared to number of columns found "
                                 f"previously ({num_columns}).")
            yield line[name_column_index], line[sequence_column_index]
                

def main(par):
    par['input'], par['output'] = Path(par['input']), Path(par['output'])
    sequence_column, name_column = par['sequence_column'], par['name_column'] 
    sequence_column_index, name_column_index = par['sequence_column_index'], par['name_column_index']
    if (sequence_column or name_column) and not par['header']:
        par["header"] = True
    if sequence_column_index and sequence_column:
        raise ValueError("Cannot specify both 'sequence_column_index' and 'sequence_column'")
    if name_column and name_column_index:
        raise ValueError("Cannot specify both 'name_column_index' and 'name_column'")
    if (sequence_column_index or name_column_index) and \
        (sequence_column_index == name_column_index):
        raise ValueError("The value specified for 'sequence_column_index' cannot be the same as "
                         "the value for 'name_column_index'.")
    with dnaio.open(par['output'], mode='w', fileformat="fasta") as writer:
        for header, sequence in csv_records(par['input'],
                                            par['delimiter'],
                                            par['quote_character'],
                                            par['header'],
                                            sequence_column,
                                            name_column,
                                            sequence_column_index,
                                            name_column_index):
            writer.write(dnaio.SequenceRecord(header, sequence))
    
if __name__ == "__main__":
    main(par)