import pytest
import re
import sys
from uuid import uuid4
from textwrap import dedent
from subprocess import CalledProcessError

## VIASH START
meta = {
    'config': 'src/sequenceformats/csv2fasta/config.vsh.yaml',
    'executable': 'target/executable/sequenceformats/csv2fasta'
}
## VIASH END

@pytest.fixture
def random_path(tmp_path):
    def wrapper(extension=None):
        extension = "" if not extension else f".{extension}"
        return tmp_path / f"{uuid4()}{extension}"
    return wrapper 

@pytest.mark.parametrize("arg,val,expected_err", [("name_column", "barcode_name",
                                                   ("sequence_column_index", "sequence_column")),
                                                  ("sequence_column", "sequence",
                                                   ("name_column", "name_column_index"))])
def test_csvtofasta_no_columns_selected_raises(run_component, random_path, arg, val, expected_err):
    csv_contents = dedent("""\
    barcode_name,some_other_column,sequence
    barcode1,foo,ACGT
    barcode2,bar,TTTA
    """)

    input_path = random_path("csv")
    with input_path.open('w') as open_input:
        open_input.write(csv_contents)
    output_path = random_path("csv")
    args = [
        "--input", input_path,
        "--output", output_path,
        "--header"
    ]
    args.extend([f"--{arg}", val])
    with pytest.raises(CalledProcessError) as err:
        run_component(args)
    assert f"ValueError: Either '{expected_err[0]}' or '{expected_err[1]}' needs to be specified." in \
        err.value.stdout.decode('utf-8')

def test_csvtofasta_column_does_not_exist_raises(run_component, random_path):
    csv_contents = dedent("""\
    barcode_name,some_other_column,sequence
    barcode1,foo,ACGT
    barcode2,bar,TTTA
    """)

    input_path = random_path("csv")
    with input_path.open('w') as open_input:
        open_input.write(csv_contents)
    output_path = random_path("csv")
    args = [
        "--input", input_path,
        "--output", output_path,
        "--sequence_column", "foo",
    ]
    with pytest.raises(CalledProcessError) as err:
        run_component(args)
    assert "ValueError: Column name 'foo' could not be found in the " + \
           "header of the CSV file." in err.value.stdout.decode('utf-8')

def test_csvtofasta_same_column_selected_raises(run_component, random_path):
    csv_contents = dedent("""\
    barcode_name,some_other_column,sequence
    barcode1,foo,ACGT
    barcode2,bar,TTTA
    """)

    input_path = random_path("csv")
    with input_path.open('w') as open_input:
        open_input.write(csv_contents)
    output_path = random_path("csv")
    args = [
        "--input", input_path,
        "--output", output_path,
        "--sequence_column_index", "1",
        "--name_column_index", "1",
    ]
    with pytest.raises(CalledProcessError) as err:
        run_component(args)
    assert "ValueError: The value specified for 'sequence_column_index' cannot " + \
           "be the same as the value for 'name_column_index'" in \
            err.value.stdout.decode('utf-8')

@pytest.mark.parametrize("arg,val,expected_err", [("sequence_column_index", "3", "sequences"),
                                                  ("name_column_index", "4", "headers")])
def test_csvtofasta_header_select_index_out_of_bounds_raises(run_component, random_path, arg, val, expected_err):
    csv_contents = dedent("""\
    barcode_name,some_other_column,sequence
    barcode1,foo,ACGT
    barcode2,bar,TTTA
    """)

    other_column_map = {
        "sequence_column_index": ["--name_column_index", "1"],
        "name_column_index": ["--sequence_column_index", "2"],
    }
    input_path = random_path("csv")
    with input_path.open('w') as open_input:
        open_input.write(csv_contents)
    output_path = random_path("csv")
    args = [
        "--input", input_path,
        "--output", output_path,
        "--header",
    ]
    args += [f"--{arg}", val]
    args += other_column_map[arg] 
    with pytest.raises(CalledProcessError) as err:
        run_component(args)
    assert f"ValueError: Requested to use column number {val} (0 based) for the FASTA " + \
           f"{expected_err}, but only 3 were found on the first line." in \
        err.value.stdout.decode('utf-8')

def test_csvtofasta_header_select_column_by_both_name_and_index(run_component, random_path):
    csv_contents = dedent("""\
    barcode_name,some_other_column,sequence
    barcode1,foo,ACGT
    barcode2,bar,TTTA
    """)

    expected= dedent("""\
    >barcode1
    ACGT
    >barcode2
    TTTA
    """)
    input_path = random_path("csv")
    with input_path.open('w') as open_input:
        open_input.write(csv_contents)
    output_path = random_path("csv")
    run_component([
        "--input", input_path,
        "--output", output_path,
        "--header", 
        "--name_column", "barcode_name",
        "--sequence_column_index", "2",
        ]
    )
    assert output_path.is_file()
    with output_path.open('r') as open_output:
        output_contents = open_output.read()
    assert output_contents == expected

def test_csvtofasta_autodetect_dialect(run_component, random_path):
    csv_contents = dedent("""\
    barcode_name\tsome_other_column\tsequence
    barcode1\tfoo\tACGT
    barcode2\tbar\tTTTA
    """)

    expected= dedent("""\
    >barcode1
    ACGT
    >barcode2
    TTTA
    """)
    input_path = random_path("csv")
    with input_path.open('w') as open_input:
        open_input.write(csv_contents)
    output_path = random_path("csv")
    run_component([
        "--input", input_path,
        "--output", output_path,
        "--header", 
        "--name_column", "barcode_name",
        "--sequence_column_index", "2",
        ]
    )
    assert output_path.is_file()
    with output_path.open('r') as open_output:
        output_contents = open_output.read()
    assert output_contents == expected

    csv_contents = dedent("""\
    "barcode_name"\t"some_other_column"\t"sequence"
    "barcode1"\t"foo"\t"ACGT"
    "barcode2"\t"bar"\t"TTTA"
    """)

    expected= dedent("""\
    >barcode1
    ACGT
    >barcode2
    TTTA
    """)
    input_path = random_path("csv")
    with input_path.open('w') as open_input:
        open_input.write(csv_contents)
    output_path = random_path("csv")
    run_component([
        "--input", input_path,
        "--output", output_path,
        "--header", 
        "--name_column", "barcode_name",
        "--sequence_column_index", "2",
        ]
    )
    assert output_path.is_file()
    with output_path.open('r') as open_output:
        output_contents = open_output.read()
    assert output_contents == expected

def test_csvtofasta_header_select_column_by_name(run_component, random_path):
    csv_contents = dedent("""\
    barcode_name,some_other_column,sequence
    barcode1,foo,ACGT
    barcode2,bar,TTTA
    """)

    expected= dedent("""\
    >barcode1
    ACGT
    >barcode2
    TTTA
    """)
    input_path = random_path("csv")
    with input_path.open('w') as open_input:
        open_input.write(csv_contents)
    output_path = random_path("csv")
    run_component([
        "--input", input_path,
        "--output", output_path,
        "--header", 
        "--name_column", "barcode_name",
        "--sequence_column", "sequence"
        ]
    )
    assert output_path.is_file()
    with output_path.open('r') as open_output:
        output_contents = open_output.read()
    assert output_contents == expected

def test_csvtofasta_header_2_columns(run_component, random_path):
    csv_contents = dedent("""\
    barcode_name,sequence
    barcode1,ACGT
    barcode2,TTTA
    """)

    expected= dedent("""\
    >barcode1
    ACGT
    >barcode2
    TTTA
    """)
    input_path = random_path("csv")
    with input_path.open('w') as open_input:
        open_input.write(csv_contents)
    output_path = random_path("csv")
    run_component([
        "--input", input_path,
        "--output", output_path,
        "--header"
        ]
    )
    assert output_path.is_file()
    with output_path.open('r') as open_output:
        output_contents = open_output.read()
    assert output_contents == expected

def test_csvtofasta_2_columns(run_component, random_path):
    csv_contents = dedent("""\
    barcode1,ACGT
    barcode2,TTTA
    """)

    expected= dedent("""\
    >barcode1
    ACGT
    >barcode2
    TTTA
    """)
    input_path = random_path("csv")
    with input_path.open('w') as open_input:
        open_input.write(csv_contents)
    output_path = random_path("csv")
    run_component([
        "--input", input_path,
        "--output", output_path]
    )
    assert output_path.is_file()
    with output_path.open('r') as open_output:
        output_contents = open_output.read()
    assert output_contents == expected

def test_csvtofasta_2_columns_but_still_swap(run_component, random_path):
    csv_contents = dedent("""\
    ACGT,barcode1
    TTTA,barcode2
    """)

    expected= dedent("""\
    >barcode1
    ACGT
    >barcode2
    TTTA
    """)
    input_path = random_path("csv")
    with input_path.open('w') as open_input:
        open_input.write(csv_contents)
    output_path = random_path("csv")
    run_component([
        "--input", input_path,
        "--output", output_path,
        "--sequence_column_index", "0",
        "--name_column_index", "1"]
    )
    assert output_path.is_file()
    with output_path.open('r') as open_output:
        output_contents = open_output.read()
    assert output_contents == expected

def test_csvtofasta_2_columns_but_not_valid_sequence(run_component, random_path):
    csv_contents = dedent("""\
    barcodes,sequences
    barcode1,ACGT
    barcode2,TTTA
    """)

    input_path = random_path("csv")
    with input_path.open('w') as open_input:
        open_input.write(csv_contents)
    output_path = random_path("csv")
    with pytest.raises(CalledProcessError) as err:
        run_component([
            "--input", input_path,
            "--output", output_path]
        )
    assert re.search(r"ValueError: The sequence \('sequences'\) found on line "
                     r"1 contains characters \(.+\) which are not valid "
                     r"IUPAC identifiers for nucleotides\.", 
                     err.value.stdout.decode('utf-8'))

    csv_contents = dedent("""\
    barcodes,sequences
    barcode1,ACGT
    barcode2,TTEA
    """)

    input_path = random_path("csv")
    with input_path.open('w') as open_input:
        open_input.write(csv_contents)
    output_path = random_path("csv")
    with pytest.raises(CalledProcessError) as err:
        run_component([
            "--input", input_path,
            "--output", output_path,
            "--header"]
        )
    assert re.search(r"ValueError: The sequence \('TTEA'\) found on line "
                     r"3 contains characters \(E\) which are not valid "
                     r"IUPAC identifiers for nucleotides\.", 
                     err.value.stdout.decode('utf-8'))



if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))