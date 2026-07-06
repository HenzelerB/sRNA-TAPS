from srnataps.three_letter import convert_fastq, convert_reference


def test_convert_reference_preserves_headers_and_coordinates(tmp_path):
    source = tmp_path / "genome.fa"
    source.write_text(">chr1\nACGTCCGA\n")
    c2t = tmp_path / "c2t.fa"
    g2a = tmp_path / "g2a.fa"

    convert_reference(source, c2t, g2a)

    assert c2t.read_text() == ">chr1\nATGTTTGA\n"
    assert g2a.read_text() == ">chr1\nACATCCAA\n"


def test_convert_fastq_preserves_quality_and_converts_only_sequence(tmp_path):
    source = tmp_path / "reads.fastq"
    source.write_text("@read1 comment\nACCG\n+\nCDEF\n")
    output = tmp_path / "converted.fastq"

    count = convert_fastq(source, output)

    assert count == 1
    assert output.read_text() == "@read1 comment\nATTG\n+\nCDEF\n"
