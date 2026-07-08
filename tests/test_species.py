import pytest

from srnataps.species import canonical_species, reference_details, supported_species


def test_supported_species_include_common_models():
    assert {"human", "mouse", "rat", "zebrafish"} <= set(supported_species())


def test_species_aliases_are_normalized():
    assert canonical_species("hg38") == "human"
    assert canonical_species("Mus musculus") == "mouse"
    assert canonical_species("fruit-fly") == "fruit_fly"


def test_unknown_species_has_clear_error():
    with pytest.raises(ValueError, match="Supported species"):
        canonical_species("capybara")


def test_mouse_reference_paths_and_urls(tmp_path):
    details = reference_details(
        {"species": "mouse", "ensembl_release": 112}, tmp_path
    )
    assert details["genome_fa"].endswith("Mus_musculus.GRCm39.dna.toplevel.fa")
    assert details["gtf"].endswith("Mus_musculus.GRCm39.112.gtf")
    assert "/release-112/fasta/mus_musculus/dna/" in details["genome_url"]
    assert details["index_stem"] == "mm39"


def test_custom_reference_paths_override_species_defaults(tmp_path):
    details = reference_details(
        {
            "species": "zebrafish",
            "genome_fa": "/refs/custom.fa",
            "gtf": "/refs/custom.gtf",
            "bowtie1_index": "/refs/custom-index",
        },
        tmp_path,
    )
    assert details["genome_fa"] == "/refs/custom.fa"
    assert details["gtf"] == "/refs/custom.gtf"
    assert details["bowtie1_index"] == "/refs/custom-index"
