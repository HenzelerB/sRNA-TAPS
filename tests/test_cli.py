# -*- coding: utf-8 -*-
"""
tests/test_cli.py — Tests for srnataps.cli

Tests the CLI entry points using click.testing.CliRunner.
Does not require a real filesystem or cluster — uses temp directories.
"""

import os
import yaml
import pytest
from pathlib import Path
from click.testing import CliRunner

from srnataps.cli import cli


@pytest.fixture
def runner():
    return CliRunner()


class TestCliVersion:

    def test_version_flag(self, runner):
        result = runner.invoke(cli, ["--version"])
        assert result.exit_code == 0
        assert "0.1.0" in result.output


class TestCliHelp:

    def test_main_help(self, runner):
        result = runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "run" in result.output
        assert "module" in result.output
        assert "init" in result.output
        assert "check" in result.output

    def test_run_help(self, runner):
        result = runner.invoke(cli, ["run", "--help"])
        assert result.exit_code == 0
        assert "--configfile" in result.output
        assert "--slurm" in result.output
        assert "--benchmark" in result.output

    def test_module_help(self, runner):
        result = runner.invoke(cli, ["module", "--help"])
        assert result.exit_code == 0
        assert "fastqc" in result.output
        assert "call" in result.output
        assert "benchmark" in result.output

    def test_check_help(self, runner):
        result = runner.invoke(cli, ["check", "--help"])
        assert result.exit_code == 0
        assert "--configfile" in result.output


class TestCliInit:

    def test_init_creates_directory(self, runner, tmp_path):
        result = runner.invoke(cli, [
            "init",
            "--outdir", str(tmp_path / "new_project"),
        ])
        assert result.exit_code == 0
        assert (tmp_path / "new_project").exists()

    def test_init_creates_config(self, runner, tmp_path):
        outdir = tmp_path / "new_project"
        runner.invoke(cli, ["init", "--outdir", str(outdir)])
        assert (outdir / "config.yaml").exists()

    def test_init_creates_samples_tsv(self, runner, tmp_path):
        outdir = tmp_path / "new_project"
        runner.invoke(cli, ["init", "--outdir", str(outdir)])
        assert (outdir / "samples.tsv").exists()

    def test_init_fills_genome_path(self, runner, tmp_path):
        outdir = tmp_path / "new_project"
        runner.invoke(cli, [
            "init",
            "--outdir", str(outdir),
            "--genome", "/fake/path/hg38.fa",
        ])
        with open(outdir / "config.yaml") as f:
            config = yaml.safe_load(f)
        assert config["reference"]["genome_fa"] == "/fake/path/hg38.fa"

    def test_init_fills_outdir_in_config(self, runner, tmp_path):
        outdir = tmp_path / "new_project"
        runner.invoke(cli, ["init", "--outdir", str(outdir)])
        with open(outdir / "config.yaml") as f:
            config = yaml.safe_load(f)
        assert config["project"]["outdir"] == str(outdir)

    def test_init_fails_if_exists_without_force(self, runner, tmp_path):
        outdir = tmp_path / "existing"
        outdir.mkdir()
        result = runner.invoke(cli, ["init", "--outdir", str(outdir)])
        assert result.exit_code != 0
        assert "already exists" in result.output.lower() or "force" in result.output.lower()

    def test_init_force_overwrites(self, runner, tmp_path):
        outdir = tmp_path / "existing"
        outdir.mkdir()
        result = runner.invoke(cli, ["init", "--outdir", str(outdir), "--force"])
        assert result.exit_code == 0


class TestCliModule:

    def test_invalid_module_fails(self, runner, tmp_path):
        # Create a minimal config so the runner doesn't fail on missing file
        config = tmp_path / "config.yaml"
        config.write_text("project:\n  outdir: /tmp\ninput:\n  samples_tsv: s.tsv\n")
        result = runner.invoke(cli, [
            "module", "nonexistent_module",
            "--configfile", str(config),
        ])
        assert result.exit_code != 0

    def test_all_module_names_accepted(self, runner, tmp_path):
        """All valid module names should be accepted by the CLI parser."""
        config = tmp_path / "config.yaml"
        config.write_text("project:\n  outdir: /tmp\n")

        valid_modules = [
            "fastqc", "trim", "index", "align",
            "biotype", "snp", "call", "benchmark", "compare",
        ]
        for mod in valid_modules:
            result = runner.invoke(cli, [
                "module", mod, "--help"
            ])
            assert result.exit_code == 0, \
                f"Module '{mod}' should be accepted, got: {result.output}"
