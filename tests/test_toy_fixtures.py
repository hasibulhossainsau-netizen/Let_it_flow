"""Guards the committed toy fixtures (a stale fixture would fail the integration test confusingly)."""
import hashlib

from utils.fasta import read_fasta

TOY = None


def _toy(repo_root):
    return repo_root / "tests" / "data" / "toy"


def test_fixture_shapes_and_embedding(repo_root):
    (_, gene), = read_fasta(_toy(repo_root) / "gene.fasta")
    (_, genome), = read_fasta(_toy(repo_root) / "genome.fasta")
    assert len(gene) == 100 and len(genome) == 1000
    assert genome.count(gene) == 1, "gene must be embedded in the genome exactly once"
    assert set(gene + genome) <= set("ACGT")


def test_fixture_checksums_match(repo_root):
    toy = _toy(repo_root)
    listed = dict(line.split("  ")[::-1] for line in
                  (toy / "checksums.md5").read_text().split("\n") if line)
    for name, dest in (("genome.fasta", "data/raw/toy_genome.fasta"),
                       ("gene.fasta", "data/raw/toy_gene.fasta")):
        assert listed[dest] == hashlib.md5((toy / name).read_bytes()).hexdigest()


def test_generator_reproduces_committed_fixtures(repo_root):
    import sys
    sys.path.insert(0, str(repo_root / "tests"))
    import make_toy_fixtures as gen
    gene, genome = gen.build(gen.DEFAULT_SEED)
    assert read_fasta(_toy(repo_root) / "gene.fasta")[0][1] == gene
    assert read_fasta(_toy(repo_root) / "genome.fasta")[0][1] == genome
