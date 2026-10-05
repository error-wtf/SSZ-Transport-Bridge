"""Provenance negative controls use explicitly synthetic repository fixtures."""
import json
import subprocess

import pytest

from transport_bridge.provenance import bind_module, load_closure_provenance, sha256_file


@pytest.fixture
def fixture_repo(tmp_path):
    repo = tmp_path / 'fixture'
    repo.mkdir()
    subprocess.run(['git', 'init', '-q', str(repo)], check=True)
    member = repo / 'member.csv'
    member.write_text('fixture-only\n')
    digest = sha256_file(member)
    files = {
        'manifest.json': {'member_hash': digest},
        'MODEL_LOCK.json': {'action_member_stream': 'member.csv',
                            'action_member_reference': 'manifest.json',
                            'action_member_sha256': digest},
        'TRUE_FULL_CLOSURE_VERDICT.json': {'git_commit': 'a' * 40, 'member_hash': digest,
                                         'verdict': 'TRUE_FULL_CLOSURE_PASS'},
        'artifacts/true_full_closure/physics_dependency_graph.json': {
            'member_hash': digest, 'nodes': []},
    }
    for name, data in files.items():
        p = repo / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(data))
    subprocess.run(['git', '-C', str(repo), 'add', '.'], check=True)
    subprocess.run(['git', '-C', str(repo), '-c', 'user.name=Fixture',
                    '-c', 'user.email=fixture@example.invalid',
                    'commit', '-qm', 'fixture'], check=True)
    return repo


def test_actual_member_and_lock_bound(fixture_repo):
    p = load_closure_provenance(fixture_repo)
    assert p['member_sha256'] == sha256_file(fixture_repo / 'member.csv')
    assert len(p['commit_sha']) == 40
    assert p['files']['MODEL_LOCK.json'] == sha256_file(fixture_repo / 'MODEL_LOCK.json')
    assert not p['historical_evidence_at_head']


def test_member_tamper_rejected(fixture_repo):
    (fixture_repo / 'member.csv').write_text('TAMPER\n')
    with pytest.raises(ValueError, match='member'):
        load_closure_provenance(fixture_repo)


def test_manifest_tamper_rejected(fixture_repo):
    (fixture_repo / 'manifest.json').write_text(json.dumps({'member_hash': '0' * 64}))
    with pytest.raises(ValueError, match='member'):
        load_closure_provenance(fixture_repo)


def test_path_escape_rejected(fixture_repo):
    p = fixture_repo / 'MODEL_LOCK.json'
    data = json.loads(p.read_text())
    data['action_member_stream'] = '../outside.csv'
    p.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='escape'):
        load_closure_provenance(fixture_repo)


def test_missing_verdict_rejected(fixture_repo):
    (fixture_repo / 'TRUE_FULL_CLOSURE_VERDICT.json').unlink()
    with pytest.raises(FileNotFoundError):
        load_closure_provenance(fixture_repo)


def test_wrong_import_origin_rejected(fixture_repo):
    import json as unrelated_module
    with pytest.raises(ValueError, match='origin'):
        bind_module(unrelated_module, fixture_repo)
