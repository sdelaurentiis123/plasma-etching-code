import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = (
    ROOT / "scripts" / "audit_reactor_global_chlorine_ionization.py")


def _load_audit():
    specification = importlib.util.spec_from_file_location(
        "audit_reactor_global_chlorine_ionization", AUDIT)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def test_evaluated_ionization_audit_preserves_unresolved_branching():
    rows, summary = _load_audit().run_audit()
    assert len(rows) == 4
    assert summary["electron_temperature_eV"] == [2.0, 3.0, 4.0, 5.0]
    assert "Lennon Eq. 6" in summary["lower_temperature_exclusion"]
    assert summary["coefficient_selection_target"] is None
    assert summary["molecular_species_branching"] == "unresolved"
    assert (
        summary["predictive_use"][
            "species_resolved_molecular_positive_ion_source"]
        == "not ready"
    )
    atomic_range = summary["atomic_cl_nist_to_lennon_ratio_range"]
    molecular_range = summary[
        "molecular_cl2_nist_to_lee_total_ratio_range"]
    assert 0.80 < atomic_range[0] < atomic_range[1] < 0.90
    assert 0.65 < molecular_range[0] < 0.75
    assert 1.10 < molecular_range[1] < 1.20
