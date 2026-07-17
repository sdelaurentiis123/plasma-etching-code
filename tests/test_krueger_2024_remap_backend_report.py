import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "krueger_2024_remap_backend_report.py"
SPEC = importlib.util.spec_from_file_location("krueger_2024_remap_backend_report", SCRIPT)
REPORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REPORT)


def test_summary_separates_global_and_local_remap_sensitivity():
    def case(depth, opening, width, thickness, residual):
        return {
            "status": "complete",
            "steps": [{
                "metrics": {
                    "etch_depth_nm": depth,
                    "mask_opening_nm": opening,
                    "top_feature_width_nm": width,
                    "remaining_mask_thickness_nm": thickness,
                },
                "maximum_remap_relative_conservation_residual": residual,
            }],
        }

    audit = {"cases": {
        "legacy_knn": case(1.0, 100.0, 80.0, 850.0, 4e-16),
        "indexed_knn": case(1.0, 100.0, 80.0, 850.0, 5e-16),
        "common_refinement": case(1.00001, 100.0002, 79.8, 849.9998, 1e-15),
        "partitioned_overlap": {
            "status": "refused",
            "exception": {"message": "nonparallel"},
        },
    }}

    summary = REPORT.build_summary(audit)

    relative = dict(summary["relative_change_ppm"])
    assert relative["Depth"] > relative["Opening"] > 0.0
    assert relative["Top width"] < 0.0
    assert summary["partitioned_overlap_status"] == "refused"
