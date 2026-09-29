import pytest

from qc import classify, infer_sex


@pytest.mark.parametrize("x_het,y_call,sex", [
    (0.001, 0.99, "male"),
    (0.30, 0.093, "female"),  # female arrays carry some noise calls on Y probes
    (0.30, 0.0, "female"),
    (0.07, 0.6, "unknown"),   # between thresholds: never guessed
    (0.001, 0.0, "unknown"),  # e.g. X0 or a chip without Y probes
])
def test_infer_sex(x_het, y_call, sex):
    assert infer_sex(x_het, y_call) == sex


@pytest.mark.parametrize("gt,cls", [("--", "no_call"), ("AA", "hom"), ("AG", "het"), ("DI", "indel"), ("II", "indel")])
def test_classify(gt, cls):
    assert classify(gt) == cls
