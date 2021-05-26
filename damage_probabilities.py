"""Calculate the required doublet and triplet DNA damage probabilities.

This module transforms experimentally derived probabilities of damages
into the requred doublet and triplet probabilities.
"""

MEAN_AP_PROBABILITY = 25 / 1_000_000  # At 1 Gy; https://doi.org/10.1021/bi9927989

PRIOR_DOSES = {"dsb": 1, "ssb": 1, "abasic": 1, "base": 1}  # In Gy.
# TODO: mean probabilities and doses
_SSB_DOUBLET_PROBABILITIES = {  # https://doi.org/10.1667/RR14886.1
    "-1 0": {
        "CC": 29 / 1_000_000,
        "GC": 29 / 1_000_000,
        "GG": 29 / 1_000_000,
        "CT": 29 / 1_000_000,
        "AG": 29 / 1_000_000,
        "TG": 29 / 1_000_000,
        "CA": 29 / 1_000_000,
        "AC": 29 / 1_000_000,
        "TC": 29 / 1_000_000,
        "GT": 29 / 1_000_000,
        "GA": 29 / 1_000_000,
        "CG": 29 / 1_000_000,
        "TT": 29 / 1_000_000,
        "AA": 29 / 1_000_000,
        "AT": 29 / 1_000_000,
        "TA": 29 / 1_000_000,
    },
    "0 +1": {
        "CC": 29 / 1_000_000,
        "GC": 29 / 1_000_000,
        "GG": 29 / 1_000_000,
        "CT": 29 / 1_000_000,
        "AG": 29 / 1_000_000,
        "TG": 29 / 1_000_000,
        "CA": 29 / 1_000_000,
        "AC": 29 / 1_000_000,
        "TC": 29 / 1_000_000,
        "GT": 29 / 1_000_000,
        "GA": 29 / 1_000_000,
        "CG": 29 / 1_000_000,
        "TT": 29 / 1_000_000,
        "AA": 29 / 1_000_000,
        "AT": 29 / 1_000_000,
        "TA": 29 / 1_000_000,
    },
}

_BASE_DOUBLET_PROBABILITIES = {  # https://doi.org/10.1080/09553002.2019.1665216
    "-1 0": {
        "GG": 29 / 1_000_000,  # G: Fapy-G or 8-oxoG
        "AG": 29 / 1_000_000,  # A: Fapy-A or 8-oxoA
        "GC": 10 / 1_000_000,  # T: Thy-Gly
        "CG": 29 / 1_000_000,
        "GT": 10 / 1_000_000,
        "GA": 29 / 1_000_000,
        "TT": 10 / 1_000_000,
        "TG": 29 / 1_000_000,
        "AC": 10 / 1_000_000,
        "AT": 10 / 1_000_000,
        "CC": 10 / 1_000_000,
        "CT": 10 / 1_000_000,
        "TC": 10 / 1_000_000,
        "AA": 29 / 1_000_000,
        "CA": 29 / 1_000_000,
        "TA": 29 / 1_000_000,
    },
    "0 +1": {
        "GG": 29 / 1_000_000,
        "GC": 29 / 1_000_000,
        "GA": 29 / 1_000_000,
        "GT": 29 / 1_000_000,
        "CT": 10 / 1_000_000,
        "TA": 10 / 1_000_000,
        "TG": 10 / 1_000_000,
        "CA": 10 / 1_000_000,
        "CG": 10 / 1_000_000,
        "TT": 10 / 1_000_000,
        "TC": 10 / 1_000_000,
        "AA": 29 / 1_000_000,
        "AT": 29 / 1_000_000,
        "AC": 29 / 1_000_000,
        "AG": 29 / 1_000_000,
        "CC": 10 / 1_000_000,
    },
}

_DSB_DOUBLET_PROBABILITIES = {  # https://doi.org/10.1007/s11033-019-04815-6
    "endo+": {
        "-1 0": {
            "CC": 3 / 1_000_000,
            "GC": 3 / 1_000_000,
            "GG": 3 / 1_000_000,
            "CT": 3 / 1_000_000,
            "AG": 3 / 1_000_000,
            "TG": 3 / 1_000_000,
            "CA": 3 / 1_000_000,
            "AC": 3 / 1_000_000,
            "TC": 3 / 1_000_000,
            "GT": 3 / 1_000_000,
            "GA": 3 / 1_000_000,
            "CG": 3 / 1_000_000,
            "TT": 3 / 1_000_000,
            "AA": 3 / 1_000_000,
            "AT": 3 / 1_000_000,
            "TA": 3 / 1_000_000,
        },
        "0 +1": {
            "CC": 3 / 1_000_000,
            "GG": 3 / 1_000_000,
            "GC": 3 / 1_000_000,
            "CA": 3 / 1_000_000,
            "CT": 3 / 1_000_000,
            "TC": 3 / 1_000_000,
            "GA": 3 / 1_000_000,
            "GT": 3 / 1_000_000,
            "TG": 3 / 1_000_000,
            "TT": 3 / 1_000_000,
            "AC": 3 / 1_000_000,
            "AG": 3 / 1_000_000,
            "AA": 3 / 1_000_000,
            "AT": 3 / 1_000_000,
            "TA": 3 / 1_000_000,
            "CG": 3 / 1_000_000,
        },
    },
    "endo-": {
        "-1 0": {
            "CC": 2.7 / 1_000_000,
            "GC": 2.7 / 1_000_000,
            "CT": 2.7 / 1_000_000,
            "GG": 2.7 / 1_000_000,
            "AG": 2.7 / 1_000_000,
            "TG": 2.7 / 1_000_000,
            "TC": 2.7 / 1_000_000,
            "AC": 2.7 / 1_000_000,
            "CA": 2.7 / 1_000_000,
            "GT": 2.7 / 1_000_000,
            "GA": 2.7 / 1_000_000,
            "TT": 2.7 / 1_000_000,
            "AA": 2.7 / 1_000_000,
            "AT": 2.7 / 1_000_000,
            "TA": 2.7 / 1_000_000,
            "CG": 2.7 / 1_000_000,
        },
        "0 +1": {
            "CC": 2.7 / 1_000_000,
            "CA": 2.7 / 1_000_000,
            "GG": 2.7 / 1_000_000,
            "TC": 2.7 / 1_000_000,
            "GC": 2.7 / 1_000_000,
            "CT": 2.7 / 1_000_000,
            "GA": 2.7 / 1_000_000,
            "TG": 2.7 / 1_000_000,
            "TT": 2.7 / 1_000_000,
            "AA": 2.7 / 1_000_000,
            "GT": 2.7 / 1_000_000,
            "AG": 2.7 / 1_000_000,
            "AC": 2.7 / 1_000_000,
            "TA": 2.7 / 1_000_000,
            "AT": 2.7 / 1_000_000,
            "CG": 2.7 / 1_000_000,
        },
    },
}

# Extract abasic probabilities from DSB data like so:
# endo+ = 1 - (1 - endo-)(1 - abasic)
# abasic = 1 - (1 - endo+) / (1 - endo-)

_ABASIC_DOUBLET_PROBABILITIES = {}
for position in ["-1 0", "0 +1"]:
    assert set(_DSB_DOUBLET_PROBABILITIES["endo+"][position].keys()) == set(
        _DSB_DOUBLET_PROBABILITIES["endo-"][position].keys()
    )
    _ABASIC_DOUBLET_PROBABILITIES[position] = {}
    for doublet, endo_plus in _DSB_DOUBLET_PROBABILITIES["endo+"][position].items():
        endo_minus = _DSB_DOUBLET_PROBABILITIES["endo-"][position][doublet]
        abasic_probability = max(0, 1 - (1 - endo_plus) / (1 - endo_minus))
        _ABASIC_DOUBLET_PROBABILITIES[position][doublet] = abasic_probability


def get_triplet_probabilities():
    """Calculate probabilities for triplets."""
    triplet_probabilities = {}
    for doublet_probabilities, type_name in zip(
        [
            _SSB_DOUBLET_PROBABILITIES,
            _BASE_DOUBLET_PROBABILITIES,
            _DSB_DOUBLET_PROBABILITIES["endo-"],
            _ABASIC_DOUBLET_PROBABILITIES,
        ],
        ["ssb", "base", "dsb", "abasic"],
    ):
        triplet_probabilities[type_name] = {}
        for left_doublet, left_prob in doublet_probabilities["-1 0"].items():
            for right_doublet, right_prob in doublet_probabilities["0 +1"].items():
                if left_doublet[-1] == right_doublet[0]:
                    triplet = left_doublet[0] + right_doublet
                    prob = (left_prob + right_prob) / 2
                    triplet_probabilities[type_name][triplet] = prob

    # Normalize abasic probabilities against mean AP probability.
    mean_ap = 1 - (1 - MEAN_AP_PROBABILITY) ** PRIOR_DOSES["abasic"]
    abasic_sum = sum(triplet_probabilities["abasic"].values())
    normalization_coeff = 4**3 * mean_ap / abasic_sum
    triplet_probabilities["abasic"] = {
        triplet: prob * normalization_coeff
        for triplet, prob in triplet_probabilities["abasic"].items()
    }
    return triplet_probabilities


def get_doublet_probabilities():
    """Return only DSB and SSB probabilities."""
    return {
        "dsb": _DSB_DOUBLET_PROBABILITIES["endo-"],
        "ssb": _SSB_DOUBLET_PROBABILITIES,
    }
