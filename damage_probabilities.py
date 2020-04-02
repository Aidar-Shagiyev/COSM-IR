"""Calculate the required doublet and triplet DNA damage probabilities.

This module transforms experimentally derived probabilities of damages
into the requred doublet and triplet probabilities.
"""

PRIOR_DOSES = {"dsb": 75, "ssb": 50, "abasic": 75, "base": 50}  # In Gy.
_SSB_DOUBLET_PROBABILITIES = {  # https://doi.org/10.1667/RR14886.1
    "-1 0": {
        "GG": 0.29,
        "GC": 0.22,
        "TC": 0.20,
        "GA": 0.18,
        "AC": 0.14,
        "TG": 0.11,
        "CC": 0.10,
        "AG": 0.09,
        "TA": 0.07,
        "CG": 0.06,
    },
    "0 +1": {
        "CC": 0.26,
        "GG": 0.25,
        "CA": 0.23,
        "GC": 0.19,
        "AG": 0.16,
        "GT": 0.10,
        "GA": 0.06,
        "AT": 0.06,
        "CT": 0.05,
        "TG": 0.05,
    },
}

_BASE_DOUBLET_PROBABILITIES = {  # https://doi.org/10.1080/09553002.2019.1665216
    "-1 0": {
        "GG": 0.82,  # G: Fapy-G or 8-oxoG
        "AG": 0.41,  # A: Fapy-A or 8-oxoA
        "GC": 0.18,  # T: Thy-Gly
        "CG": 0.16,
        "GT": 0.14,
        "GA": 0.11,
        "TT": 0.11,
        "TG": 0.10,
        "AC": 0.07,
        "AT": 0.07,
        "CC": 0.04,
        "CT": 0.03,
        "TC": 0.03,
        "AA": 0.00,
        "CA": 0.00,
        "TA": 0.00,
    },
    "0 +1": {
        "GG": 0.59,
        "GC": 0.29,
        "GA": 0.28,
        "GT": 0.27,
        "CT": 0.11,
        "TA": 0.11,
        "TG": 0.10,
        "CA": 0.07,
        "CG": 0.05,
        "TT": 0.05,
        "TC": 0.05,
        "AA": 0.04,
        "AT": 0.03,
        "AC": 0.00,
        "AG": 0.00,
        "CC": 0.00,
    },
}

_DSB_DOUBLET_PROBABILITIES = {  # https://doi.org/10.1007/s11033-019-04815-6
    "endo+": {
        "-1 0": {
            "CC": 0.158,
            "GC": 0.141,
            "GG": 0.135,
            "CT": 0.091,
            "AG": 0.082,
            "TG": 0.076,
            "CA": 0.069,
            "AC": 0.067,
            "TC": 0.063,
            "GT": 0.049,
            "GA": 0.049,
            "CG": 0.037,
            "TT": 0.032,
            "AA": 0.031,
            "AT": 0.027,
            "TA": 0.020,
        },
        "0 +1": {
            "CC": 0.141,
            "GG": 0.101,
            "GC": 0.099,
            "CA": 0.086,
            "CT": 0.085,
            "TC": 0.078,
            "GA": 0.063,
            "GT": 0.061,
            "TG": 0.058,
            "TT": 0.058,
            "AC": 0.055,
            "AG": 0.051,
            "AA": 0.044,
            "AT": 0.037,
            "TA": 0.031,
            "CG": 0.020,
        },
    },
    "endo-": {
        "-1 0": {
            "CC": 0.122,
            "GC": 0.114,
            "CT": 0.109,
            "GG": 0.100,
            "AG": 0.078,
            "TG": 0.076,
            "TC": 0.066,
            "AC": 0.064,
            "CA": 0.060,
            "GT": 0.057,
            "GA": 0.049,
            "TT": 0.047,
            "AA": 0.041,
            "AT": 0.037,
            "TA": 0.027,
            "CG": 0.022,
        },
        "0 +1": {
            "CC": 0.140,
            "CA": 0.106,
            "GG": 0.092,
            "TC": 0.088,
            "GC": 0.085,
            "CT": 0.078,
            "GA": 0.070,
            "TG": 0.062,
            "TT": 0.052,
            "AA": 0.050,
            "GT": 0.047,
            "AG": 0.046,
            "AC": 0.046,
            "TA": 0.044,
            "AT": 0.030,
            "CG": 0.019,
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
                    prob = (left_prob * right_prob) ** 0.5
                    triplet_probabilities[type_name][triplet] = prob

    # Normalize abasic probabilities against ssb.
    ssb_sum = sum(triplet_probabilities["ssb"].values())
    abasic_sum = sum(triplet_probabilities["abasic"].values())
    normalization_coeff = ssb_sum / abasic_sum
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
