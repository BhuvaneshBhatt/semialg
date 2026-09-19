"""Published TTICAD Section 8.2 benchmark dataset (29 cases)."""

from __future__ import annotations

from dataclasses import dataclass

from .cases import BenchmarkCase, BenchmarkDomain, Difficulty, Provenance

SOURCE = Provenance(
    citation="Bradford, Davenport, England, McCallum & Wilson, Truth Table Invariant Cylindrical Algebraic Decomposition, Section 8.2 dataset",
    url="https://doi.org/10.15125/BATH-00076",
    identifier="Section82-ExampleSet",
    license="CC BY-SA 4.0",
)


@dataclass(frozen=True)
class TTICADFormula:
    equational_constraints: tuple[str, ...]
    other_constraints: tuple[str, ...]


@dataclass(frozen=True)
class TTICADCase:
    number: int
    name: str
    short_name: str
    formulas: tuple[TTICADFormula, ...]
    variables: tuple[str, ...]


# EC entries are equalities to zero; other entries retain the source polynomial whose
# sign constraint is specified by the original experiment. This preserves formula grouping.
def F(ec, other=()):
    if ec is None:
        ecs = ()
    elif isinstance(ec, tuple):
        ecs = ec
    else:
        ecs = (ec,)
    return TTICADFormula(ecs, tuple(other))


_BASE_DATA = (
    TTICADCase(
        1,
        "Intersection A",
        "IntA",
        (F(("x**2-y**2/2-z**2/2", "x*z+z*y-2*x", "z**2-y")),),
        ("z", "y", "x"),
    ),
    TTICADCase(
        2,
        "Intersection B",
        "IntB",
        (F(("x**2-y**2/2-z**2/2", "x*z+z*y-2*x", "z**2-y")),),
        ("z", "x", "y"),
    ),
    TTICADCase(
        3,
        "Random A",
        "RanA",
        (F(("4*x**2+x*y**2-z+1/4", "2*x+y**2*z+1/2", "x**2*z-x/2-y**2")),),
        ("z", "y", "x"),
    ),
    TTICADCase(
        4,
        "Random B",
        "RanB",
        (F(("4*x**2+x*y**2-z+1/4", "2*x+y**2*z+1/2", "x**2*z-x/2-y**2")),),
        ("x", "y", "z"),
    ),
    TTICADCase(
        5,
        "Intersection dagger A",
        "Int†A",
        (F("x*z+z*y-2*x"), F("z**2-y", ("2*x**2-y**2-z**2",))),
        ("z", "y", "x"),
    ),
    TTICADCase(
        6,
        "Intersection dagger B",
        "Int†B",
        (F("x*z+z*y-2*x"), F("z**2-y", ("2*x**2-y**2-z**2",))),
        ("z", "x", "y"),
    ),
    TTICADCase(
        7,
        "Random dagger A",
        "Ran†A",
        (F("16*x**2+4*x*y**2-4*z+1"), F("4*x+2*y**2*z+1", ("2*x**2*z-x-2*y**2",))),
        ("z", "y", "x"),
    ),
    TTICADCase(
        8,
        "Random dagger B",
        "Ran†B",
        (F("16*x**2+4*x*y**2-4*z+1"), F("4*x+2*y**2*z+1", ("2*x**2*z-x-2*y**2",))),
        ("z", "x", "y"),
    ),
    TTICADCase(
        9,
        "Ellipse dagger A",
        "Ell†A",
        (
            F("x**2+y**2-1"),
            F(
                "b**2*x**2-2*b**2*x*c+b**2*c**2+a**2*y**2-a**2*b**2",
                ("-a", "a-1", "-b", "b-1", "-c", "c-1"),
            ),
        ),
        ("y", "x", "c", "b", "a"),
    ),
    TTICADCase(
        10,
        "Ellipse dagger B",
        "Ell†B",
        (
            F("x**2+y**2-1"),
            F(
                "b**2*x**2-2*b**2*x*c+b**2*c**2+a**2*y**2-a**2*b**2",
                ("-a", "a-1", "-b", "b-1", "-c", "c-1"),
            ),
        ),
        ("x", "y", "c", "b", "a"),
    ),
    TTICADCase(
        11,
        "Solotareff dagger A",
        "Solo†A",
        (
            F("3*x**2-2*x-a", ("-(x**3-x**2-a*x-2*b+a-2)", "-(4*a-1)", "4*a-7", "-(x+1)", "x")),
            F("3*y**2-2*y-a", ("-(y**3-y**2-a*y-a+2)", "-(4*b+3)", "4*b-3", "-y", "y-1")),
        ),
        ("y", "x", "b", "a"),
    ),
    TTICADCase(
        12,
        "Solotareff dagger B",
        "Solo†B",
        (
            F("3*x**2-2*x-a", ("-(x**3-x**2-a*x-2*b+a-2)", "-(4*a-1)", "4*a-7", "-(x+1)", "x")),
            F("3*y**2-2*y-a", ("-(y**3-y**2-a*y-a+2)", "-(4*b+3)", "4*b-3", "-y", "y-1")),
        ),
        ("y", "x", "a", "b"),
    ),
    TTICADCase(
        13,
        "Collision dagger A",
        "Coll†A",
        (
            F("x**2-2*x*t+t**2+4*y**2-80*y+396", ("-t",)),
            F("x**2-2*x*a*t+5*a**2*t**2+4*y**2-8*y*a*t-4", ("-a",)),
        ),
        ("y", "x", "t", "a"),
    ),
    TTICADCase(
        14,
        "Collision dagger B",
        "Coll†B",
        (
            F("x**2-2*x*t+t**2+4*y**2-80*y+396", ("-t",)),
            F("x**2-2*x*a*t+5*a**2*t**2+4*y**2-8*y*a*t-4", ("-a",)),
        ),
        ("t", "x", "y", "a"),
    ),
    TTICADCase(
        17,
        "ArcSin A",
        "AsinA",
        (
            F("-2*x*y", ("y**2-x**2+1",)),
            F("16*y**3*x-16*y*x**3+8*y*x", ("1-4*x**2+4*y**2+4*x**4-24*x**2*y**2+4*y**4",)),
            F("y", ("1-x",)),
            F("y", ("x+1",)),
        ),
        ("y", "x"),
    ),
    TTICADCase(
        18,
        "ArcSin B",
        "AsinB",
        (
            F("-2*x*y", ("y**2-x**2+1",)),
            F("16*y**3*x-16*y*x**3+8*y*x", ("1-4*x**2+4*y**2+4*x**4-24*x**2*y**2+4*y**4",)),
            F("y", ("1-x",)),
            F("y", ("x+1",)),
        ),
        ("x", "y"),
    ),
    TTICADCase(
        19,
        "Example Phi (JSC) A",
        "ExPhiA",
        (F("x**2+y**2-1", ("4*x*y-1",)), F("x**2-8*x+16+y**2-2*y", ("4*x*y-4*x-16*y+15",))),
        ("y", "x"),
    ),
    TTICADCase(
        20,
        "Example Phi (JSC) B",
        "ExPhiB",
        (F("x**2+y**2-1", ("4*x*y-1",)), F("x**2-8*x+16+y**2-2*y", ("4*x*y-4*x-16*y+15",))),
        ("x", "y"),
    ),
    TTICADCase(
        21,
        "Example Psi (JSC) A",
        "ExPsiA",
        (F("x**2+y**2-1", ("x*y-1/4",)), F(None, ("(x-4)**2+(y-1)**2-1", "(x-4)*(y-1)-1/4"))),
        ("y", "x"),
    ),
    TTICADCase(
        22,
        "Example Psi (JSC) B",
        "ExPsiB",
        (F("x**2+y**2-1", ("x*y-1/4",)), F(None, ("(x-4)**2+(y-1)**2-1", "(x-4)*(y-1)-1/4"))),
        ("x", "y"),
    ),
    TTICADCase(
        23,
        "3D Example A",
        "Ex25A",
        (
            F("x**2+y**2+z**2-1", ("4*x*y*z-1",)),
            F("x**2-8*x+y**2-2*y+z**2-4*z+20", ("4*y*x*z-8*y*x-4*x*z+8*x-16*y*z+32*y+16*z-33",)),
        ),
        ("z", "y", "x"),
    ),
    TTICADCase(
        24,
        "3D Example B",
        "Ex25B",
        (
            F("x**2+y**2+z**2-1", ("4*x*y*z-1",)),
            F("x**2-8*x+y**2-2*y+z**2-4*z+20", ("4*y*x*z-8*y*x-4*x*z+8*x-16*y*z+32*y+16*z-33",)),
        ),
        ("z", "x", "y"),
    ),
)
# Large Kahan cases are copied from the source with formula grouping preserved.
KAHAN_EC = "8*y**3*x+8*y*x**3+20*y**3+84*y*x**2+288*y*x+324*y"
KAHAN = (
    F(KAHAN_EC, ("-225*x**2-324*x+63*y**2-4*x**4-52*x**3+12*y**2*x+4*y**4",)),
    F("2*y", ("2*x+9",)),
    F("8*y", ("8*x**2+56*x+8*y**2+96",)),
    F("y", ("x**2+7*x+y**2+12",)),
    F(
        KAHAN_EC,
        (
            "-4*x**4-52*x**3-252*x**2+12*y**2*x-540*x+36*y**2+4*y**4-432",
            "4*x**4+52*x**3+225*x**2-12*y**2*x+324*x-63*y**2-4*y**4",
        ),
    ),
    F("2*y", ("-2*x-6", "2*x")),
    F("8*y", ("-8*x**2-56*x-8*y**2-96", "2*x**2+8*x+2*y**2")),
)
RANDOM = (
    (
        25,
        "JSC Random Example 1",
        "Rand1",
        F("-55*x-94*y+87*z-56", ("-62*y+97*z-73",)),
        F(None, ("-4*x-83*y-10*z+62", "-75-10*x**2-7*x*y-40*x*z+42*y*z-50*z**2")),
    ),
    (
        26,
        "JSC Random Example 2",
        "Rand2",
        F("-81*x-6*y-51*z-29", ("-14*x-48*y+97*z-12",)),
        F(None, ("83*x-24*y-8*z+47", "60+46*y-31*x**2-91*x*y+98*x*z+2*y**2")),
    ),
    (
        27,
        "JSC Random Example 3",
        "Rand3",
        F("-77*x+38*y+42*z+8", ("-6*x+23*y+76",)),
        F(None, ("-11*x+79*y+49*z-40", "-51-4*y+29*z-44*x*y+65*x*z+56*y**2")),
    ),
    (
        28,
        "JSC Random Example 4",
        "Rand4",
        F("27*x-95*y+84*z+48", ("14*x+64*y-88*z+18",)),
        F(None, ("56*x+96*y+77*z+54", "-94+90*x**2+85*x*y-72*x*z-2*y**2+32*z**2")),
    ),
    (
        29,
        "JSC Random Example 5",
        "Rand5",
        F("5*x-5*y-43*z+51", ("-45*x+30*y+21*z-12",)),
        F(None, ("82*x+23*y+63*z-67", "48-90*y-14*z-75*x*y+77*x*z+38*y**2")),
    ),
)
DATA = (
    _BASE_DATA[:14]
    + (
        TTICADCase(15, "Kahan A", "Ex33A", KAHAN, ("y", "x")),
        TTICADCase(16, "Kahan B", "Ex33B", KAHAN, ("x", "y")),
    )
    + _BASE_DATA[14:]
    + tuple(TTICADCase(n, name, short, (a, b), ("x", "y", "z")) for n, name, short, a, b in RANDOM)
)


# Cell counts from Table 2 of Bradford et al. (2016). These are reference
# results for the paper's named configurations, not intrinsic case properties.
_TABLE2_CELLS = {
    1: (3707, 269),
    2: (2985, 303),
    3: (2093, 435),
    4: (4097, 711),
    5: (3707, 575),
    6: (2985, 601),
    7: (2093, 663),
    8: (4097, 1075),
    9: (None, None),
    10: (None, None),
    11: (54037, None),
    12: (154527, None),
    13: (8387, 8387),
    14: (None, None),
    15: (409, 55),
    16: (1143, 39),
    17: (225, 57),
    18: (393, 25),
    19: (317, 105),
    20: (377, 153),
    21: (317, 183),
    22: (377, 233),
    23: (5453, 109),
    24: (6413, 153),
    25: (1533, 1533),
    26: (7991, 2911),
    27: (8889, 4005),
    28: (11979, 4035),
    29: (11869, 4905),
}


def tticad_cases() -> tuple[BenchmarkCase, ...]:
    out = []
    for c in DATA:
        polys = tuple(p for f in c.formulas for p in f.equational_constraints + f.other_constraints)
        out.append(
            BenchmarkCase(
                name=f"tticad_{c.number:02d}_{c.short_name.replace('†', 'dagger').lower()}",
                domain=BenchmarkDomain.CAD,
                description=c.name,
                variables=c.variables,
                polynomials=polys,
                tags=("literature", "tticad"),
                provenance=Provenance(
                    SOURCE.citation,
                    SOURCE.url,
                    f"Section82:{c.number}:{c.short_name}",
                    SOURCE.license,
                ),
                difficulty=Difficulty(len(c.variables), polynomial_count=len(polys)),
                expected={
                    "source_identifier": f"Section82:{c.number}",
                    "variable_order": c.variables,
                    "formula_count": len(c.formulas),
                    "formulae": tuple(
                        {
                            "equational_constraints": f.equational_constraints,
                            "other_constraints": f.other_constraints,
                        }
                        for f in c.formulas
                    ),
                    "reference_cell_counts": {
                        key: value
                        for key, value in zip(
                            ("published_full_cad", "published_tticad"),
                            _TABLE2_CELLS[c.number],
                            strict=True,
                        )
                        if value is not None
                    },
                },
            )
        )
    return tuple(out)
